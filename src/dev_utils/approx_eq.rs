//! Provides approximate equality for floating point types

use std::ops::{Add, Div, Mul, Sub};

pub trait ApproxEq {
    fn approx_eq(self, other: Self, epsilon: Self) -> bool;
    fn abs_rel_diff(self, other: Self, epsilon: Self) -> Self;
    fn rel_approx_eq(self, other: Self, epsilon: Self) -> bool;
    fn round_to(self, sig_decimals: u8) -> Self;
}

trait AbsPowiRound10 {
    fn abs(self) -> Self;
    fn powi(self, n: i32) -> Self;
    fn round(self) -> Self;
    fn ten() -> Self;
}

impl<T> ApproxEq for T
where
    T: Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>
        + AbsPowiRound10
        + PartialOrd
        + PartialEq
        + Copy,
{
    fn approx_eq(self, other: Self, epsilon: Self) -> bool {
        if self == other {
            // Covers infinity case
            true
        } else {
            (self - other).abs() < epsilon
        }
    }

    fn abs_rel_diff(self, other: Self, epsilon: Self) -> Self {
        let abs_diff = (self - other).abs();
        let abs_sum = self.abs() + other.abs();

        if abs_sum < epsilon {
            // case where values are too close to zero
            abs_diff // which is always <= abs_sum
        } else {
            // normal case: abs_diff / abs_mean
            (abs_diff + abs_diff) / abs_sum
        }
    }

    fn rel_approx_eq(self, other: Self, epsilon: Self) -> bool {
        // handle case when both are zero
        if self == other {
            return true;
        }
        let rel_diff = self.abs_rel_diff(other, epsilon);
        let zero = self - self;
        zero.approx_eq(rel_diff, epsilon)
    }

    fn round_to(self, sig_decimals: u8) -> Self {
        let pow = Self::ten().powi(sig_decimals as i32);
        (self * pow).round() / pow
    }
}

impl AbsPowiRound10 for f32 {
    fn abs(self) -> Self {
        f32::abs(self)
    }

    fn powi(self, n: i32) -> Self {
        f32::powi(self, n)
    }

    fn round(self) -> Self {
        f32::round(self)
    }

    fn ten() -> Self {
        10.
    }
}

impl AbsPowiRound10 for f64 {
    fn abs(self) -> Self {
        f64::abs(self)
    }

    fn powi(self, n: i32) -> Self {
        f64::powi(self, n)
    }

    fn round(self) -> Self {
        f64::round(self)
    }

    fn ten() -> Self {
        10.
    }
}

#[macro_use]
mod macros {
    #[macro_export]
    macro_rules! approx_eq {
        ($a:expr, $b:expr, $epsilon:expr $(,)?) => {
            if !$crate::dev_utils::ApproxEq::approx_eq($a, $b, $epsilon) {
                panic!(
                    "assertion for approximate equality failed: left={}, right={}, epsilon={})",
                    $a, $b, $epsilon
                );
            }
        };
    }

    #[macro_export]
    macro_rules! rel_approx_eq {
        ($a:expr, $b:expr, $epsilon:expr $(,)?) => {
            if !$crate::dev_utils::ApproxEq::rel_approx_eq($a, $b, $epsilon) {
                panic!(
                    "assertion for relative approximate equality failed: left={}, right={}, epsilon={})",
                    $a, $b, $epsilon
                );
            }
        };
    }
}

#[cfg(test)]
mod test {
    use super::ApproxEq;

    #[test]
    fn test_approx_eq() {
        {
            let w: f32 = 123.444;
            let x: f32 = 123.44444;
            let y: f32 = 123.44454;
            let z: f32 = 123.44455;
            let epsilon: f32 = 0.0001;

            assert!(x.approx_eq(y, epsilon), "x must be approx_eq to y");
            assert!(!x.approx_eq(z, epsilon), "x must not be approx_eq to z");

            assert!(f32::INFINITY.approx_eq(f32::INFINITY, epsilon));
            assert!((-f32::INFINITY).approx_eq(-f32::INFINITY, epsilon));

            assert_eq!(w, x.round_to(3), "w must equal x.round_to(4)");
            assert_ne!(w, y.round_to(3), "w must not equal y.round_to(4)");
            assert_ne!(w, z.round_to(3), "w must not equal z.round_to(4)");
        }

        {
            let x: f32 = 100_000.;
            let y: f32 = 100_009.;
            let z: f32 = 100_020.;
            let epsilon: f32 = 0.0001;

            assert!(x.rel_approx_eq(y, epsilon), "x must be rel_approx_eq to y");
            assert!(
                !x.rel_approx_eq(z, epsilon),
                "x must not be rel_approx_eq to z"
            );
        }

        {
            let w: f64 = 123.4444;
            let x: f64 = 123.444444;
            let y: f64 = 123.444454;
            let z: f64 = 123.444455;
            let epsilon: f64 = 0.00001;

            assert!(x.approx_eq(y, epsilon), "x must be approx_eq to y");
            assert!(!x.approx_eq(z, epsilon), "x must not be approx_eq to z");

            assert!(f64::INFINITY.approx_eq(f64::INFINITY, epsilon));
            assert!((-f64::INFINITY).approx_eq(-f64::INFINITY, epsilon));

            assert_eq!(w, x.round_to(4), "w must equal x.round_to(4)");
            assert_ne!(w, y.round_to(4), "w must not equal y.round_to(4)");
            assert_ne!(w, z.round_to(4), "w must not equal z.round_to(4)");

            approx_eq!(x, y, epsilon);
        }

        {
            let x: f64 = 200_000.0;
            let y: f64 = 200_001.5;
            let z: f64 = 200_003.0;
            let epsilon: f64 = 0.00001;

            assert!(x.rel_approx_eq(y, epsilon), "x must be rel_approx_eq to y");
            assert!(
                !x.rel_approx_eq(z, epsilon),
                "x must not be rel_approx_eq to z"
            );

            rel_approx_eq!(x, y, epsilon);
        }
    }
}
