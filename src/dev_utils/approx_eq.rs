//! Provides approximate equality for floating point types

pub trait ApproxEq {
    fn approx_eq(self, other: Self, epsilon: Self) -> bool;
    fn round_to(self, n: u8) -> Self;
}

impl ApproxEq for f32 {
    fn approx_eq(self, other: Self, epsilon: Self) -> bool {
        (self - other).abs() <= epsilon
    }

    fn round_to(self, sig_decimals: u8) -> f32 {
        let pow = 10.0_f32.powi(sig_decimals as i32);
        (self * pow).round() / pow
    }
}

impl ApproxEq for f64 {
    fn approx_eq(self, other: Self, epsilon: Self) -> bool {
        (self - other).abs() < epsilon
    }

    /// Rounds self to `sig_decimals` significant decimals.
    fn round_to(self, sig_decimals: u8) -> f64 {
        let pow = 10.0_f64.powi(sig_decimals as i32);
        (self * pow).round() / pow
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

            assert_eq!(w, x.round_to(3), "w must equal x.round_to(4)");
            assert_ne!(w, y.round_to(3), "w must not equal y.round_to(4)");
            assert_ne!(w, z.round_to(3), "w must not equal z.round_to(4)");
        }

        {
            let w: f64 = 123.4444;
            let x: f64 = 123.444444;
            let y: f64 = 123.444454;
            let z: f64 = 123.444455;
            let epsilon: f64 = 0.00001;

            assert!(x.approx_eq(y, epsilon), "x must be approx_eq to y");
            assert!(!x.approx_eq(z, epsilon), "x must not be approx_eq to z");

            assert_eq!(w, x.round_to(4), "w must equal x.round_to(4)");
            assert_ne!(w, y.round_to(4), "w must not equal y.round_to(4)");
            assert_ne!(w, z.round_to(4), "w must not equal z.round_to(4)");
        }
    }
}
