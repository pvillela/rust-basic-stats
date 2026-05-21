//! Supports the coercion of `Result<Value, E>` to the underlying type `Value`,
//! producing a suitable fallback output value instead of panicking in case of error.
//!
//! This module is NOT included by default. Inclusion of this module is gated by feature "**aok**".
//!
//! # Example
//!
//! Requires features **`aok`** and **`normal`**.
//! ```
#![doc = include_str!("../examples/aok.rs")]
//! ```

/// Enables coercion of `Result<T, E>` to the underlying type `T`,
/// producing a suitable fallback output value instead of panicking in case of error.
///
/// Primarily to be used for flating point numbers or types with one or more fields that are floating
/// point numbers.
pub trait Aok: Sized {
    type Value: AokValue;

    /// Returns the underlying value type of a `Result`, without panicking.
    ///
    /// If the source result is an error, this method returns a fallback value ot type `Value`.
    fn aok(self) -> Self::Value {
        Self::Value::aok_fallback()
    }
}

/// Constructs a suitable fallback value for an implementing type, to be used instead of an error [`Result`].
///
/// # Type parameter
/// `D`: dummy type used for disambiguation so that a value type V can be the target of multiple implementations
/// of Aok and AokValue.
pub trait AokValue {
    /// Returns a suitable fallback value to be used instead of an error [`Result`].
    fn aok_fallback() -> Self;

    /// Returns `true` if the receiver is a value associated with an error.
    fn is_tainted(&self) -> bool;

    /// Returns `true` if the source result was `Ok`.
    fn is_untainted(&self) -> bool {
        !self.is_tainted()
    }
}

impl<T, E> Aok for Result<T, E>
where
    T: AokValue,
{
    type Value = T;

    fn aok(self) -> T {
        self.unwrap_or_else(|_| T::aok_fallback())
    }
}

#[cfg(feature = "aok_f64")]
impl AokValue for f64 {
    fn aok_fallback() -> Self {
        f64::NAN
    }

    fn is_tainted(&self) -> bool {
        self.is_nan()
    }
}

#[cfg(feature = "aok_stats")]
mod stats {
    use super::*;
    use crate::core::{AltHyp, Ci, HypTestResult};

    impl AokValue for HypTestResult {
        fn aok_fallback() -> Self {
            HypTestResult::new(f64::NAN, f64::NAN, AltHyp::Ne)
        }

        fn is_tainted(&self) -> bool {
            self.p().is_nan() || self.alpha().is_nan()
        }
    }

    impl AokValue for Ci {
        fn aok_fallback() -> Self {
            Ci(f64::NAN, f64::NAN)
        }

        fn is_tainted(&self) -> bool {
            self.0.is_nan() || self.1.is_nan()
        }
    }

    #[cfg(feature = "wilcoxon")]
    mod wilcoxon {
        use super::*;
        use crate::wilcoxon::RankSum;

        impl AokValue for RankSum {
            fn aok_fallback() -> Self {
                RankSum {
                    n_x: 0,
                    n_y: 0,
                    w: f64::NAN,
                    ties_sum_prod: 0,
                }
            }

            fn is_tainted(&self) -> bool {
                self.w.is_nan()
            }
        }
    }
}

#[cfg(test)]
#[cfg(feature = "normal")]
mod test {
    use crate::{
        aok::{Aok, AokValue},
        core::{AltHyp, SampleMoments},
    };

    #[test]
    fn test_aok() {
        let x = [14., 15., 15., 15., 16., 18., 22., 23., 24., 25., 25.];
        let y = [
            10., 12., 14., 15., 18., 22., 24., 27., 31., 33., 34., 34., 34.,
        ];

        let moments_x = SampleMoments::from_slice(&x);
        let moments_y = SampleMoments::from_slice(&y);
        let alt_hyp = AltHyp::Gt;

        #[cfg(feature = "aok_f64")]
        {
            use crate::normal::welch_p;
            {
                println!("*** Ok scenario:");

                // Welch function calls below return Ok prior to invocation of aok().
                let p = welch_p(&moments_x, &moments_y, 0., alt_hyp).aok();
                println!("p={p}");

                assert!(p.is_untainted());
                assert!(p.is_finite());
            }

            {
                println!("*** Err scenario:");

                // Welch function call below return Err prior to invocation of aok().
                let p = welch_p(&moments_x, &SampleMoments::default(), 0., alt_hyp).aok();
                println!("p={p}");

                assert!(p.is_tainted());
                assert!(p.is_nan());
            }
        }

        #[cfg(feature = "aok_stats")]
        {
            use crate::normal::welch_alt_hyp_ci;
            {
                println!("*** Ok scenario:");

                let alpha = 0.05;

                // Welch function call below return Ok prior to invocation of aok().
                let ci = welch_alt_hyp_ci(&moments_x, &moments_y, alt_hyp, alpha).aok();
                println!("ci={ci:?}");

                assert!(ci.is_untainted());
                assert!(ci.0.is_finite());
            }

            {
                println!("*** Err scenario:");

                let alpha = 1.0;

                // Welch function call below return Err prior to invocation of aok().
                let ci = welch_alt_hyp_ci(&moments_x, &moments_y, alt_hyp, alpha).aok();
                println!("ci={ci:?}");

                assert!(ci.is_tainted());
                assert!(ci.0.is_nan());
            }
        }
    }
}
