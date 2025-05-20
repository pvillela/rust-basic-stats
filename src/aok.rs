//! Supports the coercion of `Result<Value, E>` to the underlying type `Value`,
//! producing a suitable fallback output value instead of panicking in case of error.
//!
//! Traits [`AokFloat`] and [`AokBasicStats`] are separate to enable them to be selectively
//! imported. [`AokBasicStats`] provides a blanket implementation and we want to avoid namespace pollution
//! with blanket implementations for built-in types like `f64`.
//!
//! This module is NOT included by default. Inclusion of this module is gated by feature "**aok**".
//!
//! # Example
//!
//! Requires features **`aok`** and **`normal`**.
//! ```
#![doc = include_str!("../examples/aok.rs")]
//! ```

use crate::core::{AltHyp, Ci, HypTestResult};

/// Enables coercion of `Result<T, E>` to the underlying type `T`,
/// producing a suitable fallback output value instead of panicking in case of error.
///
/// Intended to be implemented for floating point numbers.
pub trait AokFloat {
    type Value;

    /// Returns the underlying value of a `Result`, without panicking.
    ///
    /// If the source result is an error, this method returns `NaN`.
    fn aok(self) -> Self::Value;
}

impl<E> AokFloat for Result<f64, E> {
    type Value = f64;

    fn aok(self) -> Self::Value {
        self.unwrap_or(f64::NAN)
    }
}

/// Enables coercion of `Result<T, E>` to the underlying type `T`,
/// producing a suitable fallback output value instead of panicking in case of error.
///
/// Intended to be used for types constructed from floating point numbers.
pub trait AokBasicStats {
    type Value: AokBasicStatsValue;

    fn aok(self) -> Self::Value;
}

/// Constructs a suitable fallback instance for the `Value` type of [`AokBasicStats`].
///
/// For types constructed from floating point numbers, the fallback value will typically be a
/// value constructed with `NaN` fields.
pub trait AokBasicStatsValue {
    /// Returns a suitable fallback value when the source result is an error.
    fn aok_fallback() -> Self;

    /// Returns `true` if the fallback value originated from an error.
    fn is_tainted(&self) -> bool;

    /// Returns `true` if the source result was `Ok`.
    fn is_untainted(&self) -> bool {
        !self.is_tainted()
    }
}

impl<T, E> AokBasicStats for Result<T, E>
where
    T: AokBasicStatsValue,
{
    type Value = T;

    /// Returns the underlying value of a `Result`, without panicking.
    ///
    /// If the source result is an error, this method returns a suitably constructed fallback value.
    fn aok(self) -> Self::Value {
        self.unwrap_or_else(|_| T::aok_fallback())
    }
}

impl AokBasicStatsValue for HypTestResult {
    fn aok_fallback() -> Self {
        HypTestResult::new(f64::NAN, f64::NAN, AltHyp::Ne)
    }

    fn is_tainted(&self) -> bool {
        self.p().is_nan() || self.alpha().is_nan()
    }
}

impl AokBasicStatsValue for Ci {
    fn aok_fallback() -> Self {
        Ci(f64::NAN, f64::NAN)
    }

    fn is_tainted(&self) -> bool {
        self.0.is_nan() || self.1.is_nan()
    }
}

#[cfg(test)]
#[cfg(feature = "normal")]
mod test {
    //! To simulate another package that also implements `AokFloat`.

    use crate::aok::AokBasicStats;

    mod another {
        use std::{
            error::Error,
            fmt::{Debug, Display},
        };

        #[derive(Debug)]
        pub struct AnotherError;

        impl Display for AnotherError {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                Debug::fmt(&self, f)
            }
        }

        impl Error for AnotherError {}

        #[derive(Debug)]
        pub struct X {
            pub x: f64,
        }

        impl X {
            pub fn new(x: f64) -> Self {
                X { x }
            }

            pub fn div(self, y: f64) -> Result<Self, AnotherError> {
                if y != 0. {
                    Ok(X::new(self.x / y))
                } else {
                    Err(AnotherError)
                }
            }
        }

        pub trait AokFloat {
            type Output;

            fn aok(self) -> Self::Output;
        }

        impl<E> AokFloat for Result<f64, E> {
            type Output = f64;

            fn aok(self) -> Self::Output {
                self.unwrap_or(f64::NAN)
            }
        }

        pub trait AokAnother {
            type Output: AokAnotherFallback;

            fn aok(self) -> Self::Output;
        }

        pub trait AokAnotherFallback {
            fn aok_fallback() -> Self;
        }

        impl<T, E> AokAnother for Result<T, E>
        where
            T: AokAnotherFallback,
        {
            type Output = T;

            fn aok(self) -> Self::Output {
                self.unwrap_or_else(|_| T::aok_fallback())
            }
        }

        impl AokAnotherFallback for X {
            // Returns an instance constructed with `NaN`s as a fallback value.
            fn aok_fallback() -> Self {
                X { x: f64::NAN }
            }
        }
    }

    #[test]
    fn test_aok() {
        use crate::{
            core::{AltHyp, SampleMoments},
            normal::{welch_alt_hyp_ci, welch_p},
        };

        let x = [14., 15., 15., 15., 16., 18., 22., 23., 24., 25., 25.];
        let y = [
            10., 12., 14., 15., 18., 22., 24., 27., 31., 33., 34., 34., 34.,
        ];

        let moments_x = SampleMoments::from_slice(&x);
        let moments_y = SampleMoments::from_slice(&y);
        let alt_hyp = AltHyp::Gt;

        {
            use crate::aok::AokFloat;

            {
                println!("*** Ok scenario:");

                let alpha = 0.05;

                // Welch function calls below return Ok prior to invocation of noerr().

                let p = welch_p(&moments_x, &moments_y, alt_hyp).aok();
                println!("p={p}");
                let ci = welch_alt_hyp_ci(&moments_x, &moments_y, alt_hyp, alpha).aok();
                println!("ci={ci:?}");

                assert!(p.is_finite());
                assert!(ci.0.is_finite());
            }

            {
                println!("*** Err scenario:");

                let alpha = 1.0;

                // Welch function calls below return Err prior to invocation of noerr().

                let p = welch_p(&moments_x, &SampleMoments::default(), alt_hyp).aok();
                println!("p={p}");
                let ci = welch_alt_hyp_ci(&moments_x, &moments_y, alt_hyp, alpha).aok();
                println!("ci={ci:?}");

                assert!(p.is_nan());
                assert!(ci.0.is_nan());
            }
        }

        // To demonstrate the use of another implementation of `AokFloat` while sharing the same implementation
        // of `AokBasicStats`.
        {
            use another::{AokAnother, AokFloat, X};

            {
                println!("*** Ok scenario:");

                let alpha = 0.05;

                // Welch function calls below return Ok prior to invocation of noerr().

                let p = welch_p(&moments_x, &moments_y, alt_hyp).aok();
                println!("p={p}");
                let ci = welch_alt_hyp_ci(&moments_x, &moments_y, alt_hyp, alpha).aok();
                println!("ci={ci:?}");

                let x = X::new(1.);
                let y = x.div(2.).aok();
                println!("y={y:?}");

                assert!(p.is_finite());
                assert!(ci.0.is_finite());
                assert!(y.x.is_finite());
            }

            {
                println!("*** Err scenario:");

                let alpha = 1.0;

                // Welch function calls below return Err prior to invocation of noerr().

                let p = welch_p(&moments_x, &SampleMoments::default(), alt_hyp).aok();
                println!("p={p}");
                let ci = welch_alt_hyp_ci(&moments_x, &moments_y, alt_hyp, alpha).aok();
                println!("ci={ci:?}");

                let x = X::new(1.);
                let y = x.div(0.).aok();
                println!("y={y:?}");

                assert!(p.is_nan());
                assert!(ci.0.is_nan());
                assert!(y.x.is_nan());
            }
        }
    }
}
