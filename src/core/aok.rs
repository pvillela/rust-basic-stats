use super::{AltHyp, Ci, HypTestResult};

/// Enables coercion of `Result<T, E>` to the underlying type `T`,
/// producing a suitable fallback output value instead of panicking in case of error.
///
/// Intended to be implemented for floating point numbers.
pub trait AokFloat {
    type Output: AokFloatFallback;

    fn aok(self) -> Self::Output;
}

/// Constructs a suitable fallback instance for an Output type of [`Noerr`].
///
/// For floating point numbers, the fallback value will typically be `NaN`.
pub trait AokFloatFallback {
    fn aok_fallback() -> Self;
}

impl<T, E> AokFloat for Result<T, E>
where
    T: AokFloatFallback,
{
    type Output = T;

    fn aok(self) -> Self::Output {
        self.unwrap_or_else(|_| T::aok_fallback())
    }
}

impl AokFloatFallback for f64 {
    // Returns `NaN` as a fallback value.
    fn aok_fallback() -> f64 {
        f64::NAN
    }
}

/// Enables coercion of `Result<T, E>` to the underlying type `T`,
/// producing a suitable fallback output value instead of panicking in case of error.
///
/// Intended to be used for types constructed from floating point numbers.
pub trait AokBasicStats {
    type Output: AokBasicStatsFallback;

    fn aok(self) -> Self::Output;
}

/// Constructs a suitable fallback instance for an Output type of [`Noerr`].
///
/// For types constructed from floating point numbers, the fallback value will typically be a
/// value constructed with `NaN` fields.
pub trait AokBasicStatsFallback {
    fn aok_fallback() -> Self;
}

impl<T, E> AokBasicStats for Result<T, E>
where
    T: AokBasicStatsFallback,
{
    type Output = T;

    fn aok(self) -> Self::Output {
        self.unwrap_or_else(|_| T::aok_fallback())
    }
}

impl AokBasicStatsFallback for HypTestResult {
    // Returns an instance constructed with `NaN`s as a fallback value.
    fn aok_fallback() -> Self {
        HypTestResult::new(f64::NAN, f64::NAN, AltHyp::Ne)
    }
}

impl AokBasicStatsFallback for Ci {
    // Returns an instance constructed with `NaN`s as a fallback value.
    fn aok_fallback() -> Self {
        Ci(f64::NAN, f64::NAN)
    }
}

#[cfg(test)]
mod test {
    //! To simulate another package that also implements `AokFloat`.
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
            type Output: AokFloatFallback;

            fn aok(self) -> Self::Output;
        }
        pub trait AokFloatFallback {
            fn aok_fallback() -> Self;
        }

        impl<T, E> AokFloat for Result<T, E>
        where
            T: AokFloatFallback,
        {
            type Output = T;

            fn aok(self) -> Self::Output {
                self.unwrap_or_else(|_| T::aok_fallback())
            }
        }

        impl AokFloatFallback for f64 {
            // Returns `NaN` as a fallback value.
            fn aok_fallback() -> f64 {
                f64::NAN
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
            core::{AltHyp, AokBasicStats, SampleMoments},
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
            use crate::core::AokFloat;

            {
                println!("*** Ok scenario:");

                let alpha = 0.05;

                // Welch functions calls below return Ok prior to invocation of noerr().

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

                // Welch functions calls below return Err prior to invocation of noerr().

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

                // Welch functions calls below return Ok prior to invocation of noerr().

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

                // Welch functions calls below return Err prior to invocation of noerr().

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
