use super::{AltHyp, Ci, HypTestResult};

/// Enables a function that returns `Result<T, E>` to return the underlying type `T`,
/// producing a suitable fallback output value instead of panicking in case of error.
///
/// For floating point numbers and types constructed from floating point numbers, the
/// fallback value will typically be `NaN` or a value constructed with `NaN` fields.
pub trait Noerr {
    type Output: NoerrFallback;

    fn noerr(self) -> Self::Output;
}

/// Constructs a suitable fallback instance for an Output type of [`Noerr`].
///
/// For floating point numbers and types constructed from floating point numbers, the
/// fallback value will typically be `NaN` or a value constructed with `NaN` fields.
pub trait NoerrFallback {
    fn noerr_fallback() -> Self;
}

impl<T, E> Noerr for Result<T, E>
where
    T: NoerrFallback,
{
    type Output = T;

    fn noerr(self) -> Self::Output {
        self.unwrap_or_else(|_| T::noerr_fallback())
    }
}

impl NoerrFallback for f64 {
    // Returns `NaN` as a fallback value.
    fn noerr_fallback() -> f64 {
        f64::NAN
    }
}

impl NoerrFallback for HypTestResult {
    // Returns an instance constructed with `NaN`s as a fallback value.
    fn noerr_fallback() -> Self {
        HypTestResult::new(f64::NAN, f64::NAN, AltHyp::Ne)
    }
}

impl NoerrFallback for Ci {
    // Returns an instance constructed with `NaN`s as a fallback value.
    fn noerr_fallback() -> Self {
        Ci(f64::NAN, f64::NAN)
    }
}
