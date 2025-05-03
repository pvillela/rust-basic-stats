//! Errors returned by functions in this library.

use std::{
    error::Error,
    fmt::{Debug, Display},
};

use statrs::distribution::{BetaError, BinomialError};

/// Alias for results in this library.
pub(crate) type Result<V> = std::result::Result<V, StatsError>;

/// Signals an error in a function in this library.
#[derive(Debug)]
pub enum StatsError {
    /// Violation of an ordering requirement.
    Ordering,

    /// One or more arguments are outside the legal range.
    ArgOutOfRange(String),
}

impl Display for StatsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Debug::fmt(&self, f)
    }
}

impl Error for StatsError {}

impl From<BinomialError> for StatsError {
    fn from(value: BinomialError) -> Self {
        StatsError::ArgOutOfRange(value.to_string())
    }
}

impl From<BetaError> for StatsError {
    fn from(value: BetaError) -> Self {
        StatsError::ArgOutOfRange(value.to_string())
    }
}
