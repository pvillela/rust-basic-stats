//! Errors returned by functions in this library.

use std::{
    error::Error,
    fmt::{Debug, Display},
};

/// Signals a violation of an ordering requirement.
#[derive(Debug)]
pub struct OrderingError;

impl Display for OrderingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Debug::fmt(&self, f)
    }
}

impl Error for OrderingError {}
