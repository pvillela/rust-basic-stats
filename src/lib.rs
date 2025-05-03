//! A lightweight library that provides some basic parametric and non-parametric statistics and hypothesis tests.
//!
//! Some modules are gated by cargo features, as indicated in module documentation.
//! The feature **"all"** includes all features.

pub mod core;
pub mod error;
pub mod iter;
pub mod normal;

#[cfg(feature = "bernoulli")]
pub mod bernoulli;

#[cfg(feature = "wilcoxon")]
pub mod wilcoxon;

#[cfg(test)]
pub(crate) mod dev_utils;

pub mod binomial {
    pub use super::bernoulli::*;
}
