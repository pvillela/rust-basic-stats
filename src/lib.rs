//! A lightweight library that provides some basic parametric and non-parametric statistics and hypothesis tests.
//!
//! Some modules are gated by cargo features, as indicated in module documentation.
//! The feature **"all"** includes all features.

pub mod core;

#[cfg(feature = "normal")]
pub mod normal;

#[cfg(feature = "binomial")]
pub mod binomial;

#[cfg(feature = "wilcoxon")]
pub mod wilcoxon;

#[cfg(test)]
pub(crate) mod dev_utils;

// #[cfg(feature = "binomial")]
// pub mod bernoulli {
//     //! Alias for module [`crate::binomial`]. Gated by feature **bernoulli**, which includes feature **binomial**.
//     pub use super::binomial::*;
// }
