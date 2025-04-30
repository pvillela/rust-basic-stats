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
