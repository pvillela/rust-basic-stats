pub mod core;

#[cfg(test)]
pub(crate) mod dev_utils;

#[cfg(feature = "bernoulli")]
pub mod bernoulli;

pub mod normal;

#[cfg(feature = "wilcoxon")]
pub mod wilcoxon;
