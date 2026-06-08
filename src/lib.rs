#![doc = include_str!("lib.md")]

pub mod core;

#[cfg(feature = "_aok_core")]
// impls in this module have additional feature gating
pub mod aok;

#[cfg(feature = "detm_samp")]
pub mod detm_samp;

#[cfg(feature = "rand_samp")]
pub mod rand_samp;

#[cfg(feature = "normal")]
pub mod normal;

#[cfg(feature = "binomial")]
pub mod binomial;

#[cfg(feature = "wilcoxon")]
pub mod wilcoxon;

#[cfg(feature = "_dev_utils")]
pub mod dev_utils;
