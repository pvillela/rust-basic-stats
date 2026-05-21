#![doc = include_str!("lib.md")]

pub mod core;

#[cfg(feature = "_aok_core")]
// impls in this module have additional feature gating
pub mod aok;

#[cfg(feature = "normal")]
pub mod normal;

#[cfg(feature = "binomial")]
pub mod binomial;

#[cfg(feature = "wilcoxon")]
pub mod wilcoxon;

#[cfg(feature = "_dev_utils")]
pub mod dev_utils;
