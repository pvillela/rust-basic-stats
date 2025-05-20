#![doc = include_str!("lib.md")]

pub mod core;

#[cfg(feature = "aok")]
pub mod aok;

#[cfg(feature = "normal")]
pub mod normal;

#[cfg(feature = "binomial")]
pub mod binomial;

#[cfg(feature = "wilcoxon")]
pub mod wilcoxon;

#[doc(hidden)]
pub mod dev_utils;
