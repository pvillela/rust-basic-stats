#![doc = include_str!("lib.md")]

pub mod core;

#[cfg(feature = "normal")]
pub mod normal;

#[cfg(feature = "binomial")]
pub mod binomial;

#[cfg(feature = "wilcoxon")]
pub mod wilcoxon;

#[cfg(test)]
pub(crate) mod dev_utils;
