//! Functions to generate deterministic sampling for probability distributions.
//!
//! The functions return iterators that generate data using a simple pattern such that, for
//! appropriately large sample sizes, the generated sample passes the Kolmogorov-Smirnov test.
//!
//! The iterators provide a much leaner alternative to random number generators for uses that don't
//! require randomness.

mod bucket_iter;
mod finite_samp;

#[cfg(feature = "_stash")]
mod infinite_gen;

pub use finite_samp::*;

#[cfg(feature = "_stash")]
pub use infinite_gen::*;

pub(crate) use bucket_iter::*;
