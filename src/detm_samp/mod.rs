//! Functions to generate deterministic sampling for probability distributions.
//!
//! The functions return iterators that generate data using a simple pattern such that, for
//! appropriately large sample sizes (typically low 20s), the generated sample passes the
//! Kolmogorov-Smirnov test.
//!
//! The iterators provide a much leaner alternative to random number generators for uses that don't
//! require as much randomness.

mod finite_samp;
mod infinite_gen;

pub use finite_samp::*;
pub use infinite_gen::*;
