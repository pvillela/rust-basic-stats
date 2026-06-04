//! Functions to generate deterministic sampling for probability distributions.
//!
//! The functions return iterators that generate data using a simple patterns such that, for
//! appropriately large sample sizes, the generated sample passes the Kolmogorov-Smirnov test.
//!
//! The iterators provide a much leaner alternative to random number generators for uses that don't
//! require as much randomness.
//!
//! Although the infinite iterators generate samples that pass the Kolmogorov-Smirnov (KS) test, the
//! generation mechanism is such that the finite samples taken from them should be of size `2^k -1`.
//! The finite iterators provide flexibility to generate smaller sample sizes that also pass the
//! KS test.

mod finite_samp;
mod infinite_gen;

pub use finite_samp::*;
pub use infinite_gen::*;
