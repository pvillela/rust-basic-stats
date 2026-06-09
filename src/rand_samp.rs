//! For use in statistical simulations, ***not*** for cryptographic use.

use rand::{RngExt, SeedableRng, rngs::Xoshiro256PlusPlus};

/// Returns an infinite iterator that samples from the
/// probability distribution given by the inverse CDF function `inv_cdf`.
///
/// The sampling covers the output range evenly throughout the generation process.
///
/// For sufficiently large sample sizes, generated sample passes the Kolmogorov-Smirnov test.
pub fn random_gen<'a>(inv_cdf: impl Fn(f64) -> f64 + 'a) -> impl Iterator<Item = f64> + 'a {
    let unif_iter = uniform_01_rand_gen();
    unif_iter.map(inv_cdf)
}

/// Returns an infinite iterator that samples from the
/// uniform probability distribution in open interval `(0, 1)`.
///
/// The sampling covers the output range evenly throughout the generation process.
///
/// For sufficiently large sample sizes, generated sample passes the Kolmogorov-Smirnov test.
pub fn uniform_01_rand_gen() -> impl Iterator<Item = f64> {
    RandIter::new(RandIter::DEFAULT_SEED)
}

/// Returns an infinite iterator that samples from the
/// uniform probability distribution in open interval `(lo, hi)`, assuming `lo < hi`.
///
/// The sampling covers the output range evenly throughout the generation process.
///
/// For sufficiently large sample sizes, generated sample passes the Kolmogorov-Smirnov test.
pub fn uniform_rand_gen(lo: f64, hi: f64) -> impl Iterator<Item = f64> {
    uniform_01_rand_gen().map(move |v| (hi - lo) * v + lo)
}

/// Generates a deterministic sample of size `samp_size` for the
/// probability distribution given by the inverse CDF function `inv_cdf`.
///
/// The sample covers the output range evenly throughout the generation process.
///
/// For sufficiently large `samp_size`, the generated sample passes the Kolmogorov-Smirnov test
pub fn random_samp<'a>(
    inv_cdf: impl Fn(f64) -> f64 + 'a,
    samp_size: usize,
) -> impl Iterator<Item = f64> + 'a {
    let unif_iter = uniform_01_rand_samp(samp_size);
    unif_iter.map(inv_cdf)
}

/// Generates a deterministic sample of size `samp_size` for the
/// uniform probability distribution in open interval `(0, 1)`.
///
/// The sample covers the output range evenly throughout the generation process.
///
/// For sufficiently large `samp_size`, the generated sample passes the Kolmogorov-Smirnov test
pub fn uniform_01_rand_samp(samp_size: usize) -> impl Iterator<Item = f64> {
    uniform_01_rand_gen().take(samp_size)
}

/// Generates a deterministic sample of size `samp_size` for the
/// uniform probability distribution in open interval `(lo, hi)`, assuming `lo < hi`.
///
/// The sample covers the output range evenly throughout the generation process.
///
/// If `lo > hi` then the sample will be in the interval `(hi, lo)`.
/// If `lo == hi` then all samples will be equal to `lo`.
pub fn uniform_rand_samp(lo: f64, hi: f64, samp_size: usize) -> impl Iterator<Item = f64> {
    uniform_01_rand_samp(samp_size).map(move |v| (hi - lo) * v + lo)
}

struct RandIter {
    rng: Xoshiro256PlusPlus,
}

impl RandIter {
    const DEFAULT_SEED: u64 = 42;

    fn new(seed: u64) -> Self {
        Self {
            rng: Xoshiro256PlusPlus::seed_from_u64(seed),
        }
    }
}

impl Iterator for RandIter {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        // The random number generator guarantees the value is in the interval [0.0, 1.0).
        // We need to ensure 0.0 is not returned.
        let mut value: f64 = 0.0;
        while value == 0.0 {
            value = self.rng.random();
        }
        Some(value)
    }
}
