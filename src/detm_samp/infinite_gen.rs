use crate::detm_samp::{BucketIter, infinite_bucket_iterator};

/// Returns an infinite iterator that samples from the
/// probability distribution given by the inverse CDF function `inv_cdf`.
///
/// The sampling covers the output range evenly throughout the generation process.
///
/// For sufficiently large sample sizes, generated sample passes the Kolmogorov-Smirnov test.
pub fn deterministic_gen<'a>(inv_cdf: impl Fn(f64) -> f64 + 'a) -> impl Iterator<Item = f64> + 'a {
    let unif_iter = uniform_01_detm_gen();
    unif_iter.map(inv_cdf)
}

/// Returns an infinite iterator that samples from the
/// uniform probability distribution in open interval `(0, 1)`.
///
/// The sampling covers the output range evenly throughout the generation process.
///
/// For sufficiently large sample sizes, generated sample passes the Kolmogorov-Smirnov test.
pub fn uniform_01_detm_gen() -> impl Iterator<Item = f64> {
    infinite_bucket_iterator(BucketIter::DEFAULT_N_BUCKETS, 1, false)
}

/// Returns an infinite iterator that samples from the
/// uniform probability distribution in open interval `(lo, hi)`, assuming `lo < hi`.
///
/// The sampling covers the output range evenly throughout the generation process.
///
/// For sufficiently large sample sizes, generated sample passes the Kolmogorov-Smirnov test.
pub fn uniform_detm_gen(lo: f64, hi: f64) -> impl Iterator<Item = f64> {
    uniform_01_detm_gen().map(move |v| (hi - lo) * v + lo)
}

#[derive(Debug)]
enum Side {
    Left,
    Right,
}

#[cfg(test)]
// cargo test --package basic_stats --lib --all-features -- detm_samp::infinite_gen::test --nocapture
mod test {
    use super::*;
    use old_statrs::distribution::{InverseCDF, Normal, Uniform};
    use statest::ks::KSTest;

    const EPSILON: f64 = 0.005;
    const SAMPLE_SIZE: usize = 50;

    #[test]
    fn test_uniform_01() {
        let iter = uniform_01_detm_gen().take(SAMPLE_SIZE);
        let v: Vec<f64> = iter.collect();
        let dist = Uniform::new(0.0, 1.0).unwrap();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&dist);
        assert!(1. - p < EPSILON, "1.-p={}, EPSILON={EPSILON}", 1. - p);
    }

    #[test]
    fn test_uniform() {
        let iter = uniform_detm_gen(1., 4.).take(SAMPLE_SIZE);
        let v: Vec<f64> = iter.collect();
        let dist = Uniform::new(1.0, 4.0).unwrap();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&dist);
        assert!(1. - p < EPSILON, "1.-p={}, EPSILON={EPSILON}", 1. - p);
    }

    #[test]
    fn test_normal() {
        let normal = Normal::new(0., 1.).unwrap();
        let iter = deterministic_gen(|x| normal.inverse_cdf(x)).take(SAMPLE_SIZE);
        let v: Vec<f64> = iter.collect();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&normal);
        assert!(1. - p < EPSILON, "1.-p={}, EPSILON={EPSILON}", 1. - p);
    }
}
