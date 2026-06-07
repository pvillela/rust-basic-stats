use crate::detm_samp::{finite_bucket_iterator, max_sqrt_divisor_no_greater_than};

/// Generates a deterministic sample of size `samp_size` for the
/// probability distribution given by the inverse CDF function `inv_cdf`.
///
/// The sample covers the output range evenly throughout the generation process.
///
/// For sufficiently large `samp_size`, the generated sample passes the Kolmogorov-Smirnov test
pub fn deterministic_samp<'a>(
    inv_cdf: impl Fn(f64) -> f64 + 'a,
    samp_size: usize,
) -> impl Iterator<Item = f64> + 'a {
    let unif_iter = uniform_01_detm_samp(samp_size);
    unif_iter.map(inv_cdf)
}

/// Generates a deterministic sample of size `samp_size` for the
/// uniform probability distribution in open interval `(0, 1)`.
///
/// The sample covers the output range evenly throughout the generation process.
///
/// For sufficiently large `samp_size`, the generated sample passes the Kolmogorov-Smirnov test
pub fn uniform_01_detm_samp(samp_size: usize) -> impl Iterator<Item = f64> {
    let samp_size2 = (samp_size + (1 - samp_size % 2) + 1) / 2;
    let n_buckets = max_sqrt_divisor_no_greater_than(samp_size2, 10);
    let bucket_size = samp_size2 / n_buckets;
    finite_bucket_iterator(n_buckets, bucket_size).take(samp_size)
}

/// Generates a deterministic sample of size `samp_size` for the
/// uniform probability distribution in open interval `(lo, hi)`, assuming `lo < hi`.
///
/// The sample covers the output range evenly throughout the generation process.
///
/// If `lo > hi` then the sample will be in the interval `(hi, lo)`.
/// If `lo == hi` then all samples will be equal to `lo`.
pub fn uniform_detm_samp(lo: f64, hi: f64, samp_size: usize) -> impl Iterator<Item = f64> {
    uniform_01_detm_samp(samp_size).map(move |v| (hi - lo) * v + lo)
}

#[allow(unused)]
#[cfg(feature = "_stash")]
pub mod stash {
    pub(crate) struct UnifIter {
        k: usize,
        i: usize,
    }

    impl Iterator for UnifIter {
        type Item = f64;

        fn next(&mut self) -> Option<Self::Item> {
            if self.i >= 2 * self.k * self.k - 1 {
                return None;
            }
            let item = uniform_observation(self.k, self.i);
            self.i += 1;
            Some(item)
        }
    }

    /// Generates the `i`-th observation for [`uniform_01_detm_samp`].
    ///
    /// The sample covers the output range evenly throughout the generation process.
    #[inline(always)]
    fn uniform_observation(k: usize, i: usize) -> f64 {
        let side = i % 2;
        let j = i / 2;
        let bucket_idx = j % k;
        let item_idx = j / k;
        let left_idx = bucket_idx * k + item_idx + 1;
        let idx = if side == 0 {
            left_idx
        } else {
            2 * k * k - left_idx
        };
        idx as f64 / (2 * k * k) as f64
    }
}

#[cfg(test)]
// cargo test --package basic_stats --lib --all-features -- detm_samp::finite_samp::test --nocapture
mod test {
    use super::*;
    use old_statrs::distribution::{InverseCDF, Normal, Uniform};
    use statest::ks::KSTest;

    const EPSILON: f64 = 0.005;
    const SAMP_SIZE: usize = 19;

    #[test]
    fn test_uniform_01() {
        let iter = uniform_01_detm_samp(SAMP_SIZE);
        let v: Vec<f64> = iter.collect();
        let dist = Uniform::new(0.0, 1.0).unwrap();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&dist);
        assert!(1. - p < EPSILON)
    }

    #[test]
    fn test_uniform() {
        let iter = uniform_detm_samp(1., 4., SAMP_SIZE);
        let v: Vec<f64> = iter.collect();
        let dist = Uniform::new(1.0, 4.0).unwrap();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&dist);
        assert!(1. - p < EPSILON)
    }

    #[test]
    fn test_normal() {
        let normal = Normal::new(0., 1.).unwrap();
        let iter = deterministic_samp(|x| normal.inverse_cdf(x), SAMP_SIZE);
        let v: Vec<f64> = iter.collect();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&normal);
        assert!(1. - p < EPSILON)
    }

    #[test]
    fn test_sample_size() {
        let samp_size = 200;
        let normal = Normal::new(0., 1.).unwrap();
        let iter = deterministic_samp(|x| normal.inverse_cdf(x), samp_size);
        let v: Vec<f64> = iter.collect();
        assert_eq!(v.len(), samp_size);
    }
}
