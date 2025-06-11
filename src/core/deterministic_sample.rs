/// Generates a deterministic sample of size `2*k*k - 1` for the
/// probability distribution given by the inverse CDF function `inv_cdf`.
///
/// The sample covers the output range evenly throughout the generation process.
pub fn deterministic_sample<'a>(
    inv_cdf: impl Fn(f64) -> f64 + 'a,
    k: u64,
) -> impl Iterator<Item = f64> + 'a {
    let unif_iter = deterministic_uniform_01_sample(k);
    unif_iter.map(inv_cdf)
}

/// Generates a deterministic sample of size `2*k*k - 1` for the
/// uniform probability distribution in open interval `(0, 1)`.
///
/// The sample covers the output range evenly throughout the generation process.
pub fn deterministic_uniform_01_sample(k: u64) -> impl Iterator<Item = f64> {
    UnifIter { k, i: 0 }
}

/// Generates a deterministic sample of size `2*k*k - 1` for the
/// uniform probability distribution in open interval `(lo, hi)`, assuming `lo < hi`.
///
/// The sample covers the output range evenly throughout the generation process.
///
/// If `lo > hi` then the sample will be in the interval `(hi, lo)`.
/// If `lo == hi` then all samples will be equal to `lo`.
pub fn deterministic_uniform_sample(lo: f64, hi: f64, k: u64) -> impl Iterator<Item = f64> {
    deterministic_uniform_01_sample(k).map(move |v| (hi - lo) * v + lo)
}

struct UnifIter {
    k: u64,
    i: u64,
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

/// Generates the `i`-th observation for [`deterministic_uniform_sample`].
///
/// The sample covers the output range evenly throughout the generation process.
#[inline(always)]
fn uniform_observation(k: u64, i: u64) -> f64 {
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

#[cfg(test)]
mod test {
    use super::*;
    use old_statrs::distribution::{InverseCDF, Normal, Uniform};
    use statest::ks::KSTest;

    const EPSILON: f64 = 0.005;

    #[test]
    fn test_uniform_01() {
        let iter = deterministic_uniform_01_sample(10);
        let v: Vec<f64> = iter.collect();
        let dist = Uniform::new(0.0, 1.0).unwrap();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&dist);
        assert!(1. - p < EPSILON)
    }

    #[test]
    fn test_uniform() {
        let iter = deterministic_uniform_sample(1., 4., 10);
        let v: Vec<f64> = iter.collect();
        let dist = Uniform::new(1.0, 4.0).unwrap();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&dist);
        assert!(1. - p < EPSILON)
    }

    #[test]
    fn test_normal() {
        let normal = Normal::new(0., 1.).unwrap();
        let iter = deterministic_sample(|x| normal.inverse_cdf(x), 10);
        let v: Vec<f64> = iter.collect();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&normal);
        assert!(1. - p < EPSILON)
    }
}
