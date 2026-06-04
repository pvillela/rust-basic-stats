/// Returns an infinite iterator that samples from the
/// probability distribution given by the inverse CDF function `inv_cdf`.
///
/// The sampling covers the output range evenly throughout the generation process.
///
/// For sample sizes of the form `2^k - 1`, where `k` is sufficiently large, the generated sample passes
/// the Kolmogorov-Smirnov test
pub fn deterministic_gen<'a>(inv_cdf: impl Fn(f64) -> f64 + 'a) -> impl Iterator<Item = f64> + 'a {
    let unif_iter = uniform_01_detm_gen();
    unif_iter.map(inv_cdf)
}

/// Returns an infinite iterator that samples from the
/// uniform probability distribution in open interval `(0, 1)`.
///
/// The sampling covers the output range evenly throughout the generation process.
///
/// For sample sizes of the form `2^k - 1`, where `k` is sufficiently large, the generated sample passes
/// the Kolmogorov-Smirnov test
pub fn uniform_01_detm_gen() -> impl Iterator<Item = f64> {
    BucketIter::new()
}

/// Returns an infinite iterator that samples from the
/// uniform probability distribution in open interval `(lo, hi)`, assuming `lo < hi`.
///
/// The sampling covers the output range evenly throughout the generation process.
///
/// For sample sizes of the form `2^k - 1`, where `k` is sufficiently large, the generated sample passes
/// the Kolmogorov-Smirnov test
///
/// If `lo > hi` then the sample will be in the interval `(hi, lo)`.
/// If `lo == hi` then all samples will be equal to `lo`.
pub fn uniform_detm_gen(lo: f64, hi: f64) -> impl Iterator<Item = f64> {
    uniform_01_detm_gen().map(move |v| (hi - lo) * v + lo)
}

#[derive(Debug)]
enum Side {
    Left,
    Right,
}

#[derive(Debug)]
struct BucketIter {
    k: u64,
    range: u64,
    side: Side,
    granule: f64,
    last_value: f64,
    lowest_value_in_range: f64,
    items_generated: u64,
}

impl BucketIter {
    const MIDDLE: f64 = 0.5;

    fn new() -> Self {
        Self {
            k: 1,
            range: 1,
            side: Side::Left,
            granule: 0.0,
            last_value: f64::NAN,
            lowest_value_in_range: f64::NAN,
            items_generated: 0,
        }
    }

    fn increment_range(&mut self) {
        self.range *= 2;
        self.granule = 1.0 / (self.range * 2) as f64;
        self.k = 1;
        self.lowest_value_in_range = f64::INFINITY;
    }
}

impl Iterator for BucketIter {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.range == 1 {
            self.increment_range();
            self.items_generated += 1;
            return Some(Self::MIDDLE);
        }

        if self.k > self.range {
            println!("*** old struct={self:?}");
            self.increment_range();
            println!("*** new struct={self:?}");
        }

        let sign = match self.side {
            Side::Left => -1.0,
            Side::Right => 1.0,
        };

        let res = Self::MIDDLE + sign * self.k as f64 * self.granule;
        self.last_value = res;
        self.lowest_value_in_range = self.lowest_value_in_range.min(res);
        self.items_generated += 1;

        match self.side {
            Side::Left => {
                self.side = Side::Right;
            }
            Side::Right => {
                self.side = Side::Left;
                self.k += 2; // skip previously generated values
            }
        }

        Some(res)
    }
}

#[cfg(test)]
// cargo test --package basic_stats --lib --all-features -- detm_samp::test --nocapture
mod test {
    use super::*;
    use old_statrs::distribution::{InverseCDF, Normal, Uniform};
    use statest::ks::KSTest;

    const EPSILON: f64 = 0.005;
    const SAMPLE_SIZE: usize = 255; // `= 2_usize.pow(8) - 1`

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
