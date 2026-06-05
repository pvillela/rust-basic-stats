/// Returns an infinite iterator that samples from the
/// probability distribution given by the inverse CDF function `inv_cdf`.
///
/// The sampling covers the output range evenly throughout the generation process.
///
/// For sample sizes of the form `2^k - 1`, where `k` is sufficiently large, the generated sample passes
/// the Kolmogorov-Smirnov test
pub fn deterministic_gen<'a>(
    inv_cdf: impl Fn(f64) -> f64 + 'a,
    base_samp_size: u64,
) -> impl Iterator<Item = f64> + 'a {
    let unif_iter = uniform_01_detm_gen(base_samp_size);
    unif_iter.map(inv_cdf)
}

/// Returns an infinite iterator that samples from the
/// uniform probability distribution in open interval `(0, 1)`.
///
/// The sampling covers the output range evenly throughout the generation process.
///
/// For sample sizes of the form `2^k - 1`, where `k` is sufficiently large, the generated sample passes
/// the Kolmogorov-Smirnov test
pub fn uniform_01_detm_gen(base_samp_size: u64) -> impl Iterator<Item = f64> {
    BucketIter::new(true, base_samp_size)
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
pub fn uniform_detm_gen(lo: f64, hi: f64, base_samp_size: u64) -> impl Iterator<Item = f64> {
    uniform_01_detm_gen(base_samp_size).map(move |v| (hi - lo) * v + lo)
}

#[derive(Debug)]
enum Side {
    Left,
    Right,
}

#[derive(Debug)]
struct BucketIter {
    is_infinite: bool,
    samp_size: u64,
    left_base: f64,
    right_base: f64,
    k: u64,
    k_range: u64,
    side: Side,
    granule: f64,
    last_value: f64,
    lowest_value_in_range: f64,
    items_generated: u64,
}

impl BucketIter {
    fn new_empty() -> Self {
        Self {
            is_infinite: false,
            samp_size: 0,
            left_base: 0.0,
            right_base: 0.0,
            k: 0,
            k_range: 0,
            side: Side::Left,
            granule: 0.0,
            last_value: 0.0,
            lowest_value_in_range: f64::INFINITY,
            items_generated: 0,
        }
    }

    fn new(is_infinite: bool, samp_size: u64) -> Self {
        let mut it = Self::new_empty();
        it.update(is_infinite, samp_size);
        it
    }

    fn update(&mut self, is_infinite: bool, samp_size: u64) {
        let samp_sizef = samp_size as f64;
        let (left_base, right_base) = if samp_size % 2 == 1 {
            (0.5, 0.5)
        } else {
            let left_base = (1.0 / samp_sizef) * (samp_sizef / 2.0).ceil();
            let right_base = 1.0 - left_base;
            (left_base, right_base)
        };
        let k_range = samp_size / 2;
        let granule = 1.0 / samp_sizef;

        self.is_infinite = is_infinite;
        self.samp_size = samp_size;
        self.left_base = left_base;
        self.right_base = right_base;
        self.k = 1;
        self.k_range = k_range;
        self.granule = granule;
        self.last_value = f64::NAN;
    }

    fn increment_sample(&mut self) {
        self.update(self.is_infinite, self.samp_size * 2);
    }
}

impl Iterator for BucketIter {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.samp_size == 1 {
            let ret = 0.5;
            if self.is_infinite {
                self.increment_sample()
            }
            self.items_generated += 1;
            return Some(ret);
        }

        if self.k > self.k_range {
            println!("*** old struct={self:?}");
            if self.is_infinite {
                self.increment_sample();
                println!("*** new struct={self:?}");
            } else {
                return None;
            }
        }

        let (base, sign) = match self.side {
            Side::Left => (self.left_base, -1.0),
            Side::Right => (self.right_base, 1.0),
        };

        let res = base + sign * self.k as f64 * self.granule;
        self.last_value = res;
        self.lowest_value_in_range = self.lowest_value_in_range.min(res);
        self.items_generated += 1;

        match self.side {
            Side::Left => {
                self.side = Side::Right;
            }
            Side::Right => {
                self.side = Side::Left;
                if self.is_infinite {
                    // skip previously generated values
                    self.k += 2;
                } else {
                    self.k += 1
                }
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
        let iter = uniform_01_detm_gen(1).take(SAMPLE_SIZE);
        let v: Vec<f64> = iter.collect();
        let dist = Uniform::new(0.0, 1.0).unwrap();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&dist);
        assert!(1. - p < EPSILON, "1.-p={}, EPSILON={EPSILON}", 1. - p);
    }

    #[test]
    fn test_uniform() {
        let iter = uniform_detm_gen(1., 4., 1).take(SAMPLE_SIZE);
        let v: Vec<f64> = iter.collect();
        let dist = Uniform::new(1.0, 4.0).unwrap();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&dist);
        assert!(1. - p < EPSILON, "1.-p={}, EPSILON={EPSILON}", 1. - p);
    }

    #[test]
    fn test_normal() {
        let normal = Normal::new(0., 1.).unwrap();
        let iter = deterministic_gen(|x| normal.inverse_cdf(x), 1).take(SAMPLE_SIZE);
        let v: Vec<f64> = iter.collect();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&normal);
        assert!(1. - p < EPSILON, "1.-p={}, EPSILON={EPSILON}", 1. - p);
    }
}
