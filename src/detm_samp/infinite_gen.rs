/// Returns an infinite iterator that samples from the
/// probability distribution given by the inverse CDF function `inv_cdf`.
///
/// The sampling covers the output range evenly throughout the generation process.
///
/// For sample sizes of the form `2^k - 1`, where `k` is sufficiently large, the generated sample passes
/// the Kolmogorov-Smirnov test
pub fn deterministic_gen<'a>(
    inv_cdf: impl Fn(f64) -> f64 + 'a,
    samp_size2: usize,
) -> impl Iterator<Item = f64> + 'a {
    let unif_iter = uniform_01_detm_gen(samp_size2);
    unif_iter.map(inv_cdf)
}

/// Returns an infinite iterator that samples from the
/// uniform probability distribution in open interval `(0, 1)`.
///
/// The sampling covers the output range evenly throughout the generation process.
///
/// For sample sizes of the form `2^k - 1`, where `k` is sufficiently large, the generated sample passes
/// the Kolmogorov-Smirnov test
pub fn uniform_01_detm_gen(samp_size2: usize) -> impl Iterator<Item = f64> {
    BucketIter::new(true, samp_size2)
        .filter_map(|(value, flag)| if flag { None } else { Some(value) })
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
pub fn uniform_detm_gen(lo: f64, hi: f64, samp_size2: usize) -> impl Iterator<Item = f64> {
    uniform_01_detm_gen(samp_size2).map(move |v| (hi - lo) * v + lo)
}

#[derive(Debug)]
enum Side {
    Left,
    Right,
}

#[derive(Debug)]
pub(crate) struct BucketIter {
    is_initial_sample: bool,
    is_infinite: bool,
    samp_size2: usize,
    bucket_size: usize,
    n_buckets: usize,
    bucket_idx: usize,
    in_bucket_idx: usize,
    // k: usize,
    side: Side,
    granule: f64,
    last_value: f64,
    lowest_value_in_range: f64,
    items_generated: u64,
}

fn max_sqrt_divisor(n: usize) -> usize {
    if n == 0 {
        return 0; // Handle 0 case to avoid division by zero or panic
    }

    let limit = (n as f64).sqrt() as usize;

    for k in (1..=limit).rev() {
        if n % k == 0 {
            return k;
        }
    }
    1 // Fallback for prime numbers
}

impl BucketIter {
    const MIDPOINT: f64 = 0.5;

    fn new_empty() -> Self {
        Self {
            is_initial_sample: true,
            is_infinite: false,
            samp_size2: 0,
            bucket_size: 0,
            n_buckets: 0,
            bucket_idx: 0,
            in_bucket_idx: 0,
            // k: 0,
            side: Side::Left,
            granule: 0.0,
            last_value: 0.0,
            lowest_value_in_range: f64::INFINITY,
            items_generated: 0,
        }
    }

    pub(crate) fn new(is_infinite: bool, samp_size2: usize) -> Self {
        let mut it = Self::new_empty();
        it.update(samp_size2);
        it.is_infinite = is_infinite;
        it
    }

    fn update(&mut self, samp_size2: usize) {
        self.samp_size2 = samp_size2;
        self.bucket_size = max_sqrt_divisor(samp_size2);
        self.n_buckets = samp_size2 / self.bucket_size;
        let samp_size2f = samp_size2 as f64;
        self.granule = 1.0 / (samp_size2f * 2.0);
        self.bucket_idx = 1;
        self.in_bucket_idx = 1;
        self.last_value = f64::NAN;
    }

    fn increase_sample(&mut self) {
        self.update(self.samp_size2 * 2);
        self.is_infinite = true;
        self.is_initial_sample = false;
    }
}

impl Iterator for BucketIter {
    type Item = (f64, bool);

    fn next(&mut self) -> Option<Self::Item> {
        if self.samp_size2 == 0 {
            let ret = 0.5;
            if self.is_infinite {
                self.samp_size2 = 1;
                self.increase_sample()
            }
            self.items_generated += 1;
            return Some((ret, false));
        }

        if self.in_bucket_idx > self.bucket_size {
            println!("*** old struct={self:?}");
            if self.is_infinite {
                self.increase_sample();
                println!("*** new struct={self:?}");
            } else {
                return None;
            }
        }

        let sign = match self.side {
            Side::Left => -1.0,
            Side::Right => 1.0,
        };

        let k = (self.bucket_idx - 1) * self.bucket_size + self.in_bucket_idx;
        let res = Self::MIDPOINT + sign * k as f64 * self.granule;
        self.last_value = res;
        self.lowest_value_in_range = self.lowest_value_in_range.min(res);
        self.items_generated += 1;

        match self.side {
            Side::Left => {
                self.side = Side::Right;
            }
            Side::Right => {
                self.side = Side::Left;

                self.bucket_idx += 1;
                if self.bucket_idx > self.n_buckets {
                    self.bucket_idx = 1;
                    self.in_bucket_idx += 1;
                }
            }
        }

        Some((res, !self.is_initial_sample && k % 2 == 0))
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
    // cargo test --package basic_stats --lib --all-features -- detm_samp::infinite_gen::test::show_uniform_01 --exact --nocapture --include-ignored
    fn show_uniform_01() {
        // let iter = uniform_01_detm_gen(1).take(10);
        let iter = uniform_01_detm_gen(1).take(10);
        let v: Vec<f64> = iter.collect();
        println!("*** v.len()={}, v={:?}", v.len(), v);
        let dist = Uniform::new(0.0, 1.0).unwrap();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&dist);
        assert!(1. - p < EPSILON, "1.-p={}, EPSILON={EPSILON}", 1. - p);
        assert!(false);
    }
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
