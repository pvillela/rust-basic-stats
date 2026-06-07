use log::{self, trace};

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
    BucketIter::new_infinite(BucketIter::DEFAULT_N_BUCKETS, 1, false)
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

pub(crate) fn max_sqrt_divisor_no_greater_than(n: usize, upper: usize) -> usize {
    if n == 0 {
        return 0; // Handle 0 case to avoid division by zero or panic
    }

    let limit = upper.min((n as f64).sqrt() as usize);

    for k in (1..=limit).rev() {
        if n % k == 0 {
            return k;
        }
    }
    1 // Fallback for prime numbers
}

#[derive(Debug)]
enum Side {
    Left,
    Right,
}

#[derive(Debug)]
pub(crate) struct BucketIter {
    // is_initial_sample: bool,
    // samp_size2: usize,
    bucket_size: usize,
    n_buckets: usize,
    bucket_idx: usize,
    in_bucket_idx: usize,
    samp_items_generated: usize,
    // k: usize,
    side: Side,
    granule: f64,
    last_value: f64,
    lowest_value_in_range: f64,
    items_generated: usize,
}

impl BucketIter {
    const MIDPOINT: f64 = 0.5;
    pub(crate) const DEFAULT_N_BUCKETS: usize = 8;

    fn new_empty() -> Self {
        Self {
            // is_initial_sample: true,
            // samp_size2: 0,
            n_buckets: 0,
            bucket_size: 0,
            bucket_idx: 0,
            samp_items_generated: 0,
            in_bucket_idx: 0,
            // k: 0,
            side: Side::Left,
            granule: 0.0,
            last_value: f64::NAN,
            lowest_value_in_range: f64::INFINITY,
            items_generated: 0,
        }
    }

    pub(crate) fn new(n_buckets: usize, bucket_size: usize) -> Self {
        let mut it = Self::new_empty();
        it.update(n_buckets, bucket_size);
        it
    }

    pub(crate) fn new_finite(n_buckets: usize, bucket_size: usize) -> impl Iterator<Item = f64> {
        let samp_size = 2 * n_buckets * bucket_size - 1;
        Self::new(n_buckets, bucket_size)
            .map(|(value, _)| value)
            .take(samp_size)
    }

    pub(crate) fn new_infinite(
        n_buckets: usize,
        bucket_size: usize,
        allow_duplicates: bool,
    ) -> impl Iterator<Item = f64> {
        Self::new(n_buckets, bucket_size).filter_map(move |(value, flag)| {
            if !allow_duplicates && flag {
                None
            } else {
                Some(value)
            }
        })
    }

    fn update(&mut self, n_buckets: usize, bucket_size: usize) {
        // self.samp_size2 = samp_size2;
        self.n_buckets = n_buckets;
        self.bucket_size = bucket_size;
        self.granule = 0.5 / (n_buckets * bucket_size) as f64;
        self.bucket_idx = 1;
        self.samp_items_generated = 0;
        self.in_bucket_idx = 1;
        self.last_value = f64::NAN;
    }

    fn increase_sample(&mut self) {
        self.update(self.n_buckets, self.bucket_size * 2);
        // self.is_initial_sample = false;
    }

    fn is_initial_sample(&self) -> bool {
        self.items_generated == self.samp_items_generated
    }

    fn filter_flag(&self, idx: usize) -> bool {
        !self.is_initial_sample() && idx % 2 == 0
    }
}

impl Iterator for BucketIter {
    type Item = (f64, bool);

    fn next(&mut self) -> Option<Self::Item> {
        assert!(
            self.n_buckets * self.bucket_size > 0,
            "`n_buckets` and `bucket_size` must be > 0"
        );

        let samp_size = self.n_buckets * self.bucket_size * 2 - 1;
        if self.samp_items_generated >= samp_size {
            trace!("old struct={self:?}");
            self.increase_sample();
            trace!("new struct={self:?}");
        }

        if self.samp_items_generated == 0 {
            let value = 0.5;
            self.samp_items_generated += 1;
            self.items_generated += 1;
            return Some((value, self.filter_flag(0)));
        }

        let sign = match self.side {
            Side::Left => -1.0,
            Side::Right => 1.0,
        };

        let idx = (self.bucket_idx - 1) * self.bucket_size + self.in_bucket_idx;
        let value = Self::MIDPOINT + sign * idx as f64 * self.granule;
        self.last_value = value;
        self.lowest_value_in_range = self.lowest_value_in_range.min(value);
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

        self.samp_items_generated += 1;
        Some((value, self.filter_flag(idx)))
    }
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
    // cargo test --package basic_stats --lib --all-features -- detm_samp::infinite_gen::test::test_buck_iter_new_1 --exact --nocapture --include-ignored
    fn test_buck_iter_new_1() {
        _ = env_logger::try_init();

        let samp_size = 45; // fails with `samp_size = 44`
        let iter = BucketIter::new(BucketIter::DEFAULT_N_BUCKETS, 1).take(samp_size);
        let v: Vec<_> = iter.collect();
        println!("=== unfiltered v.len()={}, v={:?}", v.len(), v);
        let v: Vec<_> = BucketIter::new_infinite(BucketIter::DEFAULT_N_BUCKETS, 1, false)
            .take(samp_size)
            .collect();
        println!("=== filtered v.len()={}, v={:?}", v.len(), v);
        let dist = Uniform::new(0.0, 1.0).unwrap();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&dist);
        assert!(1. - p < EPSILON, "1.-p={}, EPSILON={EPSILON}", 1. - p);
    }

    #[test]
    // cargo test --package basic_stats --lib --all-features -- detm_samp::infinite_gen::test::test_buck_iter_new_17 --exact --nocapture --include-ignored
    fn test_buck_iter_new_17() {
        let n_buckets = 17;
        let bucket_size = 1;
        let samp_size = 2 * n_buckets * bucket_size - 1;
        let iter = BucketIter::new(n_buckets, bucket_size).take(samp_size);
        let v: Vec<_> = iter.collect();
        println!("=== unfiltered v.len()={}, v={:?}", v.len(), v);
        let v: Vec<_> = BucketIter::new_infinite(n_buckets, bucket_size, false)
            .take(samp_size)
            .collect();
        println!("=== filtered v.len()={}, v={:?}", v.len(), v);
        let dist = Uniform::new(0.0, 1.0).unwrap();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&dist);
        assert!(1. - p < EPSILON, "1.-p={}, EPSILON={EPSILON}", 1. - p);
    }

    #[test]
    // cargo test --package basic_stats --lib --all-features -- detm_samp::infinite_gen::test::test_buck_iter_new_17_delta --exact --nocapture --include-ignored
    fn test_buck_iter_new_17_delta() {
        let n_buckets = 17;
        let bucket_size = 1;
        let delta = 30;
        let samp_size = 2 * n_buckets * bucket_size - 1 + delta;
        let iter = BucketIter::new(n_buckets, bucket_size).take(samp_size);
        let v: Vec<_> = iter.collect();
        println!("=== unfiltered v.len()={}, v={:?}", v.len(), v);
        let v: Vec<_> = BucketIter::new_infinite(n_buckets, bucket_size, false)
            .take(samp_size)
            .collect();
        println!("=== filtered v.len()={}, v={:?}", v.len(), v);
        let dist = Uniform::new(0.0, 1.0).unwrap();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&dist);
        assert!(1. - p < EPSILON, "1.-p={}, EPSILON={EPSILON}", 1. - p);
    }

    #[test]
    // cargo test --package basic_stats --lib --all-features -- detm_samp::infinite_gen::test::test_skip_take --exact --nocapture --include-ignored
    fn test_skip_take() {
        let n_buckets = 17;
        let bucket_size = 1;
        let delta = 30;
        let samp_size = 2 * n_buckets * bucket_size - 1 + delta;
        let iter = BucketIter::new(n_buckets, bucket_size).take(samp_size);
        let v: Vec<_> = iter.collect();
        println!("=== unfiltered v.len()={}, v={:?}", v.len(), v);
        let v: Vec<_> = BucketIter::new_infinite(n_buckets, bucket_size, false)
            .skip(samp_size)
            .take(samp_size)
            .collect();
        println!("=== filtered v.len()={}, v={:?}", v.len(), v);
        let dist = Uniform::new(0.0, 1.0).unwrap();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&dist);
        assert!(1. - p < EPSILON, "1.-p={}, EPSILON={EPSILON}", 1. - p);
    }

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
