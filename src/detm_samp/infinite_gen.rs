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
    BucketIter::new_infinite(1)
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

#[derive(Debug)]
pub(crate) struct BucketIter {
    // is_initial_sample: bool,
    samp_size2: usize,
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
            // is_initial_sample: true,
            samp_size2: 0,
            bucket_size: 0,
            n_buckets: 0,
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

    pub(crate) fn new(samp_size2: usize) -> Self {
        let mut it = Self::new_empty();
        it.update(samp_size2);
        it
    }

    pub(crate) fn new_finite(samp_size2: usize) -> impl Iterator<Item = f64> {
        let samp_size = 2 * samp_size2 - 1;
        Self::new(samp_size2)
            .map(|(value, _)| value)
            .take(samp_size)
    }

    pub(crate) fn new_infinite(samp_size2: usize) -> impl Iterator<Item = f64> {
        Self::new(samp_size2).filter_map(|(value, flag)| if flag { None } else { Some(value) })
    }

    fn update(&mut self, samp_size2: usize) {
        self.samp_size2 = samp_size2;
        self.bucket_size = max_sqrt_divisor(samp_size2);
        self.n_buckets = samp_size2 / self.bucket_size;
        self.granule = 0.5 / samp_size2 as f64;
        self.bucket_idx = 1;
        self.samp_items_generated = 0;
        self.in_bucket_idx = 1;
        self.last_value = f64::NAN;
    }

    fn increase_sample(&mut self) {
        self.update(self.samp_size2 * 2);
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
        assert!(self.samp_size2 > 0, "samp_size2 must be > 0");

        let samp_size = self.samp_size2 * 2 - 1;
        if self.samp_items_generated >= samp_size {
            println!("*** old struct={self:?}");
            self.increase_sample();
            println!("*** new struct={self:?}");
        }

        if self.samp_items_generated == 0 {
            let ret = 0.5;
            self.samp_items_generated += 1;
            self.items_generated += 1;
            return Some((ret, self.filter_flag(0)));
        }

        let sign = match self.side {
            Side::Left => -1.0,
            Side::Right => 1.0,
        };

        let idx = (self.bucket_idx - 1) * self.bucket_size + self.in_bucket_idx;
        let res = Self::MIDPOINT + sign * idx as f64 * self.granule;
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

        self.samp_items_generated += 1;
        Some((res, self.filter_flag(idx)))
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
        let samp_size = 45; // fails with `samp_size = 44`
        // let iter = uniform_01_detm_gen(1).take(10);
        let iter = BucketIter::new(1).take(samp_size);
        let v: Vec<_> = iter.collect();
        println!("=== unfiltered v.len()={}, v={:?}", v.len(), v);
        let v: Vec<_> = BucketIter::new_infinite(1).take(samp_size).collect();
        println!("=== filtered v.len()={}, v={:?}", v.len(), v);
        let dist = Uniform::new(0.0, 1.0).unwrap();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&dist);
        assert!(1. - p < EPSILON, "1.-p={}, EPSILON={EPSILON}", 1. - p);
        // assert!(false);
    }

    #[test]
    // cargo test --package basic_stats --lib --all-features -- detm_samp::infinite_gen::test::test_buck_iter_new_17 --exact --nocapture --include-ignored
    fn test_buck_iter_new_17() {
        let samp_size2 = 17;
        let samp_size = 2 * samp_size2 - 1;
        // let iter = uniform_01_detm_gen(1).take(10);
        let iter = BucketIter::new(samp_size2).take(samp_size);
        let v: Vec<_> = iter.collect();
        println!("=== unfiltered v.len()={}, v={:?}", v.len(), v);
        let v: Vec<_> = BucketIter::new_infinite(samp_size2)
            .take(samp_size)
            .collect();
        println!("=== filtered v.len()={}, v={:?}", v.len(), v);
        let dist = Uniform::new(0.0, 1.0).unwrap();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&dist);
        assert!(1. - p < EPSILON, "1.-p={}, EPSILON={EPSILON}", 1. - p);
        // assert!(false);
    }

    #[test]
    // cargo test --package basic_stats --lib --all-features -- detm_samp::infinite_gen::test::test_buck_iter_new_17_delta --exact --nocapture --include-ignored
    fn test_buck_iter_new_17_delta() {
        let samp_size2 = 17;
        let delta = 30; // fails for delta = 29
        let samp_size = 2 * samp_size2 - 1 + delta;
        // let iter = uniform_01_detm_gen(1).take(10);
        let iter = BucketIter::new(samp_size2).take(samp_size);
        let v: Vec<_> = iter.collect();
        println!("=== unfiltered v.len()={}, v={:?}", v.len(), v);
        let v: Vec<_> = BucketIter::new_infinite(samp_size2)
            .take(samp_size)
            .collect();
        println!("=== filtered v.len()={}, v={:?}", v.len(), v);
        let dist = Uniform::new(0.0, 1.0).unwrap();
        let ks = KSTest::new(&v);
        let (p, _) = ks.ks1(&dist);
        assert!(1. - p < EPSILON, "1.-p={}, EPSILON={EPSILON}", 1. - p);
        // assert!(false);
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
