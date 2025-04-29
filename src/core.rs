use core::f64;

/// Sample mean.
///
/// Returns `NaN` if  `n == 0`.
///
/// Arguments:
/// - `n` - sample size.
/// - `sum` - sample sum.
#[inline(always)]
pub fn sample_mean(n: u64, sum: f64) -> f64 {
    if n == 0 {
        return f64::NAN;
    }
    sum / n as f64
}

/// Sample's sum of squares of deviations from sample mean.
///
/// Returns `NaN` if  `n == 0`.
///
/// Arguments:
/// - `n` - sample size.
/// - `sum` - sample sum.
/// - `sum2` - sample sum of squares.
#[inline(always)]
pub fn sample_sum2_deviations(n: u64, sum: f64, sum2: f64) -> f64 {
    if n <= 0 {
        return f64::NAN;
    }
    sum2 - sum.powi(2) / n as f64
}

/// Sample variance.
///
/// Returns `NaN` if  `n ≤ 1`.
///
/// Arguments:
/// - `n` - sample size.
/// - `sum` - sample sum.
/// - `sum2` - sample sum of squares.
#[inline(always)]
pub fn sample_var(n: u64, sum: f64, sum2: f64) -> f64 {
    if n <= 0 {
        return f64::NAN;
    }
    sample_sum2_deviations(n, sum, sum2) / (n as f64 - 1.)
}

/// Sample standard deviation.
///
/// Returns `NaN` if  `n ≤ 1`.
///
/// Arguments:
/// - `n` - sample size.
/// - `sum` - sample sum.
/// - `sum2` - sample sum of squares.
#[inline(always)]
pub fn sample_stdev(n: u64, sum: f64, sum2: f64) -> f64 {
    sample_var(n, sum, sum2).sqrt()
}

/// Contains a sample's size and first and second moments (sum and sum of squares),
/// optionally also including the sample min and max.
///
/// This is a lightweight data structure used by several functions and methods in this library.
pub struct SampleMoments {
    count: u64,
    sum: f64,
    sum2: f64,
    min: f64,
    max: f64,
}

impl SampleMoments {
    /// Instantiates `Self` with sample size, first moment, and second moment equal to `0`, and
    /// with min and max values of `NaN`.
    pub fn new_empty() -> Self {
        Self::new(0, 0., 0.)
    }

    /// Instantiates `Self` with given values for sample size, first moment, and second moment,
    /// and with min and max values defaulting to `NaN`.
    pub fn new(count: u64, sum: f64, sum2: f64) -> Self {
        Self::new_with_min_max(count, sum, sum2, f64::NAN, f64::NAN)
    }

    /// Instantiates `Self`, with given values for sample size, first and second moments, min, and max.
    pub fn new_with_min_max(count: u64, sum: f64, sum2: f64, min: f64, max: f64) -> Self {
        Self {
            count,
            sum,
            sum2,
            min,
            max,
        }
    }

    /// Updates `self` by accumulating an additional value.
    pub fn collect_value(&mut self, value: f64) {
        self.count += 1;
        self.sum += value;
        self.sum2 += value * value;
        self.min = value.min(self.min);
        self.max = value.max(self.max);
    }

    /// Instatiates `Self` from an iterator.
    pub fn from_iterator(mut dataset: impl Iterator<Item = f64>) -> Self {
        let mut moments = SampleMoments::new_empty();
        loop {
            match dataset.next() {
                Some(v) => moments.collect_value(v),
                None => break,
            }
        }
        moments
    }

    /// Instatiates `Self` from a slice.
    pub fn from_slice(dataset: &[f64]) -> Self {
        let iter = dataset.iter().cloned();
        Self::from_iterator(iter)
    }

    /// Sample size as integer.
    pub fn n(&self) -> u64 {
        self.count
    }

    /// Sample size as floating point.
    pub fn nf(&self) -> f64 {
        self.count as f64
    }

    /// Sample sum (first moment).
    pub fn sum(&self) -> f64 {
        self.sum
    }

    /// Sample mean.
    ///
    /// Returns `NaN` if  `self.n() == 0`.
    pub fn mean(&self) -> f64 {
        self.sum / self.nf()
    }

    /// Sample sum of squares (second moment).
    pub fn sum2(&self) -> f64 {
        self.sum2
    }

    /// Sample's sum of squares of deviations from sample mean.
    ///
    /// Returns `NaN` if  `self.n() == 0`.
    pub fn sum2_deviations(&self) -> f64 {
        sample_sum2_deviations(self.n(), self.sum, self.sum2)
    }

    /// Sample variance.
    ///
    /// Returns `NaN` if  `self.n() ≤ 1`.
    pub fn var(&self) -> f64 {
        sample_var(self.n(), self.sum, self.sum2)
    }

    /// Sample standard deviation.
    ///
    /// Returns `NaN` if  `self.n() ≤ 1`.
    pub fn stdev(&self) -> f64 {
        sample_stdev(self.n(), self.sum, self.sum2)
    }

    /// Sample's minimum value.
    pub fn min(&self) -> f64 {
        self.min
    }

    /// Sample's maximum value.
    pub fn max(&self) -> f64 {
        self.max
    }
}

impl Default for SampleMoments {
    fn default() -> Self {
        Self::new_empty()
    }
}

/// Alternative statistical hypothesis to the null hypothesis of equality.
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum AltHyp {
    /// Less than
    Lt,
    /// Greater than
    Gt,
    /// Not equal
    Ne,
}

/// Statistical test hypothesis.
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Hyp {
    /// Null hypothesis of equality.
    Null,
    /// Alternative hypothesis.
    Alt(AltHyp),
}

impl Hyp {
    pub fn alt_hyp(&self) -> AltHyp {
        match self {
            Self::Null => AltHyp::Ne,
            Self::Alt(h) => *h,
        }
    }
}

/// Result of a statistical hypothesis test with a null hypothesis of equality.
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct HypTestResult {
    p: f64,
    alpha: f64,
    alt_hyp: AltHyp,
    accepted: Hyp,
}

impl HypTestResult {
    /// Creates a new instance of `Self`.
    ///
    /// Arguments:
    /// - `p` - the 'p' value for the test result.
    /// - `alpha` - determines the confidence level `(1-alpha)`.
    /// - `alt_hyp` - the alternative hypothesis.
    /// - `accepted` - the accepted hypothesis (null or alternative).
    pub fn new(p: f64, alpha: f64, alt_hyp: AltHyp) -> HypTestResult {
        Self {
            p,
            alpha,
            alt_hyp,
            accepted: if p >= alpha {
                Hyp::Null
            } else {
                Hyp::Alt(alt_hyp)
            },
        }
    }

    /// The 'p' value of the test result.
    pub fn p(&self) -> f64 {
        self.p
    }

    /// The 'alpha' for the test; determines the confidence level `(1-alpha)`.
    pub fn alpha(&self) -> f64 {
        self.alpha
    }

    /// The alternative hypothesis for the test.
    pub fn alt_hyp(&self) -> AltHyp {
        self.alt_hyp
    }

    /// The hypothesis accepted by the test.
    pub fn accepted(&self) -> Hyp {
        self.accepted
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
/// Represents the position of a value with respect to a confidence interval.
pub enum PositionWrtCi {
    /// The value is lower than the low end of the confidence interval.
    Below,
    /// The value is inside the confidence interval.
    In,
    /// The value is higher than the high end of the confidence interval.
    Above,
}

#[derive(Debug, PartialEq, Clone, Copy)]
/// Confidence interval.
pub struct Ci(
    /// Low end of interval.
    pub f64,
    /// High end of interval.
    pub f64,
);

impl Ci {
    /// Returns the position of `value` with respect to `self`.
    pub fn position_of(&self, value: f64) -> PositionWrtCi {
        match value {
            _ if value <= self.0 => PositionWrtCi::Below,
            _ if self.0 < value && value < self.1 => PositionWrtCi::In,
            _ => PositionWrtCi::Above,
        }
    }
}
