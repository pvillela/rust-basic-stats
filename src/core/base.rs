//! Basic sample statistics and common types supporting inferential statistics.

use super::{StatsError, StatsResult};

/// Sample mean.
///
/// Arguments:
/// - `n` - sample size.
/// - `sum` - sample sum.
///
/// # Errors
///
/// Returns an error if `n == 0`.
#[inline(always)]
pub fn sample_mean(n: u64, sum: f64) -> StatsResult<f64> {
    if n == 0 {
        return Err(StatsError("arg `n` must be positive"));
    }
    Ok(sum / n as f64)
}

/// Sample's sum of squares of deviations from sample mean.
///
/// Arguments:
/// - `n` - sample size.
/// - `sum` - sample sum.
/// - `sum2` - sample sum of squares.
///
/// # Errors
///
/// Returns an error if `n == 0`.
#[inline(always)]
pub fn sample_sum2_deviations(n: u64, sum: f64, sum2: f64) -> StatsResult<f64> {
    if n == 0 {
        return Err(StatsError("arg `n` must be positive"));
    }
    Ok(sum2 - sum.powi(2) / n as f64)
}

/// Sample variance.
///
/// Arguments:
/// - `n` - sample size.
/// - `sum` - sample sum.
/// - `sum2` - sample sum of squares.
///
/// # Errors
///
/// Returns an error if `n <= 1`.
#[inline(always)]
pub fn sample_var(n: u64, sum: f64, sum2: f64) -> StatsResult<f64> {
    if n <= 1 {
        return Err(StatsError("arg `n` must be greater than `1`"));
    }
    Ok(sample_sum2_deviations(n, sum, sum2)? / (n as f64 - 1.))
}

/// Sample standard deviation.
///
/// Arguments:
/// - `n` - sample size.
/// - `sum` - sample sum.
/// - `sum2` - sample sum of squares.
///
/// # Errors
///
/// Returns an error if `n <= 1`.
#[inline(always)]
pub fn sample_stdev(n: u64, sum: f64, sum2: f64) -> StatsResult<f64> {
    Ok(sample_var(n, sum, sum2)?.sqrt())
}

/// Holds summary information for a sample, enabling computation of a sample's mean and variance.
///
/// Includes sample size, sum, sum of squares, and, optionally, sample min and max values.
///
/// Used by several by functions in this library.
#[derive(Debug, PartialEq)]
pub struct SampleMoments {
    n: u64,
    sum: f64,
    sum2: f64,
    min: f64,
    max: f64,
}

impl SampleMoments {
    /// Instantiates `Self` with an empty sample.
    ///
    /// Minimum and maximum sample values are defaulted to `NaN`.
    pub fn new_empty() -> Self {
        Self::new(0, 0., 0.)
    }

    /// Instantiates `Self` with given values for sample size, sample sum, and sample sum of squares.
    ///
    /// Minimum and maximum sample values are defaulted to `NaN`.
    pub fn new(n: u64, sum: f64, sum2: f64) -> Self {
        Self::new_with_min_max(n, sum, sum2, f64::NAN, f64::NAN)
    }

    /// Instantiates `Self`, with given values for sample size, sample sum, sample sum of squares, min, and max.
    pub fn new_with_min_max(n: u64, sum: f64, sum2: f64, min: f64, max: f64) -> Self {
        Self {
            n,
            sum,
            sum2,
            min,
            max,
        }
    }

    /// Updates `self` by accumulating an additional value.
    pub fn collect_value(&mut self, value: f64) {
        self.n += 1;
        self.sum += value;
        self.sum2 += value * value;
        self.min = value.min(self.min);
        self.max = value.max(self.max);
    }

    /// Instatiates `Self` from an iterator.
    pub fn from_iterator(dataset: impl Iterator<Item = f64>) -> Self {
        let mut moments = SampleMoments::new_empty();
        for v in dataset {
            moments.collect_value(v);
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
        self.n
    }

    /// Sample size as floating point.
    pub fn nf(&self) -> f64 {
        self.n as f64
    }

    /// Sample sum.
    pub fn sum(&self) -> f64 {
        self.sum
    }

    /// Sample mean.
    ///
    /// # Errors
    ///
    /// Returns an error if `n == 0`.
    pub fn mean(&self) -> StatsResult<f64> {
        sample_mean(self.n, self.sum)
    }

    /// Sample sum of squares.
    pub fn sum2(&self) -> f64 {
        self.sum2
    }

    /// Sample's sum of squares of deviations from sample mean.
    ///
    /// # Errors
    ///
    /// Returns an error if `n == 0`.
    pub fn sum2_deviations(&self) -> StatsResult<f64> {
        sample_sum2_deviations(self.n, self.sum, self.sum2)
    }

    /// Sample variance.
    ///
    /// # Errors
    ///
    /// Returns an error if `n <= 1`.
    pub fn var(&self) -> StatsResult<f64> {
        sample_var(self.n, self.sum, self.sum2)
    }

    /// Sample standard deviation.
    ///
    /// # Errors
    ///
    /// Returns an error if `n <= 1`.
    pub fn stdev(&self) -> StatsResult<f64> {
        sample_stdev(self.n, self.sum, self.sum2)
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
            accepted: if p < alpha {
                Hyp::Alt(alt_hyp)
            } else {
                Hyp::Null
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

#[cfg(test)]
mod test {
    use super::SampleMoments;
    use crate::dev_utils::ApproxEq;

    const EPSILON: f64 = 0.000005;

    #[test]
    fn test_moments() {
        let x = [14., 15., 15., 15., 16., 18., 22., 23., 24., 25., 25.];
        let n = 11;
        let nf = 11.;
        let sum = 212.;
        let sum2 = 4290.;
        let var = 20.41818;
        let sum2_dev = var * ((n - 1) as f64);
        let stdev = 4.518648;
        let min = 14.;
        let max = 25.;

        let moments = SampleMoments::new_with_min_max(n, sum, sum2, min, max);
        let moments_s = SampleMoments::from_slice(&x);

        assert_eq!(moments, moments_s);
        assert_eq!(n, moments.n());
        assert_eq!(nf, moments.nf());
        assert_eq!(sum, moments.sum());
        assert_eq!(sum2, moments.sum2());
        assert!(var.approx_eq(moments.var().unwrap(), EPSILON));
        assert!(sum2_dev.approx_eq(moments.sum2_deviations().unwrap(), EPSILON * 10.));
        assert!(stdev.approx_eq(moments.stdev().unwrap(), EPSILON));
        assert_eq!(min, moments.min());
        assert_eq!(max, moments.max());
    }
}
