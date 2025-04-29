#[inline(always)]
pub fn sample_mean(n: u64, sum: f64) -> f64 {
    sum / n as f64
}

#[inline(always)]
pub fn sample_sum2_deviations(n: u64, sum: f64, sum2: f64) -> f64 {
    sum2 - sum.powi(2) / n as f64
}

#[inline(always)]
pub fn sample_var(n: u64, sum: f64, sum2: f64) -> f64 {
    sample_sum2_deviations(n, sum, sum2) / (n as f64 - 1.)
}

#[inline(always)]
pub fn sample_stdev(n: u64, sum: f64, sum2: f64) -> f64 {
    sample_var(n, sum, sum2).sqrt()
}

pub struct SampleMoments {
    count: u64,
    sum: f64,
    sum2: f64,
    min: f64,
    max: f64,
}

impl SampleMoments {
    pub fn new(count: u64, sum: f64, sum2: f64) -> Self {
        Self {
            count,
            sum,
            sum2,
            min: f64::NAN,
            max: f64::NAN,
        }
    }

    pub fn new_with_min_max(count: u64, sum: f64, sum2: f64, min: f64, max: f64) -> Self {
        Self {
            count,
            sum,
            sum2,
            min,
            max,
        }
    }

    pub fn empty() -> Self {
        Self::new(0, 0., 0.)
    }

    pub fn collect_value(&mut self, value: f64) {
        self.count += 1;
        self.sum += value;
        self.sum2 += value * value;
        self.min = value.min(self.min);
        self.max = value.max(self.max);
    }

    pub fn from_slice(dataset: &[f64]) -> Self {
        let mut moments = SampleMoments::empty();
        for value in dataset {
            moments.collect_value(*value);
        }
        moments
    }

    pub fn n(&self) -> u64 {
        self.count
    }

    pub fn nf(&self) -> f64 {
        self.count as f64
    }

    pub fn sum(&self) -> f64 {
        self.sum
    }

    pub fn mean(&self) -> f64 {
        self.sum / self.nf()
    }

    pub fn sum2(&self) -> f64 {
        self.sum2
    }

    pub fn sum2_deviations(&self) -> f64 {
        sample_sum2_deviations(self.n(), self.sum, self.sum2)
    }

    /// Returns `Nan` if  `self.n()` == 1
    pub fn var(&self) -> f64 {
        sample_var(self.n(), self.sum, self.sum2)
    }

    /// Returns `Nan` if  `self.n()` == 1
    pub fn stdev(&self) -> f64 {
        sample_stdev(self.n(), self.sum, self.sum2)
    }

    pub fn min(&self) -> f64 {
        self.min
    }

    pub fn max(&self) -> f64 {
        self.max
    }
}

impl Default for SampleMoments {
    fn default() -> Self {
        Self {
            count: 0,
            sum: 0.,
            sum2: 0.,
            min: f64::NAN,
            max: f64::NAN,
        }
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
