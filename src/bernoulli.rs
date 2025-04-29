use super::{
    core::{AltHyp, Ci, HypTestResult},
    normal::{z_alpha, z_to_p},
};

#[allow(unused)]
/// Estimator of mean of Bernoulli distribution.
///
/// Arguments:
/// - `n`: number of trials.
/// - `successes`: number of successes (`1`s) observed.
pub fn bernoulli_p_hat(n: u64, successes: f64) -> f64 {
    successes / n as f64
}

/// Confidence interval for the probability of success of a Bernoulli distribution computed from a
/// sample of trials (Wilson score interval).
///
/// Arguments:
/// - `n`: number of trials.
/// - `p_hat`: estimate of mean of the Bernoulli distribution. See [`bernoulli_p_hat`].
/// - `alpha`: confidence level.
pub fn bernoulli_psucc_ci(n: u64, p_hat: f64, alpha: f64) -> Ci {
    let nf = n as f64;
    let p = p_hat;
    let z_alpha_2 = z_alpha(alpha / 2.);
    let mid = p + z_alpha_2.powi(2) / (2. * nf);
    let delta = z_alpha_2 * (p * (1. - p) / nf + z_alpha_2.powi(2) / (4. * nf.powi(2))).sqrt();
    let denom = 1. + z_alpha_2.powi(2) / nf;
    Ci((mid - delta) / denom, (mid + delta) / denom)
}

/// Normal approximation z-value for standardized sample mean of Bernoulli distribution under the hypothesis that
/// the probability of success is `p0`.
///
/// Arguments:
/// - `p_hat`: estimate of mean of the Bernoulli distribution. See [`bernoulli_p_hat`].
/// - `p0`: probability of success under null hypothesis.
pub fn bernoulli_normal_approx_z(n: u64, p_hat: f64, p0: f64) -> f64 {
    (p_hat - p0) / (p0 * (1. - p0) / n as f64).sqrt()
}

/// Normal approximation p-value for standardized sample mean of Bernoulli distribution under the hypothesis that
/// the probability of success is `p0`.
///
/// Arguments:
/// - `n`: number of trials.
/// - `p_hat`: estimate of mean of the Bernoulli distribution. See [`bernoulli_p_hat`].
/// - `p0`: probability of success under null hypothesis.
/// - `alt_hyp`: alternative hypothesis.
pub fn bernoulli_normal_approx_p(n: u64, p_hat: f64, p0: f64, alt_hyp: AltHyp) -> f64 {
    let z = bernoulli_normal_approx_z(n, p_hat, p0);
    z_to_p(z, alt_hyp)
}

pub fn bernoulli_test(n: u64, p_hat: f64, p0: f64, alt_hyp: AltHyp, alpha: f64) -> HypTestResult {
    let p = bernoulli_normal_approx_p(n, p_hat, p0, alt_hyp);
    HypTestResult::new(p, alpha, alt_hyp)
}
