//! Statistics related to the Bernoulli distribution. Gated by feature **bernoulli**.

use super::{
    core::{AltHyp, Ci, HypTestResult},
    normal::{z_alpha, z_to_p},
};
use crate::error::Result;
use statrs::distribution::{Beta, Binomial, ContinuousCDF, DiscreteCDF};

/// Estimator of mean of Bernoulli distribution.
///
/// Arguments:
/// - `n`: number of trials.
/// - `n_s`: number of successes (`1`s) observed.
pub fn bernoulli_p_hat(n: u64, n_s: u64) -> f64 {
    n_s as f64 / n as f64
}

/// Normal approximation z-value for the standardized sample mean of a Bernoulli distribution
/// under the hypothesis that the probability of success is `p0`.
///
/// Arguments:
/// - `n`: number of trials.
/// - `n_s`: number of successes (`1`s) observed.
/// - `p0`: probability of success under null hypothesis.
pub fn bernoulli_normal_approx_z(n: u64, n_s: u64, p0: f64) -> f64 {
    let p_hat = bernoulli_p_hat(n, n_s);
    (p_hat - p0) / (p0 * (1. - p0) / n as f64).sqrt()
}

/// Normal approximation p-value for the standardized sample mean of Bernoulli distribution
/// under the hypothesis that the probability of success is `p0`.
///
/// Arguments:
/// - `n`: number of trials.
/// - `n_s`: number of successes (`1`s) observed.
/// - `p0`: probability of success under null hypothesis.
/// - `alt_hyp`: alternative hypothesis.
pub fn bernoulli_normal_approx_p(n: u64, n_s: u64, p0: f64, alt_hyp: AltHyp) -> f64 {
    let z = bernoulli_normal_approx_z(n, n_s, p0);
    z_to_p(z, alt_hyp)
}

/// Binomial proportion confidence interval (Wilson score without continuity correction).
///
/// References:
/// [Wikipedia](https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Wilson_score_interval),
/// [Statistics How To]((https://www.statisticshowto.com/wilson-ci/),
/// [Confidence Intervals for One Proportion](https://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/PASS/Confidence_Intervals_for_One_Proportion.pdf)
///
/// Arguments:
/// - `n`: number of trials.
/// - `n_s`: number of successes (`1`s) observed.
/// - `alt_hyp`: alternative hypothesis.
/// - `alpha`: confidence level = `1 - alpha`.
pub fn binomial_ws_alt_hyp_ci(n: u64, n_s: u64, alt_hyp: AltHyp, alpha: f64) -> Ci {
    let p_hat = bernoulli_p_hat(n, n_s);
    let nr = n as f64;

    let z_alpha = if let AltHyp::Ne = alt_hyp {
        z_alpha(alpha / 2.)
    } else {
        z_alpha(alpha)
    };

    // let base = p_hat + z_alpha.powi(2) / (2. * nr);
    // let delta = z_alpha * (p_hat * (1. - p_hat) / nr + z_alpha.powi(2) / (4. * nr.powi(2))).sqrt();
    // let denom = 1. + z_alpha.powi(2) / nr;

    let base = 2. * nr * p_hat + z_alpha.powi(2);
    let delta = z_alpha * (z_alpha.powi(2) + 4. * nr * p_hat * (1. - p_hat)).sqrt();
    let denom = 2. * (nr + z_alpha.powi(2));

    let (lo, hi) = match alt_hyp {
        AltHyp::Lt => (0., (base + delta) / denom),
        AltHyp::Ne => ((base - delta) / denom, (base + delta) / denom),
        AltHyp::Gt => ((base - delta) / denom, 1.),
    };

    Ci(lo, hi)
}

/// Binomial proportion confidence interval (Wilson score without continuity correction).
/// with the alternative hypothesis of inequality (two-sided).
///
/// References:
/// [Wikipedia](https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Wilson_score_interval),
/// [Statistics How To]((https://www.statisticshowto.com/wilson-ci/),
/// [Confidence Intervals for One Proportion](https://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/PASS/Confidence_Intervals_for_One_Proportion.pdf)
///
/// Arguments:
/// - `n`: number of trials.
/// - `n_s`: number of successes (`1`s) observed.
/// - `alpha`: confidence level = `1 - alpha`.
pub fn binomial_ws_ci(n: u64, n_s: u64, alpha: f64) -> Ci {
    // let p_hat = bernoulli_p_hat(n, n_s);
    // let nr = n as f64;
    // let z_alpha_2 = z_alpha(alpha / 2.);
    // let mid = p_hat + z_alpha_2.powi(2) / (2. * nr);
    // let delta =
    //     z_alpha_2 * (p_hat * (1. - p_hat) / nr + z_alpha_2.powi(2) / (4. * nr.powi(2))).sqrt();
    // let denom = 1. + z_alpha_2.powi(2) / nr;
    // Ci((mid - delta) / denom, (mid + delta) / denom)
    binomial_ws_alt_hyp_ci(n, n_s, AltHyp::Ne, alpha)
}

/// Binomial proportion confidence interval
/// ([Clopper–Pearson](https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Clopper%E2%80%93Pearson_interval)).
pub fn binomial_cp_alt_hyp_ci(n: u64, n_s: u64, alt_hyp: AltHyp, alpha: f64) -> Result<Ci> {
    let lo_beta = Beta::new(n_s as f64, (n - n_s + 1) as f64)?;
    let hi_beta = Beta::new((n_s + 1) as f64, (n - n_s) as f64)?;
    let (lo, hi) = match alt_hyp {
        AltHyp::Lt => {
            let lo = 0.;
            let hi = hi_beta.inverse_cdf(1. - alpha);
            (lo, hi)
        }
        AltHyp::Ne => {
            let lo = lo_beta.inverse_cdf(alpha / 2.);
            let hi = hi_beta.inverse_cdf(1. - alpha / 2.);
            (lo, hi)
        }
        AltHyp::Gt => {
            let lo = lo_beta.inverse_cdf(alpha);
            let hi = 1.;
            (lo, hi)
        }
    };
    Ok(Ci(lo, hi))
}

/// Binomial proportion confidence interval
/// ([Clopper–Pearson](https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Clopper%E2%80%93Pearson_interval)),
/// with the alternative hypothesis of inequality (two-sided).
pub fn binomial_cp_ci(n: u64, n_s: u64, alpha: f64) -> Result<Ci> {
    // let lo_beta = Beta::new(n_s as f64, (n - n_s + 1) as f64)?;
    // let lo = lo_beta.inverse_cdf(alpha / 2.);
    // let hi_beta = Beta::new((n_s + 1) as f64, (n - n_s) as f64)?;
    // let hi = hi_beta.inverse_cdf(1. - alpha / 2.);
    // Ok(Ci(lo, hi))
    binomial_cp_alt_hyp_ci(n, n_s, AltHyp::Ne, alpha)
}

/// p-value for the [one-sample proportion test](exact_binomial_test) (Bernoulli distribution).
///
/// Arguments:
/// - `n`: number of trials.
/// - `n_s`: number of successes (`1`s) observed.
/// - `p0`: probability of success under null hypothesis.
/// - `alt_hyp`: alternative hypothesis.
/// - `alpha`: confidence level = `1 - alpha`.
///
/// # Errors
/// Returns an error under any of these conditions:
/// - `n == 0`.
/// - `p0` not strictly between `0` and `1`.
/// - `alpha < 0` or `alpha > 1`.
pub fn exact_binomial_p(n: u64, n_s: u64, p0: f64, alt_hyp: AltHyp) -> Result<f64> {
    let binomial = Binomial::new(p0, n)?;
    let prob_le = binomial.cdf(n_s);
    let _prob_lt = binomial.cdf(n_s - 1);
    let prob_ge = binomial.cdf(n) - _prob_lt;

    let target_s_lo = (n as f64 * p0).floor() as u64;
    let target_s_hi = (n as f64 * p0).ceil() as u64;
    let abs = n_s.abs_diff(target_s_lo).min(n_s.abs_diff(target_s_hi));

    let extreme_lo = (target_s_lo - abs).max(0);
    let extreme_hi = (target_s_hi + abs).min(n);
    let prob_extreme = (binomial.cdf(extreme_lo) + (1. - binomial.cdf(extreme_hi - 1))).min(1.);

    let p_value = match alt_hyp {
        AltHyp::Lt => prob_le,
        AltHyp::Gt => prob_ge,
        AltHyp::Ne => prob_extreme,
    };

    Ok(p_value)
}

/// One-sample proportion test (Bernoulli distribution).
///
/// Arguments:
/// - `n`: number of trials.
/// - `n_s`: number of successes (`1`s) observed.
/// - `p0`: probability of success under null hypothesis.
/// - `alt_hyp`: alternative hypothesis.
/// - `alpha`: confidence level = `1 - alpha`.
///
/// # Errors
/// Returns an error under any of these conditions:
/// - `n == 0`.
/// - `p0` not strictly between `0` and `1`.
/// - `alpha < 0` or `alpha > 1`.
pub fn exact_binomial_test(
    n: u64,
    n_s: u64,
    p0: f64,
    alt_hyp: AltHyp,
    alpha: f64,
) -> Result<HypTestResult> {
    let p_value = exact_binomial_p(n, n_s, p0, alt_hyp)?;
    let test_res = HypTestResult::new(p_value, alpha, alt_hyp);
    Ok(test_res)
}

#[cfg(test)]
mod test {
    use crate::{core::Hyp, dev_utils::ApproxEq};

    use super::*;

    const ALPHA: f64 = 0.05;
    const EPSILON: f64 = 0.0005;

    fn n_n_s() -> (u64, u64) {
        // Data generated with p0 = 0.4.
        let dat = [
            0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
            0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
            0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1,
            1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0,
        ];

        let n = dat.len() as u64;
        let n_s = dat.into_iter().sum::<u64>();
        (n, n_s)
    }
    fn check_bernoulli(
        n: u64,
        n_s: u64,
        p0: f64,
        alt_hyp: AltHyp,
        exp_p: f64,
        exp_cp_ci: Ci,
        exp_ws_ci: Ci,
        exp_accept_hyp: Hyp,
    ) {
        let cp_ci = binomial_cp_alt_hyp_ci(n, n_s, alt_hyp, ALPHA).unwrap();
        let ws_ci = binomial_ws_alt_hyp_ci(n, n_s, alt_hyp, ALPHA);
        let res = exact_binomial_test(n, n_s, p0, alt_hyp, ALPHA).unwrap();
        let p = res.p();

        assert!(
            exp_p.approx_eq(p, EPSILON),
            "alt_hyp={alt_hyp:?} -- exp_p={exp_p}, p={p}"
        );

        assert!(
            exp_cp_ci.0.approx_eq(cp_ci.0, EPSILON)
                || exp_cp_ci.0.is_infinite() && cp_ci.0.is_infinite(),
            "alt_hyp={alt_hyp:?} -- exp_ci.0={}, ci.0={}",
            exp_cp_ci.0,
            cp_ci.0
        );
        assert!(
            exp_cp_ci.1.approx_eq(cp_ci.1, EPSILON)
                || exp_cp_ci.1.is_infinite() && cp_ci.1.is_infinite(),
            "alt_hyp={alt_hyp:?} -- exp_ci.1={}, ci.1={}",
            exp_cp_ci.1,
            cp_ci.1
        );

        assert!(
            exp_ws_ci.0.approx_eq(ws_ci.0, EPSILON)
                || exp_ws_ci.0.is_infinite() && ws_ci.0.is_infinite(),
            "alt_hyp={alt_hyp:?} -- exp_ci.0={}, ci.0={}",
            exp_ws_ci.0,
            ws_ci.0
        );
        assert!(
            exp_ws_ci.1.approx_eq(ws_ci.1, EPSILON)
                || exp_ws_ci.1.is_infinite() && ws_ci.1.is_infinite(),
            "alt_hyp={alt_hyp:?} -- exp_ci.1={}, ci.1={}",
            exp_ws_ci.1,
            ws_ci.1
        );

        assert_eq!(p, res.p(), "alt_hyp={alt_hyp:?} -- res.p");
        assert_eq!(ALPHA, res.alpha(), "alt_hyp={alt_hyp:?} -- res.alpha");
        assert_eq!(alt_hyp, res.alt_hyp(), "alt_hyp={alt_hyp:?} -- res.alt_hyp");
        assert_eq!(
            exp_accept_hyp,
            res.accepted(),
            "alt_hyp={alt_hyp:?} -- res.accepted"
        );
    }

    #[test]
    fn test_bern_lt() {
        let (n, n_s) = n_n_s();
        let p0 = 0.5;
        let alt_hyp = AltHyp::Lt;
        let exp_p = 0.02844;
        let exp_cp_ci = Ci(0.0000000, 0.4870242);
        let exp_ws_ci = Ci(0.0000000, 0.4821905);
        let exp_accept_hyp = Hyp::Alt(AltHyp::Lt);

        check_bernoulli(
            n,
            n_s,
            p0,
            alt_hyp,
            exp_p,
            exp_cp_ci,
            exp_ws_ci,
            exp_accept_hyp,
        );
    }

    #[test]
    fn test_bern_eq() {
        let (n, n_s) = n_n_s();
        let p0 = 0.5;
        let alt_hyp = AltHyp::Ne;
        let exp_p = 0.05689;
        let exp_cp_ci = Ci(0.3032948, 0.5027908);
        let exp_ws_ci = Ci(0.3094013, 0.4979974);
        let exp_accept_hyp = Hyp::Null;

        check_bernoulli(
            n,
            n_s,
            p0,
            alt_hyp,
            exp_p,
            exp_cp_ci,
            exp_ws_ci,
            exp_accept_hyp,
        );
    }

    #[test]
    fn test_bern_gt() {
        let (n, n_s) = n_n_s();
        let p0 = 0.5;
        let alt_hyp = AltHyp::Gt;
        let exp_p = 0.9824;
        let exp_cp_ci = Ci(0.317526, 1.000000);
        let exp_ws_ci = Ci(0.3230781, 1.0000000);
        let exp_accept_hyp = Hyp::Null;

        check_bernoulli(
            n,
            n_s,
            p0,
            alt_hyp,
            exp_p,
            exp_cp_ci,
            exp_ws_ci,
            exp_accept_hyp,
        );
    }

    #[test]
    fn test_bern_eq_40() {
        let (n, n_s) = n_n_s();
        let p0 = 0.4;
        let alt_hyp = AltHyp::Ne;
        let exp_p = 1.;
        let exp_cp_ci = Ci(0.3032948, 0.5027908);
        let exp_ws_ci = Ci(0.3094013, 0.4979974);
        let exp_accept_hyp = Hyp::Null;

        check_bernoulli(
            n,
            n_s,
            p0,
            alt_hyp,
            exp_p,
            exp_cp_ci,
            exp_ws_ci,
            exp_accept_hyp,
        );
    }
}
