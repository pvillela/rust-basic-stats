//! Statistics related to samples of the Bernoulli distribution. The Binomial distribution with parameters `n` and `p`
//! is the distribution of the sum of `n` independent Bernoulli random variables with probability of success `p`.
//!
//! This module is included by default. However, if `default-features = false` is specified in the dependency
//! declaration for this library, then inclusion of this module is gated by feature "**binomial**".

use std::cmp::Ordering;

use super::{
    core::{AltHyp, Ci, HypTestResult},
    normal::{z_alpha, z_to_p},
};
use crate::core::{
    AsStatsResult, StatsError, StatsResult, check_alpha_in_open_0_1, check_p0_in_open_0_1,
};
use statrs::distribution::{Beta, Binomial, ContinuousCDF, Discrete, DiscreteCDF};

/// Estimator of success probability of Bernoulli distribution.
///
/// Arguments:
/// - `n`: number of trials.
/// - `n_s`: number of successes (`1`s) observed.
///
/// # Errors
///
/// Returns an error if `n == 0` or `n < n_s`.
pub fn bernoulli_p_hat(n: u64, n_s: u64) -> StatsResult<f64> {
    if n == 0 {
        return Err(StatsError("arg `n` must be positive"));
    }
    if n < n_s {
        return Err(StatsError("arg `n` must be greater than or equal to n_s"));
    }
    Ok(n_s as f64 / n as f64)
}

/// Normal approximation z-value for the standardized sample mean of a Bernoulli distribution
/// under the hypothesis that the probability of success is `p0`.
/// Without continuity correction.
///
/// Arguments:
/// - `n`: number of trials.
/// - `n_s`: number of successes (`1`s) observed.
/// - `p0`: probability of success under null hypothesis.
///
/// # Errors
///
/// Returns an error in any of these circumstances:
/// Returns an error if `n == 0` or `n < n_s`.
/// - `p0` not in `(0, 1)`.
pub fn binomial_z(n: u64, n_s: u64, p0: f64) -> StatsResult<f64> {
    check_p0_in_open_0_1(p0)?;
    let p_hat = bernoulli_p_hat(n, n_s)?;
    let ret = (p_hat - p0) / (p0 * (1. - p0) / n as f64).sqrt();
    Ok(ret)
}

/// Normal approximation p-value for the standardized sample mean of Bernoulli distribution
/// under the hypothesis that the probability of success is `p0`.
/// Without continuity correction.
///
/// Arguments:
/// - `n`: number of trials.
/// - `n_s`: number of successes (`1`s) observed.
/// - `p0`: probability of success under null hypothesis.
/// - `alt_hyp`: alternative hypothesis.
///
/// # Errors
///
/// Returns an error in any of these circumstances:
/// - `n == 0` or `n < n_s`.
/// - `p0` not in `(0, 1)`.
pub fn binomial_z_p(n: u64, n_s: u64, p0: f64, alt_hyp: AltHyp) -> StatsResult<f64> {
    let z = binomial_z(n, n_s, p0)?;
    Ok(z_to_p(z, alt_hyp))
}

/// One-sample proportion test (Bernoulli distribution) using the Binomial Normal approximation.
/// Without continuity correction.
///
/// Arguments:
/// - `n`: number of trials.
/// - `n_s`: number of successes (`1`s) observed.
/// - `p0`: probability of success under null hypothesis.
/// - `alt_hyp`: alternative hypothesis.
/// - `alpha`: confidence level = `1 - alpha`.
///
/// # Errors
///
/// Returns an error in any of these circumstances:
/// - `n == 0` or `n < n_s`.
/// - `p0` is not in `(0, 1)`.
/// - `alpha` is not in `(0, 1)`.
pub fn one_proportion_z_test(
    n: u64,
    n_s: u64,
    p0: f64,
    alt_hyp: AltHyp,
    alpha: f64,
) -> StatsResult<HypTestResult> {
    check_alpha_in_open_0_1(alpha)?;
    let p_value = binomial_z_p(n, n_s, p0, alt_hyp)?;
    let test_res = HypTestResult::new(p_value, alpha, alt_hyp);
    Ok(test_res)
}

/// Binomial proportion confidence interval (Wilson score without continuity correction).
///
/// References:
/// [Wikipedia](https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Wilson_score_interval),
/// [Statistics How To](https://www.statisticshowto.com/wilson-ci/),
/// [Confidence Intervals for One Proportion](https://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/PASS/Confidence_Intervals_for_One_Proportion.pdf)
///
/// Arguments:
/// - `n`: number of trials.
/// - `n_s`: number of successes (`1`s) observed.
/// - `alt_hyp`: alternative hypothesis.
/// - `alpha`: confidence level = `1 - alpha`.
///
/// # Errors
///
/// Returns an error in any of these circumstances:
/// - `n == 0` or `n < n_s`.
/// - `alpha` not in `(0, 1)`.
pub fn binomial_ws_alt_hyp_ci(n: u64, n_s: u64, alt_hyp: AltHyp, alpha: f64) -> StatsResult<Ci> {
    let p_hat = bernoulli_p_hat(n, n_s)?;

    let nr = n as f64;

    check_alpha_in_open_0_1(alpha)?; // need this guard because `alpha / 2.` below masks errors
    let z_alpha = if let AltHyp::Ne = alt_hyp {
        z_alpha(alpha / 2.)?
    } else {
        z_alpha(alpha)?
    };

    let base = 2. * nr * p_hat + z_alpha.powi(2);
    let delta = z_alpha * (z_alpha.powi(2) + 4. * nr * p_hat * (1. - p_hat)).sqrt();
    let denom = 2. * (nr + z_alpha.powi(2));

    let (lo, hi) = match alt_hyp {
        AltHyp::Lt => (0., (base + delta) / denom),
        AltHyp::Ne => ((base - delta) / denom, (base + delta) / denom),
        AltHyp::Gt => ((base - delta) / denom, 1.),
    };

    Ok(Ci(lo, hi))
}

/// Binomial proportion confidence interval (Wilson score without continuity correction).
/// with the alternative hypothesis of inequality (two-sided).
///
/// References:
/// [Wikipedia](https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Wilson_score_interval),
/// [Statistics How To](https://www.statisticshowto.com/wilson-ci/),
/// [Confidence Intervals for One Proportion](https://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/PASS/Confidence_Intervals_for_One_Proportion.pdf)
///
/// Arguments:
/// - `n`: number of trials.
/// - `n_s`: number of successes (`1`s) observed.
/// - `alpha`: confidence level = `1 - alpha`.
///
/// # Errors
///
/// Returns an error in any of these circumstances:
/// - `n == 0` or `n < n_s`.
/// - `alpha` not in `(0, 1)`.
pub fn binomial_ws_ci(n: u64, n_s: u64, alpha: f64) -> StatsResult<Ci> {
    binomial_ws_alt_hyp_ci(n, n_s, AltHyp::Ne, alpha)
}

/// Binomial proportion confidence interval
/// ([Clopper–Pearson](https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Clopper%E2%80%93Pearson_interval)).
///
/// See also [Confidence Intervals for One Proportion](https://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/PASS/Confidence_Intervals_for_One_Proportion.pdf).
///
/// # Errors
///
/// Returns an error in any of these circumstances:
/// - `n == 0` or `n < n_s`.
/// - `alpha` is not in `(0, 1)`.
pub fn binomial_cp_alt_hyp_ci(n: u64, n_s: u64, alt_hyp: AltHyp, alpha: f64) -> StatsResult<Ci> {
    if n == 0 {
        return Err(StatsError("arg `n` must be positive"));
    }
    if n < n_s {
        return Err(StatsError(
            "arg `n` must be greater than or equal to arg `n_s`",
        ));
    }
    check_alpha_in_open_0_1(alpha)?;

    // Include special cases not handled by Beta function.
    // Below closures based on https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/binom.test.R
    // code for lambdas p.L and p.U.

    let lo_qbeta = |alpha| -> StatsResult<f64> {
        if n_s == 0 {
            Ok(0.)
        } else {
            let lo_beta = Beta::new(n_s as f64, (n - n_s + 1) as f64)
                .stats_result("invalid arg `n` or `n_s`")?;
            Ok(lo_beta.inverse_cdf(alpha))
        }
    };

    let hi_qbeta = |alpha| -> StatsResult<f64> {
        if n_s == n {
            Ok(1.)
        } else {
            let hi_beta = Beta::new((n_s + 1) as f64, (n - n_s) as f64)
                .stats_result("invalid arg `n` or `n_s`")?;
            Ok(hi_beta.inverse_cdf(1. - alpha))
        }
    };

    let (lo, hi) = match alt_hyp {
        AltHyp::Lt => {
            let lo = 0.;
            let hi = hi_qbeta(alpha)?;
            (lo, hi)
        }
        AltHyp::Ne => {
            let lo = lo_qbeta(alpha / 2.)?;
            let hi = hi_qbeta(alpha / 2.)?;
            (lo, hi)
        }
        AltHyp::Gt => {
            let lo = lo_qbeta(alpha)?;
            let hi = 1.;
            (lo, hi)
        }
    };

    Ok(Ci(lo, hi))
}

/// Binomial proportion confidence interval
/// ([Clopper–Pearson](https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Clopper%E2%80%93Pearson_interval)),
/// with the alternative hypothesis of inequality (two-sided).
///
/// See also [Confidence Intervals for One Proportion](https://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/PASS/Confidence_Intervals_for_One_Proportion.pdf).
///
/// # Errors
///
/// Returns an error in any of these circumstances:
/// - `n == 0` or `n < n_s`.
/// - `alpha` is not in `[0, 1]`.
pub fn binomial_cp_ci(n: u64, n_s: u64, alpha: f64) -> StatsResult<Ci> {
    binomial_cp_alt_hyp_ci(n, n_s, AltHyp::Ne, alpha)
}

/// p-value for the [one-sample proportion test](exact_binomial_test) (Bernoulli distribution).
///
/// Arguments:
/// - `n`: number of trials.
/// - `n_s`: number of successes (`1`s) observed.
/// - `p0`: probability of success under null hypothesis.
/// - `alt_hyp`: alternative hypothesis.
///
/// # Errors
///
/// Returns an error in any of these circumstances:
/// - `n == 0` or `n < n_s`.
/// - `p0` is not in `[0, 1]`.
pub fn exact_binomial_p(n: u64, n_s: u64, p0: f64, alt_hyp: AltHyp) -> StatsResult<f64> {
    if n == 0 {
        return Err(StatsError("arg `n` must be positive."));
    }
    if n < n_s {
        return Err(StatsError("arg `n` must not be less than arg `n_s`."));
    }

    let binomial = Binomial::new(p0, n).stats_result("arg `p0` must be in interval [0, 1]")?;

    let prob_le = binomial.cdf(n_s);

    let prob_ge = {
        let prob_lt = if n_s == 0 { 0. } else { binomial.cdf(n_s - 1) };
        1. - prob_lt
    };

    // Sum the probabilities of all values with probability lower or equal to `n_s`'s.
    // Based on https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/binom.test.R
    // code for PVAL.
    let prob_ne = {
        if p0 == 0. {
            (n_s == 0) as u64 as f64
        } else if p0 == 1. {
            (n_s == n) as u64 as f64
        } else {
            let rel_err = 1. + 1e-7;
            let prob_n_s = binomial.pmf(n_s);
            let mode = n as f64 * p0;

            match (n_s as f64).total_cmp(&mode) {
                Ordering::Equal => {
                    // By definition of mode, all other values have prob less than `n_s`'s.
                    1.
                }

                Ordering::Less => {
                    let mut sum_prob = 0.;
                    let imode = mode.ceil() as u64;
                    for i in (imode..=n).rev() {
                        let prob_i = binomial.pmf(i);
                        if prob_i <= prob_n_s * rel_err {
                            sum_prob += prob_i;
                        } else {
                            break;
                        }
                    }
                    prob_le + sum_prob
                }

                Ordering::Greater => {
                    let mut sum_prob = 0.;
                    let imode = mode.floor() as u64;
                    for i in 0..=imode {
                        let prob_i = binomial.pmf(i);
                        if prob_i <= prob_n_s * rel_err {
                            sum_prob += prob_i;
                        } else {
                            break;
                        }
                    }
                    prob_ge + sum_prob
                }
            }
        }
    };

    let p_value = match alt_hyp {
        AltHyp::Lt => prob_le,
        AltHyp::Gt => prob_ge,
        AltHyp::Ne => prob_ne,
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
///
/// Returns an error in any of these circumstances:
/// - `n == 0` or `n < n_s`.
/// - `p0` is not in `[0, 1]`.
/// - `alpha` is not in `(0, 1)`.
pub fn exact_binomial_test(
    n: u64,
    n_s: u64,
    p0: f64,
    alt_hyp: AltHyp,
    alpha: f64,
) -> StatsResult<HypTestResult> {
    check_alpha_in_open_0_1(alpha)?;

    let p_value = exact_binomial_p(n, n_s, p0, alt_hyp)?;
    let test_res = HypTestResult::new(p_value, alpha, alt_hyp);
    Ok(test_res)
}

#[cfg(test)]
mod test {
    use crate::{core::Hyp, dev_utils::ApproxEq};

    use super::*;

    const ALPHA: f64 = 0.05;
    const EPSILON: f64 = 0.00005;

    fn check_binomial_no_z(
        n: u64,
        n_s: u64,
        p0: f64,
        alt_hyp: AltHyp,
        exp_p: f64,
        exp_cp_ci: Ci,
        exp_accept_hyp: Hyp,
    ) {
        let cp_ci = binomial_cp_alt_hyp_ci(n, n_s, alt_hyp, ALPHA).unwrap();
        let res = exact_binomial_test(n, n_s, p0, alt_hyp, ALPHA).unwrap();
        let p = res.p();

        if alt_hyp == AltHyp::Ne {
            assert_eq!(cp_ci, binomial_cp_ci(n, n_s, ALPHA).unwrap());
        }

        assert!(
            exp_p.approx_eq(p, EPSILON),
            "n={n}, n_s={n_s}, alt_hyp={alt_hyp:?} -> exp_p={exp_p}, p={p}"
        );

        assert!(
            exp_cp_ci.0.approx_eq(cp_ci.0, EPSILON)
                || exp_cp_ci.0.is_infinite() && cp_ci.0.is_infinite(),
            "n={n}, n_s={n_s}, alt_hyp={alt_hyp:?} -> exp_cp_ci.0={}, cp_ci.0={}",
            exp_cp_ci.0,
            cp_ci.0
        );
        assert!(
            exp_cp_ci.1.approx_eq(cp_ci.1, EPSILON)
                || exp_cp_ci.1.is_infinite() && cp_ci.1.is_infinite(),
            "n={n}, n_s={n_s}, alt_hyp={alt_hyp:?} -> exp_cp_ci.1={}, cp_ci.1={}",
            exp_cp_ci.1,
            cp_ci.1
        );

        assert_eq!(
            ALPHA,
            res.alpha(),
            "n={n}, n_s={n_s}, alt_hyp={alt_hyp:?} -> res.alpha"
        );
        assert_eq!(
            alt_hyp,
            res.alt_hyp(),
            "n={n}, n_s={n_s}, alt_hyp={alt_hyp:?} -> res.alt_hyp"
        );
        assert_eq!(
            exp_accept_hyp,
            res.accepted(),
            "n={n}, n_s={n_s}, alt_hyp={alt_hyp:?} -> res.accepted"
        );
    }

    fn check_binomial(
        n: u64,
        n_s: u64,
        p0: f64,
        alt_hyp: AltHyp,
        exp_p: f64,
        exp_z_p: f64,
        exp_cp_ci: Ci,
        exp_ws_ci: Ci,
        exp_accept_hyp: Hyp,
        exp_z_accept_hyp: Hyp,
    ) {
        check_binomial_no_z(n, n_s, p0, alt_hyp, exp_p, exp_cp_ci, exp_accept_hyp);

        // check z

        let ws_ci = binomial_ws_alt_hyp_ci(n, n_s, alt_hyp, ALPHA).unwrap();
        let z_res = one_proportion_z_test(n, n_s, p0, alt_hyp, ALPHA).unwrap();
        let z_p = z_res.p();

        if alt_hyp == AltHyp::Ne {
            assert_eq!(ws_ci, binomial_ws_ci(n, n_s, ALPHA).unwrap());
        }

        assert!(
            exp_z_p.approx_eq(z_p, EPSILON),
            "n={n}, n_s={n_s}, alt_hyp={alt_hyp:?} -> exp_z_p={exp_z_p}, z_p={z_p}"
        );

        assert!(
            exp_ws_ci.0.approx_eq(ws_ci.0, EPSILON)
                || exp_ws_ci.0.is_infinite() && ws_ci.0.is_infinite(),
            "n={n}, n_s={n_s}, alt_hyp={alt_hyp:?} -> exp_ws_ci.0={}, ws_ci.0={}",
            exp_ws_ci.0,
            ws_ci.0
        );
        assert!(
            exp_ws_ci.1.approx_eq(ws_ci.1, EPSILON)
                || exp_ws_ci.1.is_infinite() && ws_ci.1.is_infinite(),
            "n={n}, n_s={n_s}, alt_hyp={alt_hyp:?} -> exp_ws_ci.1={}, ws_ci.1={}",
            exp_ws_ci.1,
            ws_ci.1
        );

        assert_eq!(
            ALPHA,
            z_res.alpha(),
            "n={n}, n_s={n_s}, alt_hyp={alt_hyp:?} -> z_res.alpha"
        );
        assert_eq!(
            alt_hyp,
            z_res.alt_hyp(),
            "n={n}, n_s={n_s}, alt_hyp={alt_hyp:?} -> z_res.alt_hyp"
        );
        assert_eq!(
            exp_z_accept_hyp,
            z_res.accepted(),
            "n={n}, n_s={n_s}, alt_hyp={alt_hyp:?} -> z_res.accepted"
        );
    }

    #[test]
    fn test_binom_lt_100_40_05() {
        let (n, n_s) = (100, 40);
        let p0 = 0.5;
        let alt_hyp = AltHyp::Lt;
        let exp_p = 0.02844;
        let exp_z_p = 0.02275;
        let exp_cp_ci = Ci(0.0000000, 0.4870242);
        let exp_ws_ci = Ci(0.0000000, 0.4821905);
        let exp_accept_hyp = Hyp::Alt(AltHyp::Lt);
        let exp_z_accept_hyp = exp_accept_hyp;

        check_binomial(
            n,
            n_s,
            p0,
            alt_hyp,
            exp_p,
            exp_z_p,
            exp_cp_ci,
            exp_ws_ci,
            exp_accept_hyp,
            exp_z_accept_hyp,
        );
    }

    #[test]
    fn test_binom_eq_100_40_05() {
        let (n, n_s) = (100, 40);
        let p0 = 0.5;
        let alt_hyp = AltHyp::Ne;
        let exp_p = 0.05689;
        let exp_z_p = 0.0455;
        let exp_cp_ci = Ci(0.3032948, 0.5027908);
        let exp_ws_ci = Ci(0.3094013, 0.4979974);
        let exp_accept_hyp = Hyp::Null;
        let exp_z_accept_hyp = Hyp::Alt(AltHyp::Ne);

        check_binomial(
            n,
            n_s,
            p0,
            alt_hyp,
            exp_p,
            exp_z_p,
            exp_cp_ci,
            exp_ws_ci,
            exp_accept_hyp,
            exp_z_accept_hyp,
        );
    }

    #[test]
    fn test_binom_gt_100_40_05() {
        let (n, n_s) = (100, 40);
        let p0 = 0.5;
        let alt_hyp = AltHyp::Gt;
        let exp_p = 0.9824;
        let exp_z_p = 0.9772;
        let exp_cp_ci = Ci(0.317526, 1.000000);
        let exp_ws_ci = Ci(0.3230781, 1.0000000);
        let exp_accept_hyp = Hyp::Null;
        let exp_z_accept_hyp = exp_accept_hyp;

        check_binomial(
            n,
            n_s,
            p0,
            alt_hyp,
            exp_p,
            exp_z_p,
            exp_cp_ci,
            exp_ws_ci,
            exp_accept_hyp,
            exp_z_accept_hyp,
        );
    }

    #[test]
    fn test_binom_eq_100_40_04() {
        let (n, n_s) = (100, 40);
        let p0 = 0.4;
        let alt_hyp = AltHyp::Ne;
        let exp_p = 1.;
        let exp_z_p = 1.;
        let exp_cp_ci = Ci(0.3032948, 0.5027908);
        let exp_ws_ci = Ci(0.3094013, 0.4979974);
        let exp_accept_hyp = Hyp::Null;
        let exp_z_accept_hyp = exp_accept_hyp;

        check_binomial(
            n,
            n_s,
            p0,
            alt_hyp,
            exp_p,
            exp_z_p,
            exp_cp_ci,
            exp_ws_ci,
            exp_accept_hyp,
            exp_z_accept_hyp,
        );
    }

    #[test]
    fn test_binom_eq_100_4_0055() {
        let (n, n_s) = (100, 4);
        let p0 = 0.055;
        let alt_hyp = AltHyp::Ne;
        let exp_p = 0.6626;
        let exp_z_p = 0.5106;
        let exp_cp_ci = Ci(0.01100449, 0.09925716);
        let exp_ws_ci = Ci(0.01566330, 0.09837071);
        let exp_accept_hyp = Hyp::Null;
        let exp_z_accept_hyp = exp_accept_hyp;

        check_binomial(
            n,
            n_s,
            p0,
            alt_hyp,
            exp_p,
            exp_z_p,
            exp_cp_ci,
            exp_ws_ci,
            exp_accept_hyp,
            exp_z_accept_hyp,
        );
    }

    #[test]
    fn test_binom_eq_100_96_0945() {
        let (n, n_s) = (100, 96);
        let p0 = 0.945;
        let alt_hyp = AltHyp::Ne;
        let exp_p = 0.6626;
        let exp_z_p = 0.5106;
        let exp_cp_ci = Ci(0.9007428, 0.9889955);
        let exp_ws_ci = Ci(0.9016293, 0.9843367);
        let exp_accept_hyp = Hyp::Null;
        let exp_z_accept_hyp = exp_accept_hyp;

        check_binomial(
            n,
            n_s,
            p0,
            alt_hyp,
            exp_p,
            exp_z_p,
            exp_cp_ci,
            exp_ws_ci,
            exp_accept_hyp,
            exp_z_accept_hyp,
        );
    }

    //==================
    // Beta function corner cases.

    #[test]
    fn test_binom_eq_1_0_095() {
        let (n, n_s) = (1, 0);
        let p0 = 0.95;
        let alt_hyp = AltHyp::Ne;
        let exp_p = 0.05;
        let exp_cp_ci = Ci(0.000, 0.975);
        let exp_accept_hyp = Hyp::Null;

        check_binomial_no_z(n, n_s, p0, alt_hyp, exp_p, exp_cp_ci, exp_accept_hyp);
    }

    #[test]
    fn test_binom_eq_1_0_05() {
        let (n, n_s) = (1, 0);
        let p0 = 0.5;
        let alt_hyp = AltHyp::Ne;
        let exp_p = 1.;
        let exp_cp_ci = Ci(0.000, 0.975);
        let exp_accept_hyp = Hyp::Null;

        check_binomial_no_z(n, n_s, p0, alt_hyp, exp_p, exp_cp_ci, exp_accept_hyp);
    }

    #[test]
    fn test_binom_eq_1_1_05() {
        let (n, n_s) = (1, 1);
        let p0 = 0.5;
        let alt_hyp = AltHyp::Ne;
        let exp_p = 1.;
        let exp_cp_ci = Ci(0.025, 1.000);
        let exp_accept_hyp = Hyp::Null;

        check_binomial_no_z(n, n_s, p0, alt_hyp, exp_p, exp_cp_ci, exp_accept_hyp);
    }

    #[test]
    fn test_binom_eq_100000_0_05() {
        let (n, n_s) = (100000, 0);
        let p0 = 0.5;
        let alt_hyp = AltHyp::Ne;
        let exp_p = 2.2e-16;
        let exp_cp_ci = Ci(0.000000e+00, 3.688811e-05);
        let exp_accept_hyp = Hyp::Alt(AltHyp::Ne);

        check_binomial_no_z(n, n_s, p0, alt_hyp, exp_p, exp_cp_ci, exp_accept_hyp);
    }

    #[test]
    fn test_binom_eq_100000_1_05() {
        let (n, n_s) = (100000, 1);
        let p0 = 0.5;
        let alt_hyp = AltHyp::Ne;
        let exp_p = 2.2e-16;
        let exp_cp_ci = Ci(2.531780e-07, 5.571516e-05);
        let exp_accept_hyp = Hyp::Alt(AltHyp::Ne);

        check_binomial_no_z(n, n_s, p0, alt_hyp, exp_p, exp_cp_ci, exp_accept_hyp);
    }

    #[test]
    fn test_binom_eq_100000_99999_05() {
        let (n, n_s) = (100000, 99999);
        let p0 = 0.5;
        let alt_hyp = AltHyp::Ne;
        let exp_p = 2.2e-16;
        let exp_cp_ci = Ci(0.9999443, 0.9999997);
        let exp_accept_hyp = Hyp::Alt(AltHyp::Ne);

        check_binomial_no_z(n, n_s, p0, alt_hyp, exp_p, exp_cp_ci, exp_accept_hyp);
    }

    #[test]
    fn test_binom_eq_100000_100000_05() {
        let (n, n_s) = (100000, 100000);
        let p0 = 0.5;
        let alt_hyp = AltHyp::Ne;
        let exp_p = 2.2e-16;
        let exp_cp_ci = Ci(0.9999631, 1.0000000);
        let exp_accept_hyp = Hyp::Alt(AltHyp::Ne);

        check_binomial_no_z(n, n_s, p0, alt_hyp, exp_p, exp_cp_ci, exp_accept_hyp);
    }
}
