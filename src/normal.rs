//! Statistics related to the Normal distribution and Student's t distribution,
//! including the Student one-sample t-test and the
//! Welch two-sample t-test (for samples from distributions that may have different variances).
//!
//! This module is included by default. However, if `default-features = false` is specified in the dependency
//! declaration for this library, then inclusion of this module is gated by feature "**normal**".

use crate::core::{
    AltHyp, AsStatsResult, Ci, HypTestResult, SampleMoments, StatsError, StatsResult,
    check_alpha_in_closed_0_1, check_alpha_in_open_0_1,
};
use statrs::distribution::{ContinuousCDF, Normal, StudentsT};

/// Returns the the probability that the standard normal distribution will produce a more extreme value
/// than the argument `z`, with alternative hypothesis `alt_hyp`.
///
/// This function implements the z-test table, with `alt_hyp` defining whether the look-up is left-tailed,
/// right-tailed, or two-tailed.
pub fn z_to_p(z: f64, alt_hyp: AltHyp) -> f64 {
    let normal = Normal::standard();

    match alt_hyp {
        AltHyp::Lt => normal.cdf(z),
        AltHyp::Gt => normal.cdf(-z),
        AltHyp::Ne => normal.cdf(-z.abs()) * 2.,
    }
}

/// Returns the the probability that the Student distribution with location 0, scale 1, and `df` degrees of freedom
/// will produce a more extreme value than the argument `t`, with alternative hypothesis `alt_hyp`.
///
/// This function implements the t-test table, with `alt_hyp` defining whether the look-up is left-tailed,
/// right-tailed, or two-tailed.
///
/// # Errors
///
/// Returns an error if `df` is not `> 0`.
pub fn t_to_p(t: f64, df: f64, alt_hyp: AltHyp) -> StatsResult<f64> {
    let stud = StudentsT::new(0., 1., df).stats_result("arg `df` must be `> 0`")?;

    let value = match alt_hyp {
        AltHyp::Lt => stud.cdf(t),
        AltHyp::Gt => stud.cdf(-t),
        AltHyp::Ne => stud.cdf(-t.abs()) * 2.,
    };
    Ok(value)
}

/// Returns the value `v` for which `alpha` is the probability that the
/// standard normal distribution
/// is greater than `v`.
///
/// # Errors
///
/// Returns an error if `alpha` not in `(0, 1)`.
pub fn z_alpha(alpha: f64) -> StatsResult<f64> {
    check_alpha_in_open_0_1(alpha)?;

    let normal = Normal::standard();
    let value = normal.inverse_cdf(1. - alpha);
    Ok(value)
}

/// Returns the value `v` for which `alpha` is the probability that the
/// Student distribution with location 0, scale 1, and `df` degrees of freedom
/// is greater than `v`.
///
/// # Errors
///
/// Returns an error in any of the following conditions:
/// - `df` is not `> 0`.
/// - `alpha` not in `(0, 1)`.
pub fn t_alpha(df: f64, alpha: f64) -> StatsResult<f64> {
    check_alpha_in_open_0_1(alpha)?;

    let stud = StudentsT::new(0., 1., df).stats_result("degrees of freedom must be > 0")?;
    let value = stud.inverse_cdf(1. - alpha);
    Ok(value)
}

/// Welch's two-sample t statistic.
///
/// Arguments:
/// - `moments_x`: first sample's moments struct.
/// - `moments_y`: second sample's moments struct.
///
/// # Errors
///
/// Returns an error in any of the following conditions:
/// - `moments_x.n() <= 1`.
/// - `moments_y.n() <= 1`.
pub fn welch_t(moments_x: &SampleMoments, moments_y: &SampleMoments) -> StatsResult<f64> {
    let n_x = moments_x.nf();
    let n_y = moments_y.nf();
    let d_means = moments_x.mean()? - moments_y.mean()?;
    let s2_x = moments_x.stdev()?.powi(2);
    let s2_y = moments_y.stdev()?.powi(2);
    let s2_mean_x = s2_x / n_x;
    let s2_mean_y = s2_y / n_y;
    let s_d_means = (s2_mean_x + s2_mean_y).sqrt();
    Ok(d_means / s_d_means)
}

/// Degrees of freedom for Welch's two-sample t-test.
///
/// Arguments:
/// - `moments_x`: first sample's moments struct.
/// - `moments_y`: second sample's moments struct.
///
/// # Errors
///
/// Returns an error in any of the following conditions:
/// - `moments_x.n() <= 1`.
/// - `moments_y.n() <= 1`.
/// - `moments_x.stdev() == 0` and `moments_y.stdev() == 0`.
pub fn welch_df(moments_x: &SampleMoments, moments_y: &SampleMoments) -> StatsResult<f64> {
    if (moments_x.stdev()? + moments_y.stdev()?) == 0. {
        return Err(StatsError("sample standard deviations are zero"));
    }
    // At this point, df is guaranteed to be > 0.

    let n_x = moments_x.nf();
    let n_y = moments_y.nf();
    let s2_x = moments_x.stdev()?.powi(2);
    let s2_y = moments_y.stdev()?.powi(2);
    let s2_mean_x = s2_x / n_x;
    let s2_mean_y = s2_y / n_y;
    let numerator = (s2_mean_x + s2_mean_y).powi(2);
    let denominator = s2_mean_x.powi(2) / (n_x - 1.) + s2_mean_y.powi(2) / (n_y - 1.);
    let df = numerator / denominator;
    assert!(df > 0., "unexpected welch_df not `> 0.`");
    Ok(numerator / denominator)
}

/// p-value of Welch's two-sample t-test for equality.
///
/// Arguments:
/// - `moments_x`: first sample's moments struct.
/// - `moments_y`: second sample's moments struct.
/// - `alt_hyp`: alternative hypothesis.
///
/// # Errors
///
/// Returns an error in any of the following conditions:
/// - `moments_x.n() <= 1`.
/// - `moments_y.n() <= 1`.
/// - `moments_x.stdev() == 0` and `moments_y.stdev() == 0`.
pub fn welch_p(
    moments_x: &SampleMoments,
    moments_y: &SampleMoments,
    alt_hyp: AltHyp,
) -> StatsResult<f64> {
    let t = welch_t(moments_x, moments_y)?;
    let df = welch_df(moments_x, moments_y)?;
    t_to_p(t, df, alt_hyp)
}

/// Welch's confidence interval for the difference of means (μ(X) - μ(Y)) of two distributions.
///
/// Arguments:
/// - `moments_x`: first sample's moments struct.
/// - `moments_y`: second sample's moments struct.
/// - `alt_hyp`: alternative hypothesis.
/// - `alpha`: confidence level = `1 - alpha`.
///
/// # Errors
///
/// Returns an error in any of the following conditions:
/// - `moments_x.n() <= 1`.
/// - `moments_y.n() <= 1`.
/// - `moments_x.stdev() == 0` and `moments_y.stdev() == 0`.
/// - `alpha` not in `[0, 1]`.
pub fn welch_alt_hyp_ci(
    moments_x: &SampleMoments,
    moments_y: &SampleMoments,
    alt_hyp: AltHyp,
    alpha: f64,
) -> StatsResult<Ci> {
    let n_x = moments_x.nf();
    let n_y = moments_y.nf();
    let d_means = moments_x.mean()? - moments_y.mean()?;
    let s2_x = moments_x.stdev()?.powi(2);
    let s2_y = moments_y.stdev()?.powi(2);
    let s2_mean_x = s2_x / n_x;
    let s2_mean_y = s2_y / n_y;
    let df = welch_df(moments_x, moments_y)?;

    let stud = StudentsT::new(0., 1., df)
        .expect("StudentsT::new arg `freedom` should be guaranteed to be positive");
    let t0 = match alt_hyp {
        AltHyp::Ne => -stud.inverse_cdf(alpha / 2.),
        _ => -stud.inverse_cdf(alpha),
    };

    let mid = d_means;
    let delta = (s2_mean_x + s2_mean_y).sqrt() * t0;

    let value = match alt_hyp {
        AltHyp::Lt => Ci(-f64::INFINITY, mid + delta),
        AltHyp::Ne => Ci(mid - delta, mid + delta),
        AltHyp::Gt => Ci(mid - delta, f64::INFINITY),
    };
    Ok(value)
}

/// Welch's confidence interval for the difference of means (μ(X) - μ(Y)) of two distributions,
/// with the alternative hypothesis of inequality (two-sided).
///
/// Arguments:
/// - `moments_x`: first sample's moments struct.
/// - `moments_y`: second sample's moments struct.
/// - `alpha`: confidence level = `1 - alpha`.
///
/// # Errors
///
/// Returns an error in any of the following conditions:
/// - `moments_x.n() <= 1`.
/// - `moments_y.n() <= 1`.
/// - `moments_x.stdev() == 0` and `moments_y.stdev() == 0`.
/// - `alpha` not in `[0, 1]`.
pub fn welch_ci(
    moments_x: &SampleMoments,
    moments_y: &SampleMoments,
    alpha: f64,
) -> StatsResult<Ci> {
    welch_alt_hyp_ci(moments_x, moments_y, AltHyp::Ne, alpha)
}

/// Welch's two-sample t-test for equality of means of two distributions.
///
/// Arguments:
/// - `moments_x`: first sample's moments struct.
/// - `moments_y`: second sample's moments struct.
/// - `alt_hyp`: alternative hypothesis.
/// - `alpha`: confidence level = `1 - alpha`.
///
/// # Errors
///
/// Returns an error in any of the following conditions:
/// - `moments_x.n() <= 1`.
/// - `moments_y.n() <= 1`.
/// - `moments_x.stdev() == 0` and `moments_y.stdev() == 0`.
/// - `alpha` not in `[0, 1]`.
pub fn welch_test(
    moments_x: &SampleMoments,
    moments_y: &SampleMoments,
    alt_hyp: AltHyp,
    alpha: f64,
) -> StatsResult<HypTestResult> {
    check_alpha_in_closed_0_1(alpha)?;
    let p = welch_p(moments_x, moments_y, alt_hyp)?;
    Ok(HypTestResult::new(p, alpha, alt_hyp))
}

/// Student's one-sample t statistic.
///
/// Arguments:
/// - `moments`: sample moments struct.
/// - `mu0`: hypothesized distribution mean.
///
/// # Errors
///
/// Returns an error in any of the following conditions:
/// - `moments.n() <= 1`.
/// - `moments.stdev() == 0`.
pub fn student_one_sample_t(moments: &SampleMoments, mu0: f64) -> StatsResult<f64> {
    let n = moments.nf();
    let mean = moments.mean()?;
    let s = moments.stdev()?;
    Ok((mean - mu0) / s * n.sqrt())
}

/// Degrees of freedom for Student's one-sample t-test.
///
/// Arguments:
/// - `moments`: sample moments struct.
///
/// # Errors
///
/// Returns an error if `moments.n() <= 1`.
pub fn student_one_sample_df(moments: &SampleMoments) -> StatsResult<f64> {
    if moments.n() <= 1 {
        return Err(StatsError("`moments.n()` must be greater than 1"));
    }
    Ok(moments.nf() - 1.)
}

/// p-value of Student's one-sample t-test for equality.
///
/// Arguments:
/// - `moments`: sample moments struct.
/// - `mu0`: hypothesized distribution mean.
/// - `alt_hyp`: alternative hypothesis.
///
/// # Errors
///
/// Returns an error in any of the following conditions:
/// - `moments.n() <= 1`.
/// - `moments.stdev() == 0`.
pub fn student_one_sample_p(
    moments: &SampleMoments,
    mu0: f64,
    alt_hyp: AltHyp,
) -> StatsResult<f64> {
    let t = student_one_sample_t(moments, mu0)?;
    let df = student_one_sample_df(moments)?;
    t_to_p(t, df, alt_hyp)
}

/// Student's one-sample confidence interval for the distribution mean.
///
/// Arguments:
/// - `moments`: sample moments struct.
/// - `alt_hyp`: alternative hypothesis.
/// - `alpha`: confidence level = `1 - alpha`.
///
/// # Errors
///
/// Returns an error in any of the following conditions:
/// - `moments.n() <= 1`.
/// - `moments.stdev() == 0`.
/// - `alpha` not in `(0, 1)`.
pub fn student_one_sample_alt_hyp_ci(
    moments: &SampleMoments,
    alt_hyp: AltHyp,
    alpha: f64,
) -> StatsResult<Ci> {
    check_alpha_in_open_0_1(alpha)?;

    let df = student_one_sample_df(moments)?;

    let stud = StudentsT::new(0., 1., df)
        .expect("can't happen: degrees of freedom is always >= 3 by construction");
    let t0 = match alt_hyp {
        AltHyp::Ne => -stud.inverse_cdf(alpha / 2.),
        _ => -stud.inverse_cdf(alpha),
    };

    let mid = moments.mean()?;
    let delta = (moments.stdev()? / moments.nf().sqrt()) * t0;

    let ci = match alt_hyp {
        AltHyp::Lt => Ci(-f64::INFINITY, mid + delta),
        AltHyp::Ne => Ci(mid - delta, mid + delta),
        AltHyp::Gt => Ci(mid - delta, f64::INFINITY),
    };

    Ok(ci)
}

/// Student's one-sample confidence interval for the distribution mean,
/// with the alternative hypothesis of inequality (two-sided).
///
/// Arguments:
/// - `moments`: sample moments struct.
/// - `alpha`: confidence level = `1 - alpha`.
///
/// # Errors
///
/// Returns an error in any of the following conditions:
/// - `moments.n() <= 1`.
/// - `moments.stdev() == 0`.
/// - `alpha` not in `(0, 1)`.
pub fn student_one_sample_ci(moments: &SampleMoments, alpha: f64) -> StatsResult<Ci> {
    student_one_sample_alt_hyp_ci(moments, AltHyp::Ne, alpha)
}

/// Student's one-sample t-test for equality.
///
/// Arguments:
/// - `moments`: sample moments struct.
/// - `mu0`: hypothesized distribution mean.
/// - `alt_hyp`: alternative hypothesis.
/// - `alpha`: confidence level = `1 - alpha`.
///
/// # Errors
///
/// Returns an error in any of the following conditions:
/// - `moments.n() <= 1`.
/// - `moments.stdev() == 0`.
/// - `alpha` not in `(0, 1)`.
pub fn student_one_sample_test(
    moments: &SampleMoments,
    mu0: f64,
    alt_hyp: AltHyp,
    alpha: f64,
) -> StatsResult<HypTestResult> {
    check_alpha_in_open_0_1(alpha)?;
    let p = student_one_sample_p(moments, mu0, alt_hyp)?;
    Ok(HypTestResult::new(p, alpha, alt_hyp))
}

#[cfg(test)]
#[allow(clippy::too_many_arguments)]
mod test {
    //! Used R's t.test function to generate expected values.
    //! https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/t.test

    use super::*;
    use crate::{
        core::{AltHyp, Hyp},
        dev_utils::ApproxEq,
    };

    const ALPHA: f64 = 0.05;
    const EPSILON: f64 = 0.0005;

    fn check_welch(
        dataset_x: &[f64],
        dataset_y: &[f64],
        alt_hyp: AltHyp,
        exp_t: f64,
        exp_df: f64,
        exp_p: f64,
        exp_ci: Ci,
        exp_accept_hyp: Hyp,
    ) -> StatsResult<()> {
        let moments_x = SampleMoments::from_slice(dataset_x);
        let moments_y = SampleMoments::from_slice(dataset_y);

        let t = welch_t(&moments_x, &moments_y)?;
        let df = welch_df(&moments_x, &moments_y)?;
        let p = t_to_p(t, df, alt_hyp)?;
        let ci = welch_alt_hyp_ci(&moments_x, &moments_y, alt_hyp, ALPHA)?;
        let res = welch_test(&moments_x, &moments_y, alt_hyp, ALPHA)?;

        assert!(
            exp_t.approx_eq(t, EPSILON),
            "alt_hyp={alt_hyp:?} -- exp_t={exp_t}, t={t}"
        );
        assert!(
            exp_df.approx_eq(df, EPSILON),
            "alt_hyp={alt_hyp:?} -- exp_df={exp_df}, df={df}"
        );
        assert!(
            exp_p.approx_eq(p, EPSILON),
            "alt_hyp={alt_hyp:?} -- exp_p={exp_p}, p={p}"
        );
        assert!(
            exp_ci.0.approx_eq(ci.0, EPSILON) || exp_ci.0.is_infinite() && ci.0.is_infinite(),
            "alt_hyp={alt_hyp:?} -- exp_ci.0={}, ci.0={}",
            exp_ci.0,
            ci.0
        );
        assert!(
            exp_ci.1.approx_eq(ci.1, EPSILON) || exp_ci.1.is_infinite() && ci.1.is_infinite(),
            "alt_hyp={alt_hyp:?} -- exp_ci.1={}, ci.1={}",
            exp_ci.1,
            ci.1
        );

        assert_eq!(p, res.p(), "alt_hyp={alt_hyp:?} -- res.p");
        assert_eq!(ALPHA, res.alpha(), "alt_hyp={alt_hyp:?} -- res.alpha");
        assert_eq!(alt_hyp, res.alt_hyp(), "alt_hyp={alt_hyp:?} -- res.alt_hyp");
        assert_eq!(
            exp_accept_hyp,
            res.accepted(),
            "alt_hyp={alt_hyp:?} -- res.accepted"
        );

        Ok(())
    }

    fn check_student(
        dataset: &[f64],
        mu0: f64,
        alt_hyp: AltHyp,
        exp_t: f64,
        exp_df: f64,
        exp_p: f64,
        exp_ci: Ci,
        exp_accept_hyp: Hyp,
    ) -> StatsResult<()> {
        let moments = SampleMoments::from_slice(dataset);

        let t = student_one_sample_t(&moments, mu0)?;
        let df = student_one_sample_df(&moments)?;
        let p = t_to_p(t, df, alt_hyp)?;
        let ci = student_one_sample_alt_hyp_ci(&moments, alt_hyp, ALPHA)?;
        let res = student_one_sample_test(&moments, mu0, alt_hyp, ALPHA)?;

        assert!(
            exp_t.approx_eq(t, EPSILON),
            "alt_hyp={alt_hyp:?} -- exp_t={exp_t}, t={t}"
        );
        assert!(
            exp_df.approx_eq(df, EPSILON),
            "alt_hyp={alt_hyp:?} -- exp_df={exp_df}, df={df}"
        );
        assert!(
            exp_p.approx_eq(p, EPSILON),
            "alt_hyp={alt_hyp:?} -- exp_p={exp_p}, p={p}"
        );
        assert!(
            exp_ci.0.approx_eq(ci.0, EPSILON) || exp_ci.0.is_infinite() && ci.0.is_infinite(),
            "alt_hyp={alt_hyp:?} -- exp_ci.0={}, ci.0={}",
            exp_ci.0,
            ci.0
        );
        assert!(
            exp_ci.1.approx_eq(ci.1, EPSILON) || exp_ci.1.is_infinite() && ci.1.is_infinite(),
            "alt_hyp={alt_hyp:?} -- exp_ci.1={}, ci.1={}",
            exp_ci.1,
            ci.1
        );

        assert_eq!(p, res.p(), "alt_hyp={alt_hyp:?} -- res.p");
        assert_eq!(ALPHA, res.alpha(), "alt_hyp={alt_hyp:?} -- res.alpha");
        assert_eq!(alt_hyp, res.alt_hyp(), "alt_hyp={alt_hyp:?} -- res.alt_hyp");
        assert_eq!(
            exp_accept_hyp,
            res.accepted(),
            "alt_hyp={alt_hyp:?} -- res.accepted"
        );

        Ok(())
    }

    #[test]
    fn test_welch_eq() {
        let a = [14., 15., 15., 15., 16., 18., 22., 23., 24., 25., 25.];
        let b = [
            10., 12., 14., 15., 18., 22., 24., 27., 31., 33., 34., 34., 34.,
        ];

        let exp_t = -1.5379;
        let exp_df = 18.137;

        {
            let alt_hyp = AltHyp::Lt;
            let exp_p = 0.07067;
            let exp_ci = Ci(-f64::INFINITY, 0.5616789);
            check_welch(&a, &b, alt_hyp, exp_t, exp_df, exp_p, exp_ci, Hyp::Null).unwrap();
        }

        {
            let alt_hyp = AltHyp::Ne;
            let exp_p = 0.1413;
            let exp_ci = Ci(-10.453875, 1.614714);
            check_welch(&a, &b, alt_hyp, exp_t, exp_df, exp_p, exp_ci, Hyp::Null).unwrap();
        }

        {
            let alt_hyp = AltHyp::Gt;
            let exp_p = 0.9293;
            let exp_ci = Ci(-9.40084, f64::INFINITY);
            check_welch(&a, &b, alt_hyp, exp_t, exp_df, exp_p, exp_ci, Hyp::Null).unwrap();
        }
    }

    #[test]
    fn test_welch_gt() {
        let a = [24., 28., 32., 29., 35., 36., 30., 32., 25., 31.];
        let b = [5., 10., 25., 15., 16., 20.];

        let exp_t = 4.7857;
        let exp_df = 6.8409;

        {
            let alt_hyp = AltHyp::Lt;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.9989;
            let exp_ci = Ci(-f64::INFINITY, 21.00566);
            check_welch(
                &a,
                &b,
                alt_hyp,
                exp_t,
                exp_df,
                exp_p,
                exp_ci,
                exp_accept_hyp,
            )
            .unwrap();
        }

        {
            let alt_hyp = AltHyp::Ne;
            let exp_accept_hyp = Hyp::Alt(AltHyp::Ne);
            let exp_p = 0.00213;
            let exp_ci = Ci(7.57018, 22.49649);
            check_welch(
                &a,
                &b,
                alt_hyp,
                exp_t,
                exp_df,
                exp_p,
                exp_ci,
                exp_accept_hyp,
            )
            .unwrap();
        }

        {
            let alt_hyp = AltHyp::Gt;
            let exp_accept_hyp = Hyp::Alt(AltHyp::Gt);
            let exp_p = 0.001065;
            let exp_ci = Ci(9.061005, f64::INFINITY);
            check_welch(
                &a,
                &b,
                alt_hyp,
                exp_t,
                exp_df,
                exp_p,
                exp_ci,
                exp_accept_hyp,
            )
            .unwrap();
        }
    }

    fn student_data() -> Vec<f64> {
        vec![
            20.70, 27.46, 22.15, 19.85, 21.29, 24.75, 20.75, 22.91, 25.34, 20.33, 21.54, 21.08,
            22.14, 19.56, 21.10, 18.04, 24.12, 19.95, 19.72, 18.28, 16.26, 17.46, 20.53, 22.12,
            25.06, 22.44, 19.08, 19.88, 21.39, 22.33, 25.79,
        ]
    }
    #[test]
    fn test_student_lt() {
        let data = student_data();

        let mu0 = 23.;
        let exp_t = -3.505;
        let exp_df = 30.;

        {
            let alt_hyp = AltHyp::Lt;
            let exp_accept_hyp = Hyp::Alt(AltHyp::Lt);
            let exp_p = 0.0007288;
            let exp_ci = Ci(-f64::INFINITY, 22.17479);
            check_student(
                &data,
                mu0,
                alt_hyp,
                exp_t,
                exp_df,
                exp_p,
                exp_ci,
                exp_accept_hyp,
            )
            .unwrap();
        }

        {
            let alt_hyp = AltHyp::Ne;
            let exp_accept_hyp = Hyp::Alt(AltHyp::Ne);
            let exp_p = 0.001458;
            let exp_ci = Ci(20.46771, 22.33229);
            check_student(
                &data,
                mu0,
                alt_hyp,
                exp_t,
                exp_df,
                exp_p,
                exp_ci,
                exp_accept_hyp,
            )
            .unwrap();
        }

        {
            let alt_hyp = AltHyp::Gt;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.9993;
            let exp_ci = Ci(20.62521, f64::INFINITY);
            check_student(
                &data,
                mu0,
                alt_hyp,
                exp_t,
                exp_df,
                exp_p,
                exp_ci,
                exp_accept_hyp,
            )
            .unwrap();
        }
    }
    #[test]
    fn test_student_eq() {
        let data = student_data();

        let mu0 = 21.;
        let exp_t = 0.87624;
        let exp_df = 30.;

        {
            let alt_hyp = AltHyp::Lt;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.8061;
            let exp_ci = Ci(-f64::INFINITY, 22.17479);
            check_student(
                &data,
                mu0,
                alt_hyp,
                exp_t,
                exp_df,
                exp_p,
                exp_ci,
                exp_accept_hyp,
            )
            .unwrap();
        }

        {
            let alt_hyp = AltHyp::Ne;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.3879;
            let exp_ci = Ci(20.46771, 22.33229);
            check_student(
                &data,
                mu0,
                alt_hyp,
                exp_t,
                exp_df,
                exp_p,
                exp_ci,
                exp_accept_hyp,
            )
            .unwrap();
        }

        {
            let alt_hyp = AltHyp::Gt;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.1939;
            let exp_ci = Ci(20.62521, f64::INFINITY);
            check_student(
                &data,
                mu0,
                alt_hyp,
                exp_t,
                exp_df,
                exp_p,
                exp_ci,
                exp_accept_hyp,
            )
            .unwrap();
        }
    }

    #[test]
    fn test_student_gt() {
        let data = student_data();

        let mu0 = 20.;
        let exp_t = 3.0668;
        let exp_df = 30.;

        {
            let alt_hyp = AltHyp::Lt;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.9977;
            let exp_ci = Ci(-f64::INFINITY, 22.17479);
            check_student(
                &data,
                mu0,
                alt_hyp,
                exp_t,
                exp_df,
                exp_p,
                exp_ci,
                exp_accept_hyp,
            )
            .unwrap();
        }

        {
            let alt_hyp = AltHyp::Ne;
            let exp_accept_hyp = Hyp::Alt(AltHyp::Ne);
            let exp_p = 0.004553;
            let exp_ci = Ci(20.46771, 22.33229);
            check_student(
                &data,
                mu0,
                alt_hyp,
                exp_t,
                exp_df,
                exp_p,
                exp_ci,
                exp_accept_hyp,
            )
            .unwrap();
        }

        {
            let alt_hyp = AltHyp::Gt;
            let exp_accept_hyp = Hyp::Alt(AltHyp::Gt);
            let exp_p = 0.002276;
            let exp_ci = Ci(20.62521, f64::INFINITY);
            check_student(
                &data,
                mu0,
                alt_hyp,
                exp_t,
                exp_df,
                exp_p,
                exp_ci,
                exp_accept_hyp,
            )
            .unwrap();
        }
    }
}
