//! Used R's t.test function to generate expected values.
//! https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/t.test

#![cfg(feature = "_dev_utils")]
#![cfg(feature = "normal")]

use basic_stats::{
    aok::{AokBasicStats, AokFloat},
    core::{AcceptedHyp, AltHyp, Ci, SampleMoments},
    dev_utils::ApproxEq,
    normal::*,
};

const ALPHA: f64 = 0.05;
const EPSILON: f64 = 0.0005;

#[allow(clippy::too_many_arguments)]
fn check_welch(
    dataset_x: &[f64],
    dataset_y: &[f64],
    alt_hyp: AltHyp,
    exp_t: f64,
    exp_df: f64,
    exp_p: f64,
    exp_ci: Ci,
    exp_accept_hyp: AcceptedHyp,
) {
    let moments_x = SampleMoments::from_slice(dataset_x);
    let moments_y = SampleMoments::from_slice(dataset_y);

    let t = welch_t(&moments_x, &moments_y, 0.).aok();
    let df = welch_df(&moments_x, &moments_y).aok();
    let p = t_to_p(t, df, alt_hyp).aok();
    let ci = welch_alt_hyp_ci(&moments_x, &moments_y, alt_hyp, ALPHA).aok();
    let res = welch_test(&moments_x, &moments_y, 0., alt_hyp, ALPHA).aok();

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
        exp_ci.0.approx_eq(ci.0, EPSILON),
        "alt_hyp={alt_hyp:?} -- exp_ci.0={}, ci.0={}",
        exp_ci.0,
        ci.0
    );
    assert!(
        exp_ci.1.approx_eq(ci.1, EPSILON),
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
}

#[allow(clippy::too_many_arguments)]
fn check_student(
    dataset: &[f64],
    mu0: f64,
    alt_hyp: AltHyp,
    exp_t: f64,
    exp_df: f64,
    exp_p: f64,
    exp_ci: Ci,
    exp_accept_hyp: AcceptedHyp,
) {
    let moments = SampleMoments::from_slice(dataset);

    let t = student_1samp_t(&moments, mu0).aok();
    let df = student_1samp_df(&moments).aok();
    let p = t_to_p(t, df, alt_hyp).aok();
    let ci = student_1samp_alt_hyp_ci(&moments, alt_hyp, ALPHA).aok();
    let res = student_1samp_test(&moments, mu0, alt_hyp, ALPHA).aok();

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
        exp_ci.0.approx_eq(ci.0, EPSILON),
        "alt_hyp={alt_hyp:?} -- exp_ci.0={}, ci.0={}",
        exp_ci.0,
        ci.0
    );
    assert!(
        exp_ci.1.approx_eq(ci.1, EPSILON),
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
        check_welch(
            &a,
            &b,
            alt_hyp,
            exp_t,
            exp_df,
            exp_p,
            exp_ci,
            AcceptedHyp::Null,
        );
    }

    {
        let alt_hyp = AltHyp::Ne;
        let exp_p = 0.1413;
        let exp_ci = Ci(-10.453875, 1.614714);
        check_welch(
            &a,
            &b,
            alt_hyp,
            exp_t,
            exp_df,
            exp_p,
            exp_ci,
            AcceptedHyp::Null,
        );
    }

    {
        let alt_hyp = AltHyp::Gt;
        let exp_p = 0.9293;
        let exp_ci = Ci(-9.40084, f64::INFINITY);
        check_welch(
            &a,
            &b,
            alt_hyp,
            exp_t,
            exp_df,
            exp_p,
            exp_ci,
            AcceptedHyp::Null,
        );
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
        let exp_accept_hyp = AcceptedHyp::Null;
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
        );
    }

    {
        let alt_hyp = AltHyp::Ne;
        let exp_accept_hyp = AcceptedHyp::Alt;
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
        );
    }

    {
        let alt_hyp = AltHyp::Gt;
        let exp_accept_hyp = AcceptedHyp::Alt;
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
        );
    }
}

fn student_data() -> Vec<f64> {
    vec![
        20.70, 27.46, 22.15, 19.85, 21.29, 24.75, 20.75, 22.91, 25.34, 20.33, 21.54, 21.08, 22.14,
        19.56, 21.10, 18.04, 24.12, 19.95, 19.72, 18.28, 16.26, 17.46, 20.53, 22.12, 25.06, 22.44,
        19.08, 19.88, 21.39, 22.33, 25.79,
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
        let exp_accept_hyp = AcceptedHyp::Alt;
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
        );
    }

    {
        let alt_hyp = AltHyp::Ne;
        let exp_accept_hyp = AcceptedHyp::Alt;
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
        );
    }

    {
        let alt_hyp = AltHyp::Gt;
        let exp_accept_hyp = AcceptedHyp::Null;
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
        );
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
        let exp_accept_hyp = AcceptedHyp::Null;
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
        );
    }

    {
        let alt_hyp = AltHyp::Ne;
        let exp_accept_hyp = AcceptedHyp::Null;
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
        );
    }

    {
        let alt_hyp = AltHyp::Gt;
        let exp_accept_hyp = AcceptedHyp::Null;
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
        );
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
        let exp_accept_hyp = AcceptedHyp::Null;
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
        );
    }

    {
        let alt_hyp = AltHyp::Ne;
        let exp_accept_hyp = AcceptedHyp::Alt;
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
        );
    }

    {
        let alt_hyp = AltHyp::Gt;
        let exp_accept_hyp = AcceptedHyp::Alt;
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
        );
    }
}
