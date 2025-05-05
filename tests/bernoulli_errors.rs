#![cfg(feature = "bernoulli")]

use basic_stats::{
    bernoulli::*,
    core::{AltHyp, Ci},
};

#[test]
fn test_bernoulli_p_hat() {
    assert!(bernoulli_p_hat(0, 0).is_err(), "args: 0, 0");
    assert!(bernoulli_p_hat(1, 0).is_ok(), "args: 1, 0");
}

#[test]
fn test_bernoulli_normal_approx_z() {
    assert!(bernoulli_normal_approx_z(0, 0, 0.).is_err(),);
    assert!(bernoulli_normal_approx_z(0, 0, 0.5).is_err(),);
    assert!(bernoulli_normal_approx_z(0, 0, 1.).is_err(),);
    assert!(bernoulli_normal_approx_z(1, 0, 0.).is_err(),);
    assert!(bernoulli_normal_approx_z(1, 0, 0.5).unwrap().is_finite(),);
    assert!(bernoulli_normal_approx_z(1, 0, 1.).is_err(),);
}

#[test]
fn test_bernoulli_normal_approx_p() {
    assert!(bernoulli_normal_approx_p(0, 0, 0., AltHyp::Ne).is_err(),);
    assert!(bernoulli_normal_approx_p(0, 0, 0.5, AltHyp::Ne).is_err(),);
    assert!(bernoulli_normal_approx_p(0, 0, 1., AltHyp::Ne).is_err(),);
    assert!(bernoulli_normal_approx_p(1, 0, 0., AltHyp::Ne).is_err(),);
    assert!(
        bernoulli_normal_approx_p(1, 0, 0.5, AltHyp::Ne)
            .unwrap()
            .is_finite(),
    );
    assert!(bernoulli_normal_approx_p(1, 0, 1., AltHyp::Ne).is_err(),);
}

#[test]
fn test_binomial_ws_alt_hyp_ci() {
    // binomial_ws_alt_hyp_ci(n, n_s, alt_hyp, alpha) covered by binomial_ws_alt_hyp_ci
}

#[test]
fn test_binomial_ws_ci() {
    assert!(binomial_ws_ci(0, 0, 0.5).is_err());
    assert!(binomial_ws_ci(0, 1, 0.5).is_err());
    assert!(binomial_ws_ci(2, 1, 0.).is_err());
    assert!(binomial_ws_ci(2, 1, 1.).is_err());
    assert!(binomial_ws_ci(2, 1, f64::NAN).is_err());
    assert!(binomial_ws_ci(2, 1, f64::INFINITY).is_err());

    {
        let Ci(lo, hi) = binomial_ws_ci(2, 1, 0.5).unwrap();
        assert!(lo.is_finite());
        assert!(hi.is_finite());
    }
}

#[test]
fn test_binomial_cp_alt_hyp_ci() {
    // binomial_cp_alt_hyp_ci(n, n_s, alt_hyp, alpha) covered by binomial_cp_ci.
}

#[test]
fn test_binomial_cp_ci() {
    {
        assert!(binomial_cp_ci(0, 0, 0.5).is_err());
        assert!(binomial_cp_ci(2, 2, 0.5).is_err());
        assert!(binomial_cp_ci(3, 2, -1.).is_err());
        assert!(binomial_cp_ci(3, 2, 2.).is_err());
    }

    {
        let Ci(lo, hi) = binomial_cp_ci(3, 2, 0.).unwrap();
        assert!(lo.is_finite());
        assert!(hi.is_finite());
    }

    {
        let Ci(lo, hi) = binomial_cp_ci(3, 2, 1.).unwrap();
        assert!(lo.is_finite());
        assert!(hi.is_finite());
    }
}

#[test]
fn test_exact_binomial_p() {
    assert!(exact_binomial_p(0, 0, 0.5, AltHyp::Ne).is_err());
    assert!(exact_binomial_p(0, 1, 0.5, AltHyp::Ne).is_err());
    assert!(exact_binomial_p(2, 1, -1., AltHyp::Ne).is_err());
    assert!(exact_binomial_p(2, 1, 2., AltHyp::Ne).is_err());

    assert!(exact_binomial_p(1, 1, 0.5, AltHyp::Ne).unwrap().is_finite());
    assert!(exact_binomial_p(1, 1, 0., AltHyp::Ne).unwrap().is_finite());
    assert!(exact_binomial_p(1, 1, 1., AltHyp::Ne).unwrap().is_finite());
}

#[test]
fn test_exact_binomial_test() {
    assert!(exact_binomial_test(0, 0, 0.5, AltHyp::Ne, 0.5).is_err());
    assert!(exact_binomial_test(0, 1, 0.5, AltHyp::Ne, 0.5).is_err());
    assert!(exact_binomial_test(2, 1, -1., AltHyp::Ne, 0.5).is_err());
    assert!(exact_binomial_test(2, 1, 2., AltHyp::Ne, 0.5).is_err());
    assert!(exact_binomial_test(2, 1, 0.5, AltHyp::Ne, -1.).is_err());
    assert!(exact_binomial_test(2, 1, 0.5, AltHyp::Ne, 2.).is_err());

    assert!(exact_binomial_test(1, 1, 0.5, AltHyp::Ne, 0.5).is_ok());
    assert!(exact_binomial_test(1, 1, 0.5, AltHyp::Ne, 0.).is_ok());
    assert!(exact_binomial_test(1, 1, 0., AltHyp::Ne, 0.).is_ok());
    assert!(exact_binomial_test(1, 1, 1., AltHyp::Ne, 0.).is_ok());
    assert!(exact_binomial_test(1, 1, 0., AltHyp::Ne, 1.).is_ok());
    assert!(exact_binomial_test(1, 1, 1., AltHyp::Ne, 1.).is_ok());
}
