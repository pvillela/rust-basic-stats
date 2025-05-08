#![cfg(feature = "binomial")]

use basic_stats::{
    binomial::*,
    core::{AltHyp, Ci},
};

#[test]
fn test_bernoulli_p_hat() {
    // Returns an error if `n == 0`
    assert!(bernoulli_p_hat(0, 0).is_err());
    assert!(bernoulli_p_hat(1, 0).unwrap().is_finite());
}

#[test]
fn test_binomial_z() {
    // Returns an error in any of these circumstances:
    // - `n == 0`.
    // - `p0` not in `(0, 1)`.
    assert!(binomial_z(0, 0, 0.5).is_err(),);
    assert!(binomial_z(1, 0, 0.).is_err(),);
    assert!(binomial_z(1, 0, 0.5).unwrap().is_finite(),);
    assert!(binomial_z(1, 0, 1.).is_err(),);
}

#[test]
fn test_binomial_z_p() {
    // Returns an error in any of these circumstances:
    // - `n == 0`.
    // - `p0` not in `(0, 1)`.
    assert!(binomial_z_p(0, 0, 0.5, AltHyp::Ne).is_err(),);
    assert!(binomial_z_p(1, 0, 0., AltHyp::Ne).is_err(),);
    assert!(binomial_z_p(1, 0, 0.5, AltHyp::Ne).unwrap().is_finite(),);
    assert!(binomial_z_p(1, 0, 1., AltHyp::Ne).is_err(),);
}

#[test]
fn test_one_proportion_z_test() {
    // Returns an error in any of these circumstances:
    // - `n == 0`.
    // - `p0` not in `(0, 1)`.
    // - `alpha` not in `[0, 1]`.
    assert!(one_proportion_z_test(0, 0, 0.5, AltHyp::Ne, 0.5).is_err(),);
    assert!(one_proportion_z_test(1, 0, 0., AltHyp::Ne, 0.5).is_err(),);
    assert!(one_proportion_z_test(1, 0, 0.5, AltHyp::Ne, 0.).is_ok());
    assert!(one_proportion_z_test(1, 0, 0.5, AltHyp::Ne, 1.).is_ok());
    assert!(one_proportion_z_test(1, 0, 1., AltHyp::Ne, 0.5).is_err(),);
}

#[test]
fn test_binomial_ws_alt_hyp_ci() {
    // binomial_ws_alt_hyp_ci(n, n_s, alt_hyp, alpha) covered by binomial_ws_alt_hyp_ci
}

#[test]
fn test_binomial_ws_ci() {
    // Returns an error in any of these circumstances:
    // - `n == 0`.
    // - `alpha` not in `(0, 1)`.

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
    // Returns an error in any of these circumstances:
    // - `n == 0` or `n < n_s`.
    // - `alpha` is not in `[0, 1]`.

    {
        assert!(binomial_cp_ci(0, 0, 0.5).is_err());
        assert!(binomial_cp_ci(2, 3, 0.5).is_err());
        assert!(binomial_cp_ci(2, 1, -1.).is_err());
        assert!(binomial_cp_ci(2, 1, 2.).is_err());
    }

    {
        let Ci(lo, hi) = binomial_cp_ci(1, 0, 0.).unwrap();
        assert!(lo.is_finite());
        assert!(hi.is_finite());
    }

    {
        let Ci(lo, hi) = binomial_cp_ci(1, 0, 1.).unwrap();
        assert!(lo.is_finite());
        assert!(hi.is_finite());
    }

    {
        let Ci(lo, hi) = binomial_cp_ci(1, 1, 0.).unwrap();
        assert!(lo.is_finite());
        assert!(hi.is_finite());
    }

    {
        let Ci(lo, hi) = binomial_cp_ci(1, 1, 1.).unwrap();
        assert!(lo.is_finite());
        assert!(hi.is_finite());
    }
}

#[test]
fn test_exact_binomial_p() {
    // Returns an error in any of these circumstances:
    // - `n == 0` or `n < n_s`.
    // - `p0` is not in `[0, 1]`.

    assert!(exact_binomial_p(0, 0, 0.5, AltHyp::Ne).is_err());
    assert!(exact_binomial_p(2, 3, 0.5, AltHyp::Ne).is_err());
    assert!(exact_binomial_p(2, 1, -1., AltHyp::Ne).is_err());
    assert!(exact_binomial_p(2, 1, 2., AltHyp::Ne).is_err());

    assert!(exact_binomial_p(1, 0, 0., AltHyp::Ne).unwrap().is_finite());
    assert!(exact_binomial_p(1, 0, 1., AltHyp::Ne).unwrap().is_finite());
    assert!(exact_binomial_p(1, 1, 0., AltHyp::Ne).unwrap().is_finite());
    assert!(exact_binomial_p(1, 1, 1., AltHyp::Ne).unwrap().is_finite());
}

#[test]
fn test_exact_binomial_test() {
    // Returns an error in any of these circumstances:
    // - `n == 0` or `n < n_s`.
    // - `p0` is not in `[0, 1]`.
    // - `alpha` is not in `[0, 1]`.

    assert!(exact_binomial_test(0, 0, 0.5, AltHyp::Ne, 0.5).is_err());
    assert!(exact_binomial_test(2, 3, 0.5, AltHyp::Ne, 0.5).is_err());
    assert!(exact_binomial_test(2, 1, -1., AltHyp::Ne, 0.5).is_err());
    assert!(exact_binomial_test(2, 1, 2., AltHyp::Ne, 0.5).is_err());
    assert!(exact_binomial_test(2, 1, 0.5, AltHyp::Ne, -1.).is_err());
    assert!(exact_binomial_test(2, 1, 0.5, AltHyp::Ne, 2.).is_err());

    assert!(exact_binomial_test(1, 0, 0., AltHyp::Ne, 0.).is_ok());
    assert!(exact_binomial_test(1, 0, 1., AltHyp::Ne, 0.).is_ok());
    assert!(exact_binomial_test(1, 0, 0., AltHyp::Ne, 1.).is_ok());
    assert!(exact_binomial_test(1, 0, 1., AltHyp::Ne, 1.).is_ok());

    assert!(exact_binomial_test(1, 1, 0., AltHyp::Ne, 0.).is_ok());
    assert!(exact_binomial_test(1, 1, 1., AltHyp::Ne, 0.).is_ok());
    assert!(exact_binomial_test(1, 1, 0., AltHyp::Ne, 1.).is_ok());
    assert!(exact_binomial_test(1, 1, 1., AltHyp::Ne, 1.).is_ok());
}
