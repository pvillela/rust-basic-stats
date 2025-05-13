#![cfg(feature = "binomial")]

mod nocover;

use basic_stats::{
    binomial::*,
    core::{AltHyp, AokBasicStats, AokBasicStatsValue, AokFloat, AokFloatValue},
};
use nocover::nocover;

#[test]
fn test_bernoulli_p_hat() {
    // Returns an error if `n == 0`
    assert!(bernoulli_p_hat(0, 0).aok().is_tainted());
    if nocover() {
        assert!(bernoulli_p_hat(1, 0).aok().is_untainted());
    }
}

#[test]
fn test_binomial_z() {
    // Returns an error in any of these circumstances:
    // - `n == 0`.
    // - `p0` not in `(0, 1)`.
    assert!(binomial_z(0, 0, 0.5).aok().is_tainted());
    assert!(binomial_z(1, 0, 0.).aok().is_tainted());
    assert!(binomial_z(1, 0, 1.).aok().is_tainted());
    if nocover() {
        assert!(binomial_z(1, 0, 0.5).aok().is_untainted());
    }
}

#[test]
fn test_binomial_z_p() {
    // Returns an error in any of these circumstances:
    // - `n == 0`.
    // - `p0` not in `(0, 1)`.
    assert!(binomial_z_p(0, 0, 0.5, AltHyp::Ne).aok().is_tainted());
    assert!(binomial_z_p(1, 0, 0., AltHyp::Ne).aok().is_tainted());
    assert!(binomial_z_p(1, 0, 1., AltHyp::Ne).aok().is_tainted());
    if nocover() {
        assert!(binomial_z_p(1, 0, 0.5, AltHyp::Ne).aok().is_untainted());
    }
}

#[test]
fn test_one_proportion_z_test() {
    // Returns an error in any of these circumstances:
    // - `n == 0`.
    // - `p0` not in `(0, 1)`.
    // - `alpha` not in `(0, 1)`.
    assert!(
        one_proportion_z_test(0, 0, 0.5, AltHyp::Ne, 0.5)
            .aok()
            .is_tainted()
    );
    assert!(
        one_proportion_z_test(1, 0, 0., AltHyp::Ne, 0.5)
            .aok()
            .is_tainted()
    );
    assert!(
        one_proportion_z_test(1, 0, 1., AltHyp::Ne, 0.5)
            .aok()
            .is_tainted()
    );
    assert!(
        one_proportion_z_test(1, 0, 0.5, AltHyp::Ne, 0.)
            .aok()
            .is_tainted()
    );
    assert!(
        one_proportion_z_test(1, 0, 0.5, AltHyp::Ne, 1.)
            .aok()
            .is_tainted()
    );
    if nocover() {
        assert!(
            one_proportion_z_test(1, 0, 0.5, AltHyp::Ne, 0.5)
                .aok()
                .is_untainted()
        );
    }
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

    assert!(binomial_ws_ci(0, 0, 0.5).aok().is_tainted());
    assert!(binomial_ws_ci(0, 1, 0.5).aok().is_tainted());
    assert!(binomial_ws_ci(2, 1, 0.).aok().is_tainted());
    assert!(binomial_ws_ci(2, 1, 1.).aok().is_tainted());
    assert!(binomial_ws_ci(2, 1, f64::NAN).aok().is_tainted());
    assert!(binomial_ws_ci(2, 1, f64::INFINITY).aok().is_tainted());

    if nocover() {
        assert!(binomial_ws_ci(1, 1, 0.5).aok().is_untainted());
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
    // - `alpha` is not in `(0, 1)`.

    assert!(binomial_cp_ci(0, 0, 0.5).aok().is_tainted());
    assert!(binomial_cp_ci(2, 3, 0.5).aok().is_tainted());
    assert!(binomial_cp_ci(2, 1, 0.).aok().is_tainted());
    assert!(binomial_cp_ci(2, 1, 1.).aok().is_tainted());
}

#[test]
fn test_exact_binomial_p() {
    // Returns an error in any of these circumstances:
    // - `n == 0` or `n < n_s`.
    // - `p0` is not in `[0, 1]`.

    assert!(exact_binomial_p(0, 0, 0.5, AltHyp::Ne).aok().is_tainted());
    assert!(exact_binomial_p(2, 3, 0.5, AltHyp::Ne).aok().is_tainted());
    assert!(exact_binomial_p(2, 1, -1., AltHyp::Ne).aok().is_tainted());
    assert!(exact_binomial_p(2, 1, 2., AltHyp::Ne).aok().is_tainted());

    if nocover() {
        assert!(exact_binomial_p(1, 0, 0., AltHyp::Ne).aok().is_untainted());
        assert!(exact_binomial_p(1, 0, 1., AltHyp::Ne).aok().is_untainted());
        assert!(exact_binomial_p(1, 1, 0., AltHyp::Ne).aok().is_untainted());
        assert!(exact_binomial_p(1, 1, 1., AltHyp::Ne).aok().is_untainted());
    }
}

#[test]
fn test_exact_binomial_test() {
    // Returns an error in any of these circumstances:
    // - `n == 0` or `n < n_s`.
    // - `p0` is not in `[0, 1]`.
    // - `alpha` is not in `(0, 1)`.

    assert!(
        exact_binomial_test(0, 0, 0.5, AltHyp::Ne, 0.5)
            .aok()
            .is_tainted()
    );
    assert!(
        exact_binomial_test(2, 3, 0.5, AltHyp::Ne, 0.5)
            .aok()
            .is_tainted()
    );
    assert!(
        exact_binomial_test(2, 1, -1., AltHyp::Ne, 0.5)
            .aok()
            .is_tainted()
    );
    assert!(
        exact_binomial_test(2, 1, 2., AltHyp::Ne, 0.5)
            .aok()
            .is_tainted()
    );
    assert!(
        exact_binomial_test(2, 1, 0.5, AltHyp::Ne, 0.)
            .aok()
            .is_tainted()
    );
    assert!(
        exact_binomial_test(2, 1, 0.5, AltHyp::Ne, 1.)
            .aok()
            .is_tainted()
    );

    if nocover() {
        assert!(
            exact_binomial_test(1, 0, 0., AltHyp::Ne, 0.5)
                .aok()
                .is_untainted()
        );
        assert!(
            exact_binomial_test(1, 0, 1., AltHyp::Ne, 0.5)
                .aok()
                .is_untainted()
        );

        assert!(
            exact_binomial_test(1, 1, 0., AltHyp::Ne, 0.5)
                .aok()
                .is_untainted()
        );
        assert!(
            exact_binomial_test(1, 1, 1., AltHyp::Ne, 0.5)
                .aok()
                .is_untainted()
        );
    }
}
