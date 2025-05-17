#![cfg(feature = "wilcoxon")]

mod nocover;

use basic_stats::{core::AltHyp, wilcoxon::RankSum};
use nocover::nocover;

#[test]
fn test_from_iters_with_counts() {
    // RankSum::from_iters_with_counts covered by RankSum::from_iters.
}

#[test]
fn test_from_iters() {
    // RankSum::from_iters(it_x, it_y) covered by RankSum::from_slices
}

#[test]
fn test_from_slices() {
    // Returns an error if a slice is not sorted in non-decreasing order.

    let good = [];
    let bad = [1., 2., 1.];

    assert!(RankSum::from_slices(&bad, &bad).is_err());
    assert!(RankSum::from_slices(&bad, &good).is_err());
    assert!(RankSum::from_slices(&good, &bad).is_err());
    if nocover() {
        assert!(RankSum::from_slices(&good, &bad).is_ok());
    }
}

#[test]
fn test_z() {
    // Returns an error if `self.n_x == 0` or `self.n_y == 0`.

    let s0 = [];
    let s1 = [1.];

    let rs00 = RankSum::from_slices(&s0, &s0).unwrap();
    let rs01 = RankSum::from_slices(&s0, &s1).unwrap();
    let rs10 = RankSum::from_slices(&s1, &s0).unwrap();
    let rs11 = RankSum::from_slices(&s1, &s1).unwrap();

    assert!(rs00.z().is_err());
    assert!(rs01.z().is_err());
    assert!(rs10.z().is_err());
    if nocover() {
        assert!(rs11.z().unwrap().is_finite());
    }
}

#[test]
fn test_z_p() {
    // Returns an error if `self.n_x == 0` or `self.n_y == 0`.

    let s0 = [];
    let s1 = [1.];

    let rs00 = RankSum::from_slices(&s0, &s0).unwrap();
    let rs01 = RankSum::from_slices(&s0, &s1).unwrap();
    let rs10 = RankSum::from_slices(&s1, &s0).unwrap();
    let rs11 = RankSum::from_slices(&s1, &s1).unwrap();

    let alt_hyp = AltHyp::Ne;

    assert!(rs00.z_p(alt_hyp).is_err());
    assert!(rs01.z_p(alt_hyp).is_err());
    assert!(rs10.z_p(alt_hyp).is_err());
    if nocover() {
        assert!(rs11.z_p(alt_hyp).unwrap().is_finite());
    }
}

#[test]
fn test_z_test() {
    // Returns an error in any of these circumstances:
    // - `self.n_x == 0` or `self.n_y == 0`.
    // - `alpha` not in `(0, 1)`.

    let s0 = [];
    let s1 = [1.];

    let rs00 = RankSum::from_slices(&s0, &s0).unwrap();
    let rs01 = RankSum::from_slices(&s0, &s1).unwrap();
    let rs10 = RankSum::from_slices(&s1, &s0).unwrap();
    let rs11 = RankSum::from_slices(&s1, &s1).unwrap();

    let alt_hyp = AltHyp::Ne;

    assert!(rs00.z_test(alt_hyp, 0.5).is_err());
    assert!(rs01.z_test(alt_hyp, 0.5).is_err());
    assert!(rs10.z_test(alt_hyp, 0.5).is_err());

    assert!(rs11.z_test(alt_hyp, 0.).is_err());
    assert!(rs11.z_test(alt_hyp, 1.).is_err());

    if nocover() {
        assert!(rs11.z_test(alt_hyp, 0.5).is_ok());
    }
}
