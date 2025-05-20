#![cfg(feature = "wilcoxon")]

mod nocover;

use basic_stats::{
    core::{AltHyp, StatsError},
    wilcoxon::RankSum,
};
use nocover::nocover;

// #[test]
// fn test_from_iters_with_counts() {
//     // RankSum::from_iters_with_counts covered by RankSum::from_iters.
// }

// #[test]
// fn test_from_iters() {
//     // RankSum::from_iters(it_x, it_y) covered by RankSum::from_slices
// }

#[test]
fn test_from_slices() {
    // Returns an error if a slice is not sorted in non-decreasing order.

    let good = [];
    let bad = [1., 2., 1.];

    assert!(RankSum::from_slices(&bad, &bad).is_err());
    assert!(RankSum::from_slices(&bad, &good).is_err());
    assert!(RankSum::from_slices(&good, &bad).is_err());
    if nocover() {
        assert!(RankSum::from_slices(&good, &good).is_ok());
    }
}

#[test]
fn test_z() {
    // Returns an error in any of the following conditions:
    // - `self.n_x == 0` or `self.n_y == 0`.
    // - There are too many rank ties between the two samples (causing an intermediate `NaN` value).
    //   This is hard to quantify a priori. For example,
    //   `x = [2., 2., 2., 2.]` and `y = [2., 2., 2., 3., 3.]` is OK
    //   but `x = [2., 2., 2., 2., 2.]` and `y = [2., 2., 2., 3., 3.]` results in an error.

    let s0 = [];
    let s11 = [1.];
    let s12 = [2.];
    let s22 = [2., 2.];
    let s42 = [2., 2., 2., 2.];
    let s52 = [2., 2., 2., 2., 2.];
    let smx = [2., 2., 2., 3., 3.];

    let rs0_0 = RankSum::from_slices(&s0, &s0).unwrap();
    let rs0_11 = RankSum::from_slices(&s0, &s11).unwrap();
    let rs11_0 = RankSum::from_slices(&s11, &s0).unwrap();

    assert!(rs0_0.z().is_err());
    assert!(rs0_11.z().is_err());
    assert!(rs11_0.z().is_err());

    let get_z = |name: &str, x: &[f64], y: &[f64]| -> Result<f64, StatsError> {
        let rs = RankSum::from_slices(x, y).unwrap();
        let z = rs.z();
        println!("{name}: z={z:?}");
        z
    };

    _ = get_z("S11-s11", &s11, &s11);
    _ = get_z("S11-s12", &s11, &s12);
    _ = get_z("S12-s22", &s12, &s22);
    _ = get_z("S22-s42", &s22, &s42);
    _ = get_z("S22-s52", &s22, &s52);
    _ = get_z("S22-smx", &s22, &smx);
    _ = get_z("S42-s52", &s42, &s52);
    _ = get_z("S42-smx", &s42, &smx);
    _ = get_z("S52-smx", &s52, &smx);

    assert!(get_z("S11-s11", &s11, &s11).is_err());
    assert!(get_z("S12-s22", &s12, &s22).is_err());
    assert!(get_z("S22-s42", &s22, &s42).is_err());
    assert!(get_z("S22-s52", &s22, &s52).is_err());
    assert!(get_z("S42-s52", &s42, &s52).is_err());
    assert!(get_z("S52-smx", &s52, &smx).is_err());

    if nocover() {
        assert!(get_z("S11-s12", &s11, &s12).unwrap().is_finite());
        assert!(get_z("S22-smx", &s22, &smx).unwrap().is_finite());
        assert!(get_z("S42-smx", &s42, &smx).unwrap().is_finite());
    }
}

// #[test]
// fn test_z_p() {
//     // RankSum::z_p covered by RankSum::z test.
// }

#[test]
fn test_z_test() {
    // Returns an error in any of these conditions:
    // - `self.n_x == 0` or `self.n_y == 0`.
    // - `alpha` not in interval `(0, 1)`.
    // - There are too many rank ties between the two samples (causing an intermediate `NaN` value).
    //   This is hard to quantify a priori. For example,
    //   `x = [2., 2., 2., 2.]` and `y = [2., 2., 2., 3., 3.]` is OK
    //   but `x = [2., 2., 2., 2., 2.]` and `y = [2., 2., 2., 3., 3.]` results in an error.

    let x = [1.];
    let y = [2.];

    let rs = RankSum::from_slices(&x, &y).unwrap();

    let alt_hyp = AltHyp::Ne;

    assert!(rs.z_test(alt_hyp, 0.).is_err());
    assert!(rs.z_test(alt_hyp, 1.).is_err());
    if nocover() {
        assert!(rs.z_test(alt_hyp, 0.5).is_ok());
    }
}
