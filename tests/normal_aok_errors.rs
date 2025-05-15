mod nocover;

use basic_stats::{
    core::{AltHyp, AokBasicStats, AokBasicStatsValue, AokFloat, AokFloatValue, SampleMoments},
    normal::{
        student_1samp_ci, student_1samp_df, student_1samp_p, student_1samp_t, student_1samp_test,
        t_alpha, t_to_p, welch_ci, welch_df, welch_p, welch_t, welch_test, z_alpha,
    },
};
use nocover::nocover;

#[test]
fn test_z_to_p() {
    // `z_to_p` does not return a `Result`.
}

#[test]
fn test_t_to_p() {
    // Returns an error if `df` is not `> 0`.
    assert!(t_to_p(0., 0., AltHyp::Ne).aok().is_tainted());
    assert!(t_to_p(0., -1., AltHyp::Ne).aok().is_tainted());
    if nocover() {
        assert!(t_to_p(f64::MIN, 0.1, AltHyp::Ne).aok().is_untainted());
    }
}

#[test]
fn test_z_alpha() {
    // Returns an error if `alpha` not in `(0, 1)`.
    assert!(z_alpha(0.).aok().is_tainted());
    assert!(z_alpha(1.).aok().is_tainted());
    if nocover() {
        assert!(z_alpha(0.5).aok().is_untainted());
    }
}

#[test]
fn test_t_alpha() {
    // Returns an error in any of the following conditions:
    // - `df` is not `> 0`.
    // - `alpha` not in `(0, 1)`.
    assert!(t_alpha(0., 0.5).aok().is_tainted());
    assert!(t_alpha(-1., 0.5).aok().is_tainted());
    assert!(t_alpha(2., 0.).aok().is_tainted());
    assert!(t_alpha(2., 1.).aok().is_tainted());
    if nocover() {
        assert!(t_alpha(0.1, 0.5).aok().is_untainted());
    }
}

#[test]
fn test_welch_t() {
    // Returns an error in any of the following conditions:
    // - `moments_x.n() <= 1`.
    // - `moments_y.n() <= 1`.
    // - `moments_x.stdev() == 0` AND `moments_y.stdev() == 0`.

    let m0_0 = SampleMoments::default();
    let m1_1 = SampleMoments::new(1, 0., 1.);
    let m2_0 = SampleMoments::new(2, 0., 0.);
    let m2_1 = SampleMoments::new(2, 0., 1.);

    assert!(welch_t(&m0_0, &m2_1).aok().is_tainted());
    assert!(welch_t(&m2_1, &m0_0).aok().is_tainted());
    assert!(welch_t(&m1_1, &m2_1).aok().is_tainted());
    assert!(welch_t(&m2_1, &m1_1).aok().is_tainted());
    assert!(welch_t(&m2_0, &m2_0).aok().is_tainted());
    if nocover() {
        assert!(welch_t(&m2_0, &m2_1).aok().is_untainted());
        assert!(welch_t(&m2_1, &m2_0).aok().is_untainted());
    }
}

#[test]
fn test_welch_df() {
    // Returns an error in any of the following conditions:
    // - `moments_x.n() <= 1`.
    // - `moments_y.n() <= 1`.
    // - `moments_x.stdev() == 0` AND `moments_y.stdev() == 0`.

    let m0_0 = SampleMoments::default();
    let m1_1 = SampleMoments::new(1, 0., 1.);
    let m2_0 = SampleMoments::new(2, 0., 0.);
    let m2_1 = SampleMoments::new(2, 0., 1.);

    assert!(welch_df(&m0_0, &m2_1).aok().is_tainted());
    assert!(welch_df(&m2_1, &m0_0).aok().is_tainted());
    assert!(welch_df(&m1_1, &m2_1).aok().is_tainted());
    assert!(welch_df(&m2_1, &m1_1).aok().is_tainted());
    assert!(welch_df(&m2_0, &m2_0).aok().is_tainted());
    if nocover() {
        assert!(welch_df(&m2_0, &m2_1).aok().is_untainted());
        assert!(welch_df(&m2_1, &m2_0).aok().is_untainted());
    }
}

#[test]
fn test_welch_p() {
    // Returns an error in any of the following conditions:
    // - `moments_x.n() <= 1`.
    // - `moments_y.n() <= 1`.
    // - `moments_x.stdev() == 0` AND `moments_y.stdev() == 0`.

    let m0_0 = SampleMoments::default();
    let m1_1 = SampleMoments::new(1, 0., 1.);
    let m2_0 = SampleMoments::new(2, 0., 0.);
    let m2_1 = SampleMoments::new(2, 0., 1.);

    assert!(welch_p(&m0_0, &m2_1, AltHyp::Ne).aok().is_tainted());
    assert!(welch_p(&m2_1, &m0_0, AltHyp::Ne).aok().is_tainted());
    assert!(welch_p(&m1_1, &m2_1, AltHyp::Ne).aok().is_tainted());
    assert!(welch_p(&m2_1, &m1_1, AltHyp::Ne).aok().is_tainted());
    assert!(welch_p(&m2_0, &m2_0, AltHyp::Ne).aok().is_tainted());
    if nocover() {
        assert!(welch_p(&m2_0, &m2_1, AltHyp::Ne).aok().is_untainted());
        assert!(welch_p(&m2_1, &m2_0, AltHyp::Ne).aok().is_untainted());
    }
}

#[test]
fn test_welch_alt_hyp_ci() {
    // Covered by `test_welch_ci`.
}

#[test]
fn test_welch_ci() {
    // Returns an error in any of the following conditions:
    // - `moments_x.n() <= 1`.
    // - `moments_y.n() <= 1`.
    // - `moments_x.stdev() == 0` AND `moments_y.stdev() == 0`.
    // - `alpha` not in `(0, 1)`.

    let m0_0 = SampleMoments::default();
    let m1_1 = SampleMoments::new(1, 0., 1.);
    let m2_0 = SampleMoments::new(2, 0., 0.);
    let m2_1 = SampleMoments::new(2, 0., 1.);

    assert!(welch_ci(&m0_0, &m2_1, 0.5).aok().is_tainted());
    assert!(welch_ci(&m2_1, &m0_0, 0.5).aok().is_tainted());
    assert!(welch_ci(&m1_1, &m2_1, 0.5).aok().is_tainted());
    assert!(welch_ci(&m2_1, &m1_1, 0.5).aok().is_tainted());
    assert!(welch_ci(&m2_0, &m2_0, 0.5).aok().is_tainted());

    assert!(welch_ci(&m2_1, &m2_1, -1.).aok().is_tainted());
    assert!(welch_ci(&m2_1, &m2_1, 0.).aok().is_tainted());
    assert!(welch_ci(&m2_1, &m2_1, 1.).aok().is_tainted());
    assert!(welch_ci(&m2_1, &m2_1, 2.).aok().is_tainted());

    if nocover() {
        assert!(welch_ci(&m2_0, &m2_1, 0.5).aok().is_untainted());
        assert!(welch_ci(&m2_1, &m2_0, 0.5).aok().is_untainted());
    }
}

#[test]
fn test_welch_test() {
    // Returns an error in any of the following conditions:
    // - `moments_x.n() <= 1`.
    // - `moments_y.n() <= 1`.
    // - `moments_x.stdev() == 0` AND `moments_y.stdev() == 0`.
    // - `alpha` not in `(0, 1)`.

    let m0_0 = SampleMoments::default();
    let m1_1 = SampleMoments::new(1, 0., 1.);
    let m2_0 = SampleMoments::new(2, 0., 0.);
    let m2_1 = SampleMoments::new(2, 0., 1.);

    let alt_hyp = AltHyp::Ne;

    assert!(welch_test(&m0_0, &m2_1, alt_hyp, 0.5).aok().is_tainted());
    assert!(welch_test(&m2_1, &m0_0, alt_hyp, 0.5).aok().is_tainted());
    assert!(welch_test(&m1_1, &m2_1, alt_hyp, 0.5).aok().is_tainted());
    assert!(welch_test(&m2_1, &m1_1, alt_hyp, 0.5).aok().is_tainted());
    assert!(welch_test(&m2_0, &m2_0, alt_hyp, 0.5).aok().is_tainted());

    assert!(welch_test(&m2_1, &m2_1, alt_hyp, -1.).aok().is_tainted());
    assert!(welch_test(&m2_1, &m2_1, alt_hyp, 0.).aok().is_tainted());
    assert!(welch_test(&m2_1, &m2_1, alt_hyp, 1.).aok().is_tainted());
    assert!(welch_test(&m2_1, &m2_1, alt_hyp, 2.).aok().is_tainted());

    if nocover() {
        assert!(welch_test(&m2_0, &m2_1, alt_hyp, 0.5).aok().is_untainted());
        assert!(welch_test(&m2_1, &m2_0, alt_hyp, 0.5).aok().is_untainted());
    }
}

#[test]
fn test_student_1samp_t() {
    // Returns an error in any of the following conditions:
    // - `moments.n() <= 1`.
    // - `moments.stdev() == 0`.

    let m0_0 = SampleMoments::default();
    let m1_1 = SampleMoments::new(1, 0., 1.);
    let m2_0 = SampleMoments::new(2, 0., 0.);
    let m2_1 = SampleMoments::new(2, 0., 1.);

    assert!(student_1samp_t(&m0_0, 0.).aok().is_tainted());
    assert!(student_1samp_t(&m1_1, 0.).aok().is_tainted());
    assert!(student_1samp_t(&m2_0, 0.).aok().is_tainted());
    if nocover() {
        assert!(student_1samp_t(&m2_1, 0.).aok().is_untainted());
    }
}

#[test]
fn test_student_1samp_df() {
    // Returns an error if `moments.n() <= 1`.

    let m0_0 = SampleMoments::default();
    let m1_0 = SampleMoments::new(1, 0., 0.);
    let m2_0 = SampleMoments::new(2, 0., 0.);

    assert!(student_1samp_df(&m0_0).aok().is_tainted());
    assert!(student_1samp_df(&m1_0).aok().is_tainted());
    if nocover() {
        assert!(student_1samp_df(&m2_0).aok().is_untainted());
    }
}

#[test]
fn test_student_1samp_p() {
    // Returns an error in any of the following conditions:
    // - `moments.n() <= 1`.
    // - `moments.stdev() == 0`.

    let m0_0 = SampleMoments::default();
    let m1_1 = SampleMoments::new(1, 0., 1.);
    let m2_0 = SampleMoments::new(2, 0., 0.);
    let m2_1 = SampleMoments::new(2, 0., 1.);

    let alt_hyp = AltHyp::Ne;

    assert!(student_1samp_p(&m0_0, 0., alt_hyp).aok().is_tainted());
    assert!(student_1samp_p(&m1_1, 0., alt_hyp).aok().is_tainted());
    assert!(student_1samp_p(&m2_0, 0., alt_hyp).aok().is_tainted());
    if nocover() {
        assert!(student_1samp_p(&m2_1, 0., alt_hyp).aok().is_untainted());
    }
}

#[test]
fn test_student_1samp_alt_hyp_ci() {
    // Covered by `test_student_1samp_ci`.
}

#[test]
fn test_student_1samp_ci() {
    // Returns an error in any of the following conditions:
    // - `moments.n() <= 1`.
    // - `alpha` not in `(0, 1)`.

    let m0_0 = SampleMoments::default();
    let m1_1 = SampleMoments::new(1, 0., 1.);
    let m2_0 = SampleMoments::new(2, 0., 0.);

    assert!(student_1samp_ci(&m0_0, 0.5).aok().is_tainted());
    assert!(student_1samp_ci(&m1_1, 0.5).aok().is_tainted());

    assert!(student_1samp_ci(&m2_0, -1.).aok().is_tainted());
    assert!(student_1samp_ci(&m2_0, 0.).aok().is_tainted());
    assert!(student_1samp_ci(&m2_0, 1.).aok().is_tainted());
    assert!(student_1samp_ci(&m2_0, 2.).aok().is_tainted());

    if nocover() {
        assert!(student_1samp_ci(&m2_0, 0.5).aok().is_untainted());
    }
}

#[test]
fn test_student_1samp_test() {
    // Returns an error in any of the following conditions:
    // - `moments.n() <= 1`.
    // - `moments.stdev() == 0`.
    // - `alpha` not in `(0, 1)`.

    let m0_0 = SampleMoments::default();
    let m1_1 = SampleMoments::new(1, 0., 1.);
    let m2_0 = SampleMoments::new(2, 0., 0.);
    let m2_1 = SampleMoments::new(2, 0., 1.);

    let alt_hyp = AltHyp::Ne;

    assert!(
        student_1samp_test(&m0_0, 0., alt_hyp, 0.5)
            .aok()
            .is_tainted()
    );
    assert!(
        student_1samp_test(&m1_1, 0., alt_hyp, 0.5)
            .aok()
            .is_tainted()
    );
    assert!(
        student_1samp_test(&m2_0, 0., alt_hyp, 0.5)
            .aok()
            .is_tainted()
    );

    assert!(
        student_1samp_test(&m2_1, 0., alt_hyp, -1.)
            .aok()
            .is_tainted()
    );
    assert!(
        student_1samp_test(&m2_1, 0., alt_hyp, 0.)
            .aok()
            .is_tainted()
    );
    assert!(
        student_1samp_test(&m2_1, 0., alt_hyp, 1.)
            .aok()
            .is_tainted()
    );
    assert!(
        student_1samp_test(&m2_1, 0., alt_hyp, 2.)
            .aok()
            .is_tainted()
    );

    if nocover() {
        assert!(
            student_1samp_test(&m2_1, 0., alt_hyp, 0.5)
                .aok()
                .is_untainted()
        );
    }
}
