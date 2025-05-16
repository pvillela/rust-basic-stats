mod nocover;

use basic_stats::core::SampleMoments;
use nocover::nocover;

#[test]
pub fn test_mean() {
    let m0 = SampleMoments::new(0, 1., 2.);
    let m1 = SampleMoments::new(1, 1., 2.);

    assert!(m0.mean().is_err());
    if nocover() {
        assert!(m1.mean().unwrap().is_finite());
    }
}

#[test]
pub fn test_sum2_deviations() {
    let m0 = SampleMoments::new(0, 1., 2.);
    let m1 = SampleMoments::new(1, 1., 2.);

    assert!(m0.sum2_deviations().is_err());
    if nocover() {
        assert!(m1.sum2_deviations().unwrap().is_finite());
    }
}

#[test]
pub fn test_var() {
    let m0 = SampleMoments::new(0, 1., 2.);
    let m1 = SampleMoments::new(1, 1., 2.);
    let m2 = SampleMoments::new(2, 1., 2.);

    assert!(m0.var().is_err());
    assert!(m1.var().is_err());
    if nocover() {
        assert!(m2.var().unwrap().is_finite());
    }
}

#[test]
pub fn test_stdev() {
    let m0 = SampleMoments::new(0, 1., 2.);
    let m1 = SampleMoments::new(1, 1., 2.);
    let m2 = SampleMoments::new(2, 1., 2.);

    assert!(m0.stdev().is_err());
    assert!(m1.stdev().is_err());
    if nocover() {
        assert!(m2.stdev().unwrap().is_finite());
    }
}

#[test]
pub fn test_from_paired_slices() {
    let s0 = [];
    let s1 = [0.];

    let m00 = SampleMoments::from_paired_slices(&s0, &s0);
    let m01 = SampleMoments::from_paired_slices(&s0, &s1);
    let m10 = SampleMoments::from_paired_slices(&s1, &s0);
    let m11 = SampleMoments::from_paired_slices(&s1, &s1);

    assert!(m01.is_err());
    assert!(m10.is_err());
    if nocover() {
        assert!(m00.is_ok());
        assert!(m11.is_ok());
    }
}
