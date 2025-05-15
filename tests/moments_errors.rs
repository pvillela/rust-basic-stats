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
