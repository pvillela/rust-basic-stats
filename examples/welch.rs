use basic_stats::{
    core::{AltHyp, Hyp, SampleMoments},
    normal::{welch_ci, welch_test},
};

const ALPHA: f64 = 0.05;

fn main() {
    let dat_x = [24., 28., 32., 29., 35., 36., 30., 32., 25., 31.];
    let dat_y = [5., 10., 25., 15., 16., 20.];

    let moments_x = SampleMoments::from_slice(&dat_x);
    let moments_y = SampleMoments::from_slice(&dat_y);
    {
        let test_res_lt = welch_test(&moments_x, &moments_y, AltHyp::Lt, ALPHA).unwrap();
        assert_eq!(Hyp::Null, test_res_lt.accepted());
        println!("test result lt: {test_res_lt:?}");
        // test result lt: HypTestResult { p: 0.9989352090850637, alpha: 0.05, alt_hyp: Lt, accepted: Null }
    }

    {
        let test_res_gt = welch_test(&moments_x, &moments_y, AltHyp::Gt, ALPHA).unwrap();
        assert_eq!(Hyp::Alt(AltHyp::Gt), test_res_gt.accepted());
        println!("test result gt: {test_res_gt:?}");
        // test result gt: HypTestResult { p: 0.0010647909149362593, alpha: 0.05, alt_hyp: Gt, accepted: Alt(Gt) }
    }

    {
        let ci = welch_ci(&moments_x, &moments_y, ALPHA).unwrap();
        println!("confidence interval for difference of means: {ci:?}");
        // confidence interval for difference of means: Ci(7.570179981992128, 22.496486684674537)
    }
}
