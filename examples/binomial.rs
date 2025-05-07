use basic_stats::{
    binomial::exact_binomial_test,
    core::{AltHyp, Hyp},
};

const ALPHA: f64 = 0.05;

fn main() {
    let dat = [
        0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0,
        0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1,
        0, 1, 0, 0, 1, 1, 1, 0, 1, 0,
    ];

    let n = dat.len() as u64;
    let n_s = dat.into_iter().sum::<u64>();

    let test_res = exact_binomial_test(n, n_s, 0.5, AltHyp::Lt, ALPHA).unwrap();
    assert_eq!(Hyp::Alt(AltHyp::Lt), test_res.accepted());
    println!("test result: {test_res:?}");
    // test result: HypTestResult { p: 0.022750131947162633, alpha: 0.05, alt_hyp: Lt, accepted: Alt(Lt) }
}
