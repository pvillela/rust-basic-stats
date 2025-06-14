use basic_stats::{
    binomial::exact_binomial_test,
    core::{AcceptedHyp, AltHyp},
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
    assert_eq!(AcceptedHyp::Alt, test_res.accepted());
    println!("test result: {test_res:?}");
    // test result: HypTestResult { p: 0.028443966820498455, alpha: 0.05, alt_hyp: Lt, accepted: Alt }
}
