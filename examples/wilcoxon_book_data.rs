//! Based on https://learning.oreilly.com/library/view/nonparametric-statistical-methods/9781118553299/9781118553299c04.xhtml#c04_level1_2
//! Nonparametric Statistical Methods, 3rd Edition, by Myles Hollander, Douglas A. Wolfe, Eric Chicken
//! Example 4.1.

use basic_stats::core::{AltHyp, Hyp};
use basic_stats::wilcoxon::RankSum;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    test_book_data()
}

const ALPHA: f64 = 0.05;

fn book_data() -> (Vec<f64>, Vec<f64>) {
    let sample_x = vec![0.73, 0.80, 0.83, 1.04, 1.38, 1.45, 1.46, 1.64, 1.89, 1.91];
    let sample_y = vec![0.74, 0.88, 0.90, 1.15, 1.21];
    (sample_x, sample_y)
}

fn check_wilcoxon(rank_sum: &RankSum, alt_hyp: AltHyp, exp_accept_hyp: Hyp) {
    let w = rank_sum.w();
    let p = rank_sum.z_p(alt_hyp).unwrap();
    let res = rank_sum.z_test(alt_hyp, ALPHA).unwrap();

    println!("alt_hyp={alt_hyp:?} -- w={w}");

    assert_eq!(p, res.p(), "alt_hyp={alt_hyp:?} -- res.p");
    assert_eq!(ALPHA, res.alpha(), "alt_hyp={alt_hyp:?} -- res.alpha");
    assert_eq!(alt_hyp, res.alt_hyp(), "alt_hyp={alt_hyp:?} -- res.alt_hyp");
    assert_eq!(
        exp_accept_hyp,
        res.accepted(),
        "alt_hyp={alt_hyp:?} -- res.accepted"
    );

    let mann_whitney_u_x = rank_sum.mann_whitney_u_y();
    println!("alt_hyp={alt_hyp:?} -- mann_whitney_u_x={mann_whitney_u_x}");
    let mann_whitney_u_y = rank_sum.mann_whitney_u_x();
    println!("alt_hyp={alt_hyp:?} -- mann_whitney_u_y={mann_whitney_u_y}");
    let mann_whitney_u = rank_sum.mann_whitney_u();
    println!("alt_hyp={alt_hyp:?} -- mann_whitney_u={mann_whitney_u}");
}

fn test_book_data() -> Result<(), Box<dyn Error>> {
    let (dat_x, dat_y) = book_data();
    let rank_sum = RankSum::from_iters(dat_x.into_iter(), dat_y.into_iter())?;

    {
        let alt_hyp = AltHyp::Lt;
        let exp_accept_hyp = Hyp::Null;
        check_wilcoxon(&rank_sum, alt_hyp, exp_accept_hyp);
    }

    {
        let alt_hyp = AltHyp::Ne;
        let exp_accept_hyp = Hyp::Null;
        check_wilcoxon(&rank_sum, alt_hyp, exp_accept_hyp);
    }

    {
        let alt_hyp = AltHyp::Gt;
        let exp_accept_hyp = Hyp::Null;
        check_wilcoxon(&rank_sum, alt_hyp, exp_accept_hyp);
    }

    Ok(())
}
