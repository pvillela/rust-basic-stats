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
    let sample_a = vec![0.73, 0.80, 0.83, 1.04, 1.38, 1.45, 1.46, 1.64, 1.89, 1.91];
    let sample_b = vec![0.74, 0.88, 0.90, 1.15, 1.21];
    (sample_a, sample_b)
}

fn check_wilcoxon(rank_sum: &RankSum, alt_hyp: AltHyp, exp_accept_hyp: Hyp) {
    let w = rank_sum.w();
    let p = rank_sum.p(alt_hyp);
    let res = rank_sum.test(alt_hyp, ALPHA);

    println!("alt_hyp={alt_hyp:?} -- w={w}");

    assert_eq!(p, res.p(), "alt_hyp={alt_hyp:?} -- res.p");
    assert_eq!(ALPHA, res.alpha(), "alt_hyp={alt_hyp:?} -- res.alpha");
    assert_eq!(alt_hyp, res.alt_hyp(), "alt_hyp={alt_hyp:?} -- res.alt_hyp");
    assert_eq!(
        exp_accept_hyp,
        res.accepted(),
        "alt_hyp={alt_hyp:?} -- res.accepted"
    );

    let mann_whitney_u_a = rank_sum.mann_whitney_u_a();
    println!("alt_hyp={alt_hyp:?} -- mann_whitney_u_a={mann_whitney_u_a}");
    let mann_whitney_u_b = rank_sum.mann_whitney_u_b();
    println!("alt_hyp={alt_hyp:?} -- mann_whitney_u_b={mann_whitney_u_b}");
    let mann_whitney_u = rank_sum.mann_whitney_u();
    println!("alt_hyp={alt_hyp:?} -- mann_whitney_u={mann_whitney_u}");
}

fn test_book_data() -> Result<(), Box<dyn Error>> {
    let (dat_a, dat_b) = book_data();
    let rank_sum = RankSum::from_iter(dat_a.into_iter(), dat_b.into_iter())?;

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
