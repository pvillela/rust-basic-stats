use super::{
    core::{AltHyp, HypTestResult},
    normal::z_to_p,
};
use hdrhistogram::{
    Histogram,
    iterators::{HistogramIterator, IterationValue, recorded::Iter},
};

#[cfg(test)]
#[derive(Debug)]
#[allow(unused)] // for debugging only
struct RankedItem {
    value: u64,
    count: u64,
    rank: f64,
}

#[cfg(test)]
#[allow(unused)]
fn wilcoxon_ranked_items_ties_sum_prod(
    hist_a: &Histogram<u64>,
    hist_b: &Histogram<u64>,
) -> (Vec<RankedItem>, u64) {
    fn rank_item(
        #[cfg(test)] value: u64,
        count_i: u64,
        count_other: u64,
        iter_i: &mut HistogramIterator<u64, Iter>,
        item_opt_i: &mut Option<IterationValue<u64>>,
        prev_rank: f64,
        ties_sum_prod: &mut u64,
    ) -> (RankedItem, f64) {
        let count = count_i + count_other;
        let rank = prev_rank + (count as f64 + 1.) / 2.;
        let item = RankedItem {
            #[cfg(test)]
            value,
            count: count_i,
            rank,
        };
        let new_prev_rank = prev_rank + count as f64;
        *item_opt_i = iter_i.next();
        *ties_sum_prod += (count - 1) * count * (count + 1);
        (item, new_prev_rank)
    }

    let mut ties_sum_prod = 0;

    let ranked_items_b: Vec<RankedItem> = {
        let mut items_b = Vec::<RankedItem>::with_capacity(hist_b.distinct_values());
        let mut iter_a = hist_a.iter_recorded();
        let mut iter_b = hist_b.iter_recorded();
        let (mut item_a_opt, mut item_b_opt) = (iter_a.next(), iter_b.next());
        let mut prev_rank = 0.;

        loop {
            match (&mut item_a_opt, &mut item_b_opt) {
                (Some(item_a), None) => {
                    let count = item_a.count_at_value();
                    let (_, new_prev_rank) = rank_item(
                        #[cfg(test)]
                        item_a.value_iterated_to(),
                        count,
                        0,
                        &mut iter_a,
                        &mut item_a_opt,
                        prev_rank,
                        &mut ties_sum_prod,
                    );
                    prev_rank = new_prev_rank;
                }

                (Some(item_a), Some(item_b))
                    if item_a.value_iterated_to() < item_b.value_iterated_to() =>
                {
                    let (_, new_prev_rank) = rank_item(
                        #[cfg(test)]
                        item_a.value_iterated_to(),
                        item_a.count_at_value(),
                        0,
                        &mut iter_a,
                        &mut item_a_opt,
                        prev_rank,
                        &mut ties_sum_prod,
                    );
                    prev_rank = new_prev_rank;
                }

                (None, Some(item_b)) => {
                    let (ranked_item, new_prev_rank) = rank_item(
                        #[cfg(test)]
                        item_b.value_iterated_to(),
                        item_b.count_at_value(),
                        0,
                        &mut iter_b,
                        &mut item_b_opt,
                        prev_rank,
                        &mut ties_sum_prod,
                    );
                    items_b.push(ranked_item);
                    prev_rank = new_prev_rank;
                }

                (Some(item_a), Some(item_b))
                    if item_a.value_iterated_to() > item_b.value_iterated_to() =>
                {
                    let (ranked_item, new_prev_rank) = rank_item(
                        #[cfg(test)]
                        item_b.value_iterated_to(),
                        item_b.count_at_value(),
                        0,
                        &mut iter_b,
                        &mut item_b_opt,
                        prev_rank,
                        &mut ties_sum_prod,
                    );
                    items_b.push(ranked_item);
                    prev_rank = new_prev_rank;
                }

                // if item_a.value_iterated_to() == item_b.value_iterated_to()
                (Some(item_a), Some(item_b)) => {
                    let count_a = item_a.count_at_value();
                    let count_b = item_b.count_at_value();
                    #[cfg(test)]
                    rank_item(
                        item_a.value_iterated_to(),
                        count_a,
                        count_b,
                        &mut iter_a,
                        &mut item_a_opt,
                        prev_rank,
                        &mut ties_sum_prod,
                    );
                    let (ranked_item, new_prev_rank) = rank_item(
                        #[cfg(test)]
                        item_b.value_iterated_to(),
                        count_b,
                        count_a,
                        &mut iter_b,
                        &mut item_b_opt,
                        prev_rank,
                        &mut ties_sum_prod,
                    );
                    items_b.push(ranked_item);
                    prev_rank = new_prev_rank;
                }

                (None, None) => break,
            }
        }

        items_b
    };

    (ranked_items_b, ties_sum_prod)
}

fn wilcoxon_rank_sum_ties_sum_prod(hist_a: &Histogram<u64>, hist_b: &Histogram<u64>) -> (f64, f64) {
    fn rank_item(
        count_i: u64,
        count_other: u64,
        iter_i: &mut HistogramIterator<u64, Iter>,
        item_opt_i: &mut Option<IterationValue<u64>>,
        prev_rank: f64,
        ties_sum_prod: &mut u64,
    ) -> (f64, f64) {
        let count = count_i + count_other;
        let rank = prev_rank + (count as f64 + 1.) / 2.;
        let rank_sum = count_i as f64 * rank;
        let new_prev_rank = prev_rank + count as f64;
        *item_opt_i = iter_i.next();
        *ties_sum_prod += (count - 1) * count * (count + 1);
        (rank_sum, new_prev_rank)
    }

    let mut ties_sum_prod = 0;

    let rank_sum_b: f64 = {
        let mut rank_sum_a = 0.;
        let mut rank_sum_b = 0.;
        let mut iter_a = hist_a.iter_recorded();
        let mut iter_b = hist_b.iter_recorded();
        let (mut item_a_opt, mut item_b_opt) = (iter_a.next(), iter_b.next());
        let mut prev_rank = 0.;

        loop {
            match (&mut item_a_opt, &mut item_b_opt) {
                (Some(item_a), None) => {
                    let count = item_a.count_at_value();
                    let (rank_sum, new_prev_rank) = rank_item(
                        count,
                        0,
                        &mut iter_a,
                        &mut item_a_opt,
                        prev_rank,
                        &mut ties_sum_prod,
                    );
                    rank_sum_a += rank_sum;
                    prev_rank = new_prev_rank;
                }

                (Some(item_a), Some(item_b))
                    if item_a.value_iterated_to() < item_b.value_iterated_to() =>
                {
                    let (rank_sum, new_prev_rank) = rank_item(
                        item_a.count_at_value(),
                        0,
                        &mut iter_a,
                        &mut item_a_opt,
                        prev_rank,
                        &mut ties_sum_prod,
                    );
                    rank_sum_a += rank_sum;
                    prev_rank = new_prev_rank;
                }

                (None, Some(item_b)) => {
                    let (rank_sum, new_prev_rank) = rank_item(
                        item_b.count_at_value(),
                        0,
                        &mut iter_b,
                        &mut item_b_opt,
                        prev_rank,
                        &mut ties_sum_prod,
                    );
                    rank_sum_b += rank_sum;
                    prev_rank = new_prev_rank;
                }

                (Some(item_a), Some(item_b))
                    if item_a.value_iterated_to() > item_b.value_iterated_to() =>
                {
                    let (rank_sum, new_prev_rank) = rank_item(
                        item_b.count_at_value(),
                        0,
                        &mut iter_b,
                        &mut item_b_opt,
                        prev_rank,
                        &mut ties_sum_prod,
                    );
                    rank_sum_b += rank_sum;
                    prev_rank = new_prev_rank;
                }

                // if item_a.value_iterated_to() == item_b.value_iterated_to()
                (Some(item_a), Some(item_b)) => {
                    let count_a = item_a.count_at_value();
                    let count_b = item_b.count_at_value();

                    let (rank_sum, _) = rank_item(
                        count_a,
                        count_b,
                        &mut iter_a,
                        &mut item_a_opt,
                        prev_rank,
                        &mut ties_sum_prod,
                    );
                    rank_sum_a += rank_sum;

                    let (rank_sum, new_prev_rank) = rank_item(
                        count_b,
                        count_a,
                        &mut iter_b,
                        &mut item_b_opt,
                        prev_rank,
                        &mut ties_sum_prod,
                    );
                    rank_sum_b += rank_sum;
                    prev_rank = new_prev_rank;
                }

                (None, None) => break,
            }
        }

        // Check rank-sum calculation.
        {
            let n_a = hist_a.len() as f64;
            let n_b = hist_b.len() as f64;
            let expected_rank_sum_a = (1. + n_a + n_b) * (n_a + n_b) / 2. - rank_sum_b;
            debug_assert_eq!(expected_rank_sum_a, rank_sum_a, "rank_sum_a check");
        }

        rank_sum_b
    };

    (rank_sum_b, ties_sum_prod as f64)
}

pub fn wilcoxon_rank_sum_w(hist_a: &Histogram<u64>, hist_b: &Histogram<u64>) -> f64 {
    wilcoxon_rank_sum_ties_sum_prod(hist_a, hist_b).0
}

/// The `w` value computed by `R`'s `wilcox.test` function, which is the Mann-Whitney U for the
/// first sample (`hist_a`).
///
/// See explanation in the book Nonparametric Statistical Methods, 3rd Edition,
/// by Myles Hollander, Douglas A. Wolfe, Eric Chicken, Example 4.1.
#[cfg(test)]
pub fn wilcoxon_rank_sum_r_w(hist_a: &Histogram<u64>, hist_b: &Histogram<u64>) -> f64 {
    mann_whitney_u_a(hist_a, hist_b)
}

#[cfg(test)]
fn mann_whitney_u_b(hist_a: &Histogram<u64>, hist_b: &Histogram<u64>) -> f64 {
    let w = wilcoxon_rank_sum_w(hist_a, hist_b);
    let n_b = hist_b.len() as f64;
    w - n_b * (n_b + 1.) / 2.
}

#[cfg(test)]
fn mann_whitney_u_a(hist_a: &Histogram<u64>, hist_b: &Histogram<u64>) -> f64 {
    let n_a = hist_a.len() as f64;
    let n_b = hist_b.len() as f64;
    (n_a * n_b) - mann_whitney_u_b(hist_a, hist_b)
}

#[cfg(test)]
fn mann_whitney_u(hist_a: &Histogram<u64>, hist_b: &Histogram<u64>) -> f64 {
    mann_whitney_u_b(hist_a, hist_b).min(mann_whitney_u_a(hist_a, hist_b))
}

pub fn wilcoxon_rank_sum_z(hist_a: &Histogram<u64>, hist_b: &Histogram<u64>) -> f64 {
    let n_a = hist_a.len() as f64;
    let n_b = hist_b.len() as f64;
    let (w, ties_sum_prod) = wilcoxon_rank_sum_ties_sum_prod(hist_a, hist_b);
    let e0_w = n_b * (n_a + n_b + 1.) / 2.;
    let var0_w_base = n_a * n_b * (n_a + n_b + 1.) / 12.;
    let var0_w_ties_adjust = n_a * n_b / (12. * (n_a + n_b) * (n_a + n_b - 1.)) * ties_sum_prod;
    let var0_w = var0_w_base - var0_w_ties_adjust;
    let w_star = (w - e0_w) / var0_w.sqrt();

    -w_star
}

#[cfg(test)]
fn wilcoxon_rank_sum_z_no_ties_adjust(hist_a: &Histogram<u64>, hist_b: &Histogram<u64>) -> f64 {
    let n_a = hist_a.len() as f64;
    let n_b = hist_b.len() as f64;
    let w = wilcoxon_rank_sum_w(hist_a, hist_b);
    let e0_w = n_b * (n_a + n_b + 1.) / 2.;
    let var0_w_base = n_a * n_b * (n_a + n_b + 1.) / 12.;
    let var0_w_ties_adjust = 0.;
    let var0_w = var0_w_base - var0_w_ties_adjust;
    let w_star = (w - e0_w) / var0_w.sqrt();

    -w_star
}

pub fn wilcoxon_rank_sum_p(
    hist_a: &Histogram<u64>,
    hist_b: &Histogram<u64>,
    alt_hyp: AltHyp,
) -> f64 {
    let z = wilcoxon_rank_sum_z(hist_a, hist_b);
    z_to_p(z, alt_hyp)
}

#[cfg(test)]
fn wilcoxon_rank_sum_p_no_ties_adjust(
    hist_a: &Histogram<u64>,
    hist_b: &Histogram<u64>,
    alt_hyp: AltHyp,
) -> f64 {
    let z = wilcoxon_rank_sum_z_no_ties_adjust(hist_a, hist_b);
    z_to_p(z, alt_hyp)
}

pub fn wilcoxon_rank_sum_test(
    hist_a: &Histogram<u64>,
    hist_b: &Histogram<u64>,
    alt_hyp: AltHyp,
    alpha: f64,
) -> HypTestResult {
    let p = wilcoxon_rank_sum_p(hist_a, hist_b, alt_hyp);
    HypTestResult::new(p, alpha, alt_hyp)
}

#[cfg(test)]
pub fn wilcoxon_rank_sum_test_no_ties_adjust(
    hist_a: &Histogram<u64>,
    hist_b: &Histogram<u64>,
    alt_hyp: AltHyp,
    alpha: f64,
) -> HypTestResult {
    let p = wilcoxon_rank_sum_p_no_ties_adjust(hist_a, hist_b, alt_hyp);
    HypTestResult::new(p, alpha, alt_hyp)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod base_test {
    //! Tests other than `test_w` used R's wilcox.test function to generate expected results.
    //! https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/wilcox.test

    use super::*;
    use crate::{
        core::{AltHyp, Hyp},
        dev_utils::ApproxEq,
    };
    use hdrhistogram::Histogram;

    const ALPHA: f64 = 0.05;
    const EPSILON: f64 = 0.0005;

    fn book_data() -> (Vec<f64>, Vec<f64>) {
        let sample_a = vec![0.73, 0.80, 0.83, 1.04, 1.38, 1.45, 1.46, 1.64, 1.89, 1.91];
        let sample_b = vec![0.74, 0.88, 0.90, 1.15, 1.21];
        (sample_a, sample_b)
    }

    fn contrived_data() -> (Vec<f64>, Vec<f64>) {
        let sample_a = vec![
            85., 90., 78., 92., 88., 76., 95., 89., 91., 82., 115., 120., 108., 122., 118., 106.,
            125., 119., 121., 112., 145., 150., 138., 152., 148., 136., 155., 149., 151., 142.,
            175., 180., 168., 182., 178., 166., 185., 179., 181., 172., 205., 210., 198., 212.,
            208., 196., 215., 209., 211., 202.,
        ];
        let sample_b = vec![
            70., 85., 80., 90., 75., 88., 92., 79., 86., 81., 92., 100., 115., 110., 120., 105.,
            118., 122., 109., 116., 111., 122., 130., 145., 140., 150., 135., 148., 152., 139.,
            146., 141., 152., 160., 175., 170., 180., 165., 178., 182., 169., 176., 171., 182.,
            190., 205., 200., 210., 195., 208., 212., 199., 206., 201., 212.,
        ];
        (sample_a, sample_b)
    }

    fn shifted_contrived_data() -> (Vec<f64>, Vec<f64>) {
        let (sample_a, sample_b) = contrived_data();
        let sample_b = sample_b.into_iter().map(|v| v + 35.).collect::<Vec<_>>();
        (sample_a, sample_b)
    }

    fn data_hists(
        (data_a, data_b): (Vec<f64>, Vec<f64>),
        factor: u64,
        hist_max: u64,
        hist_sigfig: u8,
    ) -> (Histogram<u64>, Histogram<u64>) {
        let mut hist_a = Histogram::new_with_max(hist_max, hist_sigfig).unwrap();
        let mut hist_b = Histogram::new_from(&hist_a);

        for v in data_a {
            hist_a.record((v * factor as f64) as u64).unwrap();
        }

        for v in data_b {
            hist_b.record((v * factor as f64) as u64).unwrap();
        }

        (hist_a, hist_b)
    }

    #[test]
    /// Based on https://learning.oreilly.com/library/view/nonparametric-statistical-methods/9781118553299/9781118553299c04.xhtml#c04_level1_2
    /// Nonparametric Statistical Methods, 3rd Edition, by Myles Hollander, Douglas A. Wolfe, Eric Chicken
    /// Example 4.1.
    fn test_w() {
        let (hist_a, hist_b) = data_hists(book_data(), 100, 200, 3);

        let expected_w = 30.;
        let actual_w = wilcoxon_rank_sum_w(&hist_a, &hist_b);
        assert_eq!(expected_w, actual_w, "w comparison");

        let expected_r_w = 35.;
        let actual_r_w = wilcoxon_rank_sum_r_w(&hist_a, &hist_b);
        assert_eq!(expected_r_w, actual_r_w, "R w comparison");

        let expected_p_correct = 0.2544; // R: // R: wilcox.test(a, b)
        let expected_p = 0.2207; // R: wilcox.test(a, b, exact=FALSE, correct=FALSE)
        let actual_p = wilcoxon_rank_sum_p(&hist_a, &hist_b, AltHyp::Ne);
        println!(
            "expected_p_correct={expected_p_correct}, expected_p={expected_p}, actual_p={actual_p}"
        );
        assert!(expected_p.approx_eq(actual_p, EPSILON), "p comparison");
    }

    fn check_wilcoxon(
        hist_a: &Histogram<u64>,
        hist_b: &Histogram<u64>,
        alt_hyp: AltHyp,
        exp_r_w: f64,
        exp_p: f64,
        exp_accept_hyp: Hyp,
    ) {
        let w = wilcoxon_rank_sum_w(hist_a, hist_b);
        let r_w = wilcoxon_rank_sum_r_w(hist_a, hist_b);
        let p = wilcoxon_rank_sum_p(hist_a, hist_b, alt_hyp);
        let res = wilcoxon_rank_sum_test(hist_a, hist_b, alt_hyp, ALPHA);

        println!("alt_hyp={alt_hyp:?} -- w={w}");
        assert!(
            exp_r_w.approx_eq(r_w, EPSILON),
            "alt_hyp={alt_hyp:?} -- exp_r_w={exp_r_w}, r_w={r_w}"
        );
        assert!(
            exp_p.approx_eq(p, EPSILON),
            "alt_hyp={alt_hyp:?} -- exp_p={exp_p}, p={p}"
        );

        assert_eq!(p, res.p(), "alt_hyp={alt_hyp:?} -- res.p");
        assert_eq!(ALPHA, res.alpha(), "alt_hyp={alt_hyp:?} -- res.alpha");
        assert_eq!(alt_hyp, res.alt_hyp(), "alt_hyp={alt_hyp:?} -- res.alt_hyp");
        assert_eq!(
            exp_accept_hyp,
            res.accepted(),
            "alt_hyp={alt_hyp:?} -- res.accepted"
        );

        let mann_whitney_u_a = mann_whitney_u_a(hist_a, hist_b);
        println!("alt_hyp={alt_hyp:?} -- mann_whitney_u_a={mann_whitney_u_a}");
        let mann_whitney_u_b = mann_whitney_u_b(hist_a, hist_b);
        println!("alt_hyp={alt_hyp:?} -- mann_whitney_u_b={mann_whitney_u_b}");
        let mann_whitney_u = mann_whitney_u(hist_a, hist_b);
        println!("alt_hyp={alt_hyp:?} -- mann_whitney_u={mann_whitney_u}");
    }

    #[test]
    /// Based on https://learning.oreilly.com/library/view/nonparametric-statistical-methods/9781118553299/9781118553299c04.xhtml#c04_level1_2
    /// Nonparametric Statistical Methods, 3rd Edition, by Myles Hollander, Douglas A. Wolfe, Eric Chicken
    /// Example 4.1.
    fn test_book_data() {
        let (hist_a, hist_b) = data_hists(book_data(), 100, 200, 3);

        let exp_r_w = 35.;

        {
            let alt_hyp = AltHyp::Lt;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.8897;
            check_wilcoxon(&hist_a, &hist_b, alt_hyp, exp_r_w, exp_p, exp_accept_hyp);
        }

        {
            let alt_hyp = AltHyp::Ne;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.2207;
            check_wilcoxon(&hist_a, &hist_b, alt_hyp, exp_r_w, exp_p, exp_accept_hyp);
        }

        {
            let alt_hyp = AltHyp::Gt;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.1103;
            check_wilcoxon(&hist_a, &hist_b, alt_hyp, exp_r_w, exp_p, exp_accept_hyp);
        }
    }

    #[test]
    fn test_contrived_data() {
        let (hist_a, hist_b) = data_hists(contrived_data(), 1, 300, 3);

        let exp_r_w = 1442.5;

        {
            let alt_hyp = AltHyp::Lt;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.6675;
            check_wilcoxon(&hist_a, &hist_b, alt_hyp, exp_r_w, exp_p, exp_accept_hyp);
        }

        {
            let alt_hyp = AltHyp::Ne;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.6649;
            check_wilcoxon(&hist_a, &hist_b, alt_hyp, exp_r_w, exp_p, exp_accept_hyp);
        }

        {
            let alt_hyp = AltHyp::Gt;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.3325;
            check_wilcoxon(&hist_a, &hist_b, alt_hyp, exp_r_w, exp_p, exp_accept_hyp);
        }
    }

    #[test]
    fn test_shifted_contrived_data() {
        let (hist_a, hist_b) = data_hists(shifted_contrived_data(), 1, 300, 3);

        let exp_r_w = 840.;

        {
            let alt_hyp = AltHyp::Lt;
            let exp_accept_hyp = Hyp::Alt(AltHyp::Lt);
            let exp_p = 0.0002987;
            check_wilcoxon(&hist_a, &hist_b, alt_hyp, exp_r_w, exp_p, exp_accept_hyp);
        }

        {
            let alt_hyp = AltHyp::Ne;
            let exp_accept_hyp = Hyp::Alt(AltHyp::Ne);
            let exp_p = 0.0005974;
            check_wilcoxon(&hist_a, &hist_b, alt_hyp, exp_r_w, exp_p, exp_accept_hyp);
        }

        {
            let alt_hyp = AltHyp::Gt;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.9997;
            check_wilcoxon(&hist_a, &hist_b, alt_hyp, exp_r_w, exp_p, exp_accept_hyp);
        }
    }
}

#[cfg(test)]
#[cfg(feature = "_hypors")]
#[allow(clippy::unwrap_used)]
mod test_with_hypors {
    use super::*;
    use crate::{core::AltHyp, dev_utils::ApproxEq};
    use hdrhistogram::Histogram;
    use hypors::{common::TailType, mann_whitney::u_test};
    use polars::prelude::*;

    const ALPHA: f64 = 0.05;

    fn process_samples(sample_a: Vec<u64>, sample_b: Vec<u64>, hist_max: u64, hist_sigfig: u8) {
        let mut hist_a = Histogram::new_with_max(hist_max, hist_sigfig).unwrap();
        let mut hist_b = Histogram::new_from(&hist_a);

        for v in &sample_a {
            hist_a.record(*v).unwrap();
        }

        for v in &sample_b {
            hist_b.record(*v).unwrap();
        }

        let (ranked_items, _) = wilcoxon_ranked_items_ties_sum_prod(&hist_a, &hist_b);
        println!("{ranked_items:?}");

        let rank_sum_b = wilcoxon_rank_sum_w(&hist_a, &hist_b);
        println!("rank_sum_b={rank_sum_b}");

        let n_a = sample_a.len() as f64;
        let n_b = sample_b.len() as f64;
        let rank_sum_a = (1. + n_a + n_b) * (n_a + n_b) / 2. - rank_sum_b;
        println!("rank_sum_a={rank_sum_a}");

        let wilcoxon_rank_sum_a_lt_b_p = wilcoxon_rank_sum_p(&hist_a, &hist_b, AltHyp::Lt);
        println!("wilcoxon_rank_sum_a_lt_b_p={wilcoxon_rank_sum_a_lt_b_p}");
        let wilcoxon_rank_sum_a_lt_b_p_no_ties_adjust: f64 =
            wilcoxon_rank_sum_p_no_ties_adjust(&hist_a, &hist_b, AltHyp::Lt);
        println!(
            "wilcoxon_rank_sum_a_lt_b_p_no_ties_adjust={wilcoxon_rank_sum_a_lt_b_p_no_ties_adjust}"
        );
        let wilcoxon_rank_sum_a_gt_b_p = wilcoxon_rank_sum_p(&hist_a, &hist_b, AltHyp::Gt);
        println!("wilcoxon_rank_sum_a_gt_b_p={wilcoxon_rank_sum_a_gt_b_p}");
        let wilcoxon_rank_sum_a_ne_b_p: f64 = wilcoxon_rank_sum_p(&hist_a, &hist_b, AltHyp::Ne);
        println!("wilcoxon_rank_sum_a_ne_b_p={wilcoxon_rank_sum_a_ne_b_p}");
        let wilcoxon_rank_sum_a_ne_b_p_no_ties_adjust: f64 =
            wilcoxon_rank_sum_p_no_ties_adjust(&hist_a, &hist_b, AltHyp::Ne);
        println!(
            "wilcoxon_rank_sum_a_ne_b_p_no_ties_adjust={wilcoxon_rank_sum_a_ne_b_p_no_ties_adjust}"
        );

        let mann_whitney_a_lt_b_u = mann_whitney_u_b(&hist_a, &hist_b);
        println!("mann_whitney_a_lt_b_u={mann_whitney_a_lt_b_u}");
        let mann_whitney_a_gt_b_u = mann_whitney_u_a(&hist_a, &hist_b);
        println!("mann_whitney_a_gt_b_u={mann_whitney_a_gt_b_u}");
        let mann_whitney_a_ne_b_u = mann_whitney_u(&hist_a, &hist_b);
        println!("mann_whitney_a_ne_b_u={mann_whitney_a_ne_b_u}");

        {
            let series_a = Series::new(
                "a".into(),
                sample_a.iter().map(|x| *x as f64).collect::<Vec<_>>(),
            );
            let series_b = Series::new(
                "b".into(),
                sample_b.iter().map(|x| *x as f64).collect::<Vec<_>>(),
            );

            {
                let result = u_test(&series_a, &series_b, ALPHA, TailType::Two);
                println!("result(two tail)={result:?}");

                let result = result.unwrap();

                println!("U Statistic: {}", result.test_statistic);
                println!("P-value: {}", result.p_value);
                println!("Reject Null: {}", result.reject_null);

                assert_eq!(
                    result.test_statistic, mann_whitney_a_ne_b_u,
                    "comparison of U statistics"
                );

                assert_eq!(
                    result.p_value.round_to(5),
                    wilcoxon_rank_sum_a_ne_b_p_no_ties_adjust.round_to(5),
                    "comparison of p values for non-equality (no ties adjustment)"
                );
                // Below fails because `hypors` does not compute ties adjustment.
                // assert_eq!(
                //     result.p_value.round_to_sig_decimals(5),
                //     wilcoxon_rank_sum_a_ne_b_p.round_to_sig_decimals(5),
                //     "comparison of p values for non-equality"
                // );
            }
        }
    }

    fn expand_sample(sample: &[u64], delta: u64, nrepeats: usize) -> Vec<u64> {
        let mut cumulative = Vec::with_capacity(sample.len() * nrepeats);
        let mut curr = Vec::from(sample);
        for _ in 0..nrepeats {
            let next = curr.iter().map(|x| x + delta).collect::<Vec<_>>();
            cumulative.append(&mut next.clone());
            curr = next;
        }
        cumulative
    }

    #[test]
    fn test_book_data() {
        println!(
            "***** Samples from Nonparametric Statistical Methods, 3rd Edition, Example 4.1. *****"
        );
        {
            let sample_a0 = vec![0.73, 0.80, 0.83, 1.04, 1.38, 1.45, 1.46, 1.64, 1.89, 1.91];
            let sample_b0 = vec![0.74, 0.88, 0.90, 1.15, 1.21];

            let sample_a = sample_a0
                .into_iter()
                .map(|x| (x * 100.) as u64)
                .collect::<Vec<_>>();
            let sample_b = sample_b0
                .into_iter()
                .map(|x| (x * 100.) as u64)
                .collect::<Vec<_>>();

            process_samples(sample_a, sample_b, 200, 3);
        }
    }

    #[test]
    fn test_contrived_data() {
        let sample_a0 = vec![85, 90, 78, 92, 88, 76, 95, 89, 91, 82];
        let sample_b0 = vec![70, 85, 80, 90, 75, 88, 92, 79, 86, 81, 92];

        println!("***** Original samples *****");
        {
            let sample_a = sample_a0.clone();
            let sample_b = sample_b0.clone();

            {
                let mut sorted_a = sample_a.clone();
                sorted_a.sort();

                let mut sorted_b = sample_b.clone();
                sorted_b.sort();

                let mut combined = sample_a.iter().chain(sample_b.iter()).collect::<Vec<_>>();
                combined.sort();

                let exp_ranks_b = [1., 2., 5., 6., 7., 9.5, 11., 12.5, 15.5, 19., 19.];
                let exp_rank_sum_b = exp_ranks_b.iter().sum::<f64>();

                println!("sorted_a={sorted_a:?}");
                println!("sorted_b={sorted_b:?}");
                println!("combined={combined:?}");
                println!("exp_ranks_b={exp_ranks_b:?}");
                println!("exp_rank_sum_b={exp_rank_sum_b:?}");
            }
            process_samples(sample_a, sample_b, 300, 3);
        }

        println!();
        println!("***** Magnified samples *****");
        {
            let delta = 30;
            let nrepeats = 5;
            let sample_a = expand_sample(&sample_a0, delta, nrepeats);
            let sample_b = expand_sample(&sample_b0, delta, nrepeats);
            process_samples(sample_a, sample_b, 300, 3);
        }

        println!();
        println!("***** sample_a < sample_b *****");
        {
            let delta = 30;
            let nrepeats = 5;
            let sample_a = expand_sample(&sample_a0, delta, nrepeats);
            let sample_b = expand_sample(&sample_b0, delta, nrepeats);

            let shift = 35;
            let sample_b = sample_b.iter().map(|x| x + shift).collect::<Vec<_>>();

            process_samples(sample_a, sample_b, 300, 3);
        }
    }
}
