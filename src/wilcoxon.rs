use super::{
    core::{AltHyp, HypTestResult},
    normal::z_to_p,
};
use crate::iter::{IterWithCounts, OrderingError};

/// The result of the Wilcoxon rank sum computations on two data samples.
/// This struct's methods implement the Wilcoxon rank sum test and related statistics.
pub struct RankSum {
    n_a: u64,
    n_b: u64,
    w: f64,
    ties_sum_prod: u64,
}

impl RankSum {
    /// Instantiates `Self` from two samples in the form of two iterators of pairs. Each item returned
    /// by the iterators is a pair whose first component is a data value and the second component is the
    /// number of occurrences of the value in the sample.
    ///
    /// # Errors
    /// - Returns an error if the iterators do not yield data values in strictly increasing order.
    pub fn from_iter_with_counts(
        mut itc_a: impl Iterator<Item = (f64, u64)>,
        mut itc_b: impl Iterator<Item = (f64, u64)>,
    ) -> Result<RankSum, OrderingError> {
        let mut n_a = 0;
        let mut n_b = 0;
        let mut ties_sum_prod = 0;
        let mut prev_item_a: Option<(f64, u64)> = None;
        let mut prev_item_b: Option<(f64, u64)> = None;

        fn enforce_order(
            prev_item: &mut Option<(f64, u64)>,
            curr_item: &Option<(f64, u64)>,
        ) -> Result<(), OrderingError> {
            match (&prev_item, curr_item) {
                (Some((prev, _)), Some((curr, _))) if *prev < *curr => *prev_item = *curr_item,
                (Some(_), Some(_)) => return Err(OrderingError),
                (None, Some(_)) => *prev_item = *curr_item,
                (_, None) => (),
            }
            Ok(())
        }

        #[allow(clippy::too_many_arguments)]
        fn rank_item(
            count_i: u64,
            count_other: u64,
            prev_rank: f64,
            n_i: &mut u64,
            itc_i: &mut impl Iterator<Item = (f64, u64)>,
            item_opt_i: &mut Option<(f64, u64)>,
            prev_item_i: &mut Option<(f64, u64)>,
            ties_sum_prod: &mut u64,
        ) -> Result<(f64, f64), OrderingError> {
            let count = count_i + count_other;
            let rank = prev_rank + (count as f64 + 1.) / 2.;
            let rank_sum = count_i as f64 * rank;
            let new_prev_rank = prev_rank + count as f64;
            *n_i += count_i;
            *item_opt_i = itc_i.next();
            enforce_order(prev_item_i, item_opt_i)?;
            *ties_sum_prod += (count - 1) * count * (count + 1);
            Ok((rank_sum, new_prev_rank))
        }

        let rank_sum_b: f64 = {
            let mut rank_sum_a = 0.;
            let mut rank_sum_b = 0.;
            let (mut item_a_opt, mut item_b_opt) = (itc_a.next(), itc_b.next());
            enforce_order(&mut prev_item_a, &item_a_opt)?;
            enforce_order(&mut prev_item_b, &item_b_opt)?;
            let mut prev_rank = 0.;

            loop {
                match (&mut item_a_opt, &mut item_b_opt) {
                    (Some(item_a), None) => {
                        let count = item_a.1;
                        let (rank_sum, new_prev_rank) = rank_item(
                            count,
                            0,
                            prev_rank,
                            &mut n_a,
                            &mut itc_a,
                            &mut item_a_opt,
                            &mut prev_item_a,
                            &mut ties_sum_prod,
                        )?;
                        rank_sum_a += rank_sum;
                        prev_rank = new_prev_rank;
                    }

                    (Some(item_a), Some(item_b)) if item_a.0 < item_b.0 => {
                        let (rank_sum, new_prev_rank) = rank_item(
                            item_a.1,
                            0,
                            prev_rank,
                            &mut n_a,
                            &mut itc_a,
                            &mut item_a_opt,
                            &mut prev_item_a,
                            &mut ties_sum_prod,
                        )?;
                        rank_sum_a += rank_sum;
                        prev_rank = new_prev_rank;
                    }

                    (None, Some(item_b)) => {
                        let (rank_sum, new_prev_rank) = rank_item(
                            item_b.1,
                            0,
                            prev_rank,
                            &mut n_b,
                            &mut itc_b,
                            &mut item_b_opt,
                            &mut prev_item_b,
                            &mut ties_sum_prod,
                        )?;
                        rank_sum_b += rank_sum;
                        prev_rank = new_prev_rank;
                    }

                    (Some(item_a), Some(item_b)) if item_a.0 > item_b.0 => {
                        let (rank_sum, new_prev_rank) = rank_item(
                            item_b.1,
                            0,
                            prev_rank,
                            &mut n_b,
                            &mut itc_b,
                            &mut item_b_opt,
                            &mut prev_item_b,
                            &mut ties_sum_prod,
                        )?;
                        rank_sum_b += rank_sum;
                        prev_rank = new_prev_rank;
                    }

                    // if item_a.0 == item_b.0
                    (Some(item_a), Some(item_b)) => {
                        let count_a = item_a.1;
                        let count_b = item_b.1;

                        let (rank_sum, _) = rank_item(
                            count_a,
                            count_b,
                            prev_rank,
                            &mut n_a,
                            &mut itc_a,
                            &mut item_a_opt,
                            &mut prev_item_a,
                            &mut ties_sum_prod,
                        )?;
                        rank_sum_a += rank_sum;

                        let (rank_sum, new_prev_rank) = rank_item(
                            count_b,
                            count_a,
                            prev_rank,
                            &mut n_b,
                            &mut itc_b,
                            &mut item_b_opt,
                            &mut prev_item_b,
                            &mut ties_sum_prod,
                        )?;
                        rank_sum_b += rank_sum;
                        prev_rank = new_prev_rank;
                    }

                    (None, None) => break,
                }
            }

            // Check rank-sum calculation.
            {
                let nf_a = n_a as f64;
                let nf_b = n_b as f64;
                let expected_rank_sum_a = (1. + nf_a + nf_b) * (nf_a + nf_b) / 2. - rank_sum_b;
                debug_assert_eq!(expected_rank_sum_a, rank_sum_a, "rank_sum_a check");
            }

            rank_sum_b
        };

        Ok(RankSum {
            n_a,
            n_b,
            w: rank_sum_b,
            ties_sum_prod,
        })
    }

    /// Instantiates `Self` from two samples in the form of two iterators. Each item returned
    /// by the iterators is a data value.
    ///
    /// # Errors
    /// - Returns an error if the iterators do not yield data values in strictly increasing order.
    pub fn from_iter(
        it_a: impl Iterator<Item = f64>,
        it_b: impl Iterator<Item = f64>,
    ) -> Result<RankSum, OrderingError> {
        let itc_a = IterWithCounts::new(it_a);
        let itc_b = IterWithCounts::new(it_b);
        Self::from_iter_with_counts(itc_a, itc_b)
    }

    pub fn n_a(&self) -> u64 {
        self.n_a
    }

    pub fn n_b(&self) -> u64 {
        self.n_b
    }

    pub fn w(&self) -> f64 {
        self.w
    }

    /// The `w` value computed by `R`'s `wilcox.test` function, which is the Mann-Whitney U for the
    /// first sample (`hist_a`).
    ///
    /// See explanation in the book Nonparametric Statistical Methods, 3rd Edition,
    /// by Myles Hollander, Douglas A. Wolfe, Eric Chicken, Example 4.1.
    pub fn r_w(&self) -> f64 {
        self.mann_whitney_u_a()
    }

    pub fn mann_whitney_u_b(&self) -> f64 {
        let n_b = self.n_b as f64;
        self.w - n_b * (n_b + 1.) / 2.
    }

    pub fn mann_whitney_u_a(&self) -> f64 {
        let n_a = self.n_a as f64;
        let n_b = self.n_b as f64;
        (n_a * n_b) - self.mann_whitney_u_b()
    }

    pub fn mann_whitney_u(&self) -> f64 {
        self.mann_whitney_u_b().min(self.mann_whitney_u_a())
    }

    pub fn z(&self) -> f64 {
        let n_a = self.n_a as f64;
        let n_b = self.n_b as f64;
        let w = self.w;
        let ties_sum_prod = self.ties_sum_prod as f64;
        let e0_w = n_b * (n_a + n_b + 1.) / 2.;
        let var0_w_base = n_a * n_b * (n_a + n_b + 1.) / 12.;
        let var0_w_ties_adjust = n_a * n_b / (12. * (n_a + n_b) * (n_a + n_b - 1.)) * ties_sum_prod;
        let var0_w = var0_w_base - var0_w_ties_adjust;
        let w_star = (w - e0_w) / var0_w.sqrt();

        -w_star
    }

    #[cfg(test)]
    fn z_no_ties_adjust(&self) -> f64 {
        let n_a = self.n_a as f64;
        let n_b = self.n_b as f64;
        let w = self.w;
        let e0_w = n_b * (n_a + n_b + 1.) / 2.;
        let var0_w_base = n_a * n_b * (n_a + n_b + 1.) / 12.;
        let var0_w_ties_adjust = 0.;
        let var0_w = var0_w_base - var0_w_ties_adjust;
        let w_star = (w - e0_w) / var0_w.sqrt();

        -w_star
    }

    pub fn p(&self, alt_hyp: AltHyp) -> f64 {
        let z = self.z();
        z_to_p(z, alt_hyp)
    }

    #[cfg(test)]
    fn p_no_ties_adjust(&self, alt_hyp: AltHyp) -> f64 {
        let z = self.z_no_ties_adjust();
        z_to_p(z, alt_hyp)
    }

    pub fn test(&self, alt_hyp: AltHyp, alpha: f64) -> HypTestResult {
        let p = self.p(alt_hyp);
        HypTestResult::new(p, alpha, alt_hyp)
    }

    #[cfg(test)]
    pub fn test_no_ties_adjust(&self, alt_hyp: AltHyp, alpha: f64) -> HypTestResult {
        let p = self.p_no_ties_adjust(alt_hyp);
        HypTestResult::new(p, alpha, alt_hyp)
    }
}

#[cfg(test)]
fn sort_array(arr: &mut [f64]) {
    arr.sort_by(|a, b| a.partial_cmp(b).unwrap());
}

#[cfg(test)]
mod base_test {
    //! Tests other than `test_w` used R's wilcox.test function to generate expected results.
    //! https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/wilcox.test

    use std::error::Error;

    use super::*;
    use crate::{
        core::{AltHyp, Hyp},
        dev_utils::ApproxEq,
    };

    const ALPHA: f64 = 0.05;
    const EPSILON: f64 = 0.0005;

    fn book_data() -> (Vec<f64>, Vec<f64>) {
        let sample_a = vec![0.73, 0.80, 0.83, 1.04, 1.38, 1.45, 1.46, 1.64, 1.89, 1.91];
        let sample_b = vec![0.74, 0.88, 0.90, 1.15, 1.21];
        (sample_a, sample_b)
    }

    fn contrived_data() -> (Vec<f64>, Vec<f64>) {
        let mut sample_a = vec![
            85., 90., 78., 92., 88., 76., 95., 89., 91., 82., 115., 120., 108., 122., 118., 106.,
            125., 119., 121., 112., 145., 150., 138., 152., 148., 136., 155., 149., 151., 142.,
            175., 180., 168., 182., 178., 166., 185., 179., 181., 172., 205., 210., 198., 212.,
            208., 196., 215., 209., 211., 202.,
        ];
        sort_array(&mut sample_a);

        let mut sample_b = vec![
            70., 85., 80., 90., 75., 88., 92., 79., 86., 81., 92., 100., 115., 110., 120., 105.,
            118., 122., 109., 116., 111., 122., 130., 145., 140., 150., 135., 148., 152., 139.,
            146., 141., 152., 160., 175., 170., 180., 165., 178., 182., 169., 176., 171., 182.,
            190., 205., 200., 210., 195., 208., 212., 199., 206., 201., 212.,
        ];
        sort_array(&mut sample_b);

        (sample_a, sample_b)
    }

    fn shifted_contrived_data() -> (Vec<f64>, Vec<f64>) {
        let (mut sample_a, sample_b) = contrived_data();
        let mut sample_b = sample_b.into_iter().map(|v| v + 35.).collect::<Vec<_>>();
        sort_array(&mut sample_a);
        sort_array(&mut sample_b);

        (sample_a, sample_b)
    }

    #[test]
    /// Based on https://learning.oreilly.com/library/view/nonparametric-statistical-methods/9781118553299/9781118553299c04.xhtml#c04_level1_2
    /// Nonparametric Statistical Methods, 3rd Edition, by Myles Hollander, Douglas A. Wolfe, Eric Chicken
    /// Example 4.1.
    fn test_w() -> Result<(), Box<dyn Error>> {
        let (dat_a, dat_b) = book_data();

        // The iterator must yield f64 values, so use either `.iter().cloned()` or `.into_iter()`.
        // The latter is more efficient if the data set can be consumed.
        let rank_sum = RankSum::from_iter(dat_a.iter().cloned(), dat_b.into_iter())?;

        let expected_w = 30.;
        let actual_w = rank_sum.w();
        assert_eq!(expected_w, actual_w, "w comparison");

        let expected_r_w = 35.;
        let actual_r_w = rank_sum.r_w();
        assert_eq!(expected_r_w, actual_r_w, "R w comparison");

        let expected_p_correct = 0.2544; // R: // R: wilcox.test(a, b)
        let expected_p = 0.2207; // R: wilcox.test(a, b, exact=FALSE, correct=FALSE)
        let actual_p = rank_sum.p(AltHyp::Ne);
        println!(
            "expected_p_correct={expected_p_correct}, expected_p={expected_p}, actual_p={actual_p}"
        );
        assert!(expected_p.approx_eq(actual_p, EPSILON), "p comparison");
        Ok(())
    }

    fn check_wilcoxon(
        rank_sum: &RankSum,
        alt_hyp: AltHyp,
        exp_r_w: f64,
        exp_p: f64,
        exp_accept_hyp: Hyp,
    ) {
        let w = rank_sum.w();
        let r_w = rank_sum.r_w();
        let p = rank_sum.p(alt_hyp);
        let res = rank_sum.test(alt_hyp, ALPHA);

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

        let mann_whitney_u_a = rank_sum.mann_whitney_u_a();
        println!("alt_hyp={alt_hyp:?} -- mann_whitney_u_a={mann_whitney_u_a}");
        let mann_whitney_u_b = rank_sum.mann_whitney_u_b();
        println!("alt_hyp={alt_hyp:?} -- mann_whitney_u_b={mann_whitney_u_b}");
        let mann_whitney_u = rank_sum.mann_whitney_u();
        println!("alt_hyp={alt_hyp:?} -- mann_whitney_u={mann_whitney_u}");
    }

    #[test]
    /// Based on https://learning.oreilly.com/library/view/nonparametric-statistical-methods/9781118553299/9781118553299c04.xhtml#c04_level1_2
    /// Nonparametric Statistical Methods, 3rd Edition, by Myles Hollander, Douglas A. Wolfe, Eric Chicken
    /// Example 4.1.
    fn test_book_data() -> Result<(), Box<dyn Error>> {
        let (dat_a, dat_b) = book_data();
        let rank_sum = RankSum::from_iter(dat_a.into_iter(), dat_b.into_iter())?;

        let exp_r_w = 35.;

        {
            let alt_hyp = AltHyp::Lt;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.8897;
            check_wilcoxon(&rank_sum, alt_hyp, exp_r_w, exp_p, exp_accept_hyp);
        }

        {
            let alt_hyp = AltHyp::Ne;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.2207;
            check_wilcoxon(&rank_sum, alt_hyp, exp_r_w, exp_p, exp_accept_hyp);
        }

        {
            let alt_hyp = AltHyp::Gt;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.1103;
            check_wilcoxon(&rank_sum, alt_hyp, exp_r_w, exp_p, exp_accept_hyp);
        }

        Ok(())
    }

    #[test]
    fn test_contrived_data() -> Result<(), Box<dyn Error>> {
        let (dat_a, dat_b) = contrived_data();
        let rank_sum = RankSum::from_iter(dat_a.into_iter(), dat_b.into_iter())?;

        let exp_r_w = 1442.5;

        {
            let alt_hyp = AltHyp::Lt;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.6675;
            check_wilcoxon(&rank_sum, alt_hyp, exp_r_w, exp_p, exp_accept_hyp);
        }

        {
            let alt_hyp = AltHyp::Ne;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.6649;
            check_wilcoxon(&rank_sum, alt_hyp, exp_r_w, exp_p, exp_accept_hyp);
        }

        {
            let alt_hyp = AltHyp::Gt;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.3325;
            check_wilcoxon(&rank_sum, alt_hyp, exp_r_w, exp_p, exp_accept_hyp);
        }

        Ok(())
    }

    #[test]
    fn test_shifted_contrived_data() -> Result<(), Box<dyn Error>> {
        let (dat_a, dat_b) = shifted_contrived_data();
        let rank_sum = RankSum::from_iter(dat_a.into_iter(), dat_b.into_iter())?;

        let exp_r_w = 840.;

        {
            let alt_hyp = AltHyp::Lt;
            let exp_accept_hyp = Hyp::Alt(AltHyp::Lt);
            let exp_p = 0.0002987;
            check_wilcoxon(&rank_sum, alt_hyp, exp_r_w, exp_p, exp_accept_hyp);
        }

        {
            let alt_hyp = AltHyp::Ne;
            let exp_accept_hyp = Hyp::Alt(AltHyp::Ne);
            let exp_p = 0.0005974;
            check_wilcoxon(&rank_sum, alt_hyp, exp_r_w, exp_p, exp_accept_hyp);
        }

        {
            let alt_hyp = AltHyp::Gt;
            let exp_accept_hyp = Hyp::Null;
            let exp_p = 0.9997;
            check_wilcoxon(&rank_sum, alt_hyp, exp_r_w, exp_p, exp_accept_hyp);
        }

        Ok(())
    }
}

#[cfg(test)]
#[cfg(feature = "_hypors")]
#[allow(clippy::unwrap_used)]
mod test_with_hypors {
    use std::error::Error;

    use super::*;
    use crate::{core::AltHyp, dev_utils::ApproxEq};
    use hypors::{common::TailType, mann_whitney::u_test};
    use polars::prelude::*;

    const ALPHA: f64 = 0.05;

    fn process_samples(sample_a: Vec<f64>, sample_b: Vec<f64>) -> Result<(), Box<dyn Error>> {
        let rank_sum = RankSum::from_iter(sample_a.iter().cloned(), sample_b.iter().cloned())?;

        let rank_sum_b = rank_sum.w();
        println!("rank_sum_b={rank_sum_b}");

        let n_a = rank_sum.n_a() as f64;
        let n_b = rank_sum.n_b() as f64;
        let rank_sum_a = (1. + n_a + n_b) * (n_a + n_b) / 2. - rank_sum_b;
        println!("rank_sum_a={rank_sum_a}");

        let wilcoxon_rank_sum_a_lt_b_p = rank_sum.p(AltHyp::Lt);
        println!("wilcoxon_rank_sum_a_lt_b_p={wilcoxon_rank_sum_a_lt_b_p}");
        let wilcoxon_rank_sum_a_lt_b_p_no_ties_adjust: f64 = rank_sum.p_no_ties_adjust(AltHyp::Lt);
        println!(
            "wilcoxon_rank_sum_a_lt_b_p_no_ties_adjust={wilcoxon_rank_sum_a_lt_b_p_no_ties_adjust}"
        );
        let wilcoxon_rank_sum_a_gt_b_p = rank_sum.p(AltHyp::Gt);
        println!("wilcoxon_rank_sum_a_gt_b_p={wilcoxon_rank_sum_a_gt_b_p}");
        let wilcoxon_rank_sum_a_ne_b_p: f64 = rank_sum.p(AltHyp::Ne);
        println!("wilcoxon_rank_sum_a_ne_b_p={wilcoxon_rank_sum_a_ne_b_p}");
        let wilcoxon_rank_sum_a_ne_b_p_no_ties_adjust: f64 = rank_sum.p_no_ties_adjust(AltHyp::Ne);
        println!(
            "wilcoxon_rank_sum_a_ne_b_p_no_ties_adjust={wilcoxon_rank_sum_a_ne_b_p_no_ties_adjust}"
        );

        let mann_whitney_a_lt_b_u = rank_sum.mann_whitney_u_b();
        println!("mann_whitney_a_lt_b_u={mann_whitney_a_lt_b_u}");
        let mann_whitney_a_gt_b_u = rank_sum.mann_whitney_u_a();
        println!("mann_whitney_a_gt_b_u={mann_whitney_a_gt_b_u}");
        let mann_whitney_a_ne_b_u = rank_sum.mann_whitney_u();
        println!("mann_whitney_a_ne_b_u={mann_whitney_a_ne_b_u}");

        {
            let series_a = Series::new("a".into(), sample_a);
            let series_b = Series::new("b".into(), sample_b);

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

        Ok(())
    }

    fn expand_sample(sample: &[f64], delta: f64, nrepeats: usize) -> Vec<f64> {
        let mut cumulative = Vec::with_capacity(sample.len() * nrepeats);
        let mut curr = Vec::from(sample);
        for _ in 0..nrepeats {
            let next = curr.iter().map(|x| x + delta).collect::<Vec<_>>();
            cumulative.append(&mut next.clone());
            curr = next;
        }
        sort_array(&mut cumulative);
        cumulative
    }

    #[test]
    fn test_book_data() -> Result<(), Box<dyn Error>> {
        println!(
            "***** Samples from Nonparametric Statistical Methods, 3rd Edition, Example 4.1. *****"
        );
        {
            let sample_a = vec![0.73, 0.80, 0.83, 1.04, 1.38, 1.45, 1.46, 1.64, 1.89, 1.91];
            let sample_b = vec![0.74, 0.88, 0.90, 1.15, 1.21];

            process_samples(sample_a, sample_b)
        }
    }

    #[test]
    fn test_contrived_data() {
        let sample_a0 = vec![85., 90., 78., 92., 88., 76., 95., 89., 91., 82.];
        let sample_b0 = vec![70., 85., 80., 90., 75., 88., 92., 79., 86., 81., 92.];

        println!("***** Original samples *****");
        {
            let mut sorted_a = sample_a0.clone();
            sort_array(&mut sorted_a);

            let mut sorted_b = sample_b0.clone();
            sort_array(&mut sorted_b);

            {
                let mut combined = sorted_a
                    .iter()
                    .cloned()
                    .chain(sorted_b.iter().cloned())
                    .collect::<Vec<_>>();
                sort_array(&mut combined);

                let exp_ranks_b = [1., 2., 5., 6., 7., 9.5, 11., 12.5, 15.5, 19., 19.];
                let exp_rank_sum_b = exp_ranks_b.iter().sum::<f64>();

                println!("sorted_a={sorted_a:?}");
                println!("sorted_b={sorted_b:?}");
                println!("combined={combined:?}");
                println!("exp_ranks_b={exp_ranks_b:?}");
                println!("exp_rank_sum_b={exp_rank_sum_b:?}");
            }
            process_samples(sorted_a, sorted_b).unwrap();
        }

        println!();
        println!("***** Magnified samples *****");
        {
            let delta = 30.;
            let nrepeats = 5;
            let sample_a = expand_sample(&sample_a0, delta, nrepeats);
            let sample_b = expand_sample(&sample_b0, delta, nrepeats);
            process_samples(sample_a, sample_b).unwrap();
        }

        println!();
        println!("***** sample_a < sample_b *****");
        {
            let delta = 30.;
            let nrepeats = 5;
            let sample_a = expand_sample(&sample_a0, delta, nrepeats);
            let sample_b = expand_sample(&sample_b0, delta, nrepeats);

            let shift = 35.;
            let sample_b = sample_b.iter().map(|x| x + shift).collect::<Vec<_>>();

            process_samples(sample_a, sample_b).unwrap()
        }
    }
}
