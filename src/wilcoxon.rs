//! Statistics related to the Wilcoxon rank sum two-sample test, also known as the Mann-Whitney U test.
//!
//! This module is included by default. However, if `default-features = false` is specified in the dependency
//! declaration for this library, then inclusion of this module is gated by feature "**wilcoxon**".
//!
//! # Example
//!
//! ```
#![doc = include_str!("../examples/wilcoxon_ranksum.rs")]
//! ```

use crate::{
    core::{
        AltHyp, HypTestResult, StatsError, StatsResult, check_alpha_in_open_0_1, iter_with_counts,
    },
    normal::z_to_p,
};

/// Encapsulates the Wilcoxon rank sum computations on two data samples.
/// This struct's methods implement the Wilcoxon rank sum test and related statistics.
#[derive(Debug)]
pub struct RankSum {
    n_x: u64,
    n_y: u64,
    w: f64,
    ties_sum_prod: u64,
}

impl RankSum {
    /// Instantiates `Self` from two samples in the form of two iterators of pairs. Each item returned
    /// by the iterators is a pair whose first component is a data value and the second component is the
    /// number of occurrences of the value in the sample.
    ///
    /// # Errors
    ///
    /// Returns an error if an iterator does not yield data values in strictly increasing order.
    pub fn from_iters_with_counts(
        mut itc_x: impl Iterator<Item = (f64, u64)>,
        mut itc_y: impl Iterator<Item = (f64, u64)>,
    ) -> Result<RankSum, StatsError> {
        let mut n_x = 0;
        let mut n_y = 0;
        let mut ties_sum_prod = 0;
        let mut prev_item_x: Option<(f64, u64)> = None;
        let mut prev_item_y: Option<(f64, u64)> = None;

        fn enforce_order(
            prev_item: &mut Option<(f64, u64)>,
            curr_item: &Option<(f64, u64)>,
        ) -> Result<(), StatsError> {
            match (&prev_item, curr_item) {
                (Some((prev, _)), Some((curr, _))) if *prev < *curr => *prev_item = *curr_item,
                (Some(_), Some(_)) => {
                    return Err(StatsError(
                        "invalid iterator argument: items not ordered properly",
                    ));
                }
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
        ) -> Result<(f64, f64), StatsError> {
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

        let rank_sum_y: f64 = {
            let mut rank_sum_x = 0.;
            let mut rank_sum_y = 0.;
            let (mut item_x_opt, mut item_y_opt) = (itc_x.next(), itc_y.next());
            enforce_order(&mut prev_item_x, &item_x_opt)?;
            enforce_order(&mut prev_item_y, &item_y_opt)?;
            let mut prev_rank = 0.;

            loop {
                match (&mut item_x_opt, &mut item_y_opt) {
                    (Some(item_x), None) => {
                        let count = item_x.1;
                        let (rank_sum, new_prev_rank) = rank_item(
                            count,
                            0,
                            prev_rank,
                            &mut n_x,
                            &mut itc_x,
                            &mut item_x_opt,
                            &mut prev_item_x,
                            &mut ties_sum_prod,
                        )?;
                        rank_sum_x += rank_sum;
                        prev_rank = new_prev_rank;
                    }

                    (Some(item_x), Some(item_y)) if item_x.0 < item_y.0 => {
                        let (rank_sum, new_prev_rank) = rank_item(
                            item_x.1,
                            0,
                            prev_rank,
                            &mut n_x,
                            &mut itc_x,
                            &mut item_x_opt,
                            &mut prev_item_x,
                            &mut ties_sum_prod,
                        )?;
                        rank_sum_x += rank_sum;
                        prev_rank = new_prev_rank;
                    }

                    (None, Some(item_y)) => {
                        let (rank_sum, new_prev_rank) = rank_item(
                            item_y.1,
                            0,
                            prev_rank,
                            &mut n_y,
                            &mut itc_y,
                            &mut item_y_opt,
                            &mut prev_item_y,
                            &mut ties_sum_prod,
                        )?;
                        rank_sum_y += rank_sum;
                        prev_rank = new_prev_rank;
                    }

                    (Some(item_x), Some(item_y)) if item_x.0 > item_y.0 => {
                        let (rank_sum, new_prev_rank) = rank_item(
                            item_y.1,
                            0,
                            prev_rank,
                            &mut n_y,
                            &mut itc_y,
                            &mut item_y_opt,
                            &mut prev_item_y,
                            &mut ties_sum_prod,
                        )?;
                        rank_sum_y += rank_sum;
                        prev_rank = new_prev_rank;
                    }

                    // if item_x.0 == item_y.0
                    (Some(item_x), Some(item_y)) => {
                        let count_x = item_x.1;
                        let count_y = item_y.1;

                        let (rank_sum, _) = rank_item(
                            count_x,
                            count_y,
                            prev_rank,
                            &mut n_x,
                            &mut itc_x,
                            &mut item_x_opt,
                            &mut prev_item_x,
                            &mut ties_sum_prod,
                        )?;
                        rank_sum_x += rank_sum;

                        let (rank_sum, new_prev_rank) = rank_item(
                            count_y,
                            count_x,
                            prev_rank,
                            &mut n_y,
                            &mut itc_y,
                            &mut item_y_opt,
                            &mut prev_item_y,
                            &mut ties_sum_prod,
                        )?;
                        rank_sum_y += rank_sum;
                        prev_rank = new_prev_rank;
                    }

                    (None, None) => break,
                }
            }

            // Check rank-sum calculation.
            {
                let nf_x = n_x as f64;
                let nf_y = n_y as f64;
                let expected_rank_sum_x = (1. + nf_x + nf_y) * (nf_x + nf_y) / 2. - rank_sum_y;
                debug_assert_eq!(expected_rank_sum_x, rank_sum_x, "rank_sum_x check");
            }

            rank_sum_y
        };

        Ok(RankSum {
            n_x,
            n_y,
            w: rank_sum_y,
            ties_sum_prod,
        })
    }

    /// Instantiates `Self` from two samples in the form of two iterators. Each item returned
    /// by the iterators is a data value.
    ///
    /// # Errors
    ///
    /// Returns an error if an iterator does not yield data values in non-decreasing order.
    pub fn from_iters(
        it_x: impl Iterator<Item = f64>,
        it_y: impl Iterator<Item = f64>,
    ) -> Result<RankSum, StatsError> {
        let itc_x = iter_with_counts(it_x);
        let itc_y = iter_with_counts(it_y);
        Self::from_iters_with_counts(itc_x, itc_y)
    }

    /// Instantiates `Self` from two samples in the form of two slices.
    ///
    /// # Errors
    ///
    /// Returns an error if a slice is not sorted in non-decreasing order.
    pub fn from_slices(x: &[f64], y: &[f64]) -> Result<RankSum, StatsError> {
        let itc_x = x.iter().cloned();
        let itc_y = y.iter().cloned();
        Self::from_iters(itc_x, itc_y)
    }

    /// Size of first sample (X).
    pub fn n_x(&self) -> u64 {
        self.n_x
    }

    /// Size of second sample (Y).
    pub fn n_y(&self) -> u64 {
        self.n_y
    }

    /// Wilcoxon rank sum `W` statistic.
    ///
    /// As defined in the book Nonparametric Statistical Methods, 3rd Edition,
    /// by Myles Hollander, Douglas A. Wolfe, Eric Chicken, Section 4.1.
    pub fn w(&self) -> f64 {
        self.w
    }

    /// The `W` value computed by `R`'s `wilcox.test` function, which is the Mann-Whitney U for the
    /// second sample (Y).
    ///
    /// See explanation in the book Nonparametric Statistical Methods, 3rd Edition,
    /// by Myles Hollander, Douglas A. Wolfe, Eric Chicken, Example 4.1.
    pub fn r_w(&self) -> f64 {
        self.mann_whitney_u_y()
    }

    /// Mann-Whitney U statistic for the first sample (X).
    ///
    /// This is the number of X before Y predecessors.
    pub fn mann_whitney_u_x(&self) -> f64 {
        let n_y = self.n_y as f64;
        self.w - n_y * (n_y + 1.) / 2.
    }

    /// Mann-Whitney U statistic for the second sample (Y).
    ///
    /// This is the number of Y before X predecessors.
    pub fn mann_whitney_u_y(&self) -> f64 {
        let n_x = self.n_x as f64;
        let n_y = self.n_y as f64;
        (n_x * n_y) - self.mann_whitney_u_x()
    }

    /// Mann-Whitney U statistic.
    ///
    /// It is the min of [`mann_whitney_u_x`](Self::mann_whitney_u_x) and [`mann_whitney_u_y`](Self::mann_whitney_u_y).
    pub fn mann_whitney_u(&self) -> f64 {
        self.mann_whitney_u_x().min(self.mann_whitney_u_y())
    }

    /// z-value for the large sample normal approximation, without continuity correction.
    ///
    /// # Errors
    ///
    /// Returns an error in any of the following conditions:
    /// - `self.n_x == 0` or `self.n_y == 0`.
    /// - There are too many rank ties between the two samples (causing an intermediate `NaN` value).
    ///   This is hard to quantify a priori. For example,
    ///   `x = [2., 2., 2., 2.]` and `y = [2., 2., 2., 3., 3.]` are OK
    ///   but `x = [2., 2., 2., 2., 2.]` and `y = [2., 2., 2., 3., 3.]` result in an error.
    pub fn z(&self) -> StatsResult<f64> {
        // Guard against division by 0 in `var0_w_ties_adjust`.
        if self.n_x == 0 || self.n_y == 0 {
            return Err(StatsError(
                "`self.n_x` and `self.n_y` must both be positive",
            ));
        }

        let n_x = self.n_x as f64;
        let n_y = self.n_y as f64;
        let w = self.w;
        let ties_sum_prod = self.ties_sum_prod as f64;
        let e0_w = n_y * (n_x + n_y + 1.) / 2.;
        let var0_w_base = n_x * n_y * (n_x + n_y + 1.) / 12.;
        let var0_w_ties_adjust = n_x * n_y * ties_sum_prod / (12. * (n_x + n_y) * (n_x + n_y - 1.));
        let var0_w = var0_w_base - var0_w_ties_adjust;
        if var0_w <= 0. {
            return Err(StatsError("too many rank ties"));
        }
        let w_star = (w - e0_w) / var0_w.sqrt();

        Ok(-w_star)
    }

    /// p-value for the large sample normal approximation, without continuity correction.
    ///
    /// Arguments:
    /// - `alt_hyp`: alternative hypothesis.
    ///
    /// # Errors
    ///
    /// Returns an error in any of the following conditions:
    /// - `self.n_x == 0` or `self.n_y == 0`.
    /// - There are too many rank ties between the two samples (causing an intermediate `NaN` value).
    ///   This is hard to quantify a priori. For example,
    ///   `x = [2., 2., 2., 2.]` and `y = [2., 2., 2., 3., 3.]` are OK
    ///   but `x = [2., 2., 2., 2., 2.]` and `y = [2., 2., 2., 3., 3.]` result in an error.
    pub fn z_p(&self, alt_hyp: AltHyp) -> StatsResult<f64> {
        let z = self.z()?;
        Ok(z_to_p(z, alt_hyp))
    }

    /// Wilcoxon rank sum test using large sample normal approximation, without continuity correction.
    ///
    /// Arguments:
    /// - `alt_hyp`: alternative hypothesis.
    /// - `alpha`: confidence level = `1 - alpha`.
    ///
    /// # Errors
    ///
    /// Returns an error in any of these conditions:
    /// - `self.n_x == 0` or `self.n_y == 0`.
    /// - `alpha` not in interval `(0, 1)`.
    /// - There are too many rank ties between the two samples (causing an intermediate `NaN` value).
    ///   This is hard to quantify a priori. For example,
    ///   `x = [2., 2., 2., 2.]` and `y = [2., 2., 2., 3., 3.]` are OK
    ///   but `x = [2., 2., 2., 2., 2.]` and `y = [2., 2., 2., 3., 3.]` result in an error.
    pub fn z_test(&self, alt_hyp: AltHyp, alpha: f64) -> StatsResult<HypTestResult> {
        check_alpha_in_open_0_1(alpha)?;
        let p = self.z_p(alt_hyp)?;
        Ok(HypTestResult::new(p, alpha, alt_hyp))
    }
}

#[cfg(test)]
#[cfg(feature = "_dev_utils")]
mod base_test {
    //! Tests other than `test_w` used R's wilcox.test function to generate expected results.
    //! https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/wilcox.test

    use super::*;
    use crate::{
        core::{AltHyp, Hyp},
        dev_utils::ApproxEq,
    };
    use std::error::Error;

    const ALPHA: f64 = 0.05;
    const EPSILON: f64 = 0.0005;

    fn sort_array(arr: &mut [f64]) {
        arr.sort_by(|a, b| a.partial_cmp(b).unwrap());
    }

    fn book_data() -> (Vec<f64>, Vec<f64>) {
        let sample_x = vec![0.73, 0.80, 0.83, 1.04, 1.38, 1.45, 1.46, 1.64, 1.89, 1.91];
        let sample_y = vec![0.74, 0.88, 0.90, 1.15, 1.21];
        (sample_x, sample_y)
    }

    fn contrived_data() -> (Vec<f64>, Vec<f64>) {
        let mut sample_x = vec![
            85., 90., 78., 92., 88., 76., 95., 89., 91., 82., 115., 120., 108., 122., 118., 106.,
            125., 119., 121., 112., 145., 150., 138., 152., 148., 136., 155., 149., 151., 142.,
            175., 180., 168., 182., 178., 166., 185., 179., 181., 172., 205., 210., 198., 212.,
            208., 196., 215., 209., 211., 202.,
        ];
        sort_array(&mut sample_x);

        let mut sample_y = vec![
            70., 85., 80., 90., 75., 88., 92., 79., 86., 81., 92., 100., 115., 110., 120., 105.,
            118., 122., 109., 116., 111., 122., 130., 145., 140., 150., 135., 148., 152., 139.,
            146., 141., 152., 160., 175., 170., 180., 165., 178., 182., 169., 176., 171., 182.,
            190., 205., 200., 210., 195., 208., 212., 199., 206., 201., 212.,
        ];
        sort_array(&mut sample_y);

        (sample_x, sample_y)
    }

    fn shifted_contrived_data() -> (Vec<f64>, Vec<f64>) {
        let (mut sample_x, sample_y) = contrived_data();
        let mut sample_y = sample_y.into_iter().map(|v| v + 35.).collect::<Vec<_>>();
        sort_array(&mut sample_x);
        sort_array(&mut sample_y);

        (sample_x, sample_y)
    }

    #[test]
    /// Based on https://learning.oreilly.com/library/view/nonparametric-statistical-methods/9781118553299/9781118553299c04.xhtml#c04_level1_2
    /// Nonparametric Statistical Methods, 3rd Edition, by Myles Hollander, Douglas A. Wolfe, Eric Chicken
    /// Example 4.1.
    fn test_w() -> Result<(), Box<dyn Error>> {
        let (dat_x, dat_y) = book_data();
        let rank_sum = RankSum::from_slices(&dat_x, &dat_y)?;

        let expected_w = 30.;
        let actual_w = rank_sum.w();
        assert_eq!(expected_w, actual_w, "w comparison");

        let expected_r_w = 35.;
        let actual_r_w = rank_sum.r_w();
        assert_eq!(expected_r_w, actual_r_w, "R w comparison");

        let expected_p_correct = 0.2544; // R: // R: wilcox.test(a, b)
        let expected_p = 0.2207; // R: wilcox.test(a, b, exact=FALSE, correct=FALSE)
        let actual_p = rank_sum.z_p(AltHyp::Ne).unwrap();
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
        let p = rank_sum.z_p(alt_hyp).unwrap();
        let res = rank_sum.z_test(alt_hyp, ALPHA).unwrap();

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

        let mann_whitney_u_x = rank_sum.mann_whitney_u_y();
        println!("alt_hyp={alt_hyp:?} -- mann_whitney_u_x={mann_whitney_u_x}");
        let mann_whitney_u_y = rank_sum.mann_whitney_u_x();
        println!("alt_hyp={alt_hyp:?} -- mann_whitney_u_y={mann_whitney_u_y}");
        let mann_whitney_u = rank_sum.mann_whitney_u();
        println!("alt_hyp={alt_hyp:?} -- mann_whitney_u={mann_whitney_u}");
    }

    #[test]
    /// Based on https://learning.oreilly.com/library/view/nonparametric-statistical-methods/9781118553299/9781118553299c04.xhtml#c04_level1_2
    /// Nonparametric Statistical Methods, 3rd Edition, by Myles Hollander, Douglas A. Wolfe, Eric Chicken
    /// Example 4.1.
    fn test_book_data() -> Result<(), Box<dyn Error>> {
        let (dat_x, dat_y) = book_data();
        let rank_sum = RankSum::from_slices(&dat_x, &dat_y)?;

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
        let (dat_x, dat_y) = contrived_data();
        let rank_sum = RankSum::from_iters(dat_x.into_iter(), dat_y.into_iter())?;

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
    fn test_contrived_data_u() -> Result<(), Box<dyn Error>> {
        let (dat_x, dat_y) = contrived_data();
        let rank_sum = RankSum::from_iters(dat_x.into_iter(), dat_y.into_iter())?;

        // From https://www.statskingdom.com/170median_mann_whitney.html:
        let exp_u_x = 1307.5;
        let exp_u_y = 1442.5;
        let exp_u = 1307.5;

        let u_x = rank_sum.mann_whitney_u_x();
        let u_y = rank_sum.mann_whitney_u_y();
        let u = rank_sum.mann_whitney_u();

        println!("exp_u_x={exp_u_x}, u_x={u_x}");
        println!("exp_u_x={exp_u_y}, u_x={u_y}");
        println!("exp_u_x={exp_u}, u_x={u}");

        assert!(exp_u_x.approx_eq(u_x, EPSILON));
        assert!(exp_u_y.approx_eq(u_y, EPSILON));
        assert!(exp_u.approx_eq(u, EPSILON));

        Ok(())
    }

    #[test]
    fn test_shifted_contrived_data() -> Result<(), Box<dyn Error>> {
        let (dat_x, dat_y) = shifted_contrived_data();
        let rank_sum = RankSum::from_iters(dat_x.into_iter(), dat_y.into_iter())?;

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
