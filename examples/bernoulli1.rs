use basic_stats::{binomial::binomial_cp_alt_hyp_ci, core::AltHyp};

const ALPHA: f64 = 0.05;
const ALT_HYP: AltHyp = AltHyp::Ne;

fn main() {
    {
        let n = 1;
        let n_s = 0;
        let ci = binomial_cp_alt_hyp_ci(n, n_s, ALT_HYP, ALPHA);
        println!("n={n}, n_s={n_s}, ci={ci:?}");
    }

    {
        let n = 1;
        let n_s = 1;
        let ci = binomial_cp_alt_hyp_ci(n, n_s, ALT_HYP, ALPHA);
        println!("n={n}, n_s={n_s}, ci={ci:?}");
    }

    {
        let n = 100_000;
        let n_s = 0;
        let ci = binomial_cp_alt_hyp_ci(n, n_s, ALT_HYP, ALPHA);
        println!("n={n}, n_s={n_s}, ci={ci:?}");
    }

    {
        let n = 100_000;
        let n_s = 1;
        let ci = binomial_cp_alt_hyp_ci(n, n_s, ALT_HYP, ALPHA);
        println!("n={n}, n_s={n_s}, ci={ci:?}");
    }

    {
        let n = 100_000;
        let n_s = 99_999;
        let ci = binomial_cp_alt_hyp_ci(n, n_s, ALT_HYP, ALPHA);
        println!("n={n}, n_s={n_s}, ci={ci:?}");
    }

    {
        let n = 100_000;
        let n_s = 100_000;
        let ci = binomial_cp_alt_hyp_ci(n, n_s, ALT_HYP, ALPHA);
        println!("n={n}, n_s={n_s}, ci={ci:?}");
    }
}
