//! Example of use of [`Aok`] trait.

use basic_stats::{
    aok::{AokBasicStats, AokFloat},
    core::{AltHyp, SampleMoments},
    normal::{welch_alt_hyp_ci, welch_p, welch_test},
};

fn main() {
    let x = [14., 15., 15., 15., 16., 18., 22., 23., 24., 25., 25.];
    let y = [
        10., 12., 14., 15., 18., 22., 24., 27., 31., 33., 34., 34., 34.,
    ];

    let moments_x = SampleMoments::from_slice(&x);
    let moments_y = SampleMoments::from_slice(&y);
    let alt_hyp = AltHyp::Gt;

    {
        println!("*** Ok scenario:");

        let alpha = 0.05;

        // Welch function calls below return Ok prior to invocation of noerr().

        let p = welch_p(&moments_x, &moments_y, alt_hyp).aok();
        println!("p={p}");
        let ci = welch_alt_hyp_ci(&moments_x, &moments_y, alt_hyp, alpha).aok();
        println!("ci={ci:?}");
        let test_res = welch_test(&moments_x, &moments_y, alt_hyp, alpha).aok();
        println!("test_res={test_res:?}");
    }

    {
        println!("*** Err scenario:");

        let alpha = 1.0;

        // Welch function calls below return Err prior to invocation of noerr().

        let p = welch_p(&moments_x, &SampleMoments::default(), alt_hyp).aok();
        println!("p={p}");
        let ci = welch_alt_hyp_ci(&moments_x, &moments_y, alt_hyp, alpha).aok();
        println!("ci={ci:?}");
        let test_res = welch_test(&moments_x, &moments_y, alt_hyp, alpha).aok();
        println!("test_res={test_res:?}");
    }
}
