use basic_stats::core::{deterministic_sample, deterministic_uniform_sample};
use statrs::distribution::{ContinuousCDF, Normal};

fn main() {
    let iter_u = deterministic_uniform_sample(10);

    let normal = Normal::new(0., 1.).unwrap();
    let mut iter_n = deterministic_sample(|x| normal.inverse_cdf(x), 10);

    for (count, item_u) in iter_u.enumerate() {
        let item_n = iter_n.next().unwrap();
        let i = count + 1;
        println!("i={i}\tuniform={item_u}\tnormal={item_n}");
    }
}
