use basic_stats::core::SampleMoments;

fn main() {
    let x = [14., 15., 15., 15., 16., 18., 22., 23., 24., 25., 25.];
    let moments_s = SampleMoments::from_slice(&x);

    let n = 11;
    let sum = 212.;
    let sum2 = 4290.;
    let min = 14.;
    let max = 25.;
    let moments = SampleMoments::new_with_min_max(n, sum, sum2, min, max);

    assert_eq!(moments, moments_s);

    println!("n={}", moments.n());
    println!("nf={}", moments.nf());
    println!("sum={}", moments.sum());
    println!("sum2={}", moments.sum2());
    println!("min={}", moments.min());
    println!("max={}", moments.max());
    println!("mean={:?}", moments.mean());
    println!("sum2_deviations={:?}", moments.sum2_deviations());
    println!("var={:?}", moments.var());
    println!("stdev={:?}", moments.stdev());
}
