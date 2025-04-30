use basic_stats::iter::iter_with_counts;

fn main() {
    let dat = [1., 3., 9., 9., 10., 10., 10., 10., 20.];
    let dat_c = iter_with_counts(dat.into_iter()).collect::<Vec<_>>();
    let exp_dat_c = vec![(1., 1), (3., 1), (9., 2), (10., 4), (20., 1)];
    assert_eq!(exp_dat_c, dat_c);
    println!("success");
}
