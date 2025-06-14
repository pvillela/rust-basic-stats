use basic_stats::{
    core::{AltHyp, AcceptedHyp, StatsError},
    wilcoxon::RankSum,
};

fn sort_array(arr: &mut [f64]) {
    arr.sort_by(|a, b| a.partial_cmp(b).unwrap());
}

fn main() -> Result<(), StatsError> {
    let mut dat_x = vec![
        85., 90., 78., 92., 88., 76., 95., 89., 91., 82., 115., 120., 108., 122., 118., 106., 125.,
        119., 121., 112., 145., 150., 138., 152., 148., 136., 155., 149., 151., 142., 175., 180.,
        168., 182., 178., 166., 185., 179., 181., 172., 205., 210., 198., 212., 208., 196., 215.,
        209., 211., 202.,
    ];
    sort_array(&mut dat_x);

    let mut dat_y = vec![
        70., 85., 80., 90., 75., 88., 92., 79., 86., 81., 92., 100., 115., 110., 120., 105., 118.,
        122., 109., 116., 111., 122., 130., 145., 140., 150., 135., 148., 152., 139., 146., 141.,
        152., 160., 175., 170., 180., 165., 178., 182., 169., 176., 171., 182., 190., 205., 200.,
        210., 195., 208., 212., 199., 206., 201., 212.,
    ];
    sort_array(&mut dat_y);

    let rank_sum = RankSum::from_slices(&dat_x, &dat_y)?;
    let test_res = rank_sum.z_test(AltHyp::Gt, 0.05)?;
    assert_eq!(AcceptedHyp::Null, test_res.accepted());
    println!("test result: {test_res:?}");
    // test result: HypTestResult { p: 0.33244724790581637, alpha: 0.05, alt_hyp: Gt, accepted: Null }

    Ok(())
}
