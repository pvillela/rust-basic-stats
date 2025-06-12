use basic_stats::approx_eq;

fn main() {
    let x: f64 = 123.444444;
    let y: f64 = 123.444454;
    let z: f64 = 123.444455;
    let epsilon: f64 = 0.00001;

    // This assertion succeeds.
    approx_eq!(x, y, epsilon);

    // This assertion fails.
    approx_eq!(x, z, epsilon);
}
