/// Checks environment variable to include test only if coverage is not being measured.
/// Used to prevent trivial tests like `assert!(res.is_ok())` from inflating coverage numbers.
/// `NOCOVER` environment variable enables tests that are excluded from test coverage measurement.
pub fn nocover() -> bool {
    matches!(std::env::var("NOCOVER"), Ok(v) if v == "1" || v == "true")
}
