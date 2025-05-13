/// Checks environment variable to include test only if coverage is not being measured.
pub fn nocover() -> bool {
    matches!(std::env::var("NOCOVER"), Ok(v) if v == "1" || v == "true")
}
