/// Checks environment variable to include test only if coverage is not being measured.
pub fn nocover() -> bool {
    match std::env::var("NOCOVER") {
        Ok(v) if v == "1" || v == "true" => true,
        _ => false,
    }
}
