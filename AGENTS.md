# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Common commands

| Task | Command |
|------|---------|
| Check all targets/features | `./check.sh` (wrapper for `cargo check --all-targets --all-features`) |
| Check specific feature combos | `./check-features.sh` (checks default, no-default, and each feature in isolation) |
| Run all tests | `./test.sh` (uses `cargo nextest` + `cargo test --doc`) |
| Run tests for specific feature combos | `./test-features.sh` |
| Run tests with no default features | `./test-no-default.sh` |
| Run a single/named test | `./test-named.sh "<filter>"` (e.g., `./test-named.sh "welch"`) |
| Lint | `./clippy.sh` (wrapper for `cargo clippy --all-targets --all-features`) |
| Generate docs | `./gen-doc.sh` (uses `cargo makedocs` and `cargo doc`) |

All test scripts set `NOCOVER=1` which enables tests excluded from coverage measurement.

## Architecture

This is a Rust library crate (edition 2024, `basic_stats`) providing basic parametric and non-parametric statistics. The only numeric types used are `u64` and `f64`. Functions operate on primitives, iterators, or slices — there is no parallel processing or dependence on large frameworks like `polars`.

### Feature-gated module tree

```
core (always enabled)
├── base          — SampleMoments, AltHyp, AcceptedHyp, HypTestResult, Ci
├── error         — StatsError (holds &'static str or String), StatsResult<T> alias
├── iter          — iter_with_counts() grouping contiguous equal values
├── check_interval — check_alpha_in_open_0_1 (available only when `normal` is on)
└── deterministic_sample — deterministic inverse-CDF sampling

normal (default, requires statrs)
├── z-scores, t-scores, z_alpha, t_alpha
├── student_1samp_* — one-sample and paired t-tests
├── welch_* — Welch two-sample t-test
└── deterministic_lognormal_sample

binomial (default, depends on normal)
└── Bernoulli/Binomial estimators, z-tests, CI for proportions

wilcoxon (default, depends on normal)
└── RankSum struct — Wilcoxon rank sum / Mann-Whitney U test

aok (off by default, feature "aok")
└── AokFloat / AokBasicStats traits — unwrap Result without panicking (returns NaN/fallback)
```

Feature chain: `binomial` → `normal` → `statrs`; `wilcoxon` → `normal` → `statrs`.

### Private `_dev_utils` feature

Gates `dev_utils` module (`ApproxEq` trait, `approx_eq!` and `rel_approx_eq!` macros). Used by tests but also kept as a public module under this feature for use by sibling crates.

### Key design rules

- **No panics**: Every function that can fail returns `StatsResult<T>` (alias for `Result<T, StatsError>`). The `FromIterator` trait is intentionally not implemented to avoid implicit panics.
- **Natural inputs**: Functions accept iterators, slices, or scalar values directly — users don't need to pre-aggregate into custom types.
- **Module-level docs** use `#![doc = include_str!("../examples/<name>.rs")]` to embed examples.
- **Rust edition 2024** — requires a recent toolchain. Features use the edition-2024 implicit feature resolver.
- **Tests require `_dev_utils` feature**: All `#[cfg(test)]` modules are gated with `#[cfg(feature = "_dev_utils")]`.

### Test structure

- Unit tests live in `#[cfg(test)]` blocks within each source file under `src/`.
- Integration tests are in `tests/` (excluded from the published crate via `Cargo.toml` exclude).
- Tests use `cargo nextest` (not plain `cargo test`) for unit and integration tests, plus `cargo test --doc` for doc-tests.
- The `R/` directory contains reference output files (text) used to validate statistical computations against R's implementations.
