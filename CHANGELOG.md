# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2025-06-XX

### Added

- Enum `AcceptedHyp`, which replaces enum `Hyp`, removing a bit of redundancy.
- Functions to generate deterministic samples of distributions.

### Changed

- Added argument to Welch statistics to enable testing that the difference of means equals a specific value which doesn't have to be zero. To migrate from the previous version:
  - Use `0.0` as the value for the new argument.
- Changed `StatsError` so it can hold either a `&'static str` or a `String`. To migrate from the previous version:
  - Construction of a `StatsError` with a `&'static str` arguement `s` changes from `StatsError(s)` to `StatsError::new(s)`.
  - Access to the message contained in a `StatsError` instance `e` changes from `&e.0` to `e.msg()`.
- Updated `Cargo.toml` exclusions with `R` and `tests` directories.

### Removed

- Removed `Hyp` enum. It was replaced with enum `AcceptedHyp`, removing a bit of redundancy.

## [1.0.0] - 2025-05-20

Major update. Includes extensive breaking changes impacting all modules.

### Added

- Extensive error handling.
- Extensive tests, both positive and negative.
- Examples.
- Module `aok`.

### Changed

#### Breaking

- Enhanced error handling throughout the library. Many functions now return a `Result`.
- Modules `iter` and `error` are now subsumed under module `core`.
- Module `bernoulli` was renamed `binomial` and thoroughly reimplemented.
- Each module other than `core` is gated by a feature with the same name as the module.
- Introduced `default` feature which includes features `normal`, `binomial`, and `wilcoxon`.

#### Non-breaking

- Removed `Clone` constraint from iter_with_counts.
- Fixed major bugs in `z_alpha` and `t_alpha` functions.

## [0.1.1] - 2025-04-30

### Changed

- Fixed incorrect documentation for `SampleMoments` struct.

## [0.1.0] - 2025-04-30

Initial release, based on code extracted from `bench_diff` crate:
- Removed dependence on `hdrhistogram`.
- Additional types and functions.
- Additional documentation.
- Additional tests.
