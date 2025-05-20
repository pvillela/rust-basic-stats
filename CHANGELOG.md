# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-05-XX

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
