# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2025-05-XX

Includes breaking changes impacting all modules -- see below.

### Added

- Examples and tests for all modules.
- `binomial` module, which is an alias for the `bernoulli` module.

### Changed

#### Breaking

- Enhanced `error` module and implemented more appropriate error handling throughout the library.
- Thoroughly reimplemented `bernoulli` module.

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
