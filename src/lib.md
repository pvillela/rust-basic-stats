A lightweight library that provides some basic parametric and non-parametric statistics and hypothesis tests.

## Cargo features

By default, use of this library as a dependency includes modules [`core`], [`normal`], [`binomial`], and [`wilcoxon`]. The [`aok`] module is not included by default.

Each module other than [`core`] has an associated cargo feature that enables the module. To include only selected modules, specify `default-features = false` in the dependency declaration (or `--no-default-features` on the command line) and specify the desired features in the dependency declaration (or command line).

## Error handling

Functions in this library are designed not to panic. Although functions in crates this library depends on may panic, this library implements extensive error handling, validation, and testing to prevent panics. If you encounter a panic, it is a bug, so please [create an issue](https://github.com/pvillela/rust-basic-stats/issues/new) to report it.

Except in a couple of specifically documented places, functions in this library only return finite (i.e., not `NaN`, `Infinity`, or `-Infinity`) values unless a non-finite value is provided as an input.