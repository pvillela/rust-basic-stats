A lightweight library that provides some basic parametric and non-parametric statistics and hypothesis tests.

This library strives for ease of use and small size. The only numeric types used are `u64` and `f64`. Functions in this library operate on primitive types, iterators, or slices.

There is no support for parallel processing. In particular, there is no dependence on large libraries like `polars` that support parallel processing on custom data structures.

# Cargo features

By default, use of this library as a dependency includes modules [`core`], [`normal`], [`binomial`], and [`wilcoxon`]. The [`aok`] module is not included by default.

Each module other than [`core`] (which is always enabled) has an associated cargo feature that enables the module. To include only selected modules, specify `default-features = false` in the dependency declaration (or `--no-default-features` on the command line) and specify the desired features in the dependency declaration (or command line).

# Error handling

Functions in this library are designed not to panic. Although functions in crates this library depends on may panic, this library implements extensive error handling, validation, and testing to prevent panics. If you encounter a panic, it is a bug, so please [create an issue](https://github.com/pvillela/rust-basic-stats/issues/new) to report it.

Except in a couple of specifically documented places, functions in this library only return finite (i.e., not `NaN`, `Infinity`, or `-Infinity`) values unless a non-finite value is provided as an input.