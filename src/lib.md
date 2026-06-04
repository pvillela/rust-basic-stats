A lightweight library that provides some basic parametric and non-parametric statistics and hypothesis tests.

This library strives for ease of use and small size. The only numeric types used are `u64` and `f64`. Functions in this library operate on primitive types, iterators, or slices.

There is no support for parallel processing. In particular, there is no dependence on large libraries like `polars` that support parallel processing on custom data structures.

# Cargo features

By default, use of this library as a dependency includes modules [`core`], [`normal`], [`binomial`], [`wilcoxon`], and [`detm_samp`]. The `aok` module is not included by default.

Each module other than [`core`] (which is always enabled) has an associated cargo feature that enables the module. To include only selected modules, specify `default-features = false` in the dependency declaration (or `--no-default-features` on the command line) and specify the desired features in the dependency declaration (or command line).

# Error handling

Functions in this library are designed not to panic. Although functions in crates this library depends on may panic, this library implements extensive error handling, validation, and testing to prevent panics. If you encounter a panic, it is a bug, so please [create an issue](https://github.com/pvillela/rust-basic-stats/issues/new) to report it.

Except in a couple of specifically documented places, functions in this library only return finite (i.e., not `NaN`, `Infinity`, or `-Infinity`) values unless a non-finite value is provided as an input.

# Migration Guide

This section describes changes from version 1.0.0 to 2.0.0 that require updates to calling code.

## `Hyp` replaced by `AcceptedHyp`

The `Hyp` enum has been removed and replaced by `AcceptedHyp`. Update any references to `Hyp` to use `AcceptedHyp` instead.

## `StatsError` construction and message access

Construction of a `StatsError` with a `&'static str` argument `s` changes from `StatsError(s)` to `StatsError::new(s)`. Access to the message contained in a `StatsError` instance `e` changes from `&e.0` to `e.msg()`.

## Additional argument to Welch statistics

Functions that compute Welch statistics now take an additional argument specifying the value against which the difference of means is tested. Use `0.0` as that argument.

## `RankSum` constructors reject empty samples

`RankSum::from_iters_with_counts`, `from_iters`, and `from_slices` now return an error if either sample is empty. Consequently, `z()`, `z_p()`, and `z_test()` no longer return an error for empty samples.