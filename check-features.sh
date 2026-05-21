#!/bin/bash

set -e  # Stop script immediately on any error

### With default features: Externally exposed feature combinations

echo "***** (default feature)"
cargo check --lib --tests

echo "*****  --features aok"
cargo check --lib --tests  --features aok

echo "*****  --features aok_f64"
cargo check --lib --tests  --features aok_f64

echo "*****  --features aok_stats"
cargo check --lib --tests  --features aok_stats

echo "*****  --features binomial"
cargo check --lib --tests  --features binomial

echo "*****  --features normal"
cargo check --lib --tests  --features normal

echo "*****  --features wilcoxon"
cargo check --lib --tests  --features wilcoxon

echo "*****  --features _dev_utils"
cargo check --lib --tests  --features _dev_utils

### All targets and features

echo "***** --all-targets --all-features"
cargo check --all-targets --all-features

### Without default features: Externally exposed feature combinations

echo "***** --no-default-features"
cargo check --lib --tests --no-default-features

echo "***** --no-default-features --features aok"
cargo check --lib --tests --no-default-features --features aok

echo "***** --no-default-features --features aok_f64"
cargo check --lib --tests --no-default-features --features aok_f64

echo "***** --no-default-features --features aok_stats"
cargo check --lib --tests --no-default-features --features aok_stats

echo "***** --no-default-features --features binomial"
cargo check --lib --tests --no-default-features --features binomial

echo "***** --no-default-features --features normal"
cargo check --lib --tests --no-default-features --features normal

echo "***** --no-default-features --features wilcoxon"
cargo check --lib --tests --no-default-features --features wilcoxon

echo "***** --no-default-features --features _dev_utils"
cargo check --lib --tests --no-default-features --features _dev_utils

### Benches

# None

