#!/bin/bash

echo "***** --all-targets --all-features"
cargo check --all-targets --all-features

echo "***** --all-targets (default feature)"
cargo check --all-targets

echo "***** --no-default-features"
cargo check --lib --bins --tests --no-default-features

echo "***** --features normal"
cargo check --lib --bins --tests --no-default-features --features normal

echo "***** --features binomial"
cargo check --lib --bins --tests --no-default-features --features binomial

echo "***** --features wilcoxon"
cargo check --lib --bins --tests --no-default-features --features wilcoxon
