#!/bin/bash

# NOCOVER environment variable enables tests that are excluded from test coverage measurement.
export NOCOVER="1" 

echo "***** --examples --all-features"
cargo nextest run --examples  --all-features --target-dir target/test-target

echo "*****  --lib --bins --tests (default feature)"
cargo nextest run --lib --bins --tests --target-dir target/test-target

echo "***** --no-default-features"
cargo nextest run --lib --bins --tests --no-default-features --target-dir target/test-target

echo "***** --features normal"
cargo nextest run --lib --bins --tests --no-default-features --features normal --target-dir target/test-target

echo "***** --features binomial"
cargo nextest run --lib --bins --tests --no-default-features --features binomial --target-dir target/test-target

echo "***** --features wilcoxon"
cargo nextest run --lib --bins --tests --no-default-features --features wilcoxon --target-dir target/test-target

echo "***** doc"
cargo test --doc
