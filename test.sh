#!/bin/bash

export RUSTFLAGS="-Awarnings"

# NOCOVER environment variable enables tests that are excluded from test coverage measurement.
export NOCOVER="1"

./check-features.sh || { echo "Error: check-features failed"; exit 1; }

echo "***** test all except benches, all features"
cargo nextest run  --lib --bins --examples --tests --all-features --target-dir target/test-target

echo "***** test doc"
cargo test --doc --features aok --target-dir target/test-target
