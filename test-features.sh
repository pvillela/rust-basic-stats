#!/bin/bash

# NOCOVER environment variable enables tests that are excluded from test coverage measurement.
export NOCOVER="1" 

echo "*****  --lib --bins --tests (default feature)"
cargo nextest run --lib --bins --tests --features _dev_utils --target-dir target/test-target

echo "***** --no-default-features"
cargo nextest run --lib --bins --tests --no-default-features --features _dev_utils --target-dir target/test-target

echo "***** --features normal"
cargo nextest run --lib --bins --tests --no-default-features --features _dev_utils,normal --target-dir target/test-target

echo "***** --features binomial"
cargo nextest run --lib --bins --tests --no-default-features --features _dev_utils,binomial --target-dir target/test-target

echo "***** --features wilcoxon"
cargo nextest run --lib --bins --tests --no-default-features --features _dev_utils,wilcoxon --target-dir target/test-target

echo "***** --features _dev_utils"
cargo nextest run --lib --bins --tests --no-default-features --features _dev_utils --target-dir target/test-target

echo "***** doc"
cargo test --doc
