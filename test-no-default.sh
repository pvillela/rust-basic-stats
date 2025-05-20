#!/bin/bash

# NOCOVER environment variable enables tests that are excluded from test coverage measurement.
export NOCOVER="1" 

echo "***** --no-default-features"
cargo nextest run --lib --bins --tests --no-default-features --target-dir target/test-target

echo "***** doc"
cargo test --doc
