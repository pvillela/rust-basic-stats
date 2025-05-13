#!/bin/bash

# NOCOVER environment variable enables tests that are excluded from test coverage measurement.

NOCOVER="1" cargo nextest run --lib --bins --examples --tests --target-dir target/test-target
cargo test --doc
