#!/bin/bash

# NOCOVER environment variable enables tests that are disabled for test coverage measurement.

NOCOVER="1" cargo nextest run --lib --bins --examples --tests --target-dir target/test-target
cargo test --doc
