#!/bin/bash

# NOCOVER environment variable enables tests that are excluded from test coverage measurement.

NOCOVER="1" cargo nextest run --all-targets --all-features --target-dir target/test-target
cargo test --doc
