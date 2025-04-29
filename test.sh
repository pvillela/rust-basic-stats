#!/bin/bash

cargo nextest run --lib --bins --examples --tests --features all --target-dir target/test-target
cargo test --doc --features all
