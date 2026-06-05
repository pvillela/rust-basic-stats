#!/bin/bash
set -e

# Ensure llvm-tools and cargo-llvm-cov are available.
if ! command -v cargo-llvm-cov &> /dev/null; then
    echo "cargo-llvm-cov not found."
    echo "Install with: cargo install cargo-llvm-cov"
    echo "Also ensure llvm-tools is installed: rustup component add llvm-tools-preview"
    exit 1
fi

cargo llvm-cov clean

rm -rf coverage
mkdir -p coverage

export RUSTFLAGS="-Awarnings"

# Coverage instrumentation adds overhead that can cause timing-sensitive tests
# to fail. Use --no-fail-fast so all tests run and profile data is collected.
echo "***** Running tests with coverage instrumentation (nextest)"
cargo llvm-cov test --tests --all-features \
    --ignore-run-fail --no-report \
    || echo "(some tests may have failed due to instrumentation overhead; continuing)"

echo "***** Generating lcov"
cargo llvm-cov report --lcov --output-path coverage/lcov.info

echo "***** Generating html"
cargo llvm-cov report --html --output-dir coverage

echo "*****"
echo "***** Coverage reports:"
echo "*****   HTML:  coverage/html/index.html"
echo "*****   lcov:  coverage/lcov.info"
echo "*****"
