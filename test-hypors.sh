#!/bin/bash

# Below runs any tests with names containing 'hypors'.
cargo nextest run hypors --features wilcoxon --features _hypors --target-dir target/test-target
