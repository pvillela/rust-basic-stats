#!/bin/bash

rm -r target/doc

cargo makedocs -e hypors -e polars -e statrs
cargo doc -p basic_stats --no-deps --all-features
