#!/bin/bash

cargo makedocs -e log
cargo doc -p thread_map --no-deps --all-features
