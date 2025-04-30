#!/bin/bash

cargo nextest run $1 --lib --bins --tests --features all --target-dir target/test-target

# Below runs any tests with names containing 'student':
# ./test-named student
