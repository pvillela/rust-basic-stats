#!/bin/bash

cargo nextest run $1 --lib --bins --tests --target-dir target/test-target

# Below runs any tests with names containing 'student':
# ./test-named student
