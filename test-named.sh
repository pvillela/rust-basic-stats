#!/bin/bash

cargo nextest run $1 --lib --bins --tests --target-dir target/test-target

# Below runs any tests with names containing 'student' with default feature:
# ./test-named.sh student

# Below runs any tests with names containing 'eq' with normal feature only:
# ./test-named.sh "eq --no-default-features --features normal"
