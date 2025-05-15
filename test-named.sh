#!/bin/bash

NOCOVER=1 cargo nextest run $1 --target-dir target/test-target

# Below runs any tests (both unit and integration) with function names containing 'student' with default feature.
# ./test-named.sh student

# Below runs any unit tests with function names containing 'eq' with normal feature only. Notice need for '--lib' arg.
# ./test-named.sh "eq --lib --no-default-features --features normal"

# Below runs any unit and integration tests with function names containing 'eq' with normal feature only. 
# Notice need for '--tests' arg.
# ./test-named.sh "eq --tests --no-default-features --features normal"

# Below runs any integration tests with path names containing 'aok'. Notice '--test' arg and need to use '*'s.
# ./test-named.sh "--test *aok* --no-default-features --features normal"
