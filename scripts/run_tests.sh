#!/usr/bin/env bash

#Fail out on error
set -e

# We want the code to panic if there is an integer overflow
export RUSTFLAGS="-C overflow-checks=on"

cargo test --release --all -- -Zunstable-options --report-time
cargo test --no-run --features=print-trace
cargo test --no-run --no-default-features
cargo bench --no-run