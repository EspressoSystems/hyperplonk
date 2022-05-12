#!/usr/bin/env bash

# We want the code to panic if there is an integer overflow
export RUSTFLAGS="-C overflow-checks=on"

cargo test --release -- -Zunstable-options --report-time
cargo test --no-run --features=print-trace
cargo test --no-run --no-default-features
