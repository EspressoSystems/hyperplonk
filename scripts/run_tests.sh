#!/usr/bin/env bash

#Fail out on error
set -e

# We want the code to panic if there is an integer overflow
export RUSTFLAGS="-C overflow-checks=on"

cargo test --release --all
cargo test --no-run --features=print-trace
cargo bench --no-run
