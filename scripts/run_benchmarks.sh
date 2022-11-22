#!/bin/bash
cd ..
cd hyperplonk
# Run the benchmark binary
cargo bench 64 --no-default-features --features=bench
cargo bench 32 --no-default-features --features=bench
cargo bench 16 --no-default-features --features=bench
cargo bench 8 --no-default-features --features=bench
cargo bench 4 --no-default-features --features=bench
cargo bench 2 --no-default-features --features=bench
cargo bench 1 --no-default-features --features=bench