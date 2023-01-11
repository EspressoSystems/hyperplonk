#!/usr/bin/env bash

#Fail out on error
set -e

cd ..
cd hyperplonk
# Run the benchmark binary
RAYON_NUM_THREADS=64 cargo bench --no-default-features --features=bench
RAYON_NUM_THREADS=32 cargo bench --no-default-features --features=bench
RAYON_NUM_THREADS=16 cargo bench --no-default-features --features=bench
RAYON_NUM_THREADS=8 cargo bench --no-default-features --features=bench
RAYON_NUM_THREADS=4 cargo bench --no-default-features --features=bench
RAYON_NUM_THREADS=2 cargo bench --no-default-features --features=bench
RAYON_NUM_THREADS=1 cargo bench --no-default-features --features=bench