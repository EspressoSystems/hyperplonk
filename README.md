# Hyperplonk library
A linear-time FFT-free SNARK proof system (https://eprint.iacr.org/2022/1355.pdf).

## Disclaimer

**DISCLAIMER:** This software is provided "as is" and its security has not been externally audited. Use at your own risk.

## Project Overview

Hyperplonk is an implementation of a linear-time FFT-free SNARK proof system based on the research paper [Hyperplonk: Plonk with Linear-Time Prover and High-Degree Custom Gates](https://eprint.iacr.org/2022/1355.pdf). It provides efficient zero-knowledge proofs with the following features:

- Linear-time prover complexity
- No FFT requirements
- Support for high-degree custom gates
- Optimized for performance in zero-knowledge applications

## Project Structure

The repository is organized as follows:

- `subroutines/`: Core components of the Hyperplonk system
  - `src/pcs/`: Polynomial Commitment Scheme implementation using KZG
  - `src/poly_iop/`: Polynomial IOP implementations (sum checks, zero checks, product checks, permutation checks)
- `scripts/`: Utility scripts for benchmarking and testing
- `tests/`: Test suite for the library

## Development environment setup

### Install RUST

We recommend using nix for installing the correct version of rust and
additional libraries:

```bash
> curl -L https://nixos.org/nix/install | sh
```

### Compiling the project for the first time

```bash
> nix-shell
> cargo build
```

### Direnv

We recommend the following tools:

- [nix](https://nixos.org/download.html)
- [direnv](https://direnv.net/docs/installation.html)

Run `direnv allow` at the repo root. You should see dependencies (including Rust) being installed (the first time might take a while). 
Upon modification on `flake.nix`, run `direnv reload` to reflect new dependencies or environment setups.

### Tests

```
> cargo test --release --all
```

### Generate and read the documentation

#### Standard

```
> cargo doc --open
```

#### Additional Documentation Resources

For a deeper understanding of the underlying concepts and implementation details:
- Refer to the [Hyperplonk paper](https://eprint.iacr.org/2022/1355.pdf)
- Check the inline documentation in the source code
- See the readme files in the subroutines directories for component-specific information

### Code formatting

To format your code run

```
> cargo fmt
```

### Updating non-cargo dependencies

Run `nix flake update` if you would like to pin other version edit `flake.nix`
beforehand. Commit the lock file when happy.

To update only a single input specify it as an argument, for example

    nix flake update github:oxalica/rust-overlay

## Usage Examples

### Basic Usage

```rust
// Example code for creating a simple proof
use hyperplonk::{setup, prove, verify};

// Setup the proving system
let params = setup(circuit_size);

// Generate a proof
let proof = prove(&params, &circuit, &witness);

// Verify the proof
let is_valid = verify(&params, &proof, &public_inputs);
```

### Advanced Configuration

For advanced usage and configuration options, refer to the API documentation generated with `cargo doc`.

## Benchmarks

To obtain benchmarks, run the script file `scripts/run_benchmarks.sh`. 

The benchmarks measure:
- Prover time
- Verifier time
- Proof size
- Memory usage

For reference benchmarks, we refer to Table 5 and Table 6 in the [Hyperplonk paper](https://eprint.iacr.org/2022/1355.pdf). These tables provide performance comparisons with other SNARK systems across different circuit sizes.

To interpret your benchmark results:
1. Run the benchmark script
2. Compare the output metrics with the reference tables
3. Note that performance may vary based on your hardware configuration

## Contributing

Contributions to Hyperplonk are welcome! Here's how you can contribute:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

Please ensure your code follows the project's coding standards and includes appropriate tests.

## License

This project is licensed under the terms specified in the LICENSE file.

## Contact

For questions or feedback about Hyperplonk, please open an issue in the GitHub repository.
