# Hyperplonk Library

[![Rust](https://github.com/hyperplonk/hyperplonk/workflows/Rust/badge.svg)](https://github.com/hyperplonk/hyperplonk/actions)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)

Hyperplonk is an efficient implementation of a linear-time FFT-free SNARK proof system, based on the [research paper](https://eprint.iacr.org/2022/1355.pdf). The system offers improved performance and scalability compared to traditional approaches.

## Key Features

- Linear time complexity
- No FFT requirements
- Enhanced scalability
- Optimized performance

## Disclaimer

**DISCLAIMER:** This software is provided "as is" and its security has not been externally audited. Use at your own risk.

## Getting Started

### System Requirements

- Rust (latest stable version)
- Nix package manager
- Git

### Installing Rust

We recommend using nix for installing the correct version of Rust and additional libraries:

```bash
> curl -L https://nixos.org/nix/install | sh
```

### Building the Project

1. Clone the repository:
```bash
> git clone https://github.com/your-username/hyperplonk.git
> cd hyperplonk
```

2. Build the project:
```bash
> nix-shell
> cargo build
```

### Development Tools

We recommend the following tools:

- [nix](https://nixos.org/download.html) - for dependency management
- [direnv](https://direnv.net/docs/installation.html) - for environment management

After installing direnv, run:
```bash
> direnv allow
```

## Development

### Running Tests

```bash
> cargo test --release --all
```

### Documentation

Generate and view documentation:
```bash
> cargo doc --open
```

### Code Formatting

```bash
> cargo fmt
```

### Updating Dependencies

To update all dependencies:
```bash
> nix flake update
```

To update a specific dependency:
```bash
> nix flake update github:oxalica/rust-overlay
```

### Benchmarks

To run benchmarks, use the script:
```bash
> ./scripts/run_benchmarks.sh
```

Results can be compared with Tables 5 and 6 in the [original paper](https://eprint.iacr.org/2022/1355.pdf).

## Contributing

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to your fork (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Contribution Guidelines

- Follow the existing code style
- Add tests for new functionality
- Update documentation as needed
- Ensure all tests pass

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.
