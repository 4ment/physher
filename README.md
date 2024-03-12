# physher

[![CMake](https://github.com/4ment/physher/actions/workflows/cmake.yml/badge.svg)](https://github.com/4ment/physher/actions/workflows/cmake.yml)
[![License: GPL v2](https://img.shields.io/badge/License-GPLv2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0)

## About physher

``physher`` is a program for estimating evolutionary rates and divergence times from genetic, amino acid, codon, and generic data.

The current version of physher is incompatible with the first version.
Documentation for installing physher1 can be found [here](https://github.com/4ment/physher/wiki/Install) and the manual is located [here](https://github.com/4ment/physher/wiki/Usage).

## Getting Started

A C compiler such as ``gcc`` or ``clang`` is required. It is also requires the [GSL] library.
On Debian-based systems, dependencies can be installed via ``apt``:

```bash
sudo apt install gcc gsl
```

On MacOS, dependencies can be installed using a package manager such as [Homebrew](https://brew.sh).
```bash
brew install llvm gsl
```

Other package managers such as [conda](https://conda.io) and [MacPorts](https://www.macports.org) can also be used to install dependencies.

### Dependencies
 - [GSL]
 - [sse2neon] (already included in physher)

### Installation

To build ``physher`` from source you can run
```bash
git clone https://github.com/4ment/physher
cmake -S physher/ -B physher/build
cmake --build physher/build/ --target install
```

### Check install
If the installation was successful, this command should print the version of `physher`
```bash
physher
```

### Building C++ wrappers (optional)
A subset of physher's functionalities is exposed in C++ wrappers. These wrappers are used in [torchtree](https://github.com/4ment/torchtree), a python program, through bindings and [torchtree-physher](https://github.com/4ment/torchtree-physher).
A C++ compiler such as g++ or clang++ is required. Compilers can be installed using ``apt`` or ``homebrew``

```bash
git clone https://github.com/4ment/physher
cmake -S physher/ -B physher/build -DBUILD_CPP_WRAPPER=on
cmake --build physher/build/ --target install
```

### Testing (optional)

```bash
cmake -S physher/ -B physher/build -DBUILD_TESTING=on
cmake --build physher/build/ --target install
ctest --test-dir physher/build/
```

## Quick start
```bash
cd examples/fluA
physher JC69-time-ELBO.json
```

## physher in action

Some examples of projects using ``physher``
- [marginal-experiments](https://github.com/4ment/marginal-experiments): Evaluation of 19 dubious ways to compute marginal likelihood estimates. [10.1093/sysbio/syz046](https://doi.org/10.1093/sysbio/syz046).
- [phylostan](https://github.com/4ment/phylostan/tree/master/examples): Comparison of phylostan and ``physher`` using variational inference. [10.1101/702944](https://doi.org/10.1101/702944).
- [gradient-benchmark](https://github.com/4ment/gradient-benchmark): Benchmarking of automatic differentiation and analical gradients. [10.1093/gbe/evad099](https://doi.org/10.1093/gbe/evad099)
- [torchtree-physher](https://github.com/4ment/torchtree-physher): Plugin provinding fast calculation of phylogenetic functions in ``physher`` to [torchtree].

## License

Distributed under the GPLv2 License. See [LICENSE](LICENSE) for more information.

## Citing physher

Fourment M and Holmes EC. Novel non-parametric models to estimate evolutionary rates and divergence times from heterochronous sequence data. _BMC Evolutionary Biology_, 2014. doi: [10.1186/s12862-014-0163-6](https://doi.org/10.1186/s12862-014-0163-6)


[GSL]: https://www.gnu.org/software/gsl
[sse2neon]: https://github.com/DLTcollab/sse2neon
[torchtree]: https://github.com/4ment/torchtree