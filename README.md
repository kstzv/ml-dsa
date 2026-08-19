# ML-DSA

An implementation of the Module-Lattice-Based Digital Signature Algorithm (ML-DSA) in C11, based on FIPS 204.

## Development status

This project is at an early stage of development. The code is incomplete, untested, and subject to significant changes.

It must not be used for cryptographic or production purposes.

## Current work

The repository currently contains the initial ML-DSA core, including:

* algorithm constants and parameter sets;
* an initial representation of the key components;
* Barrett and Montgomery reductions;
* the forward and inverse Number Theoretic Transform (NTT).

Additional algorithm components, platform ports, tests, and the public API will be added during further development.
