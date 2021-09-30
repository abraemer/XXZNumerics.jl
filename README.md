# XXZNumerics

[![codecov](https://codecov.io/gh/abraemer/XXZNumerics.jl/branch/main/graph/badge.svg?token=XN6TT95A53)](https://codecov.io/gh/abraemer/XXZNumerics.jl)
[![CI](https://github.com/abraemer/XXZNumerics.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/abraemer/XXZNumerics.jl/actions/workflows/ci.yml)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

Just some basics for generating an XXZ hamiltonian and field term. Maybe will grow into a more sophisticated package for analyzing thermal properties of disordered spin models.

Includes:
- some spin-half basics
- different geometries to sample blockaded positions from
- a basic interface to convert particle positions to interaction strength (only `PowerLaw` is implemented)
- methods for generating an XXZ Hamiltonian and field term
- methods for computing thermal quantum ensembles
