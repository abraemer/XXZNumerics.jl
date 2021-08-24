"""
    module Symmetry

## Motivation
Oftentimes physical systems exhibit symmetries (f.e. the system behaves identically when all
spins are flipped at once.) These symmetries are manifested as operators commuting with the
systems Hamiltonian and thus there exist common eigenvectors. This can be exploited if
eigenvectors of the symmetry's operator are known. In this case the Hamiltonian can be
transformed into a block diagonal structure and each block (corresponding to a single eigen-
value of the symmetry) can be diagonalized separately, speeding up the process enormously.

## Aim of the module
Here we provide functionality to symmetrize operators and states to conform to a specific
sector of a symmetry. Note that this is not equivalent to a simple projection when combining
symmetries! Example: SpinFlip-symmetry and S_z conservation are not compatible (flipping
causes S_z -> -S_z), however in this framework you can still write down `SpinFlip(zbasis(N,k))`.
You end up with spin flip invariant states built from states with k |↑⟩ (this also means
that `SpinFlip(zbasis(N,k))` is equivalent to `SpinFlip(zbasis(N,N-k))`). If we had simply
projected on the spin flip invariant part of a magnetization block we'd end up with nothing
most of the times (special case N even and k=N/2).

## Type structure
### Abstract types:
 - [`AbstractBasis`](@!ref) - root for everything, defines the symmetrization methods
 - [`ZBasis`](@!ref) - either complete ZBasis with [`FullZBasis`](@!ref) or a single magnetization block
 with [`ZBlockBasis`](@!ref). One of these is at the bottom of every composite symmetry.
 - [`AbstractSymmetry`](@!ref)
Almost identical to `AbstractBasis`. Allows composition with other `AbstractSymmetry`s.

### Symmetries
- Select a S_z magnetization block with `ZBlockBasis` (convenience method [`zbasis`](@!ref))
- Choose symmetric/antisymmetric subspaces with `SpinFlip`
- Choose a momentum sector (in translation invariant systems) with `Shift`

## Interfaces
All bases implement:
 - [`basis_size`](@!ref): give Hilbertspace dimension
 - [`nspins`](@!ref): return number of spins
 - [`symmetrize_state`](@!ref): apply symmetry to a state
 - [`symmetrize_op`](@!ref): apply symmetry to an operator
 - the iterator protocol

All symmetries implement (additionally):
 - chaining via: `SymmetryA(SymmetryB, argsB...), argsA...)`
 - [`base`](@!ref): get the underlying symmetry

Note: When chaining symmetries the order should be irrelevant.
Note 2: `ZBasis` actually implements `base` as well and returns `nothing`.

## Implementation
The abstract supertype of all symmetries is [`AbstractBasis`](@!ref).

All transformations go from the full z-basis to the symmetrized subspace. Thus the all indices used are always in terms
of the full basis which runs from 1:2^N for N spins. The ordering of vector is such that the spin states treated as binary
binary number (with identification of |↑⟩ ≡ 0 and |↓⟩ ≡ 1) correspond to the index - 1.
So for example:
- |↑…↑⟩ ≡ |0…0⟩ ≡ |0⟩ is index 1
- |↑…↑↓⟩ ≡ |0…01⟩ ≡ |1⟩ is index 2
- |↑↓↑↓⟩ ≡ |0101⟩ ≡ |5⟩ is index 6
- |↓…↓⟩ ≡ |1…1⟩ ≡ |2ᴺ-1⟩ is index 2ᴺ

The symmetrization is done via the iteration protocol. In each step the symmetry returns the
indices of the basis vectors and the coefficient for all of them (may be a single number or
an appropriately sized array).

### Simple example
N = 2 -> Basis |↑↑⟩, |↑↓⟩, |↓↑⟩, |↓↓⟩  (in that order)
SpinFlip gives (coefficient, indices)
1. `(1/√2 , [1,2])`
2. `(1/√2 , [4,3])`
which translates to the symmetrized basis
|↑↑⟩+|↓↓⟩, |↑↓⟩+|↓↑⟩
Verify using `collect(SpinFlip(2))`.
"""
module Symmetry

import Base
import ..XXZNumerics
using ..XXZNumerics: nspins

export FullZBasis, ZBlockBasis, zbasis, SpinFlip, symmetrize_state, symmetrize_op, basis_size, parity

include("abstract.jl")

include("zbasis.jl")

include("cyclic.jl")

end #module
