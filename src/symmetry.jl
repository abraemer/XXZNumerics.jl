module Symmetry

import Base
import ..XXZNumerics
using ..XXZNumerics: nspins

export FullZBasis, ZBlockBasis, zbasis, SpinFlip, symmetrize_state, symmetrize_op, basis_size

"""
    abstract type AbstractBasis

Supertype for a bases on N spins. A subtype needs to implement:
- basis_size to yield the number of basis vectors
- Base.iterator to give pairs of (coefficient, indices) to perform the projection into the symmetrized subspace
- maybe Base.length to give the (maximal) number of vectors that need to be combined to transform from the z-basis to the subspace (f.i. 2 for spin flips)

All transformations go from the full z-basis to the symmetrized subspace. Thus the all indices used are always in terms
of the full basis which runs from 1:2^N for N spins. The ordering of vector is such that the spin states treated as binary
binary number (with identification of |↑⟩ ≡ 0 and |↓⟩ ≡ 1) correspond to the index - 1. 

So for example:
- |↑…↑⟩ ≡ |0…0⟩ ≡ |0⟩ is index 1
- |↑…↑↓⟩ ≡ |0…01⟩ ≡ |1⟩ is index 2
- |↑↓↑↓⟩ ≡ |0101⟩ ≡ |5⟩ is index 6
- |↓…↓⟩ ≡ |1…1⟩ ≡ |2ᴺ-1⟩ is index 2ᴺ
"""
abstract type AbstractBasis{N} end

"""
    basis_size(basis)

Give size of the Hilbert space for the system.
"""
function basis_size end
basis_size(b::AbstractBasis) = length(first(b)[2])

XXZNumerics.nspins(b::AbstractBasis) = nspins(typeof(b))
XXZNumerics.nspins(::Type{<:AbstractBasis{N}}) where N = N

"""
    abstract ZBasis{N} <: AbstractBasis

Simple basis for N spins. This is basically the identity transformation. One can optionally choose to consider one magnetization block only.
"""
abstract type ZBasis{N} <: AbstractBasis{N} end

Base.iterate(zb::ZBasis) = ((1, _basis_inds(zb)), nothing)
Base.iterate(::ZBasis, state) = nothing # this basis is diagonal in the z basis by definition.
Base.IteratorSize(::Type{<:ZBasis}) = Base.HasLength()
Base.length(::ZBasis) = 1

"""
    ZBasis{N} <: AbstractBasis

Simple basis for N spins. This is basically the identity transformation.
"""
struct FullZBasis{N} <: ZBasis{N}
    FullZBasis(N::Int) = new{N}()
end

basis_size(::ZBasis{N}) where N = 2^N
_basis_inds(::ZBasis{N}) where N = 1:2^N

"""
    ZBlockBasis{N, K} <: ZBasis{N}

Basis corresponding to the block of the z-basis where there are k |↓⟩ or equivalently the eigenblock of the S_z operator with eigenvalue 1/2(N-2k).
"""
struct ZBlockBasis{N, K} <: ZBasis{N}
    ZBlockBasis(N::Int, k::Int) = 0 <= k <= N ? new{N, k}() : throw(ArgumentError("block index k=$k must be between (or equal to) 0 and N=$(N)!"))
end

basis_size(::ZBlockBasis{N, k}) where {N, k} = binomial(N, k)
function _basis_inds(::ZBlockBasis{N, k}) where {N, k}
    if k == 0 || k == N
        [2^k]
    end
    inds = ones(Int, binomial(N, k))# initialize with 1 because "index = binary representation + 1"
    _zblock_inds!(inds, N, k)
    inds
end

# build indices in-place recursively
function _zblock_inds!(states, N, k)
    if k == 0
        states
    elseif N == k
        states .+= 2^N - 1
    else
        front = binomial(N-1, k)
        @views states[front+1:end] .+= 2^(N-1)
        _zblock_inds!(view(states, 1:front), N-1, k)
        _zblock_inds!(view(states, front+1:length(states)), N-1, k-1)
        states
    end
end

## convenience constructors
"""
    zbasis(N [, k])

Construct a Basis object representing the z-basis or the block thereof, where there are k |↓⟩.
"""
zbasis(N) = FullZBasis(N)
zbasis(N, k) = ZBlockBasis(N, k)


abstract type AbstractSymmetry{N, B<:AbstractBasis{N}} <: AbstractBasis{N} end

"""
    symmetrize_state(basis, state)

Project state into the given (symmetrized) basis. 

Note: This function works only full, not-yet symmetrized states!
To combine different symmetries, combine the symmetry objects and then call symmetrize_state once!
"""
function XXZNumerics.symmetrize_state(b::AbstractBasis, state)
    res = similar(state, float(eltype(state)), basis_size(b))
    fill!(res, 0)
    for (coeff, indices) in b
        res += coeff*state[indices]
    end
    res
end

"""
    symmetrize_op(basis, op)

Project operator into the given (symmetrized) basis. 

Note: This function works only full, not-yet symmetrized operators
To combine different symmetries, combine the symmetry objects and then call symmetrize_op once!
"""
function XXZNumerics.symmetrize_op(b::AbstractBasis, op)
    res = similar(op, float(eltype(op)), (basis_size(b), basis_size(b)))
    fill!(res, 0)
    for (coeff, indices) in b
        for (coeff2, indices2) in b
            res += coeff*coeff2*op[indices, indices2]
        end
    end
    res
end

"""
    SpinFlip{P, N, B} <: AbstractSymmetry

Symmetry under total spinflip (ie. σx on all spins simultaneously) with parity P. P can be +1 or -1.
"""
struct SpinFlip{P, N, B} <: AbstractSymmetry{N, B}
    basis::B
    SpinFlip(basis::AbstractBasis, parity=1) = new{Int8(sign(parity)), nspins(basis), typeof(basis)}(basis)
end
SpinFlip(N::Int, parity=1) = SpinFlip(FullZBasis(N), parity)
#SpinFlip()
parity(::SpinFlip{P,N,B}) where {P,N,B} = P
_spinflip(N, ind) = 2^N + 1 - ind

Base.length(sf::SpinFlip) = 2*length(sf.basis)

basis_size(::SpinFlip{P, N, FullZBasis{N}}) where {N,P} = 2^(N-1)
basis_size(::SpinFlip{P, N, ZBlockBasis{N, k}}) where {N,P,k} = 2k == N ? div(binomial(N,k),2) : binomial(N,k)
function Base.iterate(sf::SpinFlip)
    baseIt = Iterators.Stateful(sf.basis)
    step = 0
    Base.iterate(sf, (step, baseIt))
end

##
## TODO:
## Optimize for non-symmetric symmetric sectors (where sector ≢ 0)
## This can be done by replacing unique! with a a similar function, that simply drops
## duplicates all-together instead of keeping only one instance.
## Look at the implementation in Base._groupedunique!
##
## TODO:
## Compute symmetrized indices only once for each set of basis indices
##
function Base.iterate(sf::SpinFlip, state)
    (step, baseIt) = state
    if Iterators.isempty(baseIt)
        return nothing
    end
    (base_coeff, base_inds) = peek(baseIt)
    coeff = base_coeff*parity(sf)^step/√2
    all_inds = unique!(sort!(vcat(base_inds, _spinflip.(nspins(sf), base_inds))))
    middle = div(length(all_inds), 2)
    if step == 0
        ((coeff, all_inds[1:middle]), (1, baseIt))
    else
        Iterators.iterate(baseIt) # this vector is done - move on
        ((coeff, reverse!(all_inds[middle+1:end])), (0, baseIt))
    end
end

end #module