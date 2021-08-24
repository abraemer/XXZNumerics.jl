"""
    AbstractCyclicSymmetry{N, B}

Supertype for all cyclic symmetries.

A cyclic symmetry P has the property that it generates a cyclic group `<P>` of size L.

As such it has L sectors corresponding to the eigenvalues `λ_k = exp(2πik/L)`.

Symmetrization is easy. In sector k, we map the state
|ψ⟩ → 1/√L ∑ λ_k^l P^l |ψ⟩

All that is needed for the concrete CyclicSymmetry are 3 things:
 - `_sector_count` to determine L
 - `_sector` to determin k
 - `_shift_index` to compute the resulting index b in `|b⟩ = P |a⟩`

Note: Here we don't consider the case where the symmetry is not compatible with a z-basis representation.
If this is needed at some point, it should be implemented separately.
"""
abstract type AbstractCyclicSymmetry{N, B} <: AbstractSymmetry{N, B} end

"""
    _shift_index(cyclic_symm, i)

compute the resulting index j in `|j⟩ = P |i⟩`.
"""
function _shift_index end

"""
    _sector_count(cyclic_symm)

Compute the number L of sectors of the symmetry.
"""
function _sector_count end

"""
    _sector(cyclic_symm)

Get the sector index k. The corresponding eigenvalue is always `λ_k = exp(2πik/L)`.
Note: `k=0` is always the fully symmetric subspace.
"""
function _sector end

Base.length(csymm::AbstractCyclicSymmetry) = _sector_count(csymm)*length(base(csymm))
Base.IteratorSize(::AbstractCyclicSymmetry) = Base.HasLength()

function Base.iterate(csymm::AbstractCyclicSymmetry)
    baseIt = Iterators.Stateful(base(csymm))
    step = 0
    Base.iterate(csymm, (step, baseIt))
end

function Base.iterate(csymm::AbstractCyclicSymmetry, state)
    ## baseIt goes through the components of the underlying basis
    ## step goes through the components of this symmetry. 0 <= step < L
    (step, baseIt) = state
    if Iterators.isempty(baseIt)
        return nothing
    end
    (base_coeff, base_inds) = peek(baseIt)
    L = _sector_count(csymm)
    k = _sector(csymm)
    phase = L == 2 ? (-1)^k : exp(2π*im*k/L)
    coeff = base_coeff * phase^step / √L

    ## compute all indices
    inds = collect(base_inds)
    all_inds = [inds]
    sizehint!(all_inds, L)
    for _ in 1:L
        inds = _shift_index.(Ref(csymm), inds)
        push!(all_inds, inds)
    end

    ## only keep unique ones and then take the first 1/L as representatives
    ## indices appearing multiple times mean either
    ## 1) there are subcycles (|↑↓↑↓⟩ -> |↓↑↓↑⟩ -> |↑↓↑↓⟩ has length 2 not 4)
    ## 2) different indices below to the same cycle (|↑↓⟩ and |↓↑⟩)
    ## Need to remove the duplicates to not introduce normalization issues
    uniq_inds = unique!(vcat(all_inds...))
    repr_count = div(length(uniq_inds), L)
    representatives = uniq_inds[1:repr_count]

    ## shift representatives to correct k
    for _ in 1:step
        representatives = _shift_index.(Ref(csymm), representatives)
    end

    nextstate = (mod(step+1, L), baseIt)
    step == L-1 && Iterators.iterate(baseIt) # this vector is done - move on

    ((coeff, representatives), nextstate)
end


"""
    SpinFlip{P, N, B} <: AbstractCyclicSymmetry

Symmetry under total spinflip (ie. σx on all spins simultaneously) with parity P. P can be +1 or -1.
"""
struct SpinFlip{P, N, B} <: AbstractCyclicSymmetry{N, B}
    base::B
    SpinFlip(basis::AbstractBasis, parity=1) = new{Int8(sign(parity)), nspins(basis), typeof(basis)}(basis)
end

SpinFlip(N::Int, parity=1) = SpinFlip(FullZBasis(N), parity)

parity(::SpinFlip{P,N,B}) where {P,N,B} = P

_sector_count(::SpinFlip) = 2
_sector(::SpinFlip{P, N, B}) where {P,N,B} = P == 1 ? 0 : 1
_shift_index(::SpinFlip{P, N, B}, index) where {P,N,B} = _spinflip(N, index)
_spinflip(N, ind) = 2^N + 1 - ind

basis_size(::SpinFlip{P, N, FullZBasis{N}}) where {N,P} = 2^(N-1)
basis_size(::SpinFlip{P, N, ZBlockBasis{N, k}}) where {N,P,k} = 2k == N ? div(binomial(N,k),2) : binomial(N,k)
# function Base.iterate(sf::SpinFlip)
#     baseIt = Iterators.Stateful(sf.basis)
#     step = 0
#     Base.iterate(sf, (step, baseIt))
# end

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
# function Base.iterate(sf::SpinFlip, state)
#     (step, baseIt) = state
#     if Iterators.isempty(baseIt)
#         return nothing
#     end
#     (base_coeff, base_inds) = peek(baseIt)
#     coeff = base_coeff*parity(sf)^step/√2
#     all_inds = unique!(sort!(vcat(base_inds, _spinflip.(nspins(sf), base_inds))))
#     middle = div(length(all_inds), 2)
#     if step == 0
#         ((coeff, all_inds[1:middle]), (1, baseIt))
#     else
#         Iterators.iterate(baseIt) # this vector is done - move on
#         ((coeff, reverse!(all_inds[middle+1:end])), (0, baseIt))
#     end
# end
