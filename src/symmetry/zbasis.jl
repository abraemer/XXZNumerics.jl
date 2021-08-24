"""
    abstract ZBasis{N} <: AbstractBasis

Simple basis for N spins. This is basically the identity transformation.
One can optionally choose to consider one magnetization block only.
"""
abstract type ZBasis{N} <: AbstractBasis{N} end



## implement iterator protocol
## only thing unknown: which indices to select
Base.iterate(zb::ZBasis) = ((1, _basis_inds(zb)), nothing)
Base.iterate(::ZBasis, state) = nothing # this basis is diagonal in the z basis by definition.
Base.IteratorSize(::Type{<:ZBasis}) = Base.HasLength()
Base.length(::ZBasis) = 1

base(::ZBasis) = nothing # we are at top of the chain # is this a sensible value?

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

## build indices in-place recursively
## Note: Indices computed are from 0 to 2^N -1 -> need to +1
## Implementation based on recursive factorization:
## - {|N, k⟩} = |↑⟩ ⊗ {|N-1,k⟩} ∪ |↓⟩ ⊗ {|N-1,k-1⟩}
## - {|N, N⟩} = {|2^N -1⟩}
## - {|N, 0⟩} = {|0⟩}
## - |1⟩ ⊗ |s⟩ = |1*2^N + s⟩ (|s⟩ has N spins)
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
