"""
    abstract type AbstractBasis{N}

Supertype for bases on N spins. A subtype needs to implement:
- Base.iterator to give pairs of (coefficient, indices) to perform the projection into the symmetrized subspace

A subtype MAY want to implement (for performance reasons):
- `basis_size` to yield the number of basis vectors
- `Base.length` to give the (maximal) number of vectors that need to be combined to transform from the z-basis to the subspace (f.i. 2 for spin flips)
- `Base.IteratorSize` to compute the length (you probably know it if you implement `Base.length`)
"""
abstract type AbstractBasis{N} end

"""
    basis_size(basis)

Give size of the Hilbert space for the system.
"""
basis_size(b::AbstractBasis) = length(first(b)[2])

"""
    base(symmetry)

Get the underlying basis of the symmetry.
"""
function base end

"""
    length(::AbstractBasis)

Return number of components of the symmetry, i.e. the size of the group for cyclic operators.
"""
Base.length(b::AbstractBasis) = length(collect(b))

## iteration defaults, not clever but rather safe I think.
## subtypes should opt to overwrite this
Base.IteratorSize(::AbstractBasis) = Base.SizeUnknown()

XXZNumerics.nspins(b::AbstractBasis) = nspins(typeof(b))
XXZNumerics.nspins(::Type{<:AbstractBasis{N}}) where N = N


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
        res .+= coeff .* @view state[indices]
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
            res .+= coeff .* coeff2 .* @view op[indices, indices2]
        end
    end
    res
end

"""
    abstract type AbstractSymmetry{N, B} <: AbstractBasis{N}

Supertype for a symmetry on N spins acting on a Basis B. A subtype needs to implement:
- everything need for [`AbstractBasis`](@!ref)
- `base` to return the basis B acting upon.
Default is to use the field `base`. If you do not change this, then you do NOT need to implement
a new method.
"""
abstract type AbstractSymmetry{N, B<:AbstractBasis{N}} <: AbstractBasis{N} end

"""
    base(symmetry)

Get the underlying basis `symmetry` acts on. It may be another symmetry.
"""
base(symmetry::AbstractSymmetry) = symmetry.base
