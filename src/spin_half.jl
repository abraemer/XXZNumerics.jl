export σx, σy, σz, 𝟙, up, down, ⊗, speye, identity_op, single_spin_op, correlator, op_list, symmetrize_state, symmetrize_op, nspins


const σx = sparse([0 1; 1 0])
const σy = sparse([0 -im; im 0])
const σz = sparse([1 0; 0 -1])
const 𝟙 = spdiagm([1, 1])

const up = [1,0]
const down = [0,1]

const ⊗ = kron

speye(k::Integer) = spdiagm([1 for _ in 1:k])
identity_op(k::Integer) = speye(2^k)

"""
    nspins(obj)

Give number of spins in the system this object belongs to.
"""
function nspins(system_size::Int)
    N = Int(ceil(log2(system_size)))
    2^N == system_size || throw(ArgumentError("$(system_size) is not a power of 2!"))
    N
end
nspins(state::AbstractVector) = nspins(length(state))
nspins(op::AbstractMatrix) = size(op,1) == size(op, 2) ? nspins(size(op, 1)) : throw(ArgumentError("The operator is not a square-matrix!"))

single_spin_op(op, k::Integer, N::Integer) = identity_op(k-1) ⊗ op ⊗ identity_op(N-k)
correlator(op, i::Integer, j::Integer, N::Integer) = i > j ? correlator(op, j, i, N) : identity_op(i-1) ⊗ op ⊗ identity_op(j-i-1) ⊗ op ⊗ identity_op(N-j)

op_list(op, N::Integer) = [single_spin_op(op, k, N) for k in 1:N]
