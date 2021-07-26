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
nspins(state::AbstractVector) = Int(ceil(log2(length(state))))
nspins(op::AbstractMatrix) = Int(ceil(log2(size(op, 1))))

single_spin_op(op, k::Integer, N::Integer) = identity_op(k-1) ⊗ op ⊗ identity_op(N-k)
correlator(op, i::Integer, j::Integer, N::Integer) = i > j ? correlator(op, j, i, N) : identity_op(i-1) ⊗ op ⊗ identity_op(j-i-1) ⊗ op ⊗ identity_op(N-j)

op_list(op, N::Integer) = [single_spin_op(op, k, N) for k in 1:N]

symmetrize_state(state::AbstractVector, sign=1) = (state[1:div(length(state),2)] .+ sign*state[length(state):-1:div(length(state),2)+1])/√2
symmetrize_state(state::AbstractArray, sign=1) = mapslices(s -> symmetrize_state(s, sign), state; dims=1)

function symmetrize_op(op, sign=1)
    d = size(op, 1)
    l = div(d, 2)
    front = 1:l
    back = d:-1:l+1
    
    0.5 * (op[front, front] + op[back,back] + sign*(op[front,back]+op[back,front]))
end