export fieldterm, xxzmodel, hopping, z_field


_get_interaction(J::AbstractMatrix, i, j) = J[i,j]
_get_interaction(J::Number, i, j) = J

function fieldterm(N, op)
    res = op
    for i in 2:N
        res = 𝟙 ⊗ res + op ⊗ identity_op(i-1)
    end
    res
end

function hopping(N, i, j)
    entries = 2^(N-1)
    A,B,V = zeros(Int, entries), zeros(Int, entries), ones(Int, entries)
    _hopping_coordinates!(A, B, N, i, j)
    return sparse(A,B,V,2^N,2^N)
end

xxzmodel(J::AbstractMatrix, Δ) = xxzmodel(size(J, 1), J, Δ)

function xxzmodel(N,J,Δ)
    # 2^N diagonal entries + N(N-1)/2 * 2^(N-1) offdiagonal terms
    entries = 2^N + 2^(N-2)*N*(N-1)
    A,B,V = zeros(Int, entries), zeros(Int, entries), zeros(float(promote_type(eltype(J), typeof(Δ))), entries)

    @views _zcorrelator_coo!(A[1:2^N],B[1:2^N],V[1:2^N], N, J)
    V[1:2^N] .*= Δ

    @views _hopping!(A[2^N+1:end], B[2^N+1:end], V[2^N+1:end], N, J)
    V[2^N+1:end] .*= 2

    V ./= 4 # pauli -> spin operators
    sparse(A,B,V,2^N,2^N)
    ## In principle using Hermitian and only computing only
    ## half of the off-diagonal speeds up constructing the Hamiltonian.
    ## In practice this slows down everythin else dramatically, as Hermitian
    ## (sometimes) hides the fact that this is a sparse matrix and that leads to
    ## unfortunate interactions (f.e. Hermitian(sparse) ≈ dense takes AGES)
    #Hermitian(sparse(A,B,V,2^N,2^N), :L)
end

"""
    _zcorrelator_coo!(A, B, V, J)

Generate the COO vectors of ∑ᵢⱼJᵢⱼ σzⁱσzʲ.

# Parameters
 - `A`,`B`,`V` vectors for column (`A`) and row indices (`B`) respectively values (`V`)
 - `N` number of sites in the system
 - `J` interaction matrix of size `N`x`N` or scalar. Assumed symmetric!

`A`,`B` and `V` need to have length 2^N.
"""
function _zcorrelator_coo!(A,B,V, N, J)
    A[:] .= 1:2^N
    B[:] .= 1:2^N
    V[:] .= 0
    for i in 1:N
        for j in i+1:N
            _zcorrelator_values!(V, N, _get_interaction(J, j,i), i, j)
        end
    end
end

"""
    _zcorrelator_values!(V, N, Jij, i, j)

*ADD* to `V` the diagonal of Jᵢⱼ σzⁱσzʲ.

# Parameters
 - `V` array of length 2^`N`. Will be mutated!
 - `N` number of sites in the system
 - `Jij` interaction strength
 - `i`,`j` sites to couple
"""
function _zcorrelator_values!(V, N, Jij, i, j)
    # 𝟙(sec1) ⊗ σz ⊗ 𝟙(sec2) ⊗ σz ⊗ 𝟙(sec3)

    #             | 𝟙(sec2) |     0    |   | 𝟙(sec3) |    0     |
    # = 𝟙(sec1) ⊗ | ------- | -------- | ⊗ | ------- | -------- |
    #             |    0    | -𝟙(sec2) |   |    0    | -𝟙(sec3) |


    #             |           [ 𝟙(sec3) |     0    ] |                                  |
    #             | 𝟙(sec2) ⊗ [ ------- | -------- ] |              0                   |
    #             |           [    0    | -𝟙(sec3) ] |                                  |
    # = 𝟙(sec1) ⊗ | -------------------------------- | -------------------------------- |
    #             |                                  |           [ -𝟙(sec3) |    0    ] |
    #             |              0                   | 𝟙(sec2) ⊗ [ -------- | ------- ] |
    #             |                                  |           [    0     | 𝟙(sec3) ] |
    sec1 = 2^(i-1)
    sec2 = 2^(j-i-1)
    sec3 = 2^(N-j)

    blocksize3 = 2sec3 # σz ⊗ 𝟙(sec3)
    blocksize23 = 2sec2*blocksize3 # size of: σz ⊗ 𝟙(sec2) ⊗ σz ⊗ 𝟙(sec3)

    #sec1 loop
    for offset1 in blocksize23 .* (0:sec1-1)
        # upper block
        for offset2 in offset1 .+ blocksize3 .* (0:sec2-1)
            V[1 .+ offset2 .+ (0:sec3-1)] .+= Jij
            V[1 .+ offset2 .+ (sec3:2sec3-1)] .+= -Jij
        end
        # lower block
        for offset2 in offset1 .+ blocksize3 .* (sec2:2sec2-1)
            V[1 .+ offset2 .+ (0:sec3-1)] .+= -Jij
            V[1 .+ offset2 .+ (sec3:2sec3-1)] .+= Jij
        end
    end
end

"""
    _hopping!(A,B,V, J)

Generate the COO vectors for ∑ᵢⱼ Jᵢⱼ (σ₊ⁱσ₋ʲ + σ₋ⁱσ₊ʲ)

# Parameters
 - `A`,`B`,`V` vectors for column (`A`) and row indices (`B`) respectively values (`V`)
 - `J` interaction matrix of size `N`x`N`. Assumed symmetric!

 If `J` is `N`x`N`, then `A`,`B` and `V` need to have length (N-1)N*2^(N-2)
"""
function _hopping!(A,B,V, N, J)
    fill!(V, 1)
    length = 2^(N-1)
    at = 0
    for i in 1:N
        for j in i+1:N
            range = at .+ (1:length)
            @views _hopping_coordinates!(A[range], B[range], N, i, j)
            V[range] .*= _get_interaction(J, j, i)
            at += length
        end
    end
    ## _hopping values! does the transpose
    ## no measureable performance difference
    #copyto!(A,at+1,B,1,at)
    #copyto!(B,at+1,A,1,at)
    #copyto!(V,at+1,V,1,at)
end

"""
    _hopping_coordinates!(A, B, N, i, j)

Compute the coordinates of the matrix entries of the hopping operator between site i and J
in a spin chain. The values are 1 everywhere, so there is no input/output for these in this
function.

Hopping operator: σ₊ⁱσ₋ʲ + σ₋ⁱσ₊ʲ = 𝟙⊗…⊗𝟙⊗σ₊⊗𝟙⊗…⊗𝟙⊗σ₋⊗𝟙⊗…⊗𝟙 + transpose

# Parameters
 - `A`,`B` coordinate vectors. These should have length 2^(N-2) will be overwritten!
 - `N` number of sites in the system
 - `i`,`j` sites the excitation can hop between. 1 ≤ i,j ≤ N
"""
function _hopping_coordinates!(A,B,N,i,j)
    i > j && return _hopping_coordinates!(A,B,N,j,i)
    # 𝟙(sec1) ⊗ σ+ ⊗ 𝟙(sec2) ⊗ σ- ⊗ 𝟙(sec3) + transpose
    #             |                                  |                                  |
    #             |                0                 |                 0                |
    #             |                                  |                                  |
    # = 𝟙(sec1) ⊗ | -------------------------------- | -------------------------------- |
    #             |           [    0    |  𝟙(sec3) ] |                                  |
    #             | 𝟙(sec2) ⊗ [ ------- | -------- ] |                 0                |
    #             |           [    0    |    0     ] |                                  |
    sec1 = 2^(i-1)
    sec2 = 2^(j-i-1)
    sec3 = 2^(N-j)

    blocksize3 = 2sec3 # σz ⊗ 𝟙(sec3)
    blocksize23 = 2sec2*blocksize3 # size of: σz ⊗ 𝟙(sec2) ⊗ σz ⊗ 𝟙(sec3)

    currentIndex = 0
    #sec1 loop
    for offset1 in blocksize23 .* (0:sec1-1)
        for offset2 in blocksize3 .* (0:sec2-1)
            rowOffset2 = offset1 + offset2 # left columns # offset for the row-COORDINATE
            colOffset2 = offset1 + offset2 + blocksize3*sec2 # bottom rows # offset for the col-COORDINATE

            A[currentIndex .+ (1:sec3)] .= colOffset2 .+ (1:sec3) # upper rows
            B[currentIndex .+ (1:sec3)] .= rowOffset2 .+ (sec3+1:2sec3) # right columns
            currentIndex += sec3

            ## now transpose
            B[currentIndex .+ (1:sec3)] .= colOffset2 .+ (1:sec3) # upper rows
            A[currentIndex .+ (1:sec3)] .= rowOffset2 .+ (sec3+1:2sec3) # right columns
            currentIndex += sec3
        end
    end
end


z_field(N::Int) = z_field(ones(N))

function z_field(field_values::Vector)
    N = length(field_values)
    ret = zeros(2^N)
    for (i,hz) in enumerate(field_values)
        _z_field_values!(ret, N, hz, i)
    end
    return Diagonal(ret)
end

"""
    _z_field!(V, N, hz, i)

*ADD* to `V` the diagonal of hz σzⁱ.

# Parameters
 - `V` array of length 2^`N`. Will be mutated!
 - `N` number of sites in the system
 - `hz` field strength
 - `i` site to act on
"""
function _z_field_values!(V, N, hz, i)
    # σz^(i) = 𝟙(sec1) ⊗ σz ⊗  𝟙(sec2)
    sec1 = 2^(i-1)
    sec2 = 2^(N-i)

    blocksize2 = 2sec2 # σz ⊗ 𝟙(sec2)

    for offset1 in blocksize2 .* (0:sec1-1)
        V[offset1 .+ (1:sec2)] .+= hz
        V[offset1 .+ (sec2+1:2sec2)] .-= hz
    end
end
