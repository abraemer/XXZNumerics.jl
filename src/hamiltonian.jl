export fieldterm, xxzmodel

function xxzmodel(N, J, Δ)
    Sx = complex(σx/2); Sy = σy./2; Sz = complex(σz/2);
    res = spzeros(complex(eltype(J)), 2^N, 2^N)
    for i in 1:N
        for j in i+1:N
            res += _get_interaction(J, i, j) * (correlator(Sx, i, j, N) + correlator(Sy, i, j, N) + Δ*correlator(Sz, i, j, N))
        end
    end
    res
end

xxzmodel(J::AbstractMatrix, Δ) = xxzmodel(size(J, 1), J, Δ)

_get_interaction(J::AbstractMatrix, i, j) = J[i,j]
_get_interaction(J::Number, i, j) = J

function fieldterm(N, op)
    res = op
    for i in 2:N
        res = 𝟙 ⊗ res + op ⊗ identity_op(i-1)
    end
    res
end