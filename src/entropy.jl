export entropy, entropy!, entanglement_entropy

"""
    entropy(pvector)
    entropy(ρ)

Compute the shannon entropy: S(p) = -∑ᵢ pᵢ log pᵢ
Where pᵢ are the eigenvalues of ρ.
"""
function entropy(evals::AbstractVector)
	T = typeof(log(one(eltype(evals))))
	cutoff = convert(T, 1e-16)
	- sum(eval*log2(eval) for eval in Iterators.filter(>(cutoff), evals); init=0)
end

entropy(ρ::AbstractMatrix) = entropy(eigvals(ρ))

"""
    entropy!(ρ)

Destructive version of [`entropy`](@ref).
"""
entropy!(ρ::AbstractMatrix) = entropy(eigvals!(ρ))

"""
    entanglement_entropy(ψ[, N₁])

Compute the entanglement entropy of the left N₁ spins of the chain. Default value for N₁ is
half of the chain.
"""
function entanglement_entropy(state, N1=div(nspins(state), 2))
    return entropy(abs2.(svdvals(reshape(state, 2^N1, :))))
end
