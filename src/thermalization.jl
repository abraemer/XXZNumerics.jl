export findβ, microcanonical_ensemble, canonical_ensemble

# stupid implementation of softmax - vals should be regularized beforehand
function _softmax(vals)
    w = exp.(vals)
    w ./= sum(w)
end

"""
    _weights_canonical(β, evals)
Compute the distribution of eigenstates with given energies at temperature 1/β.
P(E) ∝ ∑ᵢexp(-β Eᵢ)

Energy values are assumed to be sorted!
"""
function _weights_canonical(evals, β)
    # assume sorted evals!
    regulator = β < 0 ? evals[end] : evals[1]
    _softmax(@. -β * (evals - regulator))
end

_ΔE(evals, E_0, β) = abs((evals ⋅ _weights_canonical(evals, β)) - E_0)

# TODO: How to optimize best? The function is pretty well behaved...
"""
    findβ(evals, E₀; β_0=0)
Compute the inverse temperature β s. t. the thermal expectation value of the energy matches E₀.
"""
findβ(evals, E_0; β_0 = 0) = Optim.optimize(β -> _ΔE(evals, E_0, β[1]), [float(β_0)], Optim.Newton(linesearch = LineSearches.HagerZhang())).minimizer[1]

"""
    canonical_ensemble(evals, E₀; β₀ = 0)
Compute the weights of the state assuming a canonical ensemble of temperature 1/β. β₀ is taken as starting value for the computation.
A state is assigned the weight ∝ exp(-β Eᵢ)

Energy values are assumed to be sorted!
"""
canonical_ensemble(evals, E_0; β_0 = 0) = _weights_canonical(evals, findβ(evals, E_0; β_0))


"""
    microcanonical_ensemble(evals, E₀[, ΔE])
Compute the weights of the state assuming a microcanonical ensemble. Default size of the window is 0.5% of the spectral width.
A state is in the ensemble if |Eᵢ-E_0| < ΔE.
If the ensemble would be empty, choose the closest state (in energy).

Energy values are assumed to be sorted!
"""
function microcanonical_ensemble(evals, E_0, ΔE)::Vector{Float64}
    micro = @. (evals < E_0+ΔE) & (evals > E_0-ΔE)
    s = sum(micro)
    if s > 0
        # enforce float conversion for type stability
        micro ./ float(s)
    else
        micro = zeros(Float64, length(micro))
        smaller = searchsortedlast(evals, E_0)
        ind = abs(evals[smaller]-E_0) < abs(evals[smaller+1]-E_0) ? smaller : smaller+1
        micro[ind] = 1
        micro
    end
end

microcanonical_ensemble(evals, E_0) = microcanonical_ensemble(evals, E_0, 0.005*(evals[end]-evals[1]))