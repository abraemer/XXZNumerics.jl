module XXZNumerics

using LinearAlgebra
using SparseArrays
import Optim, LineSearches

# Write your package code here.

include("geometry.jl")
include("interaction.jl")
include("spin_half.jl")
include("hamiltonian.jl")
include("thermalization.jl")

end
