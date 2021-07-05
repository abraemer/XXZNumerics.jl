module XXZNumerics

using LinearAlgebra
using SparseArrays

# Write your package code here.


include("geometry.jl")
include("interaction.jl")
include("spin_half.jl")
include("hamiltonian.jl")

end
