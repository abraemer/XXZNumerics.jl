module XXZNumerics

using LinearAlgebra: length, symmetric, AbstractMatrix
using Base: concurrency_violation
using LinearAlgebra
using SparseArrays

# Write your package code here.


include("geometry.jl")
include("interaction.jl")
include("spin_half.jl")

end
