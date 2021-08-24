module XXZNumerics

using LinearAlgebra
using SparseArrays
import Optim, LineSearches

# Write your package code here.

include("geometry.jl")
include("interaction.jl")
include("spin_half.jl")
include("symmetry/symmetry.jl")
include("hamiltonian.jl")
include("thermalization.jl")

using .Symmetry
export FullZBasis, ZBlockBasis, zbasis, SpinFlip, symmetrize_state, symmetrize_op, basis_size

end
