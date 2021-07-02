using Base: _sysimage_modules
using XXZNumerics
using Test, Statistics, Random, LinearAlgebra

@testset "XXZNumerics.jl" begin
    # Write your tests here.
    include("geometry.jl")

    include("interaction.jl")

    include("spin_half.jl")
end
