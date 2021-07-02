using Base: disable_library_threading, power_by_squaring
using XXZNumerics
using Test, Statistics, Random, LinearAlgebra

@testset "XXZNumerics.jl" begin
    # Write your tests here.
    include("geometry.jl")

    include("interaction.jl")
end
