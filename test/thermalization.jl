@testset "thermalization.jl" begin
    
    @testset "_weights_canonical" begin
        # check weighting
        range = collect(1:15)
        @test XXZNumerics._weights_canonical(log.(range), -1) ≈ range./sum(range)
        @test XXZNumerics._weights_canonical(log.(range), -2) ≈ range.^2 ./ sum(range .^ 2)

        # check extreme cases
        @test .!(any(isnan, XXZNumerics._weights_canonical(exp.(range), 100)))
        @test XXZNumerics._weights_canonical(exp.(range), 0) ≈ ones(length(range))/length(range)
        @test XXZNumerics._weights_canonical(exp.(range), 1) ≈ XXZNumerics._weights_canonical(sort!(-exp.(range)), -1)[end:-1:1]
    end

    @testset "findβ" begin
        @test findβ([1,2,3,4], 2.5) ≈ 0
        @test findβ([1,2,3,4], 3) < 0
        @test findβ([1,2,3,4], 2) > 0
        @test findβ([1,2,3,4], 2.5-0.7) ≈ -findβ([1,2,3,4], 2.5+0.7) # should be symmetric around the average energy as spectrum is symmetric
        # can't think of other good tests...
    end

    @testset "canonical_ensemble" begin
        l = 10
        states = sort!(rand(l)*10 .- 5)
        @test canonical_ensemble(states, states[1]) ≈ vcat(1, zeros(l-1))
        @test canonical_ensemble(states, states[end]) ≈ vcat(zeros(l-1), 1)
        @test canonical_ensemble(states, mean(states)) ≈ ones(l)/l
    end

    @testset "microcanonical_ensemble" begin
        @test microcanonical_ensemble([-2, -1, 0.5, 1, 2], 0, 1.1) ≈ [0,1,1,1,0]/3
        @test microcanonical_ensemble([-2, -1, 0.5, 1, 2], 0, 0.9) ≈ [0,0,1,0,0]
        @test microcanonical_ensemble([-2, -1, 0.5, 1, 2], 0, 0) ≈ [0,0,1,0,0] # choose closest state
        @test microcanonical_ensemble([-2, -1, 0.5, 1, 2], -0.6, 0) ≈ [0,1,0,0,0]

        @test microcanonical_ensemble([-40, -1, 0.5, 1, 40], 0.7) ≈ [0,0,1,1,0]/2 # default arg -> ΔE = 80*0.5% = 0.4
    end
end