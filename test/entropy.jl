@testset "entropy.jl" begin

    Random.seed!(6)

    @testset "entropy" begin

        @test entropy([0.5, 0.5]) ≈ 1
        @test entropy([0.25, 0.25, 0.25, 0.25]) ≈ 2

        pvec = rand(10)
        pvec ./= sum(pvec)
        ρ = Diagonal(pvec)

        @test entropy(ρ) ≈ entropy(Matrix(ρ))
        @test entropy(Hermitian(Matrix(ρ))) ≈ entropy!(Matrix(ρ))
    end


    @testset "entanglement_entropy" begin
        bell_state = normalize(up⊗down + down⊗up)
        @test entanglement_entropy(up⊗up) ≈ 0
        @test entanglement_entropy(bell_state) ≈ 1
        @test entanglement_entropy(bell_state ⊗ down) ≈ 0
        @test entanglement_entropy(bell_state ⊗ down, 2) ≈ 1
        @test entanglement_entropy(bell_state ⊗ bell_state, 1) ≈ 1
        @test isapprox(entanglement_entropy(bell_state ⊗ bell_state, 2), 0; atol=1e-15)
        @test entanglement_entropy(bell_state ⊗ bell_state, 3) ≈ 1
    end
end
