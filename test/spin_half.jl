@testset "spin_half.jl" begin
 
    @testset "single_spin_op" begin
        @test single_spin_op(σz, 1, 2) ≈ Diagonal([1,1,-1,-1])
        @test single_spin_op(σz, 2, 2) ≈ Diagonal([1,-1,1,-1])
        @test single_spin_op(σz, 2, 3) ≈ kron(Diagonal([1,-1,1,-1]), 𝟙)

        @test_throws DomainError single_spin_op(σz, 2, 1) # outside range
        @test_throws DomainError single_spin_op(σz, -1, 3) # bogus site - not in range [1, N]
        @test_throws DomainError single_spin_op(σz, 0, 3) # bogus site - not in range [1, N]
        @test_throws DomainError correlator(𝟙, 1, 2, -3) # negative length

        @test_throws MethodError single_spin_op(σz, 1.5, 3) # bogus site - not integer
        @test_throws MethodError correlator(𝟙, 1, 2, 2.5) # bogus length
    end

    @testset "symmetrize_state" begin
        # basics
        @test symmetrize_state(up ⊗ down) ≈ (up ⊗ down + down ⊗ up)[1:2]/√2
        @test symmetrize_state(up ⊗ down, -1) ≈ (up ⊗ down - down ⊗ up)[1:2]/√2
        @test symmetrize_state(up ⊗ down ⊗ down) ≈ (up ⊗ down ⊗ down + down ⊗ up ⊗ up)[1:4]/√2

        # multiple states
        @test symmetrize_state(hcat(up ⊗ down, down ⊗ up, up ⊗ up)) ≈ hcat(symmetrize_state(up ⊗ down), symmetrize_state(down ⊗ up), symmetrize_state(up ⊗ up))
    end

    @testset "symmetrize_op" begin
        @test symmetrize_op(σz ⊗ σz) ≈ σz
        @test symmetrize_op(σz ⊗ σz, -1) ≈ σz
        @test symmetrize_op(σx ⊗ σx) ≈ identity_op(1)
        @test symmetrize_op(σx ⊗ σx + σy ⊗ σy) ≈ Diagonal([0, 2]) # hopping operator
    end

    @testset "correlator" begin
        @test correlator(σx, 1, 2, 2) ≈ σx ⊗ σx
        @test correlator(σz, 1, 2, 3) ≈ σz ⊗ σz ⊗ 𝟙
        @test correlator(σz, 2, 1, 3) ≈ σz ⊗ σz ⊗ 𝟙 # symmetry

        @test_throws DomainError correlator(𝟙, 1, 1, 2) # operators on same site
        @test_throws DomainError correlator(𝟙, 1, 3, 2) # operator outside
        @test_throws DomainError correlator(𝟙, -1, 2, 3) # negative site
        @test_throws DomainError correlator(𝟙, 1, 2, -3) # negative length

        @test_throws MethodError correlator(𝟙, 1, 1.5, 3) # bogus site
        @test_throws MethodError correlator(𝟙, 1, 2, 2.5) # bogus length
    end

    @test op_list(σy, 2) ≈ [σy ⊗ 𝟙, 𝟙 ⊗ σy]
end