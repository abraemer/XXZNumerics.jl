@testset "spin_half.jl" begin
    # note 
    @test single_spin_op(σz, 1, 2) ≈ Diagonal([1,1,-1,-1])
    @test single_spin_op(σz, 2, 2) ≈ Diagonal([1,-1,1,-1])
    @test single_spin_op(σz, 2, 3) ≈ kron(Diagonal([1,-1,1,-1]), 𝟙)

    @test symmetrize_state(up ⊗ down) ≈ (up ⊗ down + down ⊗ up)[1:2]/√2
    @test symmetrize_state(up ⊗ down, -1) ≈ (up ⊗ down - down ⊗ up)[1:2]/√2
    @test symmetrize_state(up ⊗ down ⊗ down) ≈ (up ⊗ down ⊗ down + down ⊗ up ⊗ up)[1:4]/√2

    @test symmetrize_state(hcat(up ⊗ down, down ⊗ up)) ≈ hcat(symmetrize_state(up ⊗ down), symmetrize_state(down ⊗ up))

    @test symmetrize_op(σz ⊗ σz) ≈ σz
    @test symmetrize_op(σz ⊗ σz, -1) ≈ σz
    @test symmetrize_op(σx ⊗ σx) ≈ identity_op(1)
    @test symmetrize_op(σx ⊗ σx + σy ⊗ σy) ≈ Diagonal([0, 2]) # hopping operator

    @test correlator(σx, 1, 2, 2) ≈ σx ⊗ σx
    @test correlator(σz, 1, 2, 3) ≈ σz ⊗ σz ⊗ 𝟙

    @test op_list(σy, 2) ≈ [σy ⊗ 𝟙, 𝟙 ⊗ σy]
end