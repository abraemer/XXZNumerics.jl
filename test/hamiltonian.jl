@testset "hamiltonian.jl" begin

    @testset "fieldterm" begin
        @test fieldterm(2, σx) ≈ σx ⊗ 𝟙 + 𝟙 ⊗ σx
        @test fieldterm(5, 𝟙) ≈ 5*identity_op(5)
    end

    @testset "xxzmodel" begin
        @test xxzmodel(2, 1, -1)*4  ≈ [-1  0  0  0;
                                        0  1  2  0;
                                        0  2  1  0;
                                        0  0  0 -1]
        @test xxzmodel([0 1 1; 1 0 1; 1 1 0], -1) ≈ xxzmodel(3, 1, -1)
    end
end