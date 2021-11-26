@testset "hamiltonian.jl" begin

    ## Old implementation
    using SparseArrays
    function xxzmodel_old(N, J, Δ)
        Sx = complex(σx/2); Sy = σy./2; Sz = complex(σz/2);
        res = spzeros(complex(eltype(J)), 2^N, 2^N)
        for i in 1:N
            for j in i+1:N
                res += J[i, j] * (correlator(Sx, i, j, N) + correlator(Sy, i, j, N) + Δ*correlator(Sz, i, j, N))
            end
        end
        res
    end


    @testset "fieldterm" begin
        @test fieldterm(2, σx) ≈ σx ⊗ 𝟙 + 𝟙 ⊗ σx
        @test fieldterm(5, 𝟙) ≈ 5*identity_op(5)
    end

    @testset "hopping operator" begin
        @test hopping(2,1,2) ≈ [0 0 0 0;
                                0 0 1 0;
                                0 1 0 0;
                                0 0 0 0]
        @test hopping(5,4,1) == hopping(5,1,4)
    end

    @testset "xxzmodel" begin
        @test xxzmodel(2, 1, -1)*4  ≈ [-1  0  0  0;
                                        0  1  2  0;
                                        0  2  1  0;
                                        0  0  0 -1]
        @test xxzmodel([0 1 1; 1 0 1; 1 1 0], -1) ≈ xxzmodel(3, 1, -1)

        # comparison to old implementation
        for N in 2:14
            J = Symmetric(rand(N,N))
            @test xxzmodel(J,-0.7) ≈ xxzmodel_old(N, J, -0.7)
        end
    end
end
