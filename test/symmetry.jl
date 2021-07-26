@testset "symmetry.jl" begin
    
    import XXZNumerics.Symmetry
    using LinearAlgebra
    using Random

    Random.seed!(5)

    @testset "ZBasis" begin
        @testset "FullZBasis" begin
            full = Symmetry.zbasis(5)
            @test full isa Symmetry.FullZBasis
            @test Symmetry.basis_size(full) == 2^5
            @test collect(full) == [(1, 1:2^5)]
        end

        @testset "ZBlockBasis" begin
            digitsum(k) = sum(parse.(Int, [string(k-1; base=2)...]))
            for k in 0:5
                zblock = Symmetry.zbasis(5, k)
                @test zblock isa Symmetry.ZBlockBasis
                @test Symmetry.basis_size(zblock) == binomial(5, k)
                coeff, inds = first(zblock) 
                @test iterate(zblock, iterate(zblock)[2]) === nothing
                @test length(zblock) == 1
                @test all(digitsum.(inds) .== k)
            end
            @test_throws ArgumentError Symmetry.zbasis(5, 6)
            @test_throws ArgumentError Symmetry.zbasis(5, -1)
        end
        
        # TODO what state to test? Test with eigenstates of z-mag operator?
    end

    @testset "Symmetries" begin
        @testset "SpinFlip" begin
            sf = Symmetry.SpinFlip(5) # spinflip on full hilbertspace
            @test Symmetry.basis_size(sf) == 2^4
            @test length(sf) == 2
            @test first(sf) == (1/√2, 1:2^4)
            @test sf == Symmetry.SpinFlip(Symmetry.zbasis(5))
            @test Symmetry.SpinFlip(Symmetry.zbasis(5), 1) == Symmetry.SpinFlip(Symmetry.zbasis(5), 10)
            
            ## some more consistency checks
            consistent(sf, expected) = Symmetry.basis_size(sf) == length(first(sf)[2]) == expected
            @test consistent(Symmetry.SpinFlip(Symmetry.ZBlockBasis(5,2)), binomial(5,2))
            @test consistent(Symmetry.SpinFlip(Symmetry.ZBlockBasis(5,3)), binomial(5,3))
            @test consistent(Symmetry.SpinFlip(Symmetry.ZBlockBasis(6,2)), binomial(6,2))
            @test consistent(Symmetry.SpinFlip(Symmetry.ZBlockBasis(6,3)), binomial(6,3)/2)

            N = 5
            for basis in (Symmetry.FullZBasis(N), [Symmetry.ZBlockBasis(N,k) for k in 0:N]...)
                a,b = collect(Symmetry.SpinFlip(basis,  1))
                c,d = collect(Symmetry.SpinFlip(basis, -1))
                @test a == c && b[1] == -d[1] && b[2] == d[2] && a[1]^2+c[1]^2 ≈ 1
            end
        end

        @testset "Shift" begin
        
        end

        @testset "SpinFlip+Shift" begin
        
        end
    end

    @testset "symmetrize_op" begin
        N = 5
        opz = fieldterm(N, σz)
        opx = fieldterm(N, σx)

        @testset "ZBasis" begin
            @test opz == Symmetry.symmetrize_op(Symmetry.FullZBasis(N), opz) #identity
            @test opx == Symmetry.symmetrize_op(Symmetry.FullZBasis(N), opx) #identity
        end
        @testset "ZBlockBasis" for k in 0:N
            symm_op = Symmetry.symmetrize_op(Symmetry.ZBlockBasis(N, k), opz)
            @test isdiag(symm_op) && all(diag(symm_op) .≈  N-2k)
            @test all(Symmetry.symmetrize_op(Symmetry.ZBlockBasis(N, k), opx) .== 0)
        end

        @testset "SpinFlip" begin
            @test all(Symmetry.symmetrize_op(Symmetry.SpinFlip(N,  1), opz) .== 0)
            @test all(Symmetry.symmetrize_op(Symmetry.SpinFlip(N, -1), opz) .== 0)
            @test Symmetry.symmetrize_op(Symmetry.SpinFlip(N, 1), opx) + Symmetry.symmetrize_op(Symmetry.SpinFlip(N, -1), opx) ≈ 2opx[1:2^(N-1),1:2^(N-1)]
            @test Symmetry.symmetrize_op(Symmetry.SpinFlip(N, 1), opx) - Symmetry.symmetrize_op(Symmetry.SpinFlip(N, -1), opx) ≈ 2opx[2^N:-1:2^(N-1)+1,1:2^(N-1)]

            ## consistency with old methods
            H = xxzmodel(N,1,-0.7)
            @test symmetrize_op(H,  1) ≈ Symmetry.symmetrize_op(Symmetry.SpinFlip(N,  1), H)
            @test symmetrize_op(H, -1) ≈ Symmetry.symmetrize_op(Symmetry.SpinFlip(N, -1), H)
        end
    end

    @testset "symmetrize_state" begin
        N = 5
        opz = fieldterm(N, σz)
        shots = 10

        @testset "ZBasis"  begin
            for _ in 1:shots
                vec = normalize!(rand(2^N))
                @test vec == Symmetry.symmetrize_state(Symmetry.FullZBasis(N), vec)

                block_norm = 0.0
                block_mag = 0.0
                for k in 0:N
                    symm_vec = Symmetry.symmetrize_state(Symmetry.ZBlockBasis(N, k), vec)
                    block_norm += norm(symm_vec)^2
                    block_mag += dot(symm_vec, Symmetry.symmetrize_op(Symmetry.ZBlockBasis(N, k), opz), symm_vec)
                end
                @show block_norm, block_mag, dot(vec, opz, vec)
                @test block_norm ≈ 1 && block_mag ≈ dot(vec, opz, vec)
            end
        end
        ## can test whether symmetrizations added back up yield same vector
        ## test for equivalent vectors whether symmetrization is equal
        ## obv check lengths and stuff
    end

end