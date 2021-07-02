
@testset "Geometry" begin

    Random.seed!(5) # reproducability 

    @testset "generate_point" begin
        SAMPLES = 10000

        @testset "$geom" for geom in (Box, BoxPBC)
            lengths = [1,2,3]
            b = geom(lengths)
            testpoints = hcat([XXZNumerics.generate_point(b) for _ in 1:SAMPLES]...)
            
            for i in 1:length(lengths)
                @test reduce(min, testpoints[i,:]) > 0
                @test reduce(max, testpoints[i,:]) < lengths[i]
                @test abs(mean(testpoints[i,:]) - lengths[i]/2) < sqrt(SAMPLES)*lengths[i]
            end
        end

        
        @testset "$geom" for geom in (NoisyChain, NoisyChainPBC)
            L = 5
            # no noise
            chain = geom(L, 2, 0)
            @test sample_blockaded(chain, L) ≈ 2 .* collect(1:L)

            # some noise - but no overlaps
            chain = geom(L, 2, 1)
            testpoints = hcat([sample_blockaded(chain, L) for _ in 1:SAMPLES]...)
            
            for i in 1:L
                @test reduce(min, testpoints[i,:]) > 2i-1
                @test reduce(max, testpoints[i,:]) < 2i+1
                @test abs(mean(testpoints[i,:]) - 2i) < sqrt(SAMPLES)
            end
        end
    end

    @testset "distance" begin

        @testset "box" begin
            # note that box size does not matter for distances...
            box = Box([10,10,10,10])
            @test distance(box, [4,0,0,0]  , [1,0,0,0]) ≈ 3 # simple
            @test distance(box, [3,4,0,0]  , [0,0,0,0]) ≈ 5 
            @test distance(box, [1,2,3,4.0], [2,1,2,5]) ≈ 2 # more complicated (and float)
            @test distance(box, [2,1,2,5], [1,2,3,4.0]) ≈ 2 # check symmetry
        end

        @testset "boxpbc" begin
            boxpbc = BoxPBC([1,2,3])
            @test distance(boxpbc, [0.4, 0.4, 0.4], [0.5, 0.5, 0.5]) ≈ sqrt(3)*0.1  # normal distance
            @test distance(boxpbc, [0.1, 0  ,   0], [0.9, 0  ,   0]) ≈ 0.2          # wrap around
            @test distance(boxpbc, [0.4, 1.9,   0], [0.7, 0.3,   0]) ≈ 0.5          # mixed case
            @test distance(boxpbc, [1  ,   1,   2], [0,     0,   0]) ≈ sqrt(2)      # more complicated
            @test distance(boxpbc, [0,     0,   0], [1  ,   1,   2]) ≈ sqrt(2)      # more complicated (check symmetry)
        end

        @testset "NoisyChain" begin
            chain = NoisyChain(5,1,0)
            p1 = XXZNumerics.generate_point(chain, 1)
            p2 = XXZNumerics.generate_point(chain, 2)
            p3 = XXZNumerics.generate_point(chain, 3)
            p4 = XXZNumerics.generate_point(chain, 4)
            @test distance(chain, p1, p2) ≈ 1
            @test distance(chain, p2, p1) ≈ 1
            @test distance(chain, p1, p3) ≈ 2
            @test distance(chain, p1, p4) ≈ 3
        end

        @testset "NoisyChainPBC" begin
            chain = NoisyChainPBC(5,1,0)
            p1 = XXZNumerics.generate_point(chain, 1)
            p2 = XXZNumerics.generate_point(chain, 2)
            p3 = XXZNumerics.generate_point(chain, 3)
            p4 = XXZNumerics.generate_point(chain, 4)
            @test distance(chain, p1, p2) ≈ 1
            @test distance(chain, p2, p1) ≈ 1
            @test distance(chain, p1, p3) ≈ 2
            @test distance(chain, p1, p4) ≈ 2 # PBC!
            @test distance(chain, p4, p1) ≈ 2 # PBC! (check symmetry)
        end
    end

    @testset "blockade" begin

        @testset "Box" begin
            b = Box([2,])
            @test_throws ErrorException sample_blockaded(b, 3)
            points = sample_blockaded(b, 2)
            p1, p2 = points
            @test distance_matrix(b, points) == [0 distance(b, p1, p2); distance(b, p1, p2) 0]
            @test distance(b, p1, p2) > 1
        end

        @testset "BoxPBC" begin
            b = BoxPBC([2.2,])
            @test_throws ErrorException sample_blockaded(b, 3)
            points = sample_blockaded(b, 2)
            p1, p2 = points
            @test distance_matrix(b, points) == [0 distance(b, p1, p2); distance(b, p1, p2) 0]
            @test distance(b, p1, p2) > 1
        end
    end
    
end