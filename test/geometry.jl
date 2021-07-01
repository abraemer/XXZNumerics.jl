
@testset "Geometry" begin

    Random.seed!(5)
    SAMPLES = 10000
    
    @testset "Box" begin
        lengths = [1,2,3]
        b = Box(lengths)
        testpoints = hcat([XXZNumerics.generate_point(b) for _ in 1:SAMPLES]...)
        
        for i in 1:length(lengths)
            @test reduce(min, testpoints[i,:]) > 0
            @test reduce(max, testpoints[i,:]) < lengths[i]
            @test abs(mean(testpoints[i,:]) - lengths[i]/2) < sqrt(SAMPLES)*lengths[i]
        end

        # note that box size does not matter for distances...
        @test distance(b, [1,2,3,4.0], [2,1,2,5]) ≈ 2
        @test distance(b, [3,4], [0,0]) ≈ 5
    end

    @testset "BoxPBC" begin
        lengths = [1,2,3]
        b = BoxPBC(lengths)
        testpoints = hcat([XXZNumerics.generate_point(b) for _ in 1:SAMPLES]...)
        
        for i in 1:length(lengths)
            @test reduce(min, testpoints[i,:]) > 0
            @test reduce(max, testpoints[i,:]) < lengths[i]
            @test abs(mean(testpoints[i,:]) - lengths[i]/2) < sqrt(SAMPLES)*lengths[i]
        end

        @test distance(b, [0.1,0,0], [0.9,0,0]) ≈ 0.2
        @test distance(b, [0.4,1.9,0], [0.7,0.3,0]) ≈ 0.5
        @test distance(b, [1,1,2], [0,0,0]) ≈ sqrt(2)
    end

    @testset "NoisyChain" begin
        L = 5
        chain = NoisyChain(L, 2, 0)
        @test sample_blockaded(chain, L) ≈ 2 .* collect(1:L)

        chain = NoisyChain(L, 2, 1)
        testpoints = hcat([sample_blockaded(chain, L) for _ in 1:SAMPLES]...)
        
        for i in 1:L
            @test reduce(min, testpoints[i,:]) > 2i-1
            @test reduce(max, testpoints[i,:]) < 2i+1
            @test abs(mean(testpoints[i,:]) - 2i) < sqrt(SAMPLES)
        end
    end

    @testset "blockade" begin
        b = BoxPBC([2.2,])
        @test_throws ErrorException sample_blockaded(b, 3)
        points = sample_blockaded(b, 2)
        p1, p2 = points
        @test all(distance_matrix(b, points) .== [0 distance(b, p1, p2); distance(b, p1, p2) 0])
        @test distance(b, p1, p2) > 1
    end
end