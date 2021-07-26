@testset "Interaction" begin

    @testset "PowerLaw" begin
        p1 = PowerLaw(1)
        p2 = PowerLaw(2)

        @test p2 == PowerLaw(2.0) == PowerLaw(2//1)

        geom = Box([2,2])
        points = sample_blockaded(geom, 2)
        distances = distance_matrix(geom, points)

        expected = (1 ./ distances) 
        expected[1,1] = expected[2,2] = 0
        @test interaction_matrix(p1, distances) ≈ expected
        @test interaction_matrix(p1, geom, points) ≈ interaction_matrix(p1, distances)
        @test isapprox(interaction_matrix(p2, geom, points), interaction_matrix(p1, distances) .^ 2; atol=1e-10) # need a bit lower precision
    end

end