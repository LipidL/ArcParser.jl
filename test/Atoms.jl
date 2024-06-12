using ArcParser, Test

@testset "Atoms.jl" begin
    # Write your tests here.
    p1 = ArcParser.Point(1.0, 2.0)
    p2 = ArcParser.Point(3.0, 4.0)
    @test p1 + p2 ≈ ArcParser.Point(4.0, 6.0)
    p3 = ArcParser.Point(5.0, 7.0)
    @test p3 - p2 ≈ ArcParser.Point(2.0, 3.0)

    # Test unary minus
    @test -p1 ≈ ArcParser.Point(-1.0, -2.0)

    # Test multiplication by a scalar
    @test 2 * p1 ≈ ArcParser.Point(2.0, 4.0)
    @test p1 * 3 ≈ ArcParser.Point(3.0, 6.0)

    # Test division by a scalar
    @test p3 / 2 ≈ ArcParser.Point(2.5, 3.5)

    # Test indexing
    @test p1[1] == 1.0
    @test p2[2] == 4.0

    # Test approximate equality
    @test p1 ≈ ArcParser.Point(1.0, 2.0)

    # Test dist2 function
    @test ArcParser.dist2(3.0, 4.0) == 1.0
    @test ArcParser.dist2(p1, p2) == 8.0
    @test ArcParser.dist2(p1) == 5.0
end