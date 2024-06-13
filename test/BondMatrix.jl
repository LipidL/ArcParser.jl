using ArcParser
using Test

@testset "BondMatrix.jl" begin
    # Test calculate_bond_matrix
    structure = [0.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0; 1.0 0.0 0.0]
    bond_matrix = ArcParser.calculate_bond_matrix(structure, 1.1)
    @test bond_matrix == [0.0 1.0 1.0 1.0; 1.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0]
end