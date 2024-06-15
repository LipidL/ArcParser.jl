using ArcParser
using Test
using Random

@testset "similarity" begin
    # Test similarity
    ref1 = ArcParser.read_arc("ref1.arc")[1]
    ref2 = ArcParser.read_arc("ref2.arc")[1]
    ref1_matrix = ArcParser.calculate_position_matrix(ref1)
    ref2_matrix = ArcParser.calculate_position_matrix(ref2)
    # Test for Kabsch similarity
    similarity = ArcParser.calculate_similarity(ref1_matrix, ref2_matrix, ArcParser.Kabsch())
    @test similarity isa Float64
    @test similarity >= 0
    @test similarity <= 1
    # Test for SimplifiedSOAP similarity
    similarity = ArcParser.calculate_similarity(ref1.atoms, ref2.atoms, ArcParser.SimplifiedSOAP(), 1.0)
    @test similarity isa Float64
    @test similarity >= 0
end

@testset "bond_matrix" begin
    # Test bond_matrix
    structures = ArcParser.read_arc("ref1.arc")
    gp1_matrix = ArcParser.calculate_position_matrix(structures[1])
    @test ArcParser.calculate_bond_matrix(gp1_matrix, 1.0) isa Matrix
    @test size(ArcParser.calculate_bond_matrix(gp1_matrix, 1.0)) == (length(structures[1].atoms), length(structures[1].atoms))
    @test ArcParser.calculate_bond_matrix(structures[1], 1.0) ≈ ArcParser.calculate_bond_matrix(gp1_matrix, 1.0)
    structure = [0.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0; 1.0 0.0 0.0]
    bond_matrix = ArcParser.calculate_bond_matrix(structure, 1.1)
    @test bond_matrix == [0.0 1.0 1.0 1.0; 1.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0]
end