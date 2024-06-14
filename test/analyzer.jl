using ArcParser
using Test
using Random

@testset "analyzer.jl" begin
    # Test similarity
    println(pwd())
    structures = ArcParser.read_arc("test1.arc")
    atoms = structures[1].atoms
    for i in range(1, 10)
        # Randomly select 2 groups with 6 elements each, where the element field is "Fe"
        selected_groups = filter(x -> x.element == "Fe", atoms)
        selected_indices = rand(1:length(selected_groups), 12)
        selected_groups = selected_groups[selected_indices]
        group1 = selected_groups[1:6]
        group2 = selected_groups[7:12]
        # Test SimplifiedSOAP
        @test ArcParser.calculate_similarity(group1, group2, ArcParser.SimplifiedSOAP(), 0.1) isa AbstractFloat
        @test ArcParser.calculate_similarity(group1, group2, ArcParser.SimplifiedSOAP(), 0.1) >= 0.0
        @test ArcParser.calculate_similarity(group1, group1, ArcParser.SimplifiedSOAP(), 0.1) >= 1.0
        @test ArcParser.calculate_similarity(group2, group2, ArcParser.SimplifiedSOAP(), 0.1) >= 1.0
        gp1_matrix = Matrix{Float64}(undef, length(group1),3)
        gp2_matrix = Matrix{Float64}(undef, length(group2),3)
        for i in range(1, length(group1))
            gp1_matrix[i, :] = ArcParser.extract_position(group1[i])
        end
        for i in range(1, length(group2))
            gp2_matrix[i, :] = ArcParser.extract_position(group2[i])
        end
        @test ArcParser.calculate_similarity(gp1_matrix, gp2_matrix, ArcParser.Kabsch()) isa AbstractFloat
        @test ArcParser.calculate_similarity(gp1_matrix, gp2_matrix, ArcParser.Kabsch()) >= 0.0
    end
    ref1 = ArcParser.read_arc("ref.arc")[1].atoms
    ref2 = ArcParser.read_arc("ref2.arc")[1].atoms
    ref1_matrix = Matrix{Float64}(undef, length(ref1),3)
    ref2_matrix = Matrix{Float64}(undef, length(ref2),3)
    for i in range(1, length(ref1))
        ref1_matrix[i, :] = ArcParser.extract_position(ref1[i])
    end
    for i in range(1, length(ref2))
        ref2_matrix[i, :] = ArcParser.extract_position(ref2[i])
    end
    similarity = ArcParser.calculate_similarity(ref1_matrix, ref2_matrix, ArcParser.Kabsch())
    println(similarity)
end