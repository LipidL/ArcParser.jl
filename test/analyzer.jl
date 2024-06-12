using ArcParser
using Test
using Random

@testset "analyzer.jl" begin
    # Test AtomDataParser
    println(pwd())
    structures = ArcParser.read_arc("test1.arc")
    atoms = structures[1].atoms
    for i in range(1, 1)
        # Randomly select 2 groups with 6 elements each, where the element field is "Fe"
        selected_groups = filter(x -> x.element == "Fe", atoms)
        selected_indices = rand(1:length(selected_groups), 12)
        selected_groups = selected_groups[selected_indices]
        group1 = selected_groups[1:6]
        group2 = selected_groups[7:12]
        @test ArcParser.calculate_similarity(group1, group2, 1e-3) isa AbstractFloat
    end
end