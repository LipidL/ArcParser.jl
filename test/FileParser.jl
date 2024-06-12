using ArcParser
using Test

@testset "ArcParser.jl" begin
    # Test AtomDataParser
    atom_data_parser = ArcParser.AtomDataParser()
    @test ArcParser.arc_parse("H 0.0 0.0 0.0 CORE 1", atom_data_parser) == ("H", 0.0, 0.0, 0.0)
    block_header_parser = ArcParser.BlockHeaderParser()
    @test ArcParser.arc_parse("  Energy 1 -0.0 -0.0", block_header_parser) == (1, -0.0, -0.0, "C1")
    cell_header_parser = ArcParser.CellHeaderParser()
    @test ArcParser.arc_parse("CRYSTAL 10.0 10.0 10.0 90.0 90.0 90.0", cell_header_parser) == (10.0, 10.0, 10.0, 90.0, 90.0, 90.0)
end