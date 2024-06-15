struct AtomDataParser
    re::Regex
end
function AtomDataParser()
    atom_data_regex = r"^(?P<s>\w+)\s+(?P<f1>-?\d+\.\d+)\s+(?P<f2>-?\d+\.\d+)\s+(?P<f3>-?\d+\.\d+)\s+CORE\s+.*"
    return AtomDataParser(atom_data_regex)
end
struct BlockHeaderParser
    re1::Regex
    re2::Regex
end
function BlockHeaderParser()
    block_header_regex1 = r"^\s+Energy\s+(\d+)\s+(-?[0-9.]+)\s+(-?[0-9.]+)\s+(.*)$"
    block_header_regex2 = r"^\s+Energy\s+(\d+)\s+(-?[0-9.]+)\s+(-?[0-9.]+)"
    return BlockHeaderParser(block_header_regex1, block_header_regex2)
end
struct CellHeaderParser
    re::Regex
end
function CellHeaderParser()
    cell_header_regex = r"^\w+\s+(?P<x>\d+\.\d+)\s+(?P<y>\d+\.\d+)\s+(?P<z>\d+\.\d+)\s+(?P<alpha>\d+.\d+)\s+(?P<beta>\d+.\d+)\s+(?P<gamma>\d+.\d+)"
    return CellHeaderParser(cell_header_regex)
end
function arc_parse(line::String, parser::AtomDataParser)
    m = match(parser.re, line)
    if m === nothing
        return nothing
    end
    s::String = m.captures[1]
    f1 = parse(Float64, m.captures[2])
    f2 = parse(Float64, m.captures[3])
    f3 = parse(Float64, m.captures[4])
    return (s, f1, f2, f3)
end
function arc_parse(line::String, parser::BlockHeaderParser)
    m1 = match(parser.re1, line)
    m2 = match(parser.re2, line)
    if m1 === nothing
        if m2 === nothing
            return nothing
        end
        number = parse(Int64, m2.captures[1])
        float1 = parse(Float64, m2.captures[2])
        energy = parse(Float64, m2.captures[3])
        symmetry = "C1"
        return (number, float1, energy, symmetry)
    end
    n = parse(Int64, m1.captures[1])
    f1 = parse(Float64, m1.captures[2])
    f2 = parse(Float64, m1.captures[3])
    s = m1.captures[4]
    return (n, f1, f2, s)
end
function arc_parse(line::String, parser::CellHeaderParser)
    m = match(parser.re, line)
    if m === nothing
        return nothing
    end
    x = parse(Float64, m.captures[1])
    y = parse(Float64, m.captures[2])
    z = parse(Float64, m.captures[3])
    α = parse(Float64, m.captures[4])
    β = parse(Float64, m.captures[5])
    γ = parse(Float64, m.captures[6])
    cell_info = (x, y, z, α, β, γ)
    return cell_info
end
function add_atom(structure_block::StructureBlock, new_atom::Atom)
    return StructureBlock(structure_block.n_atoms + 1, structure_block.cell, structure_block.energy, push!(copy(structure_block.atoms), new_atom), structure_block.symmetry, structure_block.additional_information)
end
function change_cell(structure_block::StructureBlock, new_cell::Cell)
    return StructureBlock(structure_block.n_atoms, new_cell, structure_block.energy, structure_block.atoms, structure_block.symmetry, structure_block.additional_information)
end

function read_arc(filename::String)
    f = open(filename)
    blocks = ArcParser.StructureBlock{Float64}[]
    current_block = nothing
    atom_data_parser = AtomDataParser()
    block_header_parser = BlockHeaderParser()
    cell_header_parser = CellHeaderParser()
    cell_info = ArcParser.Cell{Float64}(0,0,0,0,0,0)
    while !eof(f)
        line = readline(f)
        atom_data_parse_result = arc_parse(line, atom_data_parser)
        if atom_data_parse_result !== nothing
            s, f1, f2, f3 = atom_data_parse_result
            new_point = ArcParser.Point(f1, f2, f3)
            new_atom = ArcParser.Atom{Float64}(s, 0, 0, 0, new_point)
            current_block = add_atom(current_block, new_atom)
            continue
        end
        header_parse_result = arc_parse(line, block_header_parser)
        if header_parse_result !== nothing
            n, f1, energy, s = header_parse_result
            if current_block !== nothing
                push!(blocks, current_block)
            end
            current_block = ArcParser.StructureBlock{Float64}(0, cell_info, energy, [], s, [])
            continue
        end
        cell_parse_result = arc_parse(line, cell_header_parser)
        if cell_parse_result !== nothing
            current_block = ArcParser.change_cell(current_block, ArcParser.Cell{Float64}(cell_parse_result...))
            continue
        end
    end
    if current_block !== nothing
        push!(blocks, current_block)
    end
    return blocks
end

function write_to_arc(structure::StructureBlock, filename::String)
    f = open(filename, "w")
    write(f, "!BIOSYM archive 2\n")
    write(f, "PBC=ON\n")
    println(f, lpad("Energy", 28), lpad(0, 10), lpad("0.0", 16), lpad(structure.energy, 18), lpad(structure.symmetry, 10))
    println(f, "!DATE")
    println(f, "PBC ", lpad(structure.cell.a, 14), lpad(structure.cell.b, 14), lpad(structure.cell.c, 14), lpad(structure.cell.α, 14), lpad(structure.cell.β, 14), lpad(structure.cell.γ, 14))
    for (i, atom) in enumerate(structure.atoms)
        println(f, lpad(atom.element, 5), lpad(atom.position[1], 15), lpad(atom.position[2], 15), lpad(atom.position[3], 15), " CORE ", lpad(i+1, 5), " ", lpad(atom.element, 3), lpad(atom.element, 5), lpad(0.0, 6), lpad(i+1, 5))
    end
    println(f, "end")
    println(f, "end")
    println(f, "\n")
    close(f)
end