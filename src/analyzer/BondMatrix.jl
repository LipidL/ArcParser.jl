using ArcParser
using LinearAlgebra



function calculate_bond_matrix(structure::Matrix{T}, shreshold::T) where T <: Real
    if size(structure, 2) != 3
        error("Input matrix should be a n*3 matrix")
    end
    n_atoms = size(structure, 1)
    bond_matrix = zeros(T, n_atoms, n_atoms)
    for i in range(1, n_atoms)
        for j in range(i+1, n_atoms)
            r_i = structure[i, :]
            r_j = structure[j, :]
            r_ij = r_i - r_j
            r_ij_norm = LinearAlgebra.norm(r_ij)
            if r_ij_norm < shreshold
                bond_matrix[i, j] = 1
                bond_matrix[j, i] = 1
            end
        end
    end
    return bond_matrix
end

function calculate_bond_matrix(structure::StructureBlock{T}, shreshold::T) where T <: Real
    return calculate_bond_matrix(calculate_position_matrix(structure), shreshold)   
end

function calculate_position_matrix(structure::StructureBlock{T}) where T <: Real
    transposed_mat = reduce(hcat, ([ArcParser.extract_position(atom) for atom in structure.atoms]))
    retmat = zeros(T, size(transposed_mat, 2), size(transposed_mat, 1))
    adjoint!(retmat, transposed_mat)
    return retmat
end