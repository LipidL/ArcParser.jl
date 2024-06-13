using ArcParser
using LinearAlgebra
function calculate_similarity(structure1::Vector{Atom{T}}, structure2::Vector{Atom{T}}, θ::T) where T
    n_atoms = length(structure1)
    if n_atoms != length(structure2)
        return -Inf64
    end
    similarity = 0.0
    for atomi in structure1
        for atomj in structure2
            r_i = ArcParser.extract_position(atomi)
            r_j = ArcParser.extract_position(atomj)
            tmp = -LinearAlgebra.norm(r_i - r_j)/(4*θ^2)
            similarity += exp(tmp)
        end 
    end
    similarity /= n_atoms
    return similarity
end