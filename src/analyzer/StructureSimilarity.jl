using ArcParser
using LinearAlgebra
abstract type AbstractSimilarityMeasurement end
struct SimplifiedSOAP <: AbstractSimilarityMeasurement end
struct Kabsch <: AbstractSimilarityMeasurement end
function calculate_similarity(structure1::Vector{Atom{T}}, structure2::Vector{Atom{T}}, ::SimplifiedSOAP, θ::T) where T
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

function calculate_similarity(struct1::Matrix{T}, struct2::Matrix{T}, ::Kabsch ) where T
    if size(struct1, 2) != 3 || size(struct2, 2) != 3
        error("Input matrix should be a n*3 matrix")
    end
    if size(struct1, 1) != size(struct2, 1)
        error("Two structures should have the same number of atoms")
    end
    num_atoms = size(struct1, 1)
    center1 = sum(struct1, dims=1) / size(struct1, 1)
    center2 = sum(struct2, dims=1) / size(struct2, 1)
    moved1 = struct1 .- center1
    moved2 = struct2 .- center2
    H = moved1' * moved2
    U, S, V = svd(H)
    d = sign(det(U) * det(V))
    R = U * Diagonal([1, 1, d]) * V'
    moved1 = moved1 * R
    rmsd = LinearAlgebra.norm(moved1 - moved2) / sqrt(num_atoms)
    return rmsd
end