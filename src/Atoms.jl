struct Point{N,T}
    coo::NTuple{N,T}
end
Point(arg::T, args::T...) where T<:Number = Point((arg, args...))
Base.:(+)(x::Point{N,T}, y::Point{N,T}) where {N, T} = Point(x.coo .+ y.coo)
Base.:(-)(x::Point{N,T}, y::Point{N,T}) where {N, T} = Point(x.coo .- y.coo)
Base.:(-)(x::Point{N,T}) where {N, T} = Point(Base.:(-).(x.coo))
Base.adjoint(x::Point) = x
Base.:(*)(x::Number, y::Point) = Point(y.coo .* x)
Base.:(*)(y::Point, x::Number) = Point(y.coo .* x)
Base.:(/)(y::Point, x::Number) = Point(y.coo ./ x)
Base.iterate(x::Point, args...) = Base.iterate(x.coo, args...)
Base.getindex(x::Point, i::Int) = x.coo[i]
Base.isapprox(x::Point, y::Point; kwargs...) = all(isapprox.(x.coo, y.coo; kwargs...))

dist2(x::Number, y::Number) = abs2(x - y)
dist2(x::Point, y::Point) = sum(abs2, x - y)
dist2(x::Point) = sum(abs2, x)

struct Atom{T}
    element::String
    type::Int
    mass::T
    charge::T
    position::Point{3,T}
end
Atom(;type::Int = 1, mass::T = 1.0, charge::T = 0.0) where T<:Number = Atom{T}(type, mass, charge)
Base.isapprox(x::Atom, y::Atom; kwargs...) = x.element == y.element && isapprox(x.position, y.position; kwargs...)
function extract_position(atom::Atom{T})::Vector{T} where T
    return collect(atom.position.coo)
end
struct Cell{T}
    a::T
    b::T
    c::T
    α::T
    β::T
    γ::T
end
Base.isapprox(x::Cell, y::Cell; kwargs...) = isapprox(x.a, y.a; kwargs...) && isapprox(x.b, y.b; kwargs...) && isapprox(x.c, y.c; kwargs...) && isapprox(x.α, y.α; kwargs...) && isapprox(x.β, y.β; kwargs...) && isapprox(x.γ, y.γ; kwargs...)
struct StructureBlock{T}
    n_atoms::Int64
    cell::Cell{T}
    energy::T
    atoms::Vector{Atom{T}}
    symmetry::String
    additional_information::Vector{String}
end
Base.isapprox(x::StructureBlock, y::StructureBlock; kwargs...) = x.n_atoms == y.n_atoms && isapprox(x.cell, y.cell; kwargs...) && all(isapprox.(x.atoms, y.atoms; kwargs...)) && x.symmetry == y.symmetry && all(isapprox.(x.additional_information, y.additional_information; kwargs...))