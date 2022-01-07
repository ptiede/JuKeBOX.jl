abstract type Chart{M} end

abstract type Manifold end


struct Coordinate{A<:Chart, T} <: FieldVector{4, T}
    x0::N
    x1::N
    x2::N
    x3::N
end

const Coordinates{A} = Coordinates{A, T} where {N, A}

include(joinpath(@__DIR__, "kerr.jl"))
include(joinpath(@__DIR__, "minkowski.jl"))
