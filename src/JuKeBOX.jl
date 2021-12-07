module JuKeBOX

# Write your package code here.
using Elliptic
using DocStringExtensions
using LinearAlgebra
using StaticArrays


abstract type AccretionModel end

abstract type AbstractEDF end

abstract type AbstractStopCondition end



include("model_helpers.jl")
include("traits.jl")


struct ReverseNMax <: AbstractStopCondition
    nmax::Int
end

propogationdirection(::Type{NMax}) = Reverse()


include("observer.jl")
include("profile.jl")
include("manifold.jl")
include("kerr.jl")
include("bam.jl")
end
