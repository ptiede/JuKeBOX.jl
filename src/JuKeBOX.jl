module JuKeBOX

# Write your package code here.
using Elliptic
using DocStringExtensions
using LinearAlgebra
using StaticArrays
using ROSE


include("model_helpers.jl")
include("kerr.jl")
include("bam.jl")
end
