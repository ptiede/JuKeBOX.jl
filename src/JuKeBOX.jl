module JuKeBOX

# Write your package code here.
using Elliptic
using DocStringExtensions
using LinearAlgebra
using Polyester
using StaticArrays
using StructArrays
using ROSE


include("model_helpers.jl")
include("kerr.jl")
include("bam.jl")
end
