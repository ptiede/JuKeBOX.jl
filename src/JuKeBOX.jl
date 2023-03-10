module JuKeBOX

# Write your package code here.
using DocStringExtensions
using LinearAlgebra
#using Polyester
using StaticArrays
#using StructArrays

export Observer, BAM, bam, bamDblPower, raytrace, intensity_point

include("geodesics.jl")
include("metric/kerr.jl")
include("bam.jl")
end
