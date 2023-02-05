module JuKeBOX

# Write your package code here.
using ComradeBase
using DocStringExtensions
using LinearAlgebra
using StaticArrays

export raytrace, intensity_point, profile, Kerr

#include("geodesics.jl")
#nclude("bam.jl")
include("metrics.jl")
include("profiles.jl")
include("raytracers/kerr_analytic.jl")
include("accretionModel.jl")
end
