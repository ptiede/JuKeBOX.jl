module JuKeBOX

# Write your package code here.
using ComradeBase
using ConcreteStructs: @concrete
using Elliptic
using DocStringExtensions
using Elliptic
using LinearAlgebra
using Polyester
using StaticArrays
using StructArrays

export Observer, BAM, bam, raytrace, intensity_point

abstract type AccretionModel <: ComradeBase.AbstractPolarizedModel end
ComradeBase.visanalytic(::Type{<:AccretionModel}) = ComradeBase.NotAnalytic()
ComradeBase.imanalytic(::Type{<:AccretionModel}) = ComradeBase.IsAnalytic()

struct SimpleModel{B,G,O} <: AccretionModel
    acc::B
    g::G
    o::O
end
ComradeBase.isprimitive(::Type{<:SimpleModel}) = ComradeBase.IsPrimitive()



ComradeBase.intensity_point(s::AccretionModel, α, β) = raytrace(α, β, s.g, s.o, s.acc)


include("model_helpers.jl")
include("metric/kerr.jl")
include("bam.jl")
end
