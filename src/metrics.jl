abstract type AbstractMetric end

struct Schwarzschild <: AbstractMetric
    mass::Float64
end

struct Kerr <: AbstractMetric
    spin::Float64
    mass::Float64
end

