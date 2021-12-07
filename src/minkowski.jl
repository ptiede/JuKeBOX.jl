struct Minkowski{T} <: Manifold end
Minkowski() = Minkowski{Float64}()

struct Cartesian{M<:Minkowski} <: Chart{M} end
coordinate(m::Cartesian{M<:Minkowski}, t::T, x::T, y::T, z::T) = Coordinate{typeof(m), T}(t,x,y,z)

metric(m::Minkowski, t, r, θ, ϕ) = MetricTensor(m, (t,r,θ,ϕ), (;))
components(::MetricTensor{<:Minkowski{T}, P,C}) where {T,P,C} = Diagonal(SVector(-one(T), one(T), one(T), one(T)))
