"""
    $(TYPEDEF)
Kerr spacetime with spin `spin`. We assume that we are in unitless scales where
the black hole mass defines the length, so `M=1`.
"""
struct Kerr{S}
    spin::S
end

"""
    $(TYPEDEF)
For a spacetime `g` construct metric tensor at position `p`.
There is an optional cache that can be used to cache intermediate results for increased
speed.
"""
struct MetricTensor{M, P, C}
    g::M
    p::P
    cache::C
end

"""
    $(SIGNATURES)
Construct the metric tensor for spacetime `g` at the location `t`, `r`, `θ`, `ϕ`.

Returns a `MetricTensor` object. To access the components of the tensor please
the components function.
"""
function metric(g::Kerr, t, r, θ, ϕ)
    Δ = _Δ(g, r)
    Σ = _Σ(g, r, θ)
    Ξ = _Ξ(g, r, θ,  Δ)
    Ξsin2 = Ξ*sin(θ)^2
    ω = _ωZamo(g, r, Ξ,)
    cache = (;Δ, Σ, Ξ, Ξsin2, ω)
    return MetricTensor(g, (t,r,θ,ϕ), cache)
end

"""
    $(SIGNATURES)
Returns the components of the metric tensor `g` in SMatrix.
"""
function components(g::MetricTensor{K, C}) where {K<:Kerr, C}
    Δ, Σ, Ξ, Ξsin2, ω = g.cache
    gtt = -Δ*Σ/Ξ + ω^2*Ξsin2/Σ
    gtϕ = -Ξsin2/Σ*ω
    grr = Σ/Δ
    gθθ = Σ
    gϕϕ = Ξsin2/Σ
    T = eltype(gtt)
    return @SMatrix [ gtt     zero(T) zero(T) gtϕ    ;
                      zero(T) grr     zero(T) zero(T);
                      zero(T) zero(T) gθθ     zero(T);
                      gtϕ     zero(T) zero(T) gϕϕ
                    ]
end



"""
    $(SIGNATURES)
Computes the radial potential of the metric tensor assuming the energy scaled angular momentum
`λ`` and carter constant `η``
"""
function ℛ(g::MetricTensor, η, λ)
    r = g.p[2]
    Δ = g.cache.Δ
    return _ℛ(g.g, r, λ, η, Δ)
end


struct Minkowski{T} end
Minkowski() = Minkowski{Float64}()
metric(m::Minkowski, t, r, θ, ϕ) = MetricTensor(m, (t,r,θ,ϕ), (;))
components(::MetricTensor{<:Minkowski{T}, P,C}) where {T,P,C} = Diagonal(SVector(-one(T), one(T), one(T), one(T)))
