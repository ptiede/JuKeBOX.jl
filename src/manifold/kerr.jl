"""
    $(TYPEDEF)
Kerr spacetime with spin `spin`. We assume that we are in unitless scales where
the black hole mass defines the length, so `M=1`.
"""
struct Kerr{S} <: Manifold
    spin::S
end

struct BoyerLindquivst{M} <: Chart{M}
    m::M
end


@inline function coordinate(b::BoyerLindquivst{K}, t::T, r::T, θ::T, ϕ::T) where {T}
    Coordinate{T,  typeof(b)}(t, r, θ, ϕ)
end


"""
    $(TYPEDEF)
For a spacetime `g` construct metric tensor at position `p`.
There is an optional cache that can be used to cache intermediate results for increased
speed.
"""
struct MetricTensor{M<:Chart, P, C}
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
function metric(g::Kerr, coordinate::BoyerLindquivst)
    t, r, θ, ϕ = coordinate
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

@inline function get_λ_η(α, β,θ, g::Kerr)
    a = g.spin
    λ = -α * sin(θ)
    η = (α^2 - a^2)*cos(θ)^2+β^2
    return λ, η
end

@inline function get_up_um(λ, η, g::Kerr)
    a = g.spin
    Δ_θ  = 1/2*(1-(η+λ^2)/a^2)
    up = Δ_θ + √(Δ_θ^2 +η/a^2)
    um = Δ_θ - √(Δ_θ^2 +η/a^2)
    return up, um
end

function radialroots(λ::Real, η::Real, g::Kerr)
    a = g.spin
    T = promote_type(eltype(λ), eltype(η))
    return radialroots(Complex{T}(λ), Complex{T}(η), a)
end


"""
    $(SIGNATURES)
Finds the radial roots using the of the geodesic with
energy-scaled angular momentum `λ` and carter constant `η`.
"""
function radialroots(λ::Complex, η::Complex, g::Kerr)
    a = g.spin
    A = a^2 - η - λ^2
    B = 2*(η+(λ-a)^2)
    C = -a^2 * η
    P = -A^2 / 12 - C
    Q = -A/3 * ((A/6)^2-C)-B^2/8
    H = -9*Q+√(12*P^3+81*Q^2)
    z = √((-2*(3^(1/3)*P)+2^(1/3)*H^(2/3))/(2*6^(2/3)*H^(1/3) - A/6))
    r1 = -z-√(-A/2-z^2+B/(4*z))
    r2 = -z+√(-A/2-z^2+B/(4*z))
    r3 = z-√(-A/2-z^2+B/(4*z))
    r4 = z+√(-A/2-z^2+B/(4*z))
    return r1, r2, r3, r4
end

@inline function _Δ(g::Kerr, r)
    return r^2 - 2*r + g.spin^2
end

@inline function _Σ(g::Kerr, r, θ)
    return r^2 + (g.spin*cos(θ))^2
end

@inline function _Ξ(g::Kerr, r, θ, Δ)
    return (r^2 + g.spin^2)^2 - Δ*(g.spin*sin(θ))^2
end

@inline function _ωZamo(g, r, Ξ)
    return 2*g.spin*r/Ξ
end

@inline function _ℛ(g, r, λ, η, Δ)
    a = g.spin
    return (r^2 + a^2 - a*λ)^2 + Δ*(η + (a-λ)^2)
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
