abstract type AbstractProfile end

profilesymmetry(::Type{<:AbstractProfile}) = NoSymmetry()

"""
    Returns the profile of the fluid in the local fluid frame
"""
function profile end
struct GaussianRing{R,W} <: AbstractProfile
    r0::R
    width::W
end

profilesymmetry(::Type{<:GaussianRing}) = Axisymmetric()

@inline function profile(a::GaussianRing, t, r, θ, ϕ)
    return exp(-4*log(2)*abs2( (r-a.r0)/a.width) )
end
