abstract type AbstractProfile end

struct FluidVelocity{T}
    β::T
    χ::T
    βr::T
    βϕ::T
    function FluidVelocity(β, χ)
        T = promote_type(eltype(β), eltype(χ))
        βT, χT = promote(β, χ)
        s,c = sincos(χT)
        return new{T}(βT, χT, βT*c, βT*s)
    end
end

struct MagneticField{T}
    ι::T
    η::T
    beq::T
    bz::T
    br::T
    bϕ::T
    function MagneticField(ι, η)
        T = promote_type(typeof(ι), typeof(η))
        ιT, ηT = promote(ι, η)
        beq,bz = sincos(ι)
        s,c = sincos(η)
        br = beq*c
        bϕ = beq*s
        return new{T}(ιT, ηT, beq, bz, br, bϕ)
    end
end
"""
    GaussianRing
A gaussian ring profile
"""
struct GaussianRing <: AbstractProfile
    """ FluidVelocity struct"""
    fluid_velocity::FluidVelocity{Float64}
    """ MagneticField struct"""
    magnetic_field::MagneticField{Float64}
    """emission spectral index"""
    α::Float64
    """emission cross product spectral index"""
    αζ::Float64
    """Function that computes the profile as radius `r` and redshift `z`"""
    rpeak::Float64
    width::Float64
end
@inline profile(p::GaussianRing, r) = exp(-4.0*log(2.0)*abs2( (r-p.rpeak)/p.width ))

"""
    DoublePower
A double power law for a profile
    ``\\frac{r}{r_0}^a\\left(1+\\frac{r}{r_0}^{-(a+b)}\\right)^{-1}``
"""
struct DoublePower <: AbstractProfile
    """ FluidVelocity struct"""
    fluid_velocity::FluidVelocity{Float64}
    """ MagneticField struct"""
    magnetic_field::MagneticField{Float64}
    """emission spectral index"""
    α::Float64
    """emission cross product spectral index"""
    αζ::Float64
    """Function that computes the profile as radius `r` and redshift `z`"""
    r0::Float64
    p1::Float64
    p2::Float64
end
@inline profile(p::DoublePower, r) = (r/p.r0)^p.p1/(1 + (r/p.r0)^(p.p1+p.p2))

