

"""
    $(TYPDEF)
Local EDF in the frame of the fluid. We include cutoffs to ensure the distribution
is integrable.
"""
struct SimplePowerLaw{T,B} <: AbstractEDF
    α::T
    νmin::B
    νmax::B
end

function SimplePowerLaw(α)
    return SimplePowerLaw(α, 1e-30, 1e30)
end

struct BAM{M, P, E, V, B, S} <: AccretionModel
    manifold::M
    """Function that computes the profile as radius `r`"""
    profile::P
    """local EDF function relative to fluid frame"""
    edf::E
    """velocity field"""
    β::V
    """magnetic field relative to fluid frame"""
    b::B
    """ Stop Condition """
    stop::S
end

"""
    $(SIGNATURES)
A convinience helper function used to construct the canonical BAM problem from
Palumbo (2022) using the BAMVelocity profile, BAMMagneticField, and the ReverseNMax
stop condition.

# Arguments

- `nmax`: The maximum numbe of turning points considered before stopping ray-Tracing
- `profile`: The profile function we will be using
- `α`: The local electron **spectral** index assumed in the flow of the fluid
- `βv`: The fluid flow velocity magnitude.
- `χ`: The angle of the fluid velocity relative to the ZAMO reference frame
- `ι`: ``B_{eq} = \\sin(\\iota)`` which defines the relative strength of the B-field
- `η`: The angle of the magnetic field relative to outward radial unit vector
"""
function BAM(nmax::Int, profile::AbstractProfile, α, βv, χ, ι, η=χ+π)
    edf = SimplePowerLaw(α)
    β = BAMVelocity(βv, χ)
    b = BAMMagneticField(ι, η)
    stop = ReverseNMax(nmax)
    return BAM(profile, edf, β, b, stop)
end

@inline velocityfield(b::BAM, p::Coordinate) = b.β
@inline magneticfield(b::BAM, p::Coo) = b.b

struct BAMVelocity{T}
    β::T
    χ::T
    βr::T
    βϕ::T
    function FluidVelocity(β, χ)
        T = promote_type(eltype(β), eltype(χ))
        χT, βT = promote(β, χ)
        s,c = sincos(χT)
        return new{T}(βT, χT, βT*c, βT*s)
    end
end

struct BAMMagneticField{T}
    ι::T
    η::T
    beq::T
    bz::T
    br::T
    bϕ::T
    function MagneticField(ι, η)
        T = promote_type(ι, η)
        ιT, ηT = promote(ι, η)
        beq,bz = sincos(ι)
        s,c = sincos(η)
        br = beq*c
        bϕ = beq*s
        return new{T}(ιT, ηT, beq, bz, br, bϕ)
    end
end

@inline function MagneticField(ι, v::FluidVelocity)
    return MagneticField(ι, v.χ + π)
end


@inline function magnetic_vector(b::MagneticField)
    return SVector(b.br, b.bϕ, b.bz)
end

@inline function profile(b::BAM, r, z)
    b.profile(r, z)
end



abstract type RayCache end

struct InteriorRayCache{T,R} <: RayCache
    sβ::T
    λ::T
    η::T
    u₋a²::T
    roots::R
    rootdiffs::R
    F₀::T
    K₀::T
    ffac::T
    # This is the specific stuff for interior rays
    rp::T
    Agl::T
    Bgl::T
    k::T
    I3r::T
    It_total::T
end

struct ExteriorRayCache{T,R} <: RayCache
    sβ::T
    λ::T
    η::T
    u₋a²::T
    roots::R
    rootdiffs::R
    F₀::T
    K₀::T
    ffac::T
    #This is the specific stuff for exterior ray
    k::T
    f₀::T
    Ir_turn::T
    It_total::T
end



@inline function _I3(Agl, Bgl, k)
    I3r_angle = acos((Agl-Bgl)/(Agl+Bgl))
    I3r = F(I3r_angle, k)/sqrt(Agl*Bgl)
    return I3r
end

@inline function _I3total(I3r, Agl, Bgl, k, rp, r1::Real, r2::Real)
    I3rp_angle = acos((Agl*(rp-r1)-Bgl*(rp-r2))/(Agl*(rp-r1)+Bgl*(rp-r2)))
    I3rp = ellipkinc(I3rp_angle, k)/sqrt(Agl*Bgl)
    return I3r - I3rp
end

function InteriorRayCache(sβ, λ, η, u₋a², roots, rootdiffs, F0, K0, ffac, a)
    r1,r2,_,_ = roots
    r31, r32, r42, r41 = rootdiffs
    rp = 1+sqrt(1-a^2)
    Agl = real(sqrt(r32*r42))
    Bgl = real(sqrt(r31*r41))
    rr1 = real(r1)
    rr2 = real(r2)
    k = ((Agl+Bgl)^2 - (rr2-rr1)^2)/sqrt(4*Agl*Bgl)
    I3r = _I3(Agl, Bgl, k)
    I3total = _I3total(I3r, Agl, Bgl, k, rp, rr1, rr2)

    return InteriorRayCache(sβ, λ, η, u₋a², roots, rootdiffs, F0, K0, ffac,
                            rp, Agl, Bgl, k, I3r, I3total)

end

function ExteriorRayCache(sβ, λ, η, u₋a², roots, rootdiffs, F0, K0, ffac)
    r31, r32, r42, r41 = rootdiffs
    k = real(r32*r41/(r31*r42))
    f₀ = F(real(asin(sqrt(r31/r41))), k)
    Ir_turn = real(2/sqrt(r31*r42)*f₀)
    Ir_total = 2*Ir_turn

    return ExteriorRayCache(sβ, λ, η, u₋a², roots, rootdiffs, F0, K0, ffac,
                            k, f₀, Ir_turn, Ir_total)

end



function raycache(sβ, λ, η, g::Kerr, o::Observer)
    a = g.spin
    up, um = get_up_um(λ, η, a)
    urat = up/um
    u₋a² = um*a^2

    roots = radialroots(λ, η, a)
    r1,r2,r3,r4 = roots
    rootdiffs = (r3-r1, r3-r2, r4-r2, r4-r1)
    r31, _, r42, _ = rootdiffs

    F₀_sin = cos(o.inclination)/sqrt(up)
    F₀ = F(asin(F₀_sin), urat)
    K₀ = K(urat)
    ffac = 1/2*sqrt(real(r31*r42))

    # case 3 interior rays
    if (imag(r3) > eps(typeof(um)))
        return InteriorRayCache(sβ, λ, η, u₋a², roots, rootdiffs, F₀, K₀, ffac, a)
    else
        return ExteriorRayCache(sβ, λ, η, u₋a², roots, rootdiffs, F₀, K₀, ffac)
    end

end


function raytrace(α, β, o::Observer, bam::BAM)
    g = bam.manifold
    # Find the geodesic charges
    λ, η = get_λ_η(α, β, o.inclination, g)
    # Kill the vortical geodesics
    η < 0 && return (
                     zero(eltype(η)), #i
                     zero(eltype(η)), #q
                     zero(eltype(η)), #u
                     zero(eltype(η))  #v
                     )

    _raytrace(symmetry(bam), α, β, o, bam)
    # precompute some stuff that will be constant over ray
    # TODO: Make sure this isn't causing dynamic dispatch...
    cache = raycache(sign(β), λ, η, g, o)

    # create the stokes i,q,u allocators
    T = eltype(λ)
    stokesi = zero(T); stokesq = zero(T); stokesu = zero(T); stokesv = zero(T)
    for n in 0:bam.nmax
        _, i, q, u, v = trace_nring(n, α, β, cache, g, o, bam)
        stokesi += i
        stokesq += q
        stokesu += u
        stokesv += v

    end

    return stokesi, stokesq, stokesu, stokesv
end

function trace_nring(n::Int, α, β, cache::RayCache, g::Kerr, o::Observer, bam::BAM)
    r, spr = _emission_radius(n, cache)
    T = eltype(r)
    # bail out
    isnan(r) && return zero(T), zero(T), zero(T), zero(T), zero(T)
    r > cache.Ir_total && return zero(T), zero(T), zero(T), zero(T), zero(T)

    return _emission(n, α, β, cache.λ, cache.η, r, spr, g, o, bam)
end

function momentum_1form(λ, η, spr, Δs, ℛs)
    pt = -1
    pr = clamp(spr*sqrt(ℛs)/Δs, -10, 10)
    pθ = λ
    pθ = (-1)^m*sβ*sqrt(η)
    return SVector(pt, pr, pθ, pϕ)
end

function zamo_tetrad(gs::MetricTensor{M,C}) where {M<:Kerr, C}
    Δs, Σs, Ξs, _, ωs = gs.cache
    rtΞ = sqrt(Ξ)
    rtΔ = sqrt(Δ)
    ett, etϕ = 1/r*rtΞ/rtΔ.*(1, ω)
    err = 1/r*rtΔ
    eϕϕ = r/rtΞ
    eθθ = -1/r
    T = zero(Δs)
    return @SMatrix [ ett     zero(T) zero(T) etϕ    ;
                      zero(T) err     zero(T) zero(T);
                      zero(T) zero(T) zero(T) eϕϕ    ;
                      zero(T) zero(T) eθθ     zero(T)
                    ]
end

@inline function penrose_walker(pcon, k, spin, rs)
    pt = pcon[1]; pr = pcon[2]; pθ = pcon[2]; pϕ = pcon[3]
    kt = k[1]; kr = k[2]; kθ = k[3]; kϕ = k[4]

    AA = (pt * kfr - pr * kft) + spin * (pr * kfphi - pphi * kfr)
    BB = (rs^2 + spin^2) * (pphi * kftheta - ptheta * kfphi) - spin * (pt * kftheta - ptheta * kft)
    κ1 = rs * AA
    κ2 = -rs * BB

    return κ1, κ2

end

function _emission(n, α, β, λ, η, r, spr, g, o, bam)
    sβ = sign(β)
    m = get_nturns(n, sβ)
    # Get the metric stuff
    gs = metric(g, 0, r, π/2, 0)
    Δs, Σs, Ξs, _, ωs = gs.cache
    ℛs = ℛ(g, λ, η)

    # For a 4x4
    gij = components(gs)
    invgij = inv(gij)
    ηij = components(MetricComponents(Minkowski()))
    # Construct the photon momentum 1-form and vector
    pform = momentum_1form(λ, η, spr, Δs, ℛs)
    pcon = invgij*pform

    # Construct the ZAMO tetrad
    ezamo = zamo_tetrad(gs)

    # Transform to fluid frame
    vfield = velocity(bam)
    Λ = get_Λ(vfield.β, vfield.χ)
    efluid = ηij*Λ*ηij*ezamo
    pfluid = ηij*efluid*pform # Computes and raises the index

    z = 1/pfluid[1]
    lp = abs(pfluid[1]/pfluid[3])

    #Now lets get the polarization vector
    b = magneticfield(bam)
    p3 = @view pfluid[2:end]
    f3 = cross(p3, magnetic_vector(b))
    f = SVector(zero(eltype(f3)), f3[1], f3[2], f3[3])

    #Now move back to coordinate basis this needs a transpose right?
    fkerr = efluid'*f

    # Constuct PW constant
    pt = pcon
    κ1, κ2 = penrose_walker(pcon, fkerr, g.spin, r)

    # screen appearance
    ν = -(α + g.spin*sin(o.inclination))
    eα = (β*κ2 - ν*κ1)/(ν^2 + β^2)
    eβ = (β*κ1 + ν*κ2)/(ν^2 + β^2)

    # Get the profile value at the emission radius
    prof = profile(bam, r, 0.0)*z^(3+bam.α)
    q = -(eα^2 - eβ^2)*lp*prof
    u = -2*eα*eβ*lp*prof
    i = hypot(q,u)

    return z, i, q, u
end

@inline function get_nturns(n, sβ)
    return sβ+1+n
end

@inline function _Ir(u₋a², m, K₀, sβ, F₀)
    return 1/sqrt(-u₋a²)*(2*m*K₀ - sβ*F₀)
end


function _emission_radius(n::Int, cache::ExteriorRayCache)
    # get the correct number of turning point
    m = get_nturns(n, cache.sβ)
    r1,r2,r3,r4 = cache.roots
    r31,r32,r42,r41 = cache.rootdiffs

    Ir = _Ir(cache.u₋a², m, cache.K₀, cache.sβ, cache.F₀)
    sn = Jacobi.sn(cache.ffac*Ir - cache.f₀, cache.k)
    snsqr = sn^2
    #abs(snsqr) > 1.1 && return NaN*one(eltype(real(sn))), one(eltype(real(sn)))
    r = real((r4*r31 - r3*r41*snsqr)/(r31-r41*snsqr))
    return r, one(eltype(r))
end

function _emission_radius(n::Int, cache::InteriorRayCache)
    # get the correct number of turning point
    m = get_nturns(n, cache.sβ)
    r1,r2,r3,r4 = cache.roots
    rr1 = real(r1)
    rr2 = real(r2)

    Ir = _Ir(-cache.u₋a², m, cache.K₀, cache.sβ, cache.F₀)
    spr = sign(cache.Ir_turn - Ir)
    X = sqrt(Agl*Bgl)*(Ir - spr*cache.I3r)
    cn = Jacobi.cn(X, cache.k)
    Agl, Bgl = cache.Agl, cache.Bgl
    return ((Bgl*rr2 - Agl*rr1) + (Bgl*rr2+Agl*rr1)*cn) / ((Bgl-Agl)+(Bgl+Agl)*cn), spr
end
