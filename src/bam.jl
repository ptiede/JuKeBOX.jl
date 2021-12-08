"""
    Differentiable accelerated Ray Tracing

This is a differentiable general relativistic ray-tracer. I think that's it
"""




struct GaussianRing{R,W}
    rpeak::R
    width::W
end

@inline profile(p::GaussianRing, r) = exp(-4*log(2)*abs2( (r-p.rpeak)/p.width ))


struct Observer{D,O}
    distance::D
    inclination::O
end


struct BAM{S, FD, V, B}
    nmax::Int
    """spectral index"""
    α::S
    """Function that computes the profile as radius `r` and redshift `z`"""
    profile::FD
    """ FluidVelocity struct"""
    β::V
    """ MagneticField struct"""
    b::B
end

function bam(nmax, spin, α, rpeak, width, βv, χ, ι, η=χ+π)
    g = Kerr(spin)
    v = FluidVelocity(βv, χ)
    b = MagneticField(ι, η)
    profile = GaussianRing(rpeak, width)
    return g, BAM{typeof(α), typeof(profile), typeof(v), typeof(b)}(nmax, α, profile, v, b)
end


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

@inline function MagneticField(ι, v::FluidVelocity)
    return MagneticField(ι, v.χ + π)
end


@inline function magnetic_vector(b::MagneticField)
    return SVector(b.br, b.bϕ, b.bz)
end

@inline function profile(b::BAM, r)
    profile(b.profile, r)
end




struct RayCache{T,R}
    interiorflag::Bool
    sβ::T
    λ::T
    η::T
    u₋a²::T
    roots::R
    rootdiffs::R
    F₀::T
    K₀::T
    ffac::T
    k::T
    Ir_total::T
    # This is the specific stuff for interior rays
    rp::T
    Agl::T
    Bgl::T
    I3r::T
    # This is specific stuff for exterior rays
    f₀::T
    Ir_turn::T
end





@inline function _I3(Agl, Bgl, k)
    I3r_angle = acos((Agl-Bgl)/(Agl+Bgl))
    I3r = F(I3r_angle, k)/sqrt(Agl*Bgl)
    return I3r
end

@inline function _I3total(I3r, Agl, Bgl, k, rp, r1::Real, r2::Real)
    I3rp_angle = acos( (Agl*(rp - r1) - Bgl*(rp - r2))/(Agl*(rp - r1) + Bgl*(rp - r2)) )
    I3rp = F(I3rp_angle, k)/sqrt(Agl*Bgl)
    return I3r - I3rp
end

function interior_ray_cache(sβ, λ, η, u₋a², roots, rootdiffs, F0, K0, ffac, a)
    r1,r2,_,_ = roots
    r31, r32, r42, r41 = rootdiffs
    rp = 1+sqrt(1-a^2)
    Agl = real(sqrt(r32*r42))
    Bgl = real(sqrt(r31*r41))
    rr1 = real(r1)
    rr2 = real(r2)
    k = ((Agl+Bgl)^2 - (rr2-rr1)^2)/(4*Agl*Bgl)
    I3r = _I3(Agl, Bgl, k)
    I3total = _I3total(I3r, Agl, Bgl, k, rp, rr1, rr2)

    return RayCache(true, sβ, λ, η, u₋a², roots, rootdiffs, F0, K0, ffac,
                            k, I3total,
                            rp, Agl, Bgl, I3r,
                            zero(typeof(λ)), zero(typeof(λ)) #zero this out
                    )

end

function exterior_ray_cache(sβ, λ, η, u₋a², roots, rootdiffs, F0, K0, ffac)
    r31, r32, r42, r41 = rootdiffs
    k = real(r32*r41/(r31*r42))
    f₀ = F(real(asin(sqrt(r31/r41))), k)
    Ir_turn = real(2/sqrt(r31*r42)*f₀)
    Ir_total = 2*Ir_turn
    T = typeof(λ)
    return RayCache(false, sβ, λ, η, u₋a², roots, rootdiffs, F0, K0, ffac,
                            k, Ir_total,
                            zero(T), zero(T), zero(T), zero(T),
                            f₀, Ir_turn)

end


function raycache(sβ, λ, η, g::Kerr, o::Observer)
    a = g.spin
    up, um = get_up_um(λ, η, a)
    urat = up/um
    u₋a² = um*a^2

    roots = radialroots(λ, η, a)
    r1,r2,r3,r4 = roots
    rootdiffs = (r3-r1, r3-r2, r4-r2, r4-r1)
    r31, r32, r42, r41 = rootdiffs

    F₀_sin = clamp(cos(o.inclination)/sqrt(up), -1.0, 1.0)
    F₀ = F(asin(F₀_sin), urat)
    K₀ = K(urat)
    ffac = 1/2*sqrt(real(r31*r42))

    # case 3 interior rays
    if abs(imag(r3)) > 1e-14
        return interior_ray_cache(sβ, λ, η, u₋a², roots, rootdiffs, F₀, K₀, ffac, a)
    else
        return exterior_ray_cache(sβ, λ, η, u₋a², roots, rootdiffs, F₀, K₀, ffac)
    end

end

function traceimg!(polim, alpha, beta, g, o, bam)
    stokesi = polim.I
    stokesq = polim.Q
    stokesu = polim.U
    @inbounds for C in CartesianIndices(stokesi)
        iy,ix = Tuple(C)
        x = alpha[ix]
        y = beta[iy]
        i,q,u = raytrace(x, y, g, o, bam)
        stokesi[C] = i
        stokesq[C] = q
        stokesu[C] = u
    end
end



function raytrace(α, β, g::Kerr, o::Observer, bam::BAM)

    # Find the geodesic charges
    λ, η = get_λ_η(α, β, o.inclination, g.spin)
    # Kill the vortical geodesics
    η < 0 && return (
                     zero(eltype(η)), #i
                     zero(eltype(η)), #q
                     zero(eltype(η)), #u
                     )


    # precompute some stuff that will be constant over ray
    # TODO: Make sure this isn't causing dynamic dispatch...
    cache = raycache(sign(β), λ, η, g, o)

    # create the stokes i,q,u allocators
    T = eltype(λ)
    stokesi = zero(T); stokesq = zero(T); stokesu = zero(T)
    for n in 0:bam.nmax
        _, i, q, u = trace_nring(n, α, β, cache, g, o, bam)
        stokesi += i
        stokesq += q
        stokesu += u
    end

    return stokesi, stokesq, stokesu
end

function trace_nring(n::Int, α, β, cache::RayCache, g::Kerr, o::Observer, bam::BAM)
    r, spr = _emission_radius(n, cache)
    T = eltype(r)
    # bail out
    r < cache.rp*1.01 && return zero(T), zero(T), zero(T), zero(T)
    r < zero(T) && return zero(T), zero(T), zero(T), zero(T)
    return _emission(n, α, β, cache.λ, cache.η, r, spr, g, o, bam)
end

function momentum_1form(m, sβ, λ, η, spr, Δs, ℛs)
    pt = -1
    pr = clamp(spr*sqrt(ℛs)/Δs, -10, 10)
    pϕ = λ
    pθ = (-1)^m*sβ*sqrt(η)
    return SVector(pt, pr, pθ, pϕ)
end

function zamo_tetrad(gs::MetricTensor{M,C}) where {M<:Kerr, C}
    Δs, Σs, Ξs, _, ωs = gs.cache
    _,r,_,_ = gs.p
    rtΞ = sqrt(Ξs)
    rtΔ = sqrt(Δs)
    ett, etϕ = 1/r*rtΞ/rtΔ, 1/r*rtΞ/rtΔ*ωs
    err = 1/r*rtΔ
    eϕϕ = r/rtΞ
    eθθ = -1/r
    T = zero(Δs)
    return @SMatrix [ ett     zero(T) zero(T) etϕ    ;
                      zero(T) err     zero(T) zero(T);
                      zero(T) zero(T) eθθ     zero(T);
                      zero(T) zero(T) zero(T) eϕϕ
                    ]
end

function penrose_walker(pcon, k, spin, rs)
    pt = pcon[1]; pr = pcon[2]; pθ = pcon[3]; pϕ = pcon[4]
    kt = k[1];    kr = k[2];    kθ = k[3];    kϕ = k[4]

    AA = (pt*kr - pr*kt) + spin*(pr*kϕ - pϕ*kr)
    BB = (rs^2 + spin^2) * (pϕ*kθ - pθ*kϕ) - spin * (pt*kθ - pθ*kt)
    κ1 = rs * AA
    κ2 = -rs * BB

    return κ1, κ2

end

# function momentum_vec(m, sβ, λ, η, spr, Δs, ℛs, r, a)
#     pt = -1/r^2*(-a*(a-λ) + (r^2 + a^2)*(r^2 + a^2 - a*λ)/Δs)
#     pr = spr*1/r^2*sqrt(ℛs)
#     pϕ = 1/r^2*(-(a-λ) + a/Δs*(r^2 + a^2 - a*λ))
#     pθ = (-1)^m*sβ*sqrt(η)/r^2
#     return SVector(pt, pr, pθ, pϕ)
# end



function _emission(n, α, β, λ, η, r, spr, g, o, bam)
    sβ = sign(β)
    m = get_nturns(n, sβ)
    # Get the metric stuff
    gs = metric(g, 0, r, π/2, 0)
    Δs, _, _, _, _ = gs.cache
    ℛs = ℛ(gs, λ, η)

    # For a 4x4
    gij = components(gs)
    # Don't worry this uses a special inverse for 4x4 matrices that should be stable
    invgij = inv(gij)
    ηij = components(metric(Minkowski(), 0, 0, 0, 0))

    # # Construct the photon momentum 1-form and vector
    pform = momentum_1form(m, sβ, λ, η, spr, Δs, ℛs)
    pcon = invgij*pform

    # # Construct the ZAMO tetrad
    ezamo = zamo_tetrad(gs)

    # # Transform to fluid frame
     vfield = bam.β
     #Λ = get_Λ(vfield.β, vfield.χ)
     #efluid = ηij*Λ*ηij*ezamo
     # -β is inverse and vectors are covariant so apply inverse transform
     Λ = get_Λ(-vfield.β, vfield.χ)
     efluid = Λ*ezamo
     pfluid = ηij*efluid*pform

    @inbounds begin
        z = 1/pfluid[1]
        lp = abs(pfluid[1]/pfluid[3])
        # #Now lets get the polarization vector
        b = bam.b
        p3 = SVector(pfluid[2], pfluid[3], pfluid[4])
    end
    bvec = magnetic_vector(b)
    f3 = cross(p3, bvec)
    f = @inbounds SVector(zero(eltype(f3)), f3[1], f3[2], f3[3])

    #Now move back to coordinate basis
    fkerr = efluid'*f
    # Constuct PW constant
    κ1, κ2 = penrose_walker(pcon, fkerr, g.spin, r)

    # screen appearance
    ν = -(α + g.spin * sin(o.inclination))
    eα = (β*κ2 - ν*κ1) / (ν^2 + β^2)
    eβ = (β*κ1 + ν*κ2) / (ν^2 + β^2)
    # Get the profile value at the emission radius
    prof = profile(bam, r)*z^(3+bam.α)
    q = -(eα^2 - eβ^2)*lp*prof
    u = -2*eα*eβ*lp*prof
    i = hypot(q,u)
    return z,i,q,u
end

@inline function get_nturns(n, sβ)
    if sβ > 0
        return sβ+n
    else
        return sβ+n+1
    end
end

@inline function minotime(m, cache)
    _Ir(cache.u₋a², m, cache.K₀, cache.sβ, cache.F₀)
end

@inline function _Ir(u₋a², m, K₀, sβ, F₀)
    return (2*m*K₀ - sβ*F₀)/sqrt(-u₋a²)
end


function _emission_radius(n::Int, cache::RayCache)
    if cache.interiorflag
        return _interior_emission_radius(n, cache)
    else
        return _exterior_emission_radius(n, cache)
    end
end

function _exterior_emission_radius(n::Int, cache)
    # get the correct number of turning point
    m = get_nturns(n, cache.sβ)
    _,_,r3,r4 = cache.roots
    r31,_,_,r41 = cache.rootdiffs

    Ir = minotime(m, cache)
    Ir > cache.Ir_total && return -one(typeof(Ir)), one(eltype(Ir))
    spr = sign(cache.Ir_turn - Ir)
    sn = Jacobi.sn(cache.ffac*Ir - cache.f₀, cache.k)
    snsqr = sn^2
    #abs(snsqr) > 1.1 && return NaN*one(eltype(real(sn))), one(eltype(real(sn)))
    r = real((r4*r31 - r3*r41*snsqr)/(r31-r41*snsqr))
    return r, spr
end

function _interior_emission_radius(n::Int, cache)
    # get the correct number of turning point
    m = get_nturns(n, cache.sβ)
    r1,r2,r3,r4 = cache.roots
    rr1 = real(r1)
    rr2 = real(r2)
    Agl = cache.Agl; Bgl = cache.Bgl
    Ir = minotime(m, cache)
    Ir > cache.Ir_total && return -one(typeof(Ir)), one(eltype(Ir))
    X = sqrt(Agl*Bgl)*(Ir - cache.I3r)
    cn = Jacobi.cn(X, cache.k)
    return ((Bgl*rr2 - Agl*rr1) + (Bgl*rr2 + Agl*rr1)*cn) / ((Bgl-Agl)+(Bgl+Agl)*cn), one(typeof(X))
end
