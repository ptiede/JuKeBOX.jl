"""
    Differentiable accelerated Ray Tracing

This is a differentiable general relativistic ray-tracer. I think that's it
"""


struct GaussianRing{R,W}
    rpeak::R
    width::W
end

@inline profile(p::GaussianRing, r) = exp(-4*log(2)*abs2( (r-p.rpeak)/p.width ))

"""
    DblPower
A double power law for a profile
"""
struct DblPower{R,A,B}
    r0::R
    a::A
    b::B
end

@inline profile(p::DblPower, r) = (r/p.r0)^p.a/(1 + (r/p.r0)^(p.a+p.b))

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
    return SVector(b.br, b.bz, b.bϕ)
end

@inline function profile(b::BAM, r)
    profile(b.profile, r)
end




@concrete struct RayCache
    interiorflag::Bool
    sβ
    λ
    η
    u₋a²
    roots
    rootdiffs
    F₀
    K₀
    k
    Ir_total
    # This is the specific stuff for interior rays
    rp
    Agl
    Bgl
    I3r
    # This is specific stuff for exterior rays
    f₀
    ffac
    Ir_turn
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

 function interior_ray_cache(sβ, λ, η, u₋a², roots, rootdiffs, F0, K0, a)
    r1,r2,_,_ = roots
    r31, r32, r42, r41 = rootdiffs
    Agl = real(sqrt(r32*r42))
    Bgl = real(sqrt(r31*r41))
    rr1 = real(r1)
    rr2 = real(r2)
    k = ((Agl+Bgl)^2 - (rr2-rr1)^2)/(4*Agl*Bgl)
    I3r = _I3(Agl, Bgl, k)
    rp = 1+sqrt(1-a^2)
    I3total = _I3total(I3r, Agl, Bgl, k, rp, rr1, rr2)

    return RayCache(true, sβ, λ, η, u₋a², roots, rootdiffs, F0, K0,
                            k, I3total,
                            rp, Agl, Bgl, I3r,
                            zero(typeof(λ)), zero(typeof(λ)), zero(typeof(λ)) #zero this out
                    )

end

 function exterior_ray_cache(sβ, λ, η, u₋a², roots, rootdiffs, F0, K0, a)
    r31, r32, r42, r41 = rootdiffs
    k = real(r32*r41/(r31*r42))
    f₀ =1# F(real(asin(sqrt(r31/r41))), k)
    rp = 1+sqrt(1-a^2)
    Ir_turn = 1#real(2/sqrt(r31*r42)*f₀)
    Ir_total = !#2*Ir_turn
    ffac = 1# 1/2*sqrt(real(r31*r42))
    T = typeof(λ)
    return RayCache(false, sβ, λ, η, u₋a², roots, rootdiffs, F0, K0,
                            k, Ir_total,
                            rp, zero(T), zero(T), zero(T),
                            f₀, ffac, Ir_turn)

end





 function raycache(sβ, λ, η, g::Kerr, o::Observer)
    a = g.spin
    up, um = get_up_um(λ, η, a)
    urat = up/um
    u₋a² = um*a^2

    roots = radialroots(λ, η, a)
    r1,r2,r3,r4 = roots
    rootdiffs = (r3-r1, r3-r2, r4-r2, r4-r1)

    F₀_sin = 1#clamp(cos(o.inclination)/sqrt(up), -1.0, 1.0)
    F₀ = 1#F(asin(F₀_sin), urat)
    K₀ = 1#K(urat)
    # case 3 interior rays
    return exterior_ray_cache(sβ, λ, η, u₋a², roots, rootdiffs, F₀, K₀, a)

end

function traceimg!(polim, alpha, beta, g, θs, o, bam)
    stokesi = polim.I
    stokesq = polim.Q
    stokesu = polim.U
    @batch for C in CartesianIndices(stokesi)
        iy,ix = Tuple(C)
        x = alpha[ix]
        y = beta[iy]
        i,q,u = raytrace(x, y, g, θs, o, bam, true) + raytrace(x, y, g, θs, o, bam, false)
        stokesi[C] = i
        stokesq[C] = q
        stokesu[C] = u
    end
end




function tracepolmap(alpha, beta, g, θs, o, bam)
    nx,ny = length(alpha), length(beta)
    polim = StructArray((I=zeros(ny, nx), Q = zeros(ny, nx), U = zeros(ny, nx)))
    traceimg!(polim, alpha, beta, g, θs, o, bam)
    return polim
end




function raytrace(α, β, g::Kerr, θs, o::Observer, bam::BAM, isindir)

    # Find the geodesic charges
    # λ, η = get_λ_η(α, β, o.inclination, g.spin)
    # Kill the vortical geodesics
    #η < 0 && return StokesVector{typeof(λ)}(0,0,0,0)


    # precompute some stuff that will be constant over ray
    #cache = raycache(sign(β), λ, η, g, o)
    # create the stokes i,q,u allocators
    
    stokesi = 0.; stokesq = 0.; stokesu = 0.
    for n in 0:bam.nmax
        _, i, q, u = trace_nring(n, α, β, g, θs, o, bam, isindir)
        stokesi += i
        stokesq += q
        stokesu += u
    end

    return StokesVector(stokesi, stokesq, stokesu, zero(typeof(stokesi)))
end

function trace_nring(n::Int, α, β, g::Kerr, θs, o::Observer, bam::BAM, isindir)
    #r, spr = _emission_radius(n, cache)
    r, νr, _ = rs(α, β, θs, o.inclination, g.spin, isindir, n) 
    T = eltype(r)

    νθ =  cos(θs)< abs(cos(o.inclination)) ? (o.inclination>θs) ⊻ (n%2==1) : !isindir
    # bail out
    #cache.k > 1 && return zero(T), zero(T), zero(T), zero(T)
    r < 1+ √(1-g.spin^2) && return zero(T), zero(T), zero(T), zero(T)
    #r < zero(T) && return zero(T), zero(T), zero(T), zero(T)
    #r > 100 && return zero(T), zero(T), zero(T), zero(T)
    return _emission(α, β, r, νr, νθ, g, o, bam, θs)
end

function mpow(m)
    m == 0 && return one(eltype(m))
    res = -1
    for _ in 1:m
        res *= -1
    end
    return res
end

function momentum_1form(m, sβ, λ, η, spr, Δs, ℛs)
    pt = -1
    pr = clamp(spr*sqrt(abs(ℛs))/Δs, -15, 15)
    pϕ = λ
    pθ = mpow(m)*sβ*sqrt(η)
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

#function penrose_walker(pcon, k, spin, rs)
#    pt = pcon[1]; pr = pcon[2]; pθ = pcon[3]; pϕ = pcon[4]
#    kt = k[1];    kr = k[2];    kθ = k[3];    kϕ = k[4]
#
#    AA = (pt*kr - pr*kt) + spin*(pr*kϕ - pϕ*kr)
#    BB = (rs^2 + spin^2) * (pϕ*kθ - pθ*kϕ) - spin * (pt*kθ - pθ*kt)
#    κ1 = rs * AA
#    κ2 = -rs * BB
#
#    return κ1, κ2
#
#end

# function momentum_vec(m, sβ, λ, η, spr, Δs, ℛs, r, a)
#     pt = -1/r^2*(-a*(a-λ) + (r^2 + a^2)*(r^2 + a^2 - a*λ)/Δs)
#     pr = spr*1/r^2*sqrt(ℛs)
#     pϕ = 1/r^2*(-(a-λ) + a/Δs*(r^2 + a^2 - a*λ))
#     pθ = (-1)^m*sβ*sqrt(η)/r^2
#     return SVector(pt, pr, pθ, pϕ)
# end



function _emission(α, β, r, νr, νθ, g, o, bam, θs)
    mag_field = bam.b
    fluid_vel = bam.β
    B = @SVector [mag_field.br, mag_field.bϕ, mag_field.bz]
    fluidβ = @SVector[fluid_vel.β, 0. ,fluid_vel. χ]
    κ1, κ2, redshift, lp = calcPol(α, β, r,  θs, o.inclination, g.spin, B, fluidβ , νr, νθ)

    # screen appearance
    ν = -(α + g.spin * sin(o.inclination))
    enorm = (ν^2 + β^2)*sqrt(κ1^2+κ2^2) + eps()
    eα = (β*κ2 - ν*κ1) / enorm
    eβ = (β*κ1 + ν*κ2) / enorm
    # Get the profile value at the emission radius
    prof = profile(bam, r)*redshift^(3+bam.α)
    # We add a small perturbation to q and u. This prevent taking a
    # derivative of hypot at (0,0) which gives  a NaN. This is silly because this
    # occurs when there is no emission so the derivative should be zero there.
    q = -(eα^2 - eβ^2)*lp*prof + eps()
    u = -2*eα*eβ*lp*prof + eps()
    i = hypot(q, u)
    return redshift,i,q,u
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


