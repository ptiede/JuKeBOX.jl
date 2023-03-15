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
    """emission cross product spectral index"""
    αζ::S
    """Function that computes the profile as radius `r` and redshift `z`"""
    profile::FD
    """ FluidVelocity struct"""
    β::V
    """ MagneticField struct"""
    b::B
    #BAM(nmax, α, αζ, profile, β, b) = new{typeof(α), typeof(profile), typeof(β), typeof(b)}(nmax, Vector{Tuple{Float64,Bool,Int}}(undef,nmax+1), α, αζ, profile, β, b)
end

function bam(nmax, spin, α, αζ, rpeak, width, βv, χ, ι, η=χ+π)
    g = Kerr(spin)
    v = FluidVelocity(βv, χ)
    b = MagneticField(ι, η)
    profile = GaussianRing(rpeak, width)
    return g, BAM{typeof(α), typeof(profile), typeof(v), typeof(b)}(nmax, α, αζ, profile, v, b)
    #return g, BAM(nmax, α, αζ, profile, v, b)

end
bam(nmax, spin, α, rpeak, width, βv, χ, ι, η=χ+π) = bam(nmax, spin, α, α, rpeak, width, βv, χ, ι, η)

function bamDblPower(nmax, spin, α, αζ, rpeak, p1, p2, βv, χ, ι, η=χ+π)
    g = Kerr(spin)
    v = FluidVelocity(βv, χ)
    b = MagneticField(ι, η)
    profile = DblPower(rpeak, p1, p2)
    return g, BAM{typeof(α), typeof(profile), typeof(v), typeof(b)}(nmax, α, αζ, profile, v, b)
    #return g, BAM(nmax, α, αζ, profile, v, b)
end
bamDblPower(nmax, spin, α, rpeak, p1, p2, βv, χ, ι, η=χ+π) = bamDblPower(nmax, spin, α, α, rpeak, p1, p2, βv, χ, ι, η)

bamGauss(nmax, spin, α, rpeak, width, βv, χ, ι, η=χ+π) = bam(nmax, spin, α, rpeak, width, βv, χ, ι, η)


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

function raytrace(α, β, g::Kerr, θs, o::Observer, bam::BAM, isindir)
    stokesi = 0.; stokesq = 0.; stokesu = 0.
    for n in 0:bam.nmax
        _, i, q, u = trace_nring(n, α, β, g, θs, o, bam, isindir)
        stokesi += i
        stokesq += q
        stokesu += u
    end

    return @SVector[stokesi, stokesq, stokesu, zero(typeof(stokesi))]
end

function raytrace_n(n, α, β, g::Kerr, θs, o::Observer, bam::BAM, isindir)
    _, stokesi, stokesq, stokesu = trace_nring(n, α, β, g, θs, o, bam, isindir)

    return @SVector[stokesi, stokesq, stokesu, zero(typeof(stokesi))]
end

function raytrace_and_get_mask(n, α, β, g::Kerr, θs, o::Observer, bam::BAM, isindir)
   (_, stokesi, stokesq, stokesu), mask = trace_nring_and_get_mask(n, α, β, g, θs, o, bam, isindir)#, radBuffer)

    return @SVector[stokesi, stokesq, stokesu, zero(typeof(stokesi))], mask
end

function trace_nring(n::Int, α, β, g::Kerr, θs, o::Observer, bam::BAM, isindir)
    r, νr, _ = rs(α, β, θs, o.inclination, g.spin, isindir, n) 
    T = eltype(r)

    νθ =  cos(θs)< abs(cos(o.inclination)) ? (o.inclination>θs) ⊻ (n%2==1) : !isindir

    # bail out if r is within horizon
    (r == Inf || r < 1+ √(1-g.spin^2) + eps()) && return zero(T), zero(T), zero(T), zero(T)

    return _emission(α, β, r, νr, νθ, g, o, bam, θs)
end

function trace_nring_and_get_mask(n, α, β, g::Kerr, θs, o::Observer, bam::BAM, isindir)#, radBuffer)
    (r, νr, _), mask = rs_mask(n, α, β, θs, o.inclination, g.spin, isindir) 

    T = eltype(r)
    νθ =  cos(θs)< abs(cos(o.inclination)) ? (o.inclination>θs) ⊻ (n%2==1) : !isindir
    if (r == Inf || r < 1+ √(1-g.spin^2) + eps())
        return (zero(T), zero(T), zero(T), zero(T)), mask
    else
        return _emission(α, β, r, νr, νθ, g, o, bam, θs), mask
    end
end

function mpow(m)
    m == 0 && return one(eltype(m))
    res = -1
    for _ in 1:m
        res *= -1
    end
    return res
end

function _emission(α, β, r, νr, νθ, g, o, bam, θs)
    fluid_vel = bam.β
    fluidβ = @SVector[fluid_vel.β, π/2. , fluid_vel.χ]
    eα, eβ, redshift, lp = calcPol(α, β, r, θs, o.inclination, g.spin, bam.αζ,  magnetic_vector(bam.b), fluidβ , νr, νθ)

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