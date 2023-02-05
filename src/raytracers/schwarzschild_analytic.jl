import FastElliptic
export p, ϕ, ψforward, φ
export get_roots, rt, rs, findb, ψt, ψbackwards_in, ψbackwards_out_back, ψbackwards_out_front, magnification, magnification_numeric, dψ_db
##
# Follows the Formalism of Schwarzschild Magnification Notes
##

##----------------------------------------------------------------------------------------------------------------------
# Useful functions
##----------------------------------------------------------------------------------------------------------------------

function βboundary(α, θo, θs) 
  cosθs2 = cos(θs)^2
  √max((cos(θo)^2-cosθs2)*(α^2)/(cosθs2 -1), 0.0)
end

"""
  r_potential(met::Schwarzschild, r, η, λ)

Radial potential of a kerr blackhole

  `met` : Schwarzschild Metric

  `η`   : Reduced Carter constant

  `λ`   : Reduced azimuthal agular momentum

  `r`   : Boyer Lindquist radius

"""
function r_potential(met::Schwarzschild, η, λ, r) 
  λ2 = λ^2
  r*(2*η+2*λ2+r*(-η-λ2+r^2)) # Eq 7 PhysRevD.101.044032
end

"""
  θ_potential(met::Schwarzschild, r, η, λ)

Theta potential of a kerr blackhole

  `met` : Schwarzschild Metric

  `η`   : Reduced Carter constant

  `λ`   : Reduced azimuthal agular momentum

  `θ`   : Boyer Lindquist inclination

"""
θ_potential(met::Schwarzschild, η, λ, θ) = η - λ^2*cot(θ)^2

"""
  get_roots(b)

Returns double the number `x` plus `1`.
"""
function get_roots(b)
  q = 2*b^2 + 0im
  p = -b^2 + 0im
  u = (-q/2 + (q^2/4 + p^3/27)^(1/2))^(1/3) 
  C = ((-1+0im)^(2/3), (-1+0im)^(4/3), 1) .* u
  v = -p ./ ((3+0im) .* C)
  return sort([i for i in C .+ v], lt=(x,y)->real(x)<real(y))
end

"""
  ϕ(θs, φ, θo, n)

Screen arival angle of emission.

  `θs`: Emitter Inclination

  `ϕ` : Boyer Lindquist ϕ

  `θo`: Observer Inclination

  'n' : Index of nth image
"""
function φ(θs::Real, ϕ::Real, θo::Real, n::Real)#φ screen coordinate
  p1 = cos(ϕ)sin(θs)
  p2 = cos(θs)sin(θo) + cos(θo)sin(θs)sin(ϕ)
  return (sign(p2)*acos(p1 / √(p1^2 + p2^2))) + n*π 
end

"""
  ψforward(θs, ϕ, θo)

Winding angle of direct emission in the forward raytracing problem.

  `θs`: The half opening angle of the cone that the point lies on.

  `ϕ`: The parameterized angle around the point.

  `θo`: The inclination angle of the point.
"""
function ψforward(θs::Real, ϕ::Real, θo::Real)#ψ angle swept out 
  #p3 = sin(α)cos(θ) - cos(α)sin(θ)sin(φ)
  #p3 = (cos(α)cos(θ)^2+(cos(α)cos(φ)^2+sin(φ)^2)*sin(θ)^2)
  #p3 = cos(α)*cos(θ)^2 + (cos(φ)^2 + cos(α)*sin(φ)^2)*sin(θ)^2
  p3 = cos(θs)*cos(θo) - sin(θs)*sin(θo)*sin(ϕ)
  return acos(min(p3, 1))
end


η(α, β, θo) = (α^2)*cos(θo)^2 + β^2
λ(α, θo) = -α*sin(θo)
##----------------------------------------------------------------------------------------------------------------------
# Radial Stuff
##----------------------------------------------------------------------------------------------------------------------

function rs(α, β, θs, θo, isindir, n) 
    if cos(θs) > abs(cos(θo))
        βbound = (abs(α) >= eps() ? βboundary(α, θo, θs) : 0.)
        if abs(β) + eps() < βbound
            return 0., true, 4
        end
    end
    ηtemp = η(α, β, θo)
    λtemp = λ(α, θo)
    b = √(α^2 + β^2)
    τ, _, _ = _Gθ(sign(β), θs, θo, a, isindir, n, ηtemp, λtemp)
    if τ != Inf
      return _rs(b, τ)
    else 
      return (0., true, 4) 
    end
end

"""
  _rs(b, τ)

Emission radius for emission that lies outside the photon ring and whose ray intersects the equatorial plane

  `λ` : Impact Parameter (Total Normalized Angular Momentum)

  `τ` : Mino Time
"""
function _rs(b, τ)
  ans = 0.0
  νr = true

  roots = get_roots(b)
  rh = 2
  numreals = (abs(imag(roots[1])) > 1e-10 ? 0 : 2) + (abs(imag(roots[3])) > 1e-10 ? 0 : 2)

  if numreals == 4 #case 1 & 2
    ans, νr = _rs_case1(real.(roots), rh, τ)
  elseif numreals == 2 #case3
    if abs(imag(roots[4])) < 1e-10
     roots= (roots[1], roots[4],roots[2],roots[3])
    end

    ans, νr = _rs_case3(roots, τ)
  else #case 4
    ans, νr = _rs_case4(roots, rh, τ)
    end
  return ans, νr, numreals
end
"""
  rsout(b::Real, ψ::Real)::Real

Emission radius for emission that lies outside the photon ring
"""
function rsout(b::Real, ψ::Real)::Real
  r1, r3, r4 = get_roots(b)
  r32 = r3
  r41 = r4 - r1
  r31 = r3 - r1
  r42 = r4
  k = abs((r32*r41) / (r31*r42))
  fo = Elliptic.F(asin(real(√(r31/r41))), min(k,1))
  sn = r41 * Elliptic.Jacobi.sn(fo- real(1/2 * √(r31*r42) * ψ/b) , k)^2
  return real((r31*r4 - r3*sn) / (r31 - sn))
end

"""
  rsin(b::Real, ψ::Real)::Real

Emission radius for emission that lies inside the photon ring

  `b` : Impact parameter

  `ψ` : Winding angle
"""
function rsin(b::Real, ψ::Real)::Real
  r1, r3, r4 = get_roots(b)
  r2 = 0
  r32 = r3 - r2
  r21 = r2 - r1
  r41 = r4 - r1
  r31 = r3 - r1
  r42 = r4
  A = √(r32*r42)
  B = √(r31*r41)
  k =  min(real(((A + B)^2 - r21^2)/(4*A*B)),1)
  fo = Elliptic.F(real(acos((A-B)/(A+B))), k)
  cn = Elliptic.Jacobi.cn(real(fo - √(A*B)*ψ/b), k)
  num = -A*r1 + B*r2 + (A*r1+B*r2)*cn
  den = -A + B + (A+B)*cn
  return real(num/den)
end

"""
  rs(b::Real, ψ::Real)

Emission radius of source.

  `b`: Impact parameter

  `ψ`: Winding angle
"""
rs(b::Real, ψ::Real) = b^2 > 27 ? rsout(b,ψ) : rsin(b,ψ)

##----------------------------------------------------------------------------------------------------------------------
# θ Stuff
##----------------------------------------------------------------------------------------------------------------------
"""
  τ(α, β, θs, θo, isindir::Bool, n::Int)

Mino time of trajectory between two inclinations for a given screen coordinate

  `α` : Bardeen α

  `β` : Bardeen β 

  `θs` : Emission inclination

  `θo` : Observer inclination

  `isindir` : Is the path direct or indirect?

  `n` : nth image in orde of amount of minotime traversed
"""
τ(α, β, a, θs, θo, isindir, n) = _τ(sign(β), θs, θo, a, isindir, n, η(α, β, θo), λ(α, θo))

function _τ(signβ, θs, θo, a, isindir, n, η, λ)
  return minotime, Ghat, isvortical
end 



"""
  p(α::Real, φ::Real, θ::Real, χ::Real)

Parameterized point in space bulk spacetime. The point is assumed to lie on a cone whose apex coincides with the 
center of the blackhole.

  `α`: The half opening angle of the cone that the point lies on.

  `φ`: The parameterized angle around the point.

  `θ`: The inclination angle of the point.

  `χ`: The azimuthal rotation of the point.
"""
function p(α::Real, φ::Real, θ::Real, χ::Real )

    return (
    cos(χ)cos(φ)sin(α) + sin(χ)(cos(α)sin(θ) - cos(θ)sin(α)sin(φ)),
    sin(χ)cos(φ)sin(α) - cos(χ)(cos(α)sin(θ) - cos(θ)sin(α)sin(φ)),
    cos(α)cos(θ) + sin(α)sin(θ)sin(φ)
  ) 
end


"""
  g(rs::real)

Redshift factor
  `rs` : The emission radius of the photon
"""
g(rs::Real) = √(1- 2/rs) 

dψ_drdb(b::Real, r::Real)::Real = real(r^4 / ((r^4 - b^2*r*(r-2+0im))^(3/2)))

function ψt(b::Real)
  if b^2 < 27
    return Inf
  else
    r1,r3,r4 = get_roots(b)
    r2 = 0
    r32 = r3 - r2
    r41 = r4 - r1
    r31 = r3 - r1
    r42 = r4 - r2

    k = min(real((r32*r41 / (r42*r31))),1)
    g = 2 / √(r42 * r31)
    return real(b * g * Elliptic.F(real(asin(√(r31/r41))), k))
  end
end

function dψ_du(b::Real, u::Real)
  """
      u^2 = r-r4
  """
  r1, r3, r4 = get_roots(b)
  r2 = 0

  return 2*b / √((u^2 + r4 - r1) * (u^2 + r4 - r2) * (u^2 + r4 - r3))
end

function droots_db(b::Real)
  p = -b^2
  q = 2 * b^2
  disc = (q^2 / (4 + 0im) + p^3 / (27 + 0im))^(1/2)
  C = (-q/(2+0im) + disc)^(1/3)

  temp_desc = √(81+0im-3*b^2)
  dC_db = -(3)^(1/3)*(2*b + b*(b^2 - 18)/temp_desc) / (b^2 * (temp_desc - 9))^(2/3)
  dp_db = -2*b
  dp_3C_db = (dp_db/C - p*dC_db/C^2) / 3

  dr1_db = dC_db*(-1+√(3)*1im)/2 - dp_3C_db/((-1+√(3)*1im)/2)
  dr3_db = dC_db*(-1-√(3)*1im)/2 - dp_3C_db/((-1-√(3)*1im)/2)
  dr4_db = dC_db - dp_3C_db
  return dr1_db, dr3_db, dr4_db
end

"""
  ψbackwards_out_front(b::Real, rs::Real)::Real

Winding angle of the backwards raytracing problem on the front side of a ring

  `b` : Impact parameter

  `rs` : Emission radius
"""
function ψbackwards_out_front(b::Real, rs::Real)::Real
  r1,r3,r4 = get_roots(b)
  r2 = 0
  r32 = r3 - r2
  r41 = r4 - r1
  r31 = r3 - r1
  r42 = r4 - r2

  k = min(real((r32*r41) / (r31*r42)),1)
  g = 2 / √(r42 * r31)
  fo = Elliptic.F(real(asin(√(r31/r41))), k)
  return real(b*real(g*(fo - Elliptic.F(real(asin(√(r31*(rs-r4)/(r41*(rs-r3))))), k))))
end

"""
  ψbackwards_out_bank(b::Real, rs::Real)::Real

Winding angle of the backwards raytracing problem on the back side of a ring

  `b` : Impact parameter

  `rs` : Emission radius
"""
function ψbackwards_out_back(b::Real, rs::Real)::Real
  return real(2ψt(b) - ψbackwards_out_front(b,rs))
end

"""
  ψbackwards_in(b::Real, rs::Real)::Real

Winding angle of the backwards raytracing problem for rays that fall into the blackhole

  `b` : Impact parameter

  `rs` : Emission radius
"""
function ψbackwards_in(b::Real, rs::Real)::Real
  r1,r3,r4 = get_roots(b)
  r2 = 0
  r32 = r3 - r2
  r41 = r4 - r1
  r31 = r3 - r1
  r42 = r4 - r2
  r21 = r2 - r1

  A = √(r32*r42)
  B = √(r31*r41)
  k =  min(real(((A + B)^2 - r21^2)/(4*A*B)),1)
  fo = Elliptic.F(real(acos((A-B)/(A+B))), k)
  x = real((A*(rs-r1)-B*(rs-r2))/(A*(rs-r1)+B*(rs-r2)))
  return real(b*(fo - Elliptic.F(acos(x), k)) / √(A * B))
end


dsnInv_dx(x::Real, k::Real) = 1 / √((1-x^2)*(1-k*x^2))
dsnInv_dk(x::Real, k::Real)  = -(2Elliptic.E(asin(x),k) + 2(k-1)Elliptic.F(asin(x),k) - k*sin(2asin(x))/√(1-k*x^2))/(4(k-1)*k)
dcnInv_dx(x::Real, k::Real) = -1 / √((1-x^2)*(1-k*(1-x^2)))
dcnInv_dk(x::Real, k::Real)  = -(2Elliptic.E(acos(x),k) + 2(k-1)Elliptic.F(acos(x),k) - k*sin(2acos(x))/√(1+k*(x^2-1)))/(4(k-1)*k)

function dgfo_out_db(b::Real)::Real
  r1,r3,r4 = get_roots(b)
  r2 = 0
  r32 = r3 - r2
  r41 = r4 - r1
  r31 = r3 - r1
  r42 = r4 - r2

  dr1, dr3, dr4 = droots_db(b)

  k = min(real(r32*r41 / (r42*r31)),1)
  dk = (dr3/r32 + (dr4-dr1)/r41 - dr4/r42 - (dr3-dr1)/r31)*k
  g = 2 / √(r42 * r31)
  dg = -(dr4/r42 + (dr3-dr1)/r31)*g/2
  x = real(√(r31/r41))
  dx = ((dr3-dr1)/r31 - (dr4-dr1)/r41)*x/2

  ans = dg * Elliptic.F(real(asin(x)), k)
  ans += g * dsnInv_dx(x, k)*dx
  ans += g * dsnInv_dk(x, k)*dk
  return real(ans)
end

function dΔτout_db(b::Real, rs::Real)::Real
  r1,r3,r4 = get_roots(b)
  r2 = 0
  r32 = r3 - r2
  r41 = r4 - r1
  r31 = r3 - r1
  r42 = r4 - r2

  dr1, dr3, dr4 = droots_db(b)

  
  k = min(real(r32*r41 / (r42*r31)),1)
  dk_db = (dr3/r32 + (dr4-dr1)/r41 - dr4/r42 - (dr3-dr1)/r31)*k
  g = 2 / √(r42 * r31)
  dg = -(dr4/r42 + (dr3-dr1)/r31)*(g/2)
  x = real(√(r31*(rs-r4)/(r41*(rs-r3))))
  dx_db = ((dr3-dr1)/r31 - dr4/(rs-r4) - (dr4-dr1)/r41 + dr3/(rs-r3))*x/2

  ans = dg * Elliptic.F(asin(x), k)
  ans += g * dsnInv_dx(x, k)*dx_db
  ans += g * dsnInv_dk(x, k)*dk_db

  return real(ans)
end

function dgfo_in_db(b::Real)::Real
  r1, r3, r4 = get_roots(b)
  r2 = 0
  r21 = r2 - r1
  r32 = r3 - r2
  r41 = r4 - r1
  r31 = r3 - r1
  r42 = r4 - r2

  dr1, dr3, dr4 = droots_db(b)
  dr2 = 0

  A = √(r32*r42)
  dA = (dr3/r32 + dr4/r42)*(A/2)
  B = √(r31*r41)
  dB = ((dr3-dr1)/r31 + (dr4-dr1)/r41)*(B/2)
  g = 1 / √(A*B)
  dg = (-dA/A-dB/B)*(g/2)
  x = real((A-B)/(A+B))
  dx = ((dA-dB) - (A-B)*(dA+dB)/(A+B))/(A+B)
  k = min(real(((A+B)^2 - r21^2)/(4A*B)),1)
  dk = (2((dA+dB)*(A+B) - (dr2-dr1)*(r21)) + ((A+B)^2 - r21^2)*(-dA/A-dB/B))/(4A*B)

  ans = dg * Elliptic.F(acos(x), k)
  ans += g * dcnInv_dx(x, k)*dx
  ans += g * dcnInv_dk(x, k)*dk
  return real(ans)
end


function dΔτin_db(b::Real, rs::Real)
  r1, r3, r4 = get_roots(b)
  r2 = 0
  r32 = r3 - r2
  r41 = r4 - r1
  r31 = r3 - r1
  r42 = r4 - r2
  r21 = r2 - r1

  dr1, dr3, dr4 = droots_db(b)
  dr2 = 0

  A = √(r32*r42)
  dA = (dr3/r32 + dr4/r42)*(A/2)
  B = √(r31*r41)
  dB = ((dr3-dr1)/r31 + (dr4-dr1)/r41)*(B/2)
  k = min(real(((A+B)^2 - r21^2)/(4A*B)),1)
  dk = (2((dA+dB)*(A+B) - (dr2-dr1)*(r21)) + ((A+B)^2 - r21^2)*(-dA/A-dB/B))/(4A*B)
  g = 1 / √(A*B)
  dg = (-dA/A-dB/B)*(g/2)
  x = real((A*(rs-r1)-B*(rs-r2))/(A*(rs-r1)+B*(rs-r2)))
  dx = ((dA*(rs-r1)-A*dr1-dB*(rs-r2))-(A*(rs-r1)-B*(rs-r2))*(dA*(rs-r1)-A*dr1+dB*(rs-r2))/(A*(rs-r1)+B*(rs-r2)))/(A*(rs-r1)+B*(rs-r2))

  ans = dg * Elliptic.F(acos(x), k)
  ans += g * dcnInv_dx(x, k)*dx
  ans += g * dcnInv_dk(x, k)*dk
  return real(ans)
end

dψout_front_db(b::Real, rs::Real) = b*(dgfo_out_db(b) - dΔτout_db(b, rs)) + ψbackwards_out_front(b,rs)/b
dψout_back_db(b::Real, rs::Real) = b*(dgfo_out_db(b) + dΔτout_db(b, rs)) + ψbackwards_out_back(b,rs)/b
dψin_db(b::Real, rs::Real) = b*(dgfo_in_db(b) - dΔτin_db(b, rs)) + ψbackwards_in(b,rs)/b


function dψ_dudb(b::Real, u::Real)
  """
      u^2 = r-r4
  """
  r1, r3, r4 = get_roots(b) 
  r2 = 0
  dr1_db, dr3_db, dr4_db = droots_db(b)
  dr2_db = 0
  n1 = -1/2 * (dr4_db - dr1_db) / (u^2 + r4 - r1)
  n2 = -1/2 * (dr4_db - dr2_db) / (u^2 + r4 - r2)
  n3 = -1/2 * (dr4_db - dr3_db) / (u^2 + r4 - r3)

  d = 2 / √((u^2 + r4 - r1) * (u^2 + r4 - r2) * (u^2 + r4 - r3))

  return d * (b * (n1 + n2 + n3) + 1)
end

"""
  magnification(b::Real, ψ::Real)::Real

Analytic value for the magnifications of a point lensed by a blackhole

  `b` : Real

  `ψ` : Real
"""
function magnification(b::Real, ψ::Real)::Real
  ri = rs(b,ψ)
  dψ_db = 0
  if b^2 >= 27
    if ψ > ψt(b)
      dψ_db = dψout_back_db(b,ri)
    else
      dψ_db = dψout_front_db(b,ri)
    end
  else
    dψ_db=dψin_db(b, ri)
  end
  return abs(b/√radial_potential(b, ri)/(sin(ψ)*dψ_db))
end

"""
  magnification_numeric(b::Real, ψ::Real)::Real

Numeric value for the magnifications of a point lensed by a blackhole

  `b` : Real

  `ψ` : Real
"""
function magnification_numeric(b::Real, ψ::Real)::Real
  ri = rs(b,ψ)
  dψ_db = quadgk(r->dψ_drdb(b, r), ri, 1000)[1]
  if b^2 > 27 && ψ > ψt(b)
    p = -b^2
    q = 2b^2
    disc = √(q^2/4 + p^3/27+0im)
    C = (-q/2 + disc)^(1/3)

    temp_desc = √(81-3b^2+0im)
    dC_db = -3^(1/3)*(2b+b*(b^2-18)/temp_desc)/(b^2*(temp_desc-9))^(2/3)
    dp_db = -2*b
    dp_3C_db = (dp_db/C - p*dC_db/C^2) / 3

    rT = C - p / (3 * C)
    us = real(√(ri-rT))
    dψ_db += 2*quadgk(u->dψ_dudb(b, real(u)), 0, us)[1]
    dψ_db -= dψ_du(b, us)/us * (dC_db - dp_3C_db)
  end
  return abs(b/√radial_potential(b, ri)/(sin(ψ)*dψ_db))
end

ψequatorial(φ::Real, θ::Real) =  acos(-sin(θ)*sin(φ))
ϕequatorial(φ::Real, θ::Real) = sign(sin(φ))*acos(cos(φ) / (√((cos(θ)*sin(φ))^2 +cos(φ)^2)))

