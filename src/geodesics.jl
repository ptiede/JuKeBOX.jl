import Elliptic
using LinearAlgebra, StaticArrays
export get_roots, Gθ, rs, calcPol, η, λ, r_potential, θ_potential, λcrit, ηcrit, ϕ, Ipm, Δ, Σ, A
##
# Follows the Formalism of Gralla & Lupsasca (https://arxiv.org/pdf/1910.12881.pdf)
##

##----------------------------------------------------------------------------------------------------------------------
# Useful functions
##----------------------------------------------------------------------------------------------------------------------

αboundary(a::Real, θs::Real) = a*sin(θs)

function βboundary(α::Real, θo::Real, a::Real, θs::Real) 
    √max((cos(θo)^2-cos(θs)^2)*(α^2-a^2+a^2*cos(θs)^2)/(cos(θs)^2 -1), 0.0)
end



"""
  r_potential(r::Real, η::Real, λ::Real, a::Real)::Real

Radial potential of a kerr blackhole

  `η` : Reduced Carter constant

  `λ` : Reduced azimuthal agular momentum

  `a` : Blackhole spin

  `r` : Boyer Lindquist radius

"""
r_potential(η, λ, a::Real, r::Real) = (r^2+a^2-a*λ)^2-(r^2 -2r + a^2)*(η+(λ-a)^2) # Eq 7 PhysRevD.101.044032

"""
  θ_potential(r::Real, η::Real, λ::Real, a::Real)::Real

Theta potential of a kerr blackhole

  `η` : Reduced Carter constant

  `λ` : Reduced azimuthal agular momentum

  `a` : Blackhole spin

  `θ` : Boyer Lindquist inclination

"""
θ_potential(η::Real, λ::Real, a::Real, θ::Real) = η + a^2*cos(θ)^2 - λ^2*cot(θ)^2

"""
  get_roots(η::Real, λ::Real, a::Real)

Returns roots of r⁴ + (a²-η-λ²)r² + 2(η+(a-λ)²)r - a²η

  `η` : Normalized Carter constant

  `λ` : Normalized angular momentum

  `a` : Blackhole spin
"""
function get_roots(η::Real, λ::Real, a::Real)
    A = a^2 - η - λ^2
    B = 2(η + (λ-a)^2)
    C = -a^2*η   

    P = -A^2 / 12 - C
    Q = -A/3*((A/6 +0im)^2 - C) - B^2/8

    Δ3 = -4*P^3 - 27*Q^2
    ωp = (-Q/2 + (-Δ3/108)^(1/2) + 0im)^(1/3)

    C = ((-1+0im)^(2/3), (-1+0im)^(4/3), 1) .* ωp
    v = -P ./ ((3+0im) .* C)
    
    ξ = sort([i for i in C .+ v], lt=(x,y)->real(x)<real(y))
    ξ1 = ξ[1]
    ξ2 = ξ[2]
    ξ0 = ξ[3]

    ξ0 = sort([ξ0,ξ1,ξ2], lt=(x,y)->real(x)<real(y))[end]-A/3

    r1 = (-(2ξ0)^(1/2) - (-(2A + 2ξ0 - (√2B)/(√ξ0)))^(1/2))/2
    r2 = (-(2ξ0)^(1/2) + (-(2A + 2ξ0 - (√2B)/(√ξ0)))^(1/2))/2
    r3 = ((2ξ0)^(1/2) - (-(2A + 2ξ0 + (√2B)/(√ξ0)))^(1/2))/2
    r4 = ((2ξ0)^(1/2) + (-(2A + 2ξ0 + (√2B)/(√ξ0)))^(1/2))/2

    return [r1, r2, r3, r4]
end

Δ(r::Real, a::Real) = r^2 -2r + a^2
Σ(r::Real, θ::Real, a::Real) = r^2 + a^2*cos(θ)^2
A(r::Real, θ::Real, a::Real) = (r^2 + a^2)^2 - a^2*Δ(r, a)*sin(θ)^2
Ξ(r::Real, θ::Real, a::Real) = (r^2+a^2)^2-Δ(r, a)*a^2*sin(θ)^2
ω(r::Real, θ::Real, a::Real) = 2*a*r/ Ξ(r, θ, a)

η(α::Real, β::Real, θo::Real, a::Real) = (α^2 - a^2)*cos(θo)^2 + β^2
λ(α::Real, θo::Real) = -α*sin(θo)

rtildep(a::Real) = 2*(1+Cos(2/3*acos(a)))
rtilden(a::Real) = 2*(1+Cos(2/3*acos(-a)))


"""
    λcrit(r::Complex, a::Real)

Returns λ values on the critical curve associated with a given r.

  `r` : Radius of orbit associated with critical η value

  `a` : Blackhole spin
"""
λcrit(r::Real, a::Real) = a + r/a*(r- 2Δ(r, a)/(r-1))
"""
    ηcrit(r::Complex, a::Real)

Returns η values on the critical curve associated with a given r.

  `r` : Radius of orbit associated with critical η value

  `a` : Blackhole spin
"""
ηcrit(r::Real, a::Real) = (r^3/a^2)*(4*Δ(r,a)/(r-1)^2 - r)


##----------------------------------------------------------------------------------------------------------------------
# Radial Stuff
##----------------------------------------------------------------------------------------------------------------------

function get_root_diffs(r1, r2, r3, r4)
  r21 = r2 - r1
  r31 = r3 - r1
  r32 = r3 - r2
  r41 = r4 - r1
  r42 = r4 - r2

  return r21, r31, r32, r41, r42
end

function rs(α, β, θs, θo, a, isindir, n) 

    if cos(θs) > abs(cos(θo))
        αmin = αboundary(a, θs)
        βbound = (abs(α) >= αmin ? βboundary(α, θo, a, θs) : 0.)
        if abs(β) < βbound
            return 0., true, 4
        end
    end
    #βtemp = ((θo > π/2) ⊻ (n%2 == 1)) ? -β : β
    ηtemp = η(α, β, θo, a)
    λtemp = λ(α, θo)
    τ = Gθ(α, β, a, θs, θo, isindir, n)[1]
    if τ != Inf; return _rs(ηtemp, λtemp, a, τ)  else (0., true, 4) end
end

"""
  _rs(η::Real, λ::Real, a::Real, τ::Real)::Real

Emission radius for emission that lies outside the photon ring and whose ray intersects the equatorial plane

  `η` : Normalized Carter constant

  `λ` : Normalized angular momentum

  `a` : Angular Momentum

  `τ` : Mino Time
"""
function _rs(η::Real, λ::Real, a::Real, τ::Real)
  ans = 0
  νr = true

  r1, r2, r3, r4 = get_roots(η, λ, a)
  rh = 1 + √(1-a^2)
  numreals = (abs(imag(r1)) > 1e-10 ? 0 : 2) + (abs(imag(r3)) > 1e-10 ? 0 : 2)

  if numreals == 4
    r1, r2, r3, r4 = real(r1), real(r2), real(r3), real(r4)
  elseif abs(imag(r4)) < 1e-10
     r2,r3,r4 = r4,r2,r3
  end

  roots = [r1, r2, r3, r4]
  r21, r31, r32, r41, r42 = get_root_diffs(r1, r2, r3, r4)
  root_diffs = [r21, r31, r32, r41, r42]

  if numreals == 4. #case 1 & 2
    if r4 >= rh && τ > 2I2r_turn(root_diffs) # invalid case1
      ans = 0
    elseif r4 < rh && τ > I2r(roots, root_diffs, rh, true) # invalid case2
      ans = 0
    else
      k = (r32*r41) / (r31*r42)

      fo = Elliptic.F(asin(√(r31/r41)), k)
      X2 = fo - √(r31*r42) * τ / 2
      sn = r41 * Elliptic.Jacobi.sn(X2, k)^2

      ans = (r31*r4 - r3*sn) / (r31 - sn)
      νr = X2 > 0
    end
  elseif numreals == 2. #case3

    if τ > I3r_full(root_diffs)#(roots, root_diffs, rh)
      ans = 0
    else
      A = √real(r32*r42)
      B = √real(r31*r41)
      k =  real(((A + B)^2 - r21^2)/(4*A*B))

      fo = Elliptic.F(acos((A-B)/(A+B)), k)
      X3  = real(fo - √(A*B)*τ)
      cn = Elliptic.Jacobi.cn(X3, k)
      num = -A*r1 + B*r2 + (A*r1+B*r2)*cn
      den = -A + B + (A+B)*cn

      ans = real(num/den)
      νr = X3 > 0 
    end
  else #case 4
    if τ > I4r(roots, root_diffs, rh)
      ans = 0
    else
      a2 = abs(imag(r1))
      b1 = real(r4)
      C = √real(r31*r42)
      D = √real(r32*r41)
      k4 = 4*C*D/(C+D)^2
      
      go = √(4a2^2 - (C-D)^2) / ((C+D)^2 - 4a2^2)
      fo = 2/(C+D)*Elliptic.F(π/2 + atan(go), k4) 
      X4 =  (C+D)/2*(fo - τ)
      num = go - Elliptic.Jacobi.sc(X4, k4)
      den = 1 + go*Elliptic.Jacobi.sc(X4, k4)

      ans = -(a2*num/den + b1)
      νr = X4 > 0 
    end
  end
  return ans, νr, numreals
end

function I2r_turn(root_diffs::Vector{Float64})
  _, r31, r32, r41, r42 = root_diffs
  k = r32*r41/(r31*r42)
  return 2/√real(r31*r42)*Elliptic.F(asin(√(r31/r41)), k)
end

function I2r(roots::Vector{Float64}, root_diffs::Vector{Float64}, rs, isindir)
  _, _, r3, r4 = roots
  _, r31, r32, r41, r42 = root_diffs

  k = r32*r41/(r31*r42)
  x2_s = √((rs-r4)/(rs-r3)*r31/r41)
  if !(-1 < x2_s < 1); return 0.; end
  Ir_s = 2/√real(r31*r42)*Elliptic.F(asin(x2_s), k)
  Ir_turn = I2r_turn(root_diffs)

  return Ir_turn + (isindir ? Ir_s : -Ir_s)
end

function I3r_full(root_diffs)
  r21, r31, r32, r41, r42 = map(abs, root_diffs)
  A2 = r32*r42
  B2 = r31*r41
  if A2 < 0. || B2 < 0; return Inf; end

  A, B = √A2, √B2
  k3 = ((A+B)^2 - r21^2)/(4A*B)

  temprat = B/A
  x3_turn = real(√((1+0im - temprat)/(1 + temprat)))
  return 2/√real(A*B)*Elliptic.F(acos(x3_turn), k3)
end

function I3r(roots, root_diffs, rs)
  r1, r2, _, _ = roots
  r21, r31, r32, r41, r42 = root_diffs

  A2 = real(r32*r42)
  B2 = real(r31*r41)
  if A2 < 0. || B2 < 0; return 0; end

  A, B = √A2, √B2


  k3 = real(((A+B)^2 - r21^2)/(4A*B))
  temprat = B*(rs-r2)/(A*(rs-r1))
  x3_s = real(√((1+0im - temprat)/(1 + temprat)))
  Ir_s = 2/√real(A*B)*Elliptic.F(real(acos(x3_s)), k3)
  Ir_full = I3r_full(root_diffs)

  return abs(Ir_full - Ir_s)
end

function I4r_full(roots, root_diffs)
  r1, _, _, r4 = roots
  _, r31, r32, r41, r42 = root_diffs

  try
    C = √real(r31*r42)
    D = √real(r32*r41)
    k4 = 4C*D/(C+D)^2
    a2 = abs(imag(r1))

    k4 = 4*C*D/(C+D)^2
    
    go = √max((4a2^2 - (C-D)^2) / ((C+D)^2 - 4a2^2), 0.)
    return 2/(C+D)*Elliptic.F(π/2 + atan(go), k4) 
  catch e
    return 0
  end
end

function I4r(roots, root_diffs, rs)
  r1, _, _, r4 = roots
  _, r31, r32, r41, r42 = root_diffs


  C = √real(r31*r42)
  D = √real(r32*r41)
  k4 = 4C*D/(C+D)^2
  a2 = abs(imag(r1))
  b1 = real(r4)

  k4 = 4*C*D/(C+D)^2
  
  go = √max((4a2^2 - (C-D)^2) / ((C+D)^2 - 4a2^2), 0.)
  x4_s = (rs + b1)/a2
  Ir_s = 2/(C+D)*Elliptic.F(atan(x4_s) + atan(go), k4) 
  Ir_full = I4r_full(roots, root_diffs)

  return Ir_full - Ir_s
end


#
##----------------------------------------------------------------------------------------------------------------------
# θ Stuff
##----------------------------------------------------------------------------------------------------------------------
"""
  Gθ(η, λ, a::Real, θs::Real, θo::Real, isindir::Bool, n::Int64)::Real

Mino time of trajectory between two inclinations for a given screen coordinate

  `α` : Bardeen α

  `β` : Bardeen β 

  `a` : Blackhole angular Momentum

  `θs` : Emission inclination

  `θo` : Observer inclination

  `isindir` : Is the path direct or indirect?

  `n` : nth image in orde of amount of minotime traversed
"""
function Gθ(α::Real, β::Real, a::Real, θs::Real, θo::Real, isindir::Bool, n::Int64)
  Go, Gs, Ghat, minotime, isvortical = 0, 0, 0, 0, false
  ηtemp = η(α, β, θo, a)
  λtemp = λ(α, θo)

  Δθ = 1/2*(1 - (ηtemp + λtemp^2)/a^2)
  up = Δθ + √(Δθ^2 + ηtemp/a^2)
  um = Δθ - √(Δθ^2 + ηtemp/a^2)
  m = up/um

  #if um == 0; return Inf, isvortical ;end


  if ηtemp > 0. #Ordinary motion
    k = m
    # Skip sign flip because another sign flip occurs when stitching together the mino time
    args =  cos(θs)/√(up)
    argo =  cos(θo)/√(up)

    if !(-1 < args < 1) || !(-1 < argo < 1);return Inf, isvortical;end

    Gs = 1/√(-um*a^2 +0im)*Elliptic.F(asin(args), k)
    Go = 1/√(-um*a^2 +0im)*Elliptic.F(asin(argo), k)
    Ghat = 2/√(-um*a^2 +0im)*(Elliptic.K(k))
  
  else #Vortical motion
    isvortical = true 

    argo = ((cos(θo)^2 - um)/(up-um))
    args = ((cos(θs)^2 - um)/(up-um))

    if !(0. < argo < 1.) ||  !(0. < args <  1.);return Inf, isvortical;end

    Υo = asin(√argo)
    Υs = asin(√args)    
    k = 1. - m
    Go = 1/√(um*a^2)*Elliptic.F(Υo, k)
    Gs = 1/√(um*a^2)*Elliptic.F(Υs, k)
    Ghat = 2/√(um*a^2)*Elliptic.K(k)

  end

  if cos(θs) < abs(cos(θo))
        if isindir != ((β > 0) ⊻ (θo > π/2))
          return Inf, isvortical
        end
  end

  if cos(θs) < abs(cos(θo))
    #if θo < π/2
      #minotime = real(isindir ? (n+1)*Ghat - Go - Gs : n*Ghat + Go - Gs ) #Sign of Go indicates whether the ray is from the forward cone or the rear cone
    #else
    νθ =  θo < π/2 ? -1 : 1
    minotime = real(isindir ? (n+1)*Ghat -sign(β)*Go + (n%2==1 ? -1 : 1)*νθ*Gs : n*Ghat - sign(β)*Go + (n%2==1 ? -1 : 1)*νθ*Gs ) #Sign of Go indicates whether the ray is from the forward cone or the rear cone
    #end
  else
    #minotime = real((isindir ? (-cos(θo) > cos(θs) ? (n+1)*Ghat+(Gs + Go) : (n+1)*Ghat-(Gs + Go)) : (β < 0 ? n*Ghat + Gs + Go : n*Ghat + Gs - Go) ))
    minotime = real((isindir ? (n+1)*Ghat-(Gs + sign(β)*Go) : (β < 0 ? n*Ghat + Gs - (isvortical ? -1 : 1)*sign(β)*Go : n*Ghat + Gs - sign(β)*Go) ))
  end

  if (((β < 0) ⊻ (n%2==1)) && cos(θs) > abs(cos(θo)) && !isvortical) || (isvortical && θo >= π/2)
    return Inf, isvortical
  end
  return minotime, isvortical
end 

##----------------------------------------------------------------------------------------------------------------------
#Polarization stuff
##----------------------------------------------------------------------------------------------------------------------
MinkowskiMet() = [-1. 0. 0. 0.; 0. 1. 0. 0.; 0. 0. 1. 0.; 0. 0. 0. 1.]

p_boyer_lindquist_d(r::Real, θ::Real, a::Real, η::Real, λ::Real, νr::Bool, νθ::Bool) = [-1,(νr ? 1 : -1)*√(r_potential(η, λ, a, r))/Δ(r, a), λ, (νθ ? 1 : -1)*√θ_potential(η, λ, a, θ)]

"""
    kerr_met_uu(r::Real, θ::Real, a::Real)

Inverse Kerr Metric in Boyer Lindquist (BL) coordinates.

    `r` : Radius
    
    `θ` : Inclination 
    
    `a` : Blackhole spin
"""
kerr_met_uu(r::Real, θ::Real, a::Real) = [ #Eq 1 2105.09440
        -Ξ(r,θ,a)/(Σ(r,θ,a)*Δ(r,a))             0.              -Ξ(r,θ,a)*ω(r,θ,a)/(Σ(r,θ,a)*Δ(r,a))                                0.;
        0.                                      Δ(r,a)/Σ(r,θ,a) 0.                                                                  0.;
        -Ξ(r,θ,a)*ω(r,θ,a)/(Σ(r,θ,a)*Δ(r,a))    0.              Σ(r,θ,a)*csc(θ)^2/Ξ(r,θ,a)-Ξ(r,θ,a)*ω(r,θ,a)^2/(Σ(r,θ,a)*Δ(r,a))    0.;
        0.                                      0.              0.                                                                  1/Σ(r,θ,a)
    ]
"""
    jac_bl2zamo_du(r::Real, θ::Real, a::Real)

Jacobian which converts Boyer-Lindquist (BL) covector on the right to a ZAMO covector

    `r` : Radius
    
    `θ` : Inclination 
    
    `a` : Blackhole spin
"""
jac_bl2zamo_du(r::Real, θ::Real, a::Real) = [# Eq 3.1 1972ApJ...178..347B
    # coords = {t, r, ϕ, θ}
        √(A(r,θ,a)/(Σ(r,θ,a)*Δ(r,a)))   0.                  2*a*r/√(A(r,θ,a)*Σ(r,θ,a)*Δ(r,a))   0.;
        0.                              √(Δ(r,a)/Σ(r,θ,a))  0.                                  0.;
        0.                              0.                  √(Σ(r,θ,a)/A(r,θ,a))*csc(θ)         0.;
        0.                              0.                  0.                                  -1/√Σ(r,θ,a)
    ]
"""
    jac_zamo2zbl_du(r::Real, θ::Real, a::Real)

Jacobian which converts ZAMO covector on the right to a Boyer-Lindquist (BL) covector

    `r` : Radius
    
    `θ` : Inclination 
    
    `a` : Blackhole spin
"""
jac_zamo2bl_du(r::Real, θ::Real, a::Real) = @SMatrix [
    # coords = {t, r, ϕ, θ}
        √((Σ(r,θ,a)*Δ(r,a))/A(r,θ,a))   0.                      -2*a*r*sin(θ)/√(A(r,θ,a)*Σ(r,θ,a))  0.;
        0.                              √(Σ(r, θ, a)/Δ(r, a))   0.                                  0.;
        0.                              0.                      √(A(r,θ,a)/Σ(r,θ,a))*sin(θ)         0.;
        0.                              0.                      0.                                  -√Σ(r,θ,a)
    ]


jac_bl2zamo_ud(r::Real, θ::Real, a::Real) = [#  Eq 3.2 1972ApJ...178..347B
    # coords = {t, r, ϕ, θ}
    √((Σ(r,θ,a)*Δ(r,a))/A(r,θ,a))       0.                      0.                          0.;
    0.                                  √(Σ(r, θ, a)/Δ(r, a))   0.                          0.;
    -(2a*r*sin(θ))/√(A(r,θ,a)*Σ(r,θ,a)) 0.                      √(A(r,θ,a)/Σ(r,θ,a))*sin(θ) 0.;
    0.                                  0.                      0.                          -√Σ(r,θ,a)
    ]

jac_zamo2bl_ud(r::Real, θ::Real, a::Real) = [
    # coords = {t, r, ϕ, θ}
    √(A(r,θ,a)/(Σ(r,θ,a)*Δ(r,a)))       0.                  0.                          0.;
    0.                                  √(Δ(r,a)/Σ(r,θ,a))  0.                          0.;
    2a*r/√(Σ(r,θ,a)*Δ(r,a)*A(r,θ,a))    0.                  √(Σ(r,θ,a)/A(r,θ,a))*csc(θ) 0.;
    0.                                  0.                  0.                          -1/√Σ(r,θ,a)
]

function jac_zamo2fluid_ud(β::Real, θ::Real, φ::Real)
    γ = 1 / √(1 - β^2)

    return [
        γ                   -β*γ*cos(φ)*sin(θ)                              -β*γ*sin(φ)*sin(θ)                          -β*γ*cos(θ);
        -β*γ*cos(φ)*sin(θ)  cos(θ)^2*cos(φ)^2+γ*cos(φ)^2*sin(θ)^2+sin(φ)^2  (γ-1)*cos(φ)*sin(θ)^2*sin(φ)                (γ-1)*cos(θ)*cos(φ)*sin(θ);
        -β*γ*sin(θ)*sin(φ)  (γ-1)*cos(φ)*sin(θ)^2*sin(φ)                    cos(φ)^2+(cos(θ)^2+γ*sin(θ)^2)*sin(φ)^2     (γ-1)*cos(θ)*sin(θ)*sin(φ);
        -β*γ*cos(θ)         (γ-1)*cos(θ)*cos(φ)*sin(θ)                      (γ-1)*cos(θ)*sin(θ)*sin(φ)                  γ*cos(θ)^2+sin(θ)^2         
    ]
end

function penrose_walker(r::Real, θ::Real, a::Real, p_u, f_u)# Eq 6 arXiv:2001.08750v1
    pt, pr, pϕ, pθ = p_u
    ft, fr, fϕ, fθ = f_u

    A = pt*fr - pr*ft + a*sin(θ)^2(pr*fϕ - pϕ*fr)
    B = ((r^2 + a^2)*(pϕ*fθ-pθ*fϕ) - a*(pt*fθ - pθ*ft))*sin(θ)
    κ = (A - B*im)*(r - a*cos(θ)*im)
    return real(κ), imag(κ)
end

function screen_polarisation(κ::Complex, θ::Real, a::Real, α::Real, β::Real)# Eq 31 10.1103/PhysRevD.104.044060
    #TODO: Check which is real and which is imaginary
    κ1 = real(κ)
    κ2 = imag(κ)

    μ = -(α+a*sin(θ))
    fα = (β*κ2 - μ*κ1)/√((κ1^2+κ2^2)*(μ^2+β^2))
    fβ = (β*κ1 + μ*κ2)/√((κ1^2+κ2^2)*(μ^2+β^2))


    return fα, fβ
end

evpa(fα,fβ) = atan(-fα, fβ)


function calcPol(α::Real, β::Real, ri::Real, θs::Real, θo::Real, a::Real, B::Vector{Float64}, βfluid::Vector{Float64}, νr::Bool, θsign::Bool)
    βv::Real = βfluid[1]
    θz::Real = βfluid[2]
    ϕz::Real = βfluid[3]

    ηtemp::Real = η(α, β, θo, a)
    λtemp::Real = λ(α, θo)
    p_bl_d = p_boyer_lindquist_d(ri, θs, a, ηtemp, λtemp, νr, θsign)

    
    p_bl_u = kerr_met_uu(ri, θs, a) * p_bl_d
    p_zamo_u = jac_bl2zamo_ud(ri, θs, a) * p_bl_u
    p_fluid_u = jac_zamo2fluid_ud(βv, θz, ϕz) *  p_zamo_u
    f_fluid_u = similar(p_fluid_u)
    vec = cross(normalize(p_fluid_u[begin+1:end]), B)
    norm = √dot(vec, vec)
    #f_fluid_u = cat([0], cross(normalize(p_fluid_u[begin+1:end]), B), dims=1)
    f_fluid_u = cat([0], vec / norm, dims=1)
    pt = p_zamo_u[1]
    pz = p_zamo_u[4]
    mag = √abs(norm^2 / (pz*pt))
    f_zamo_u = jac_zamo2fluid_ud(-βv, θz, ϕz) * f_fluid_u
    f_bl_u = jac_zamo2bl_ud(ri, θs, a) * f_zamo_u
    κ = penrose_walker(ri, θs, a, p_bl_u, f_bl_u)
    f_screen = screen_polarisation(κ, θo, a, α, β)

    fα = f_screen[1] 
    fβ = f_screen[2] 

    evpatemp = atan(fα, fβ)
    sinϕ = sin(evpatemp) #* mag
    cosϕ = cos(evpatemp) #* mag

    return  sinϕ, cosϕ
end




