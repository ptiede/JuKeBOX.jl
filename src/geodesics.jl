import Elliptic
using LinearAlgebra, StaticArrays
export get_roots, Gθ, rs, calcPol, η, λ, r_potential, θ_potential, λcrit, ηcrit, ϕ, Ipm, Δ, Σ, A
##
# Follows the Formalism of Gralla & Lupsasca (https://arxiv.org/pdf/1910.12881.pdf)
##

##----------------------------------------------------------------------------------------------------------------------
# Useful functions
##----------------------------------------------------------------------------------------------------------------------

αboundary(a, θs) = a*sin(θs)

function βboundary(α, θo, a, θs) 
    √max((cos(θo)^2-cos(θs)^2)*(α^2-a^2+a^2*cos(θs)^2)/(cos(θs)^2 -1), 0.0)
end



"""
  r_potential(r, η, λ, a)

Radial potential of a kerr blackhole

  `η` : Reduced Carter constant

  `λ` : Reduced azimuthal agular momentum

  `a` : Blackhole spin

  `r` : Boyer Lindquist radius

"""
r_potential(η, λ, a, r) = (r^2+a^2-a*λ)^2-(r^2 -2r + a^2)*(η+(λ-a)^2) # Eq 7 PhysRevD.101.044032

"""
  θ_potential(r, η, λ, a)

Theta potential of a kerr blackhole

  `η` : Reduced Carter constant

  `λ` : Reduced azimuthal agular momentum

  `a` : Blackhole spin

  `θ` : Boyer Lindquist inclination

"""
θ_potential(η, λ, a, θ) = η + a^2*cos(θ)^2 - λ^2*cot(θ)^2

sort_roots(roots) = NTuple{4,ComplexF64}(sort(sort([roots...],by=x->abs(imag(x)),rev=true),by=y->real(y)))

"""
  get_roots(η::Float64, λ::Float64, a::Float64)

Returns roots of r⁴ + (a²-η-λ²)r² + 2(η+(a-λ)²)r - a²η

  `η` : Normalized Carter constant

  `λ` : Normalized angular momentum

  `a` : Blackhole spin
"""
function get_roots(η::Float64, λ::Float64, a::Float64)
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

    return r1, r2, r3, r4
end

Δ(r, a) = r^2 -2r + a^2
Σ(r, θ, a) = r^2 + a^2*cos(θ)^2
A(r, θ, a) = (r^2 + a^2)^2 - a^2*Δ(r, a)*sin(θ)^2
Ξ(r, θ, a) = (r^2+a^2)^2-Δ(r, a)*a^2*sin(θ)^2
ω(r, θ, a) = 2*a*r/ Ξ(r, θ, a)

η(α, β, θo, a) = (α^2 - a^2)*cos(θo)^2 + β^2
λ(α, θo) = -α*sin(θo)

rtildep(a) = 2*(1+Cos(2/3*acos(a)))
rtilden(a) = 2*(1+Cos(2/3*acos(-a)))


"""
    λcrit(r::Complex, a)

Returns λ values on the critical curve associated with a given r.

  `r` : Radius of orbit associated with critical η value

  `a` : Blackhole spin
"""
λcrit(r, a) = a + r/a*(r- 2Δ(r, a)/(r-1))
"""
    ηcrit(r::Complex, a)

Returns η values on the critical curve associated with a given r.

  `r` : Radius of orbit associated with critical η value

  `a` : Blackhole spin
"""
ηcrit(r, a) = (r^3/a^2)*(4*Δ(r,a)/(r-1)^2 - r)


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
  _rs(η, λ, a, τ)

Emission radius for emission that lies outside the photon ring and whose ray intersects the equatorial plane

  `η` : Normalized Carter constant

  `λ` : Normalized angular momentum

  `a` : Angular Momentum

  `τ` : Mino Time
"""
function _rs(η, λ, a, τ)
  ans = 0.0
  num = 0.0
  fo = 0.0
  k = 0.0
  νr = true

  roots = get_roots(η, λ, a)
  rh = 1 + √(1-a^2)
  numreals::Int = (abs(imag(roots[1])) > 1e-10 ? 0 : 2) + (abs(imag(roots[3])) > 1e-10 ? 0 : 2)
  num::ComplexF64 = 0.0
  den::ComplexF64 = 0.0


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

function _rs_case1(roots, rh, τ)
  ans = 0.0
  νr = true
  _, _, r3, r4 = roots
  root_diffs = get_root_diffs(roots...)
  _, r31, r32, r41, r42 = root_diffs

  if r4 >= rh && τ > 2I2r_turn(root_diffs) # invalid case1
    ans = 0.
  elseif r4 < rh && τ > I2r(roots, root_diffs, rh, true) # invalid case2
    ans = 0.
  else
    k = (r32*r41) / (r31*r42)

    fo = Elliptic.F(asin(√(r31/r41)), k)
    X2 = fo - √(r31*r42) * τ / 2
    sn = r41 * Elliptic.Jacobi.sn(X2, k)^2

    ans = (r31*r4 - r3*sn) / (r31 - sn)
    νr = X2 > 0
  end
  return ans::Float64, νr::Bool
end
function _rs_case3(roots, τ)
  ans = 0.0
  νr = true
  r1, r2, _, _ = roots
  root_diffs = get_root_diffs(roots...)
  r21, r31, r32, r41, r42 = root_diffs

  if τ > I3r_full(root_diffs)#(roots, root_diffs, rh)
    ans = 0.
  else
    A = √abs(r32*r42)
    B = √abs(r31*r41)
    k =  real(((A + B)^2 - r21^2)/(4*A*B))

    fo = Elliptic.F(acos((A-B)/(A+B)), k)
    X3  = real(fo - √(A*B)*τ)
    cn = Elliptic.Jacobi.cn(X3, k)
    num = -A*r1 + B*r2 + (A*r1+B*r2)*cn
    den = -A + B + (A+B)*cn

    ans = real(num/den)
    νr = X3 > 0 
  end

    return ans, νr
end
function _rs_case4(roots, rh, τ)
  ans = 0.0
  νr = true

  r1, _, _, r4 = roots
  root_diffs = get_root_diffs(roots...)
  _, r31, r32, r41, r42 = root_diffs


  if τ > I4r(roots, root_diffs, rh)
    ans = 0.
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
  return ans, νr
end


function I2r_turn(root_diffs::NTuple{5})
  _, r31, r32, r41, r42 = root_diffs
  k = r32*r41/(r31*r42)
  return 2/√real(r31*r42)*Elliptic.F(asin(√(r31/r41)), k)
end

function I2r(roots::NTuple{4}, root_diffs::NTuple{5}, rs, isindir)
  _, _, r3, r4 = roots
  _, r31, r32, r41, r42 = root_diffs

  k = r32*r41/(r31*r42)
  x2_s = √((rs-r4)/(rs-r3)*r31/r41)
  if !(-1.0 < x2_s < 1.0); return 0.0; end
  Ir_s = 2.0/√real(r31*r42)*Elliptic.F(asin(x2_s), k)
  Ir_turn = I2r_turn(root_diffs)

  return Ir_turn + (isindir ? Ir_s : -Ir_s)
end

function I3r_full(root_diffs::NTuple{5})
  r21, r31, r32, r41, r42 = map(abs, root_diffs)
  A2 = r32*r42
  B2 = r31*r41
  if A2 < 0.0 || B2 < 0.0; return Inf; end

  A, B = √A2, √B2
  k3 = ((A+B)^2.0 - r21^2.0)/(4A*B)

  temprat = B/A
  x3_turn = real(√((1.0+0.0im - temprat)/(1.0 + temprat)))
  return 2.0/√real(A*B)*Elliptic.F(acos(x3_turn), k3)
end

function I3r(roots::NTuple{4}, root_diffs::NTuple{5}, rs)
  r1, r2, _, _ = roots
  r21, r31, r32, r41, r42 = root_diffs

  A2 = real(r32*r42)
  B2 = real(r31*r41)
  if A2 < 0.0 || B2 < 0.0; return 0.0; end

  A, B = √A2, √B2


  k3 = real(((A+B)^2.0 - r21^2.0)/(4.0A*B))
  temprat = B*(rs-r2)/(A*(rs-r1))
  x3_s = real(√((1.0+0.0im - temprat)/(1.0 + temprat)))
  Ir_s = 2.0/√real(A*B)*Elliptic.F(real(acos(x3_s)), k3)
  Ir_full = I3r_full(root_diffs)

  return abs(Ir_full - Ir_s)
end

function I4r_full(roots::NTuple{4}, root_diffs::NTuple{5})
  r1, _, _, r4 = roots
  _, r31, r32, r41, r42 = root_diffs

  try
    C = √real(r31*r42)
    D = √real(r32*r41)
    k4 = 4C*D/(C+D)^2.0
    a2 = abs(imag(r1))

    k4 = 4*C*D/(C+D)^2.0
    
    go = √max((4a2^2.0 - (C-D)^2.0) / ((C+D)^2.0 - 4a2^2.0), 0.0)
    return 2.0/(C+D)*Elliptic.F(π/2.0 + atan(go), k4) 
  catch e
    return 0.0
  end
end

function I4r(roots::NTuple{4}, root_diffs::NTuple{5}, rs)
  r1, _, _, r4 = roots
  _, r31, r32, r41, r42 = root_diffs

  if real(r32*r41) < 0 || real(r31*r42) < 0
    return 0
  end
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
  Gθ(η, λ, a, θs, θo, isindir::Bool, n::Int64)

Mino time of trajectory between two inclinations for a given screen coordinate

  `α` : Bardeen α

  `β` : Bardeen β 

  `a` : Blackhole angular Momentum

  `θs` : Emission inclination

  `θo` : Observer inclination

  `isindir` : Is the path direct or indirect?

  `n` : nth image in orde of amount of minotime traversed
"""
function Gθ(α, β, a, θs, θo, isindir::Bool, n::Int64)
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

  νθ =  cos(θs) < abs(cos(θo)) ? (θo > θs) ⊻ (n%2==1) : !isindir
  minotime = real(isindir ? (n+1)*Ghat -sign(β)*Go + (νθ ? 1 : -1)*Gs : n*Ghat - sign(β)*Go + (νθ ? 1 : -1)*Gs ) #Sign of Go indicates whether the ray is from the forward cone or the rear cone

  #if cos(θs) < abs(cos(θo))
  #  #if θo < π/2
  #    #minotime = real(isindir ? (n+1)*Ghat - Go - Gs : n*Ghat + Go - Gs ) #Sign of Go indicates whether the ray is from the forward cone or the rear cone
  #  #else
  #  #νθ =  θo < π/2 ? -1 : 1
  #  minotime = real(isindir ? (n+1)*Ghat -sign(β)*Go + (n%2==1 ? -1 : 1)*νθ*Gs : n*Ghat - sign(β)*Go + (n%2==1 ? -1 : 1)*νθ*Gs ) #Sign of Go indicates whether the ray is from the forward cone or the rear cone
  #  #end
  #else
  #  #minotime = real((isindir ? (-cos(θo) > cos(θs) ? (n+1)*Ghat+(Gs + Go) : (n+1)*Ghat-(Gs + Go)) : (β < 0 ? n*Ghat + Gs + Go : n*Ghat + Gs - Go) ))
  #  minotime = real((isindir ? (n+1)*Ghat-(Gs + sign(β)*Go) : n*Ghat + Gs - sign(β)*Go))
  #end

  if (((β < 0) ⊻ (n%2==1)) && cos(θs) > abs(cos(θo)) && !isvortical) || (isvortical && θo >= π/2)
    return Inf, isvortical
  end
  return minotime, isvortical
end 

##----------------------------------------------------------------------------------------------------------------------
#Polarization stuff
##----------------------------------------------------------------------------------------------------------------------
MinkowskiMet() = @SMatrix [-1. 0. 0. 0.; 0. 1. 0. 0.; 0. 0. 1. 0.; 0. 0. 0. 1.]

p_boyer_lindquist_d(r, θ, a, η, λ, νr::Bool, νθ::Bool) = @SVector [-1,(νr ? 1 : -1)*√(r_potential(η, λ, a, r))/Δ(r, a), λ, (νθ ? 1 : -1)*√θ_potential(η, λ, a, θ)]

"""
    kerr_met_uu(r, θ, a)

Inverse Kerr Metric in Boyer Lindquist (BL) coordinates.

    `r` : Radius
    
    `θ` : Inclination 
    
    `a` : Blackhole spin
"""
kerr_met_uu(r, θ, a) = @SMatrix [ #Eq 1 2105.09440
        -Ξ(r,θ,a)/(Σ(r,θ,a)*Δ(r,a))             0.              -Ξ(r,θ,a)*ω(r,θ,a)/(Σ(r,θ,a)*Δ(r,a))                                0.;
        0.                                      Δ(r,a)/Σ(r,θ,a) 0.                                                                  0.;
        -Ξ(r,θ,a)*ω(r,θ,a)/(Σ(r,θ,a)*Δ(r,a))    0.              Σ(r,θ,a)*csc(θ)^2/Ξ(r,θ,a)-Ξ(r,θ,a)*ω(r,θ,a)^2/(Σ(r,θ,a)*Δ(r,a))    0.;
        0.                                      0.              0.                                                                  1/Σ(r,θ,a)
    ]
"""
    jac_bl2zamo_du(r, θ, a)

Jacobian which converts Boyer-Lindquist (BL) covector on the right to a ZAMO covector

    `r` : Radius
    
    `θ` : Inclination 
    
    `a` : Blackhole spin
"""
jac_bl2zamo_du(r, θ, a) = @SMatrix [# Eq 3.1 1972ApJ...178..347B
    # coords = {t, r, ϕ, θ}
        √(A(r,θ,a)/(Σ(r,θ,a)*Δ(r,a)))   0.                  2*a*r/√(A(r,θ,a)*Σ(r,θ,a)*Δ(r,a))   0.;
        0.                              √(Δ(r,a)/Σ(r,θ,a))  0.                                  0.;
        0.                              0.                  √(Σ(r,θ,a)/A(r,θ,a))*csc(θ)         0.;
        0.                              0.                  0.                                  -1/√Σ(r,θ,a)
    ]
"""
    jac_zamo2zbl_du(r, θ, a)

Jacobian which converts ZAMO covector on the right to a Boyer-Lindquist (BL) covector

    `r` : Radius
    
    `θ` : Inclination 
    
    `a` : Blackhole spin
"""
jac_zamo2bl_du(r, θ, a) = @SMatrix [
    # coords = {t, r, ϕ, θ}
        √((Σ(r,θ,a)*Δ(r,a))/A(r,θ,a))   0.                      -2*a*r*sin(θ)/√(A(r,θ,a)*Σ(r,θ,a))  0.;
        0.                              √(Σ(r, θ, a)/Δ(r, a))   0.                                  0.;
        0.                              0.                      √(A(r,θ,a)/Σ(r,θ,a))*sin(θ)         0.;
        0.                              0.                      0.                                  -√Σ(r,θ,a)
    ]


jac_bl2zamo_ud(r, θ, a) = @SMatrix [#  Eq 3.2 1972ApJ...178..347B
    # coords = {t, r, ϕ, θ}
    √((Σ(r,θ,a)*Δ(r,a))/A(r,θ,a))       0.                      0.                          0.;
    0.                                  √(Σ(r, θ, a)/Δ(r, a))   0.                          0.;
    -(2a*r*sin(θ))/√(A(r,θ,a)*Σ(r,θ,a)) 0.                      √(A(r,θ,a)/Σ(r,θ,a))*sin(θ) 0.;
    0.                                  0.                      0.                          -√Σ(r,θ,a)
    ]

jac_zamo2bl_ud(r, θ, a) = @SMatrix [
    # coords = {t, r, ϕ, θ}
    √(A(r,θ,a)/(Σ(r,θ,a)*Δ(r,a)))       0.                  0.                          0.;
    0.                                  √(Δ(r,a)/Σ(r,θ,a))  0.                          0.;
    2a*r/√(Σ(r,θ,a)*Δ(r,a)*A(r,θ,a))    0.                  √(Σ(r,θ,a)/A(r,θ,a))*csc(θ) 0.;
    0.                                  0.                  0.                          -1/√Σ(r,θ,a)
]

function jac_zamo2fluid_ud(β, θ, φ)
    γ = 1 / √(1 - β^2)

    return @SMatrix  [
        γ                   -β*γ*cos(φ)*sin(θ)                              -β*γ*sin(φ)*sin(θ)                          -β*γ*cos(θ);
        -β*γ*cos(φ)*sin(θ)  cos(θ)^2*cos(φ)^2+γ*cos(φ)^2*sin(θ)^2+sin(φ)^2  (γ-1)*cos(φ)*sin(θ)^2*sin(φ)                (γ-1)*cos(θ)*cos(φ)*sin(θ);
        -β*γ*sin(θ)*sin(φ)  (γ-1)*cos(φ)*sin(θ)^2*sin(φ)                    cos(φ)^2+(cos(θ)^2+γ*sin(θ)^2)*sin(φ)^2     (γ-1)*cos(θ)*sin(θ)*sin(φ);
        -β*γ*cos(θ)         (γ-1)*cos(θ)*cos(φ)*sin(θ)                      (γ-1)*cos(θ)*sin(θ)*sin(φ)                  γ*cos(θ)^2+sin(θ)^2         
    ]
end

function penrose_walker(r, θ, a, p_u::AbstractVector, f_u::AbstractVector)# Eq 6 arXiv:2001.08750v1
    pt, pr, pϕ, pθ = p_u
    ft, fr, fϕ, fθ = f_u

    A = pt*fr - pr*ft + a*sin(θ)^2(pr*fϕ - pϕ*fr)
    B = ((r^2 + a^2)*(pϕ*fθ-pθ*fϕ) - a*(pt*fθ - pθ*ft))*sin(θ)
    return (A - B*im)*(r - a*cos(θ)*im)
    #return A*r - B*a*cos(θ), -(A*a*cos(θ) - B*r)

end

function screen_polarisation(κ::Complex, θ, a, α, β)# Eq 31 10.1103/PhysRevD.104.044060
    #TODO: Check which is real and which is imaginary
    κ1 = real(κ)
    κ2 = imag(κ)

    μ = -(α+a*sin(θ))
    fα = (β*κ2 - μ*κ1)/((μ^2+β^2))
    fβ = (β*κ1 + μ*κ2)/((μ^2+β^2))


    return fα, fβ
end

evpa(fα,fβ) = atan(-fα, fβ)


function calcPol(α, β, ri, θs, θo, a, spec_index, B::AbstractArray{Float64}, βfluid::AbstractArray{Float64}, νr::Bool, θsign::Bool)
    βv = βfluid[1]
    θz = βfluid[2]
    ϕz = βfluid[3]

    ηtemp = η(α, β, θo, a)
    λtemp = λ(α, θo)
    p_bl_d = p_boyer_lindquist_d(ri, θs, a, ηtemp, λtemp, νr, θsign)

    
    p_bl_u = kerr_met_uu(ri, θs, a) * p_bl_d
    p_zamo_u = jac_bl2zamo_ud(ri, θs, a) * p_bl_u
    p_fluid_u = jac_zamo2fluid_ud(βv, θz, ϕz) *  p_zamo_u
    vec = cross( (@view p_fluid_u[begin+1:end]) / p_fluid_u[1], B)
    norm = √dot(vec, vec) + eps()
    f_fluid_u = zeros(4)
    f_fluid_u[2:end] .= vec 
    f_zamo_u = jac_zamo2fluid_ud(-βv, θz, ϕz) * f_fluid_u
    f_bl_u = jac_zamo2bl_ud(ri, θs, a) * f_zamo_u
    #κ = penrose_walker(ri, θs, a, p_bl_u, f_bl_u) 
    A = @SMatrix [
      0.0 1.0 0.0 0.0;
      -1.0 0.0 a*sin(θs)^2 0.;
      0. -a*sin(θs)^2 0. 0.;
      0. 0. 0. 0.
    ]
    B = @SMatrix [
        0. 0. 0. -a*sin(θs);
        0. 0. 0. 0;
        0. 0. 0. (ri^2 + a^2);
        a*sin(θs) 0. -(ri^2 + a^2) 0
    ]
    f_temp_d = ((A - B*im)*(ri - a*cos(θs)*im))* (f_bl_u)
    κ = dot(p_bl_u, f_temp_d)
    κ = κ / √(conj(κ)*κ)


    eα, eβ = screen_polarisation(κ, θo, a, α, β) .* (norm^((spec_index+1.)/2))

    #evpatemp = atan(f_screen...)
    #return  sin(evpatemp)/ p_fluid_u[1],  cos(evpatemp) / p_fluid_u[1]
    return eα, eβ, 1/p_fluid_u[1], abs(p_fluid_u[1]/p_fluid_u[4])
end



