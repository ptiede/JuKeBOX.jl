import FastElliptic
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
  cosθs2 = cos(θs)^2
  √max((cos(θo)^2-cosθs2)*(α^2-a^2(1-cosθs2))/(cosθs2 -1), 0.0)
end

"""
  r_potential(r, η, λ, a)

Radial potential of a kerr blackhole

  `η` : Reduced Carter constant

  `λ` : Reduced azimuthal agular momentum

  `a` : Blackhole spin

  `r` : Boyer Lindquist radius

"""
function r_potential(η, λ, a, r) 
  λ2 = λ^2
  a*(a*(r*(r+2)-η)-4*λ*r)+r*(2*η+2*λ2+r*(-η-λ2+r^2)) # Eq 7 PhysRevD.101.044032
end

"""
  θ_potential(r, η, λ, a)

Theta potential of a kerr blackhole

  `η` : Reduced Carter constant

  `λ` : Reduced azimuthal agular momentum

  `a` : Blackhole spin

  `θ` : Boyer Lindquist inclination

"""
θ_potential(η, λ, a, θ) = η + a^2*cos(θ)^2 - λ^2*cot(θ)^2

"""
  get_roots(η, λ, a)

Returns roots of r⁴ + (a²-η-λ²)r² + 2(η+(a-λ)²)r - a²η

  `η` : Normalized Carter constant

  `λ` : Normalized angular momentum

  `a` : Blackhole spin
"""
function get_roots(η, λ, a)
    a2 = a*a
    A = a2 - η - λ*λ
    B = 2(η + (λ-a)^2)
    C = -a2*η   

    P = -A*A / 12 - C
    Q = -A/3*(A*A/36 +0im - C) - B*B/8

    Δ3 = -4*P^3 - 27*Q^2
    ωp = (-Q/2 + sqrt(-Δ3/108)+ 0im)^(1/3)

    #C = ((-1+0im)^(2/3), (-1+0im)^(4/3), 1) .* ωp
    C = (-0.4999999999999999 + 0.8660254037844387im, -0.5000000000000002 - 0.8660254037844385im, 1) .* ωp

    v = -P ./ ((3+0im) .* C)
    
    ξ0 = argmax(real, (C .+ v))  - A/3.

    r1 = (-sqrt(2ξ0) - sqrt(-(2A + 2ξ0 - (√2B)/(√ξ0))))/2
    r2 = (-sqrt(2ξ0) + sqrt(-(2A + 2ξ0 - (√2B)/(√ξ0))))/2
    r3 = (sqrt(2ξ0) - sqrt(-(2A + 2ξ0 + (√2B)/(√ξ0))))/2
    r4 = (sqrt(2ξ0) + sqrt(-(2A + 2ξ0 + (√2B)/(√ξ0))))/2

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

get_root_diffs(r1, r2, r3, r4) = r2 - r1, r3 - r1, r3 - r2, r4 - r1, r4 - r2


function rs(α, β, θs, θo, a, isindir, n) 

    if cos(θs) > abs(cos(θo))
        αmin = αboundary(a, θs)
        βbound = (abs(α) >= αmin - eps() ? βboundary(α, θo, a, θs) : 0.)
        if abs(β) + eps() < βbound
            return 0., true, 4
        end
    end
    ηtemp = η(α, β, θo, a)
    λtemp = λ(α, θo)
    τ, _, _ = _Gθ(sign(β), θs, θo, a, isindir, n, ηtemp, λtemp)
    if τ != Inf
      return _rs(ηtemp, λtemp, a, τ)
    else 
      return (0., true, 4) 
    end
end

function rs_mask(n_init, α, β, θs, θo, a, isindir) 
    if abs(cos(θs)) > abs(cos(θo))
        αmin = αboundary(a, θs)
        βbound = (abs(α) >= αmin ? βboundary(α, θo, a, θs) : 0.)
        abs(β) < βbound && return (0., true, 4), true
    end
    ηtemp = η(α, β, θo, a)
    λtemp = λ(α, θo)
    #τ = Gθ(α, β, a, θs, θo, isindir, n)[1]
    τ0, τhat, _ = _Gθ(sign(β), θs, θo, a, isindir, n_init, ηtemp, λtemp)
    if τ0 != Inf
      τ = τ0 + τhat
      return _rsmask(ηtemp, λtemp, a, τ0, τ)
    else
      return (0., true, 4), true
    end
end


function rs_cumulative!(rsvals, n_init, α, β, θs, θo, a, isindir) 
    n = length(rsvals)
    if abs(cos(θs)) > abs(cos(θo))
        αmin = αboundary(a, θs)
        βbound = (abs(α) >= αmin ? βboundary(α, θo, a, θs) : 0.)
        if abs(β) < βbound
          for i in 1:n
            rsvals[i] = (0., true, 4)
          end
        end
    end
    ηtemp = η(α, β, θo, a)
    λtemp = λ(α, θo)
    #τ = Gθ(α, β, a, θs, θo, isindir, n)[1]
    τ0, τhat, _ = _Gθ(sign(β), θs, θo, a, isindir, 0, ηtemp, λtemp)
    if τ0 != Inf
      τ = τ0 + n_init*τhat
      currrs = (1000., true, 4)
      for i in 1:n
        if currrs[1] != 0.0
          currrs = _rs(ηtemp, λtemp, a, τ)
          rsvals[i] = currrs
          τ += τhat
        else
          rsvals[i] = (0., true, 4)
        end
      end
    else
      for i in 1:n
        rsvals[i] = (0., true, 4)
      end
    end
    return nothing
end

function _rsmask(η, λ, a, τ0, τ)
  ans = 0.0
  ansmask = 0.0
  νr = true
  νrmask = true

  roots = get_roots(η, λ, a)
  rh = 1 + √(1-a*a)
  numreals = (abs(imag(roots[1])) > 1e-10 ? 0 : 2) + (abs(imag(roots[3])) > 1e-10 ? 0 : 2)

  if numreals == 4 #case 1 & 2
    ans, νr = _rs_case1(real.(roots), rh, τ0)
    ansmask, νrmask = _rs_case1(real.(roots), rh, τ)

  elseif numreals == 2 #case3
    if abs(imag(roots[4])) < 1e-10
     roots= (roots[1], roots[4],roots[2],roots[3])
    end

    ans, νr = _rs_case3(roots, τ0)
    ansmask, νrmask = _rs_case3(roots, τ)
  else #case 4
    ans, νr = _rs_case4(roots, rh, τ0)
    ansmask, νrmask = _rs_case4(roots, rh, τ)
    end
  return (ans, νr, numreals), (ansmask <= rh)
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
  νr = true

  roots = get_roots(η, λ, a)
  rh = 1 + √(1-a^2)
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

function _rs_case1(roots, rh, τ)
  _, _, r3, r4 = roots
  root_diffs = get_root_diffs(roots...)
  _, r31, r32, r41, r42 = root_diffs

  r4 >= rh && τ > 2I2r_turn(root_diffs) && return 0.0, true # invalid case1
  r4 < rh && τ > I2r(roots, root_diffs, rh, true) && return 0.0, true# invalid case2

  k = (r32*r41) / (r31*r42)
  fo = FastElliptic.F(asin(√(r31/r41)), k)
  X2 = fo - √(r31*r42) * τ / 2
  sn = r41 * FastElliptic.sn(X2, k)^2
  return (r31*r4 - r3*sn) / (r31 - sn), X2 > 0
end
function _rs_case3(roots, τ)
  r1, r2, _, _ = roots
  root_diffs = get_root_diffs(roots...)
  r21, r31, r32, r41, r42 = root_diffs

  τ > I3r_full(root_diffs) && return 0.0, true

  A = √abs(r32*r42)
  B = √abs(r31*r41)
  k =  real(((A + B)^2 - r21^2)/(4*A*B))

  fo = FastElliptic.F(acos((A-B)/(A+B)), k)
  X3  = real(fo - √(A*B)*τ)
  cn = FastElliptic.cn(X3, k)
  num = -A*r1 + B*r2 + (A*r1+B*r2)*cn
  den = -A + B + (A+B)*cn

  return real(num/den), X3 > 0 
end
function _rs_case4(roots, rh, τ)
  r1, _, _, r4 = roots
  root_diffs = get_root_diffs(roots...)
  _, r31, r32, r41, r42 = root_diffs

  τ > I4r(roots, root_diffs, rh) && return 0.0, true
  a2 = abs(imag(r1))
  b1 = real(r4)
  C = √real(r31*r42)
  D = √real(r32*r41)
  k4 = 4*C*D/(C+D)^2
  
  go = √(4a2^2 - (C-D)^2) / ((C+D)^2 - 4a2^2)
  fo = 2/(C+D)*FastElliptic.F(π/2 + atan(go), k4) 
  X4 =  (C+D)/2*(fo - τ)
  num = go - FastElliptic.sc(X4, k4)
  den = 1 + go*FastElliptic.sc(X4, k4)

  return -(a2*num/den + b1), X4 > 0 
end

function I2r_turn(root_diffs::NTuple{5})
  _, r31, r32, r41, r42 = root_diffs
  k = r32*r41/(r31*r42)
  return 2/√real(r31*r42)*FastElliptic.F(asin(√(r31/r41)), k)
end

function I2r(roots::NTuple{4}, root_diffs::NTuple{5}, rs, isindir)
  _, _, r3, r4 = roots
  _, r31, r32, r41, r42 = root_diffs

  k = r32*r41/(r31*r42)
  x2_s = √((rs-r4)/(rs-r3)*r31/r41)
  !(-1 < x2_s < 1) && return 0.

  Ir_s = 2/√real(r31*r42)*FastElliptic.F(asin(x2_s), k)
  Ir_turn = I2r_turn(root_diffs)

  if isindir 
    return Ir_turn + Ir_s
  else
    return Ir_turn - Ir_s
  end
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
  return 2/√real(A*B)*FastElliptic.F(acos(x3_turn), k3)
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
  Ir_s = 2/√real(A*B)*FastElliptic.F(real(acos(x3_s)), k3)
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
    return 2/(C+D)*FastElliptic.F(π/2 + atan(go), k4) 
  catch e
    return 0
  end
end

function I4r(roots, root_diffs, rs)
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
  Ir_s = 2/(C+D)*FastElliptic.F(atan(x4_s) + atan(go), k4) 
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
Gθ(α, β, a, θs, θo, isindir, n) = _Gθ(sign(β), θs, θo, a, isindir, n, η(α, β, θo, a), λ(α, θo))

function _Gθ(signβ, θs, θo, a, isindir, n, η, λ)
  Go, Gs, Ghat, minotime, isvortical = 0., 0., 0., 0., η < 0.


  isincone = abs(cos(θs)) < abs(cos(θo))
  if isincone && (isindir != ((signβ > 0) ⊻ (θo > π/2) ))
    return Inf, Inf, isvortical
  end
  if ((((signβ < 0) ⊻ (θs > π/2)) ⊻ (n%2==1)) && !isincone && !isvortical) || (isvortical && ((θo >= π/2) ⊻ (θs > π/2)))
    return Inf, Inf, isvortical
  end

  Δθ = 1/2*(1 - (η + λ^2)/a^2)
  up = Δθ + √(Δθ^2 + η/a^2)
  um = Δθ - √(Δθ^2 + η/a^2)
  m = up/um
  k = m

  #isvortical = η < 0.
  args = 0 
  argo = 0 
  k = 0
  if isvortical 
    args = (cos(θs)^2 - um)/(up-um)
    argo = (cos(θo)^2 - um)/(up-um)
    k = 1. - m 
    if (!(0. < argo < 1.) ||  !(0. < args <  1.))
      return Inf, Ghat, isvortical
    end
    tempfac = 1/√abs(um*a^2)
    Go = ((θs > π/2) ? -1 : 1)*tempfac*FastElliptic.F(asin(√argo), k)
    Gs = ((θs > π/2) ? -1 : 1)*tempfac*FastElliptic.F(asin(√args), k)
    Ghat = 2tempfac*FastElliptic.K(k)
  else 
    args = cos(θs)/√(up)
    argo = cos(θo)/√(up)
    k = m
    if !(-1 < args < 1) || !(-1 < argo < 1)
     return Inf, Ghat, isvortical
    end
    tempfac = 1/√abs(um*a^2)
    Go = tempfac*FastElliptic.F(asin(argo), k)
    Gs = tempfac*FastElliptic.F(asin(args), k)
    Ghat = 2tempfac*FastElliptic.K(k)

  end
 
  νθ =  isincone ? (n%2==1) ⊻ (θo > θs) : !isindir ⊻ (θs > π/2) 
  minotime = real(isindir ? (n+1)*Ghat -signβ*Go + (νθ ? 1 : -1)*Gs : n*Ghat - signβ*Go + (νθ ? 1 : -1)*Gs ) #Sign of Go indicates whether the ray is from the forward cone or the rear cone
  return minotime, Ghat, isvortical
end 

##----------------------------------------------------------------------------------------------------------------------
#Polarization stuff
##----------------------------------------------------------------------------------------------------------------------
MinkowskiMet() = @SMatrix [-1. 0. 0. 0.; 0. 1. 0. 0.; 0. 0. 1. 0.; 0. 0. 0. 1.]

p_boyer_lindquist_d(r, θ, a, η, λ, νr::Bool, νθ::Bool) = @SVector [-1,(νr ? 1 : -1)*√abs(r_potential(η, λ, a, r))/Δ(r, a), λ, (νθ ? 1 : -1)*√abs(θ_potential(η, λ, a, θ))]

"""
    kerr_met_uu(r, θ, a)

Inverse Kerr Metric in Boyer Lindquist (BL) coordinates.

    `r` : Radius
    
    `θ` : Inclination 
    
    `a` : Blackhole spin
"""
function kerr_met_uu(r, θ, a) 
  Ξt = Ξ(r,θ,a)
  Σt = Σ(r,θ,a)
  Δt = Δ(r,a)
  ωt = ω(r,θ,a)

  return @SMatrix [ #Eq 1 2105.09440
    -Ξt/(Σt*Δt)       0.      -Ξt*ωt/(Σt*Δt)                    0.;
    0.                Δt/Σt   0.                                0.;
    -Ξt*ωt/(Σt*Δt)    0.      Σt*csc(θ)^2/Ξt-Ξt*ωt^2/(Σt*Δt)    0.;
    0.                0.      0.                                1/Σt
  ]
end

"""
    jac_bl2zamo_du(r, θ, a)

Jacobian which converts Boyer-Lindquist (BL) covector on the right to a ZAMO covector

    `r` : Radius
    
    `θ` : Inclination 
    
    `a` : Blackhole spin
"""
function jac_bl2zamo_du(r, θ, a)
  # coords = {t, r, ϕ, θ}
  Σt = Σ(r,θ,a)
  Δt = Δ(r,a)
  At = A(r,θ,a)
  return @SMatrix [# Eq 3.1 1972ApJ...178..347B
    √(At/(Σt*Δt))   0.          2*a*r/√(At*Σt*Δt)   0.;
    0.              √(Δt/Σt)    0.                  0.;
    0.              0.          √(Σt/At)*csc(θ)     0.;
    0.              0.          0.                  -1/√Σt
  ]
end

"""
    jac_zamo2zbl_du(r, θ, a)

Jacobian which converts ZAMO covector on the right to a Boyer-Lindquist (BL) covector

    `r` : Radius
    
    `θ` : Inclination 
    
    `a` : Blackhole spin
"""
function jac_zamo2bl_du(r, θ, a)    
  Σt = Σ(r,θ,a)
  Δt = Δ(r,a)
  At = A(r,θ,a)

  return @SMatrix [
    # coords = {t, r, ϕ, θ}
      √((Σt*Δt)/At)   0.                      -2*a*r*sin(θ)/√(At*Σt)  0.;
      0.                  √(Σ(r, θ, a)/Δ(r, a))   0.                      0.;
      0.                  0.                      √(At/Σt)*sin(θ)         0.;
      0.                  0.                      0.                      -√Σt
    ]
end


function jac_bl2zamo_ud(r, θ, a)
  Σt = Σ(r,θ,a)
  Δt = Δ(r,a)
  At = A(r,θ,a)

  return @SMatrix [#  Eq 3.2 1972ApJ...178..347B
    # coords = {t, r, ϕ, θ}
    √((Σt*Δt)/At)           0.                      0.               0.;
    0.                      √(Σt/Δt)   0.               0.;
    -(2a*r*sin(θ))/√(At*Σt) 0.                      √(At/Σt)*sin(θ)  0.;
    0.                      0.                      0.               -√Σt
  ]
end

function jac_zamo2bl_ud(r, θ, a) 
  Σt = Σ(r,θ,a)
  Δt = Δ(r,a)
  At = A(r,θ,a)


  return @SMatrix [
    # coords = {t, r, ϕ, θ}
    √(At/(Σt*Δt))       0.                  0.                          0.;
    0.                  √(Δt/Σt)            0.                          0.;
    2a*r/√(Σt*Δt*At)    0.                  √(Σt/At)*csc(θ)             0.;
    0.                  0.                  0.                          -1/√Σt
]
end

function jac_zamo2fluid_ud(β, θ, φ)
    γ = 1 / √(1 - β^2)
    sinφ = sin(φ)
    cosφ = cos(φ) 
    sinθ = sin(θ)
    cosθ = cos(θ) 

    return @SMatrix  [
        γ               -β*γ*cosφ*sinθ                         -β*γ*sinφ*sinθ                        -β*γ*cosθ;
        -β*γ*cosφ*sinθ  cosθ^2*cosφ^2+γ*cosφ^2*sinθ^2+sinφ^2   (γ-1)*cosφ*sinθ^2*sinφ                (γ-1)*cosθ*cosφ*sinθ;
        -β*γ*sinθ*sinφ  (γ-1)*cosφ*sinθ^2*sinφ                 cosφ^2+(cosθ^2+γ*sinθ^2)*sinφ^2     (γ-1)*cosθ*sinθ*sinφ;
        -β*γ*cosθ       (γ-1)*cosθ*cosφ*sinθ                   (γ-1)*cosθ*sinθ*sinφ                  γ*cosθ^2+sinθ^2         
    ]
end

function penrose_walker(r, θ, a, p_u::AbstractVector, f_u::AbstractVector)# Eq 6 arXiv:2001.08750v1
    pt, pr, pϕ, pθ = p_u
    ft, fr, fϕ, fθ = f_u
    sinθ = sin(θ)
    cosθ = cos(θ)

    A = pt*fr - pr*ft + a*sinθ^2(pr*fϕ - pϕ*fr)
    B = ((r^2 + a^2)*(pϕ*fθ-pθ*fϕ) - a*(pt*fθ - pθ*ft))*sinθ
    return A*r - B*a*cosθ, -(A*a*cosθ - B*r)
end

function screen_polarisation(κ::Complex, θ, a, α, β)# Eq 31 10.1103/PhysRevD.104.044060
    #TODO: Check which is real and which is imaginary
    κ1 = real(κ)
    κ2 = imag(κ)

    μ = -(α+a*sin(θ))
    fα = (β*κ2 - μ*κ1)/(μ^2+β^2)
    fβ = (β*κ1 + μ*κ2)/(μ^2+β^2)


    return fα, fβ
end

evpa(fα,fβ) = atan(-fα, fβ)



function calcPol(α, β, ri, θs, θo, a, cross_spec_index, magfield::AbstractArray{Float64}, βfluid::AbstractArray{Float64}, νr::Bool, θsign::Bool)
    βv = βfluid[1]
    θz = βfluid[2]
    ϕz = βfluid[3]

    ηtemp = η(α, β, θo, a)
    λtemp = λ(α, θo)
    p_bl_d = p_boyer_lindquist_d(ri, θs, a, ηtemp, λtemp, νr, θsign)

    p_bl_u = kerr_met_uu(ri, θs, a) * p_bl_d
    p_zamo_u = jac_bl2zamo_ud(ri, θs, a) * p_bl_u
    p_fluid_u = jac_zamo2fluid_ud(βv, θz, ϕz) *  p_zamo_u
    magfieldx, magfieldy, magfieldz = magfield
    _, p_fluid_ux, p_fluid_uy, p_fluid_uz = p_fluid_u ./ p_fluid_u[1]
    vec = @SVector[p_fluid_uy*magfieldz - p_fluid_uz*magfieldy, p_fluid_uz*magfieldx - p_fluid_ux*magfieldz, p_fluid_ux*magfieldy - p_fluid_uy*magfieldx]  
    norm = √dot(vec, vec) + eps()
    f_fluid_u = SVector(zero(eltype(vec)), vec[1], vec[2], vec[3])
    #f_fluid_u[2:end] .= vec 
    f_zamo_u = jac_zamo2fluid_ud(-βv, θz, ϕz) * f_fluid_u
    f_bl_u = jac_zamo2bl_ud(ri, θs, a) * f_zamo_u
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


    eα, eβ = screen_polarisation(κ, θo, a, α, β) .* (norm^((cross_spec_index+1.)/2))

    return eα, eβ, 1/p_fluid_u[1], abs(p_fluid_u[1]/p_fluid_u[4])
end



