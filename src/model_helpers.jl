using LinearAlgebra
using StaticArrays

@inline function get_λ_η(α, β,θ, a)
    λ = -α * sin(θ)
    η = (α^2 - a^2)*cos(θ)^2+β^2
    return λ, η
end

@inline function get_up_um(λ, η, a)
    Δ_θ  = 1/2*(1-(η+λ^2)/a^2)
    up = Δ_θ + √(Δ_θ^2 +η/a^2)
    um = Δ_θ - √(Δ_θ^2 +η/a^2)
    return up, um
end

function radialroots(λ::Real, η::Real, a)
    T = promote_type(typeof(λ), typeof(η))
    return radialroots(Complex{T}(λ), Complex{T}(η), a)
end


"""
    $(SIGNATURES)
Finds the radial roots using the of the geodesic with
energy-scaled angular momentum `λ` and carter constant `η`.
"""
@fastmath function radialroots(λ::Complex, η::Complex, a)
    A = a^2 - η - λ^2
    B = 2*(η+(λ-a)^2)
    C = -a^2 * η
    P = -A^2 / 12 - C
    Q = -A/3 * ((A/6)^2-C)-B^2/8
    H = -9*Q + √(12*P^3 + 81*Q^2)
    z = √((-2*(3^(1/3)*P)+2^(1/3)*H^(2/3))/(2*6^(2/3)*H^(1/3)) - A/6)
    discp = -A/2 - z^2 + B/(4*z)
    dp = sqrt(discp)
    discm = discp - B/(2*z)
    dm = sqrt(discm)
    r1 = -z - dp
    r2 = -z + dp
    r3 =  z - dm
    r4 =  z + dm
    #roots = SVector(r1,r2,r3,r4)
    #ind = sortperm(real.(roots))
    #println(real(λ), " ", real(η))
    return r1,r2,r3,r4#Tuple(roots[ind])
end

function get_Λ(β, χ)
    s,c = sincos(χ)
    γ = 1/sqrt(1-β^2)

    Λ00 = γ
    Λ01 = -γ*β*c
    Λ03 = -γ*β*s
    Λ11 = (γ-1)*c^2 +1
    Λ13 = (γ-1)*s*c
    Λ33 = (γ-1)*s^2 + 1
    z = zero(typeof(γ))
    o =  one(typeof(γ))
    return @SMatrix[Λ00 Λ01  z  Λ03;
                    Λ01 Λ11  z  Λ13;
                    z   z    o  z
                    Λ03 Λ13  z  Λ33
                   ]
end
