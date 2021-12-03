using LinearAlgebra

function get_λ_η(α, β,θ, a)
    λ = -α * sin(θ)
    η = (α^2 - a^2)*cos(θ)^2+β^2
    return λ, η
end

function get_up_um(λ, η, a)
    Δ_θ  = 1/2*(1-(η+λ^2)/a^2)
    up = Δ_θ + √(Δ_θ^2 +η/a^2)
    um = Δ_θ - √(Δ_θ^2 +η/a^2)
    return up, um
end

function get_radial_roots(λ, η, a)
    A = a^2 - η - λ^2
    B = 2*(η+(λ-a)^2)
    C = -a^2 * η
    P = -A^2 / 12 - C
    Q = -A/3 * ((A/6)^2-C)-B^2/8
    H = -9*Q+√(12*P^3+81*Q^2)
    z = √((-2*(3^(1/3)*P)+2^(1/3)*H^(2/3))/(2*6^(2/3)*H^(1/3) - A/6))
    r1 = -z-√(-A/2-z^2+B/(4*z))
    r2 = -z+√(-A/2-z^2+B/(4*z))
    r3 = z-√(-A/2-z^2+B/(4*z))
    r4 = z+√(-A/2-z^2+B/(4*z))
    return r1, r2, r3, r4
end

@inline function Δ(r,a)
    return r^2 - 2^r +a^2
end

@inline function Ξ(r, a, θ)
    return (r^2+a^2)^2 - Δ(r,a)*a^2*sin(θ)^2
end

@inline function ω(r, a, θ)
    return 2*a*r/Ξ(r,a,θ)
end

@inline function R(r, a, λ, η)
    return (r^2+a^2-a*λ)^2 - Δ(r,a)*(η+(a-λ)^2)
end

function get_Λ(β, χ)
    γ = 1/√(1-β^2)
    Λ = [γ -γ*β*cos(χ) -γ*β*sin(χ) 0;
         -γ*β*cos(χ) (γ-1)*cos(χ)^2+1 (γ-1)*sin(χ)*cos(χ) 0;
         -γ*β*sin(χ) (γ-1)*sin(χ)*cos(χ) (γ-1)*sin(χ)^2+1 0;
         0 0 0 1]
    return Λ
end


