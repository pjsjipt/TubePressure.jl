module TubePressure

import SpecialFunctions: besselj0, besselj

export PressureLine, pressurecorr, phaseangle

struct PressureLine{T}
    "Radius of each section"
    R::Vector{T}
    "Length of each section"
    L::Vector{T}
    "Volume of each section"
    V::Vector{T}
    "Tube volume of each section"
    Vt::Vector{T}
    "Diaphragm deformation"
    σ::Vector{T}
    "Polytropic exponent for volumes"
    k::Vector{T}
    "Undisturbed density of air"
    ρ::T
    "Isentropic coefficient"
    γ::T
    "Speed of sound"
    a₀::T
    "Prandtl Number"
    Pr::T
    "Viscosity"
    μ::T
    α::Vector{Complex{T}}
    n::Vector{Complex{T}}
    ϕ::Vector{Complex{T}}
end

function PressureLine(D::AbstractVector{T}, L::AbstractVector{T},
                      V::AbstractVector{T}, 
function PressureLine(D::T, L::T, V::T;
                      σ=0.0, Ta=293.15, γ=1.4, k=1.4, Ra=287.05,
                      Pa=101325.0,Pr=0.707, mu=1.82e-5) where {T}
    ρ = Pa / (Ra * Ta)
    R = D/2
    a₀ = sqrt(γ*Ra*Ta)
    Vt = π*R^2 * L
    α = zeros(Complex{T},1)
    n = zeros(Complex{T},1)
    ϕ = zeros(Complex{T},1)
    return PressureLine{T}([R],[L],[V], [Vt], [σ], [k], ρ, γ, a₀, Pr, mu, α, n, ϕ)
end

function pressurecorr(press::PressureLine{T}, f) where {T}

    ω = 2π*f
    s = sqrt(press.Pr)
    N = length(press.R)
    i32 = one(T)*im * sqrt(one(T)*im)
    α = press.α
    ϕ = press.ϕ
    n = press.n
    L = press.L
    R = press.R
    Vt = press.Vt
    V = press.V
    σ = press.σ
    
    for j in 1:N
        α[j] = i32 * R[j] * sqrt(ω*press.ρ/press.μ)
        n[j] = one(T) / ( one(T) + (press.γ - one(T))/press.γ *
            besselj(2,α[j]*s) /besselj0(α[j]*s))
        ϕ[j] = ω/press.a₀ * sqrt( besselj0(α[j]) / besselj(2,α[j]) ) *
            sqrt(press.γ/n[j])
    end
    p = one(T) + zero(T)*im
    rⱼ₊₁ = p
    for j in N:-1:1
        ϕL = ϕ[j]*L[j]
        term1 = cosh(ϕL)
        term2 = V[j] / Vt[j] * (σ[j] + one(T)/press.k[j]) * n[j] * ϕ[j] * L[j] *
            sinh(ϕL)
        if j == N
            term3 = zero(T)*im
        else
            x1 = (R[j+1]/R[j])^2
            x2 = (ϕ[j+1]/ϕ[j])
            x3 = besselj0(α[j]) / besselj0(α[j+1])
            x4 = besselj(2,α[j+1]) / besselj(2, α[j])
            x5 = sinh(ϕL) / sinh(ϕ[j+1]*L[j+1])
            x6 = cosh(ϕ[j+1]*L[j+1]) - one(T)/rⱼ₊₁
            term3 = x1*x2*x3*x4*x5*x6
        end
        rⱼ₊₁ = term1 + term2 + term3
        p *= 1/ rⱼ₊₁
    end
    return p
            
end

function phaseangle(p)
    ϕ = angle.(p)
    np = length(p)
    s = 0.0
    for i = 2:np
        δϕ = (ϕ[i]+s) - ϕ[i-1]
        if δϕ > 6
            s -= 2π
        elseif δϕ < -6
            s += 2π
        end
        ϕ[i] += s
    end
    return ϕ
end

end


