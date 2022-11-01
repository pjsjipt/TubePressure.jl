module TubePressure

using LinearAlgebra
using FFTW
import SpecialFunctions: besselj0, besselj

export PressureLine, correctcoef, phaseangle, presscorrect

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
    "Fluid gas constant"
    Ra::T
    "Undisturbed temperature (K)"
    Ta::T
    "Undisturbed pressure (Pa)"
    Pa
    "Undisturbed density of air"
    ρ::T
    "Isentropic coefficient"
    γ::T
    "Speed of sound"
    a₀::T
    "Prandtl Number"
    Pr::T
    "Viscosity (Pa.s)"
    μ::T
    "Temporary storage"
    α::Vector{Complex{T}}
    n::Vector{Complex{T}}
    ϕ::Vector{Complex{T}}
end

Base.Broadcast.broadcastable(p::PressureLine) = Ref(p)

function PressureLine(::Type{T}, D, L, V, sigma=0.0, k=1.4;
                      gamma=1.4, Ta=288.15, Pa=101325.0,
                      Ra=287.05, Pr=0.707, mu=1.82e-5) where {T}
    nD = length(D)
    nL = length(L)
    nV = length(V)
    nσ = length(sigma)
    nk = length(k)
    
    nsec = max(nD, nL, nV, nσ, nk)

    nD ∉ (1,nsec) && error("Number of tube sections inconsistent")
    nL ∉ (1,nsec) && error("Number of tube sections inconsistent")
    nV ∉ (1,nsec) && error("Number of tube sections inconsistent")
    nσ ∉ (1,nsec) && error("Number of tube sections inconsistent")
    nk ∉ (1,nsec) && error("Number of tube sections inconsistent")

    if nD==1
        R = fill(T(0.5*D[1]), nsec)
    else
        R = T.(D) ./ 2
    end

    if nL==1
        L = fill(T(L[1]), nsec)
    end
    
    if nV==1
        V = fill(T(V[1]), nsec)
    end

    if nσ==1
        sigma = fill(sigma[1], nsec)
    end

    if nk==1
        k = fill(k[1], nsec)
    end

    α = zeros(T, nsec)
    n = zeros(T, nsec)
    ϕ = zeros(T, nsec)

    ρ = T(Pa / (Ra * Ta))
    a₀ = T(sqrt(gamma*Ra*Ta))
    Vt = T.(π .* R .^ 2 .* L)
    return PressureLine{T}(R, T.(L), T.(V), Vt, T.(sigma), T.(k),
                                 T(Ra), T(Ta), T(Pa), T(ρ), T(gamma),
                                 T(a₀), T(Pr), T(mu), α, n, ϕ)
end



    
    
PressureLine(D, L, V, sigma=0.0, k=1.4;
             gamma=1.4, Ta=288.15, Pa=101325.0,
             Ra=287.05, Pr=0.707, mu=1.82e-5)  = PressureLine(Float64, D, L, V,
                                                              sigma, k; gamma=gamma,
                                                              Ta=Ta,Pa=Pa,Ra=Ra,
                                                              Pr=Pr,mu=mu)



function correctcoef(press::PressureLine{T}, f) where {T}

    if f==0
        return one(T)+zero(T)*im
    end
    
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

(p::PressureLine)(f) = correctcoef(p, f)

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


function presscorrect(line::PressureLine{T}, fs, p::AbstractVector{T}) where {T}
    N = length(p)
    f = rfftfreq(N, fs)
    return irfft(rfft(p) ./ line.(f), N)
end

function presscorrect(line::PressureLine{T},fs,p::AbstractMatrix{T}; dim=2) where {T}
    N = size(p,dim)
    f = rfftfreq(N,fs)
    r = line.(f)
    P = rfft(p, dim)
  
    if dim==1
        return irfft(Diagonal(1 ./ r) * P, N, 1)
    elseif dim==2
        return irfft(P * Diagonal(1 ./ r), N, 2)
    end
    
end


struct PressureLineSet{T}
    lines::Vector{PressureLine{T}}
    chans::Vector{Int}
end

end


