module TubePressure

using LinearAlgebra
using FFTW
import SpecialFunctions: besselj0, besselj

export PressureLine, correctcoef, phaseangle, presscorrect, PressureLineSet

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

"""
`PressureLine(::Type{T}, D, L, V, sigma=0.0, k=1.4;
                      gamma=1.4, Ta=288.15, Pa=101325.0,
                      Ra=287.05, Pr=0.707, mu=1.82e-5)`

Creates a model of a pressure line.

A pressure line consists of a sequence of tubes and volumes.
The tubes are assumed to be rigid. The volumes might deform.

The calculation procedure is carried out according to the
technical report

 * Berg, H. and Tijdeman, H., "Theoretical and experimental results for the dynamic response of pressure measuring systems", National Aero- and Astronautical Research Institute, NLR-TR F.238, Amsterdam, January 1965.

This function calculates the pressure correction for systems as shown in the schematic
figure below. `P₀` is the input pressure, there are `n` sections where each section has
a length of `Lᵢ`, diameter `Dᵢ` and volume `Vᵢ`. Each volume can expand becaus of a diaphgragm and this is taken into account by the parameter `σᵢ`.

```
P₀ ----L₁,D₁-------|V₁,P₁|----L₂,D₂----|V₂,P₂|----L₃,D₃---|V₃,P₃|...----Lₙ,Dₙ----|Vₙ,Pₙ|
```
## Arguments

 * `::Type{T}` number type to be used. If not provided `Float64` is assumed.
 * `D` Diameter in meters of each section of the pressure line.
 * `L` Length in meters of each section of the pressure line.
 * `V` Volume at the end of each section.
 * `sigma` Dimensionless increase in transducer for each volume due to diaphragm deflection
 * `k` Polytropic constant for each volume

Each of the arguments above can be a scalar or a vector. If any of them are vectors,
the longest one is assumed to be the number of sections. Every other parameter above
should have this length or 1. If it is 1, this value is assumed for every section.

There are other keyword arguments that specify the properties of the fluid:

 * `gamma`  specific heat ratio
 * `Ta` Undisturbed fluid temperature (K)
 * `Pa` Undisturbed fluid absolute pressure (Pa)
 * `Ra` Ideal gas constant (J/kg.K)
 * `Pr` Prandtl number of the undisturbed fluid
 * `mu` Viscosity of the fluid (Pa.s)

This method returns a `PressureLine` object. This object can be used to calculate
the pressure correction for each frequency. The function [`correctcoef`](@ref)
 calculates the coefficient for a frequency which is a complex number.
The `PressureLine` object can be used as a function to compute the pressure
correction coefficient. The modulus of this number represents the amplification
of the pressure. The method [`phaseangle`](@ref) computes the phase angle in
radians of a sequence of correction coefficients. It is used to plot a
continuous phase angle.

To apply the pressure correction to a pressure time series,
use method [`presscorrect`](@ref). Again the `PressureLine` object can be used to apply the pressure correction directly. If instead of an abstract vector a matrix is
used, the same correction is applied to every pressure, where the dimension along
which time varies is given by keyword argument `dim` (default value 2 - each
time series is a row of the matrix).


## Examples

```julia-repl
julia> # Create a pressure line with a single section 1mm in diameter and 0.7 m in length with a transducer volume of 10mm³.

julia> pline = PressureLine(1e-3, 0.7, 10e-9, Ta=293.15, Pa=93e3);

julia> f = 1.0:300.0;

julia> r = pline.(f); # Equivalent to `correctcoef.(pline, f)`

julia> ϕ = phaseangle(r);

julia> # Let's plot this using Makie

julia> fig=Figure(); ax1=Axis(fig[1,1]); lines!(ax1, f, abs.(r)); ax2=Axis(fig[1,2]); lines!(ax2, ϕ, phaseangle(r)*180/π); fig

julia> # Apply the same pressure correction to every pressure line
       # in measurement given by each row of matrix `P1`

julia> P1 = presscorrect(pline, 1000.0, P1); # Alternatively `P1 = pline(1000.0,P1)`


```

"""
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


"""
`correctcoef(press, f)`

Compute the pressure correction coefficient for frequency `f`.
 This correction coefficient is a complex number.

For more details on the calculations, see [`PressureLine`](@ref).
To calculate the phase angle, use the [`Base.angle`](@ref) and use
[`Base.abs`](@ref) to compute the amplitude correction.

A common operation is to plot the gain and phase angle for a range of
frequencies. The phase might present discontinuities because of the way
method `angle` is defined. The function [`phaseangle`](@ref) tries to compute a continuous phase angle.

## Example
```julia-repl
julia> pline = PressureLine(1e-3, 0.7, 10e-9, Ta=293.15, Pa=93e3);

julia> correctcoef(pline, 300.0)
-0.08575010622823653 + 0.9362983308678333im

julia> pline(300.0)
-0.08575010622823653 + 0.9362983308678333im


```
"""
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

"""
`phaseangle(p)`

Computes a continuous phase angle.


This is just an utility to help plot phase angles.


## Examples

When computing the pressure correction for a sequence of frequencies,
it is commont that at certain frequencies, the quadrant of the correction
coefficient changes and then he have a large change in phase angle:


```julia-repl
julia> pline = PressureLine(1e-3, 0.7, 10e-9, Ta=293.15, Pa=93e3);

julia> r = pline(180); (r, angle(r))
(-0.7672228900344953 - 0.12088663931898555im, -2.985313571575921)

julia> r = pline(190); (r, angle(r))
(-0.7573088763918315 - 0.045963425507916295im, -3.0809739105853797)

julia> r = pline(200); (r, angle(r))
(-0.7500980418456674 + 0.026766917836382593im, 3.1059232297793367)

```

Using the function `phaseangle` we can get a smoother variation of the phase angle:
```julia-repl
julia> pline = PressureLine(1e-3, 0.7, 10e-9, Ta=293.15, Pa=93e3);

julia> f = [180, 190, 200];

julia> r = pline.(f);

julia> ϕ = phaseangle(r)
3-element Vector{Float64}:
 -2.985313571575921
 -3.0809739105853797
 -3.1772620774002496
```


"""
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

(pline::PressureLine{T})(fs, p::AbstractVector{T}) where {T} =
    presscorrect(pline, fs, p)
(pline::PressureLine{T})(fs, p::AbstractMatrix{T}; dim=2) where {T} =
    presscorrect(pline, fs, p, dim=dim)


"""
`presscorrect(line, fs, p)`
`presscorrect(line, fs, p; dim=2)`

Computes the corrected pressure due to the effects of a pressure line.

Uses [`correctcoef`](@ref) to compute the pressure correction coefficient and
applies this to a pressure measurement.

The [`PressureLine`](@ref) and [`PressureLineSet`](@ref) can be used directly
to compute the pressure correction.


## Arguments

 * `line` a `PressureLine` or `PressureLineSet` object characterizing pressure tubing
 * `fs` sampling frequency in Hz
 * `p` Time series - vector or matrix
 * `dim` Dimension of the matrix along which the pressure correction is computed.

## Examples

For examples, see [`PressureLine`](@ref) or [`PressureLineSet`](@ref)
"""
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

"""
`PressureLineSet(lines, chans)`

When measuring pressure, it is very common to use a pressure scanner that measures
several pressures simultaneously. If every pressure line connected to the scanner
is identical, [`PressureLine`](@ref) can be used directly to correct every
scanner channel.

Often, however, different scanner channels have different tubes (different lengths
usually). In this case, you can specify different pressure lines characteristics
to different channels.

## Arguments

 * `lines` Different pressure lines used
 * `chans` Vector specifying which channel is connected to which pressure line

## Example
```julia-repl
julia> pline1 = PressureLine(1e-3, 0.2, 10e-9, Ta=293.15, Pa=93e3);

julia> pline2 = PressureLine(1e-3, 1.0, 10e-9, Ta=293.15, Pa=93e3);

julia> chans = [fill(1,32); fill(2,32)];

julia> pset = PressureLineSet([pline1, pline2], chans)
PressureLineSet{Float64}(PressureLine{Float64}[PressureLine{Float64}([0.0005], [0.2], [1.0e-8], [1.5707963267948966e-7], [0.0], [1.4], 287.05, 293.15, 93000.0, 1.1051863155473898, 1.4, 343.23197767690584, 0.707, 1.82e-5, ComplexF64[0.0 + 0.0im], ComplexF64[0.0 + 0.0im], ComplexF64[0.0 + 0.0im]), PressureLine{Float64}([0.0005], [1.0], [1.0e-8], [7.853981633974482e-7], [0.0], [1.4], 287.05, 293.15, 93000.0, 1.1051863155473898, 1.4, 343.23197767690584, 0.707, 1.82e-5, ComplexF64[0.0 + 0.0im], ComplexF64[0.0 + 0.0im], ComplexF64[0.0 + 0.0im])], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1  …  2, 2, 2, 2, 2, 2, 2, 2, 2, 2])

julia> P1 = pset(1000.0, P0, dim=2);
```



"""
struct PressureLineSet{T}
    lines::Vector{PressureLine{T}}
    chans::Vector{Int}
end

function presscorrect(lines::PressureLineSet{T},fs,p::AbstractMatrix{T};dim=2) where {T}

    N = size(p,dim)
    f = rfftfreq(N,fs)
    N₁ = length(f)
    r = [line.(f) for line in lines.lines]
    P = rfft(p, dim)
    
    for (i,ch) in enumerate(lines.chans)
        rr  = r[ch]
        if dim==1
            for k in 1:N₁
                P[k,i] /= rr[k]
            end
        elseif dim==2
            for k in 1:N₁
                P[i,k] /= rr[k]
            end
        end
    end

    return irfft(P, N, dim)
end

(lines::PressureLineSet{T})(fs, p::AbstractMatrix{T}; dim=2) where {T} =
    presscorrect(lines, fs, p; dim=dim)


    
end

