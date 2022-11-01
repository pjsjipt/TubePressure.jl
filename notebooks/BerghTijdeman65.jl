### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 5a31aa00-f652-40b7-989b-0238c8067b87
begin
	import Pkg
	Pkg.activate(".")
end

# ╔═╡ a1b1204a-e892-419d-80bf-901e086d3330
using Revise

# ╔═╡ 49a1a599-4f49-4b39-a5cc-3a8c1f35de53
begin
	using CairoMakie
	using DelimitedFiles
	using TubePressure
end

# ╔═╡ 95d2c466-5910-11ed-2815-0d5b5cb4134b
md"""
# Pressure correction
"""

# ╔═╡ 98f0215d-0427-4bfb-a1e3-aadfe08cdccf
md"""
## Effects of radius with L=500 mm
"""

# ╔═╡ ca214731-2c87-4dcd-98d1-3c1885673f63
fig9 = [("bergh65-fig9-R0_50mm-L500mm-Vv300mm3-k1.4-s0-ampl.txt", "bergh65-fig9-R0_50mm-L500mm-Vv300mm3-k1.4-s0-phase.txt"), ("bergh65-fig9-R0_75mm-L500mm-Vv300mm3-k1.4-s0-ampl.txt", "bergh65-fig9-R0_75mm-L500mm-Vv300mm3-k1.4-s0-phase.txt"), ("bergh65-fig9-R1_25mm-L500mm-Vv300mm3-k1.4-s0-ampl.txt", "bergh65-fig9-R1_25mm-L500mm-Vv300mm3-k1.4-s0-phase.txt")]

# ╔═╡ 8fe1cd47-6717-4391-856a-58e430da800a
function loadpoints(flst, path="points", sep=',')
	f = joinpath.(path, flst)
	return readdlm.(f, sep, Float64)
end

# ╔═╡ a636902f-87cc-4432-9678-1aaba7efdc28
pbergh9 = loadpoints.(fig9, "points")

# ╔═╡ 3062ca4c-c3f6-4a92-8162-c6b40c634a2d
freq = 1.0:0.2:300.0

# ╔═╡ 4ce7dc80-b1d7-4e30-b434-5388836bdc4b
R = [0.5, 0.75, 1.25]

# ╔═╡ 895014f8-a512-4028-be26-c2283f13ec5f
press_fig9 = let r=R, 
			 L=500e-3, Vv=300e-9, 
             σ=0.0, γ=1.4, Pa=101325, Ta=288.15, 
			 Pr=0.707, μ=1.82e-5
	D = r .* 2e-3
	PressureLine.(D, L, Vv, σ, γ; Ta=Ta, gamma=γ, Pa=Pa, Pr=Pr, mu=μ)
end

# ╔═╡ 3a733a31-90ce-40e3-a497-a5b7aa5add84
let fig=Figure(), press=press_fig9
	pL = [correctcoef.(pL, freq) for pL in press]
	ax1 = Axis(fig[1,1], title="Bergh and Tijdeman (1965) figure 9", xlabel="Frequency (Hz)", 
		ylabel="Amplitude")
	ax2 = Axis(fig[1,2], xlabel="Frequency (Hz)", 
		ylabel="Phase")
	
	for i in 1:3
		p = pL[i]
		bA = pbergh9[i][1]
		bϕ = pbergh9[i][2]
		
		A = abs.(p)
		ϕ = phaseangle(p) * 180/π
		lines!(ax1, freq, A, label="R = $(R[i])")
		scatter!(ax1, bA[:,1], bA[:,2])
		lines!(ax2, freq, ϕ)
		scatter!(ax2, bϕ[:,1], bϕ[:,2])
	end
	axislegend(ax1)
	fig
end


# ╔═╡ 2a17be96-a5ac-4457-b631-35a61c0b583b
md"""
## Effects of radius with L=1000 mm
"""

# ╔═╡ e3e0ea67-cb2f-4fdd-a111-fdcafdfe878a
fig10 = [("bergh65-fig9-R0_50mm-L1000mm-Vv300mm3-k1.4-s0-ampl.txt", "bergh65-fig9-R0_50mm-L1000mm-Vv300mm3-k1.4-s0-phase.txt"), ("bergh65-fig9-R0_75mm-L1000mm-Vv300mm3-k1.4-s0-ampl.txt", "bergh65-fig9-R0_75mm-L1000mm-Vv300mm3-k1.4-s0-phase.txt"), ("bergh65-fig9-R1_25mm-L1000mm-Vv300mm3-k1.4-s0-ampl.txt", "bergh65-fig9-R1_25mm-L1000mm-Vv300mm3-k1.4-s0-phase.txt")]

# ╔═╡ ae71115a-567d-4300-a88c-5ea3f3004cf3
pbergh10 = loadpoints.(fig10, "points")

# ╔═╡ 60098f06-055b-413d-9ff3-919497727050
press_fig10 = let r=R, 
			 L=1000e-3, Vv=300e-9, 
             σ=0.0, γ=1.4, Pa=101325, Ta=288.15, 
			 Pr=0.707, μ=1.82e-5
	D = r .* 2e-3
	PressureLine.(D, L, Vv, σ, γ; Ta=Ta, gamma=γ, Pa=Pa, Pr=Pr, mu=μ)
end
	

# ╔═╡ 445099a5-e96c-41ac-85b7-b29645222fd6
let fig=Figure(), press=press_fig10
	pL = [correctcoef.(pL, freq) for pL in press]
	ax1 = Axis(fig[1,1], title="Bergh and Tijdeman (1965) figure 10", xlabel="Frequency (Hz)", 
		ylabel="Amplitude")
	ax2 = Axis(fig[1,2], xlabel="Frequency (Hz)", 
		ylabel="Phase")
	
	for i in 1:3
		p = pL[i]
		bA = pbergh10[i][1]
		bϕ = pbergh10[i][2]
		
		A = abs.(p)
		ϕ = phaseangle(p) * 180/π
		lines!(ax1, freq, A, label="R = $(R[i])")
		scatter!(ax1, bA[:,1], bA[:,2])
		lines!(ax2, freq, ϕ)
		scatter!(ax2, bϕ[:,1], bϕ[:,2])
	end
	axislegend(ax1)
	fig
end

# ╔═╡ 2ea662c3-5eb3-4d29-8172-98df45721b47
md"""
## Effects of radius with L=3000 mm
"""

# ╔═╡ 493f38d9-e876-4436-9ded-4f3c4507a498
fig11 = [("bergh65-fig11-R0_50mm-L3000mm-Vv300mm3-k1.4-s0-ampl.txt", "bergh65-fig11-R0_50mm-L3000mm-Vv300mm3-k1.4-s0-phase.txt"), ("bergh65-fig11-R0_75mm-L3000mm-Vv300mm3-k1.4-s0-ampl.txt", "bergh65-fig11-R0_75mm-L3000mm-Vv300mm3-k1.4-s0-phase.txt"), ("bergh65-fig11-R1_25mm-L3000mm-Vv300mm3-k1.4-s0-ampl.txt", "bergh65-fig11-R1_25mm-L3000mm-Vv300mm3-k1.4-s0-phase.txt")]

# ╔═╡ e2a208d8-c5c3-4ab7-bc83-061a06dc8b43
pbergh11 = loadpoints.(fig11, "points")

# ╔═╡ 40eebd93-5097-4f6d-b984-d738e5240fd4
press_fig11 = let r=R, 
			 L=3000e-3, Vv=300e-9, 
             σ=0.0, γ=1.4, Pa=101325, Ta=288.15, 
			 Pr=0.707, μ=1.82e-5
	D = r .* 2e-3
	PressureLine.(D, L, Vv, σ, γ; Ta=Ta, gamma=γ, Pa=Pa, Pr=Pr, mu=μ)
end
	

# ╔═╡ e881f434-9895-400e-bb49-29b689e37df3
let fig=Figure(), press=press_fig11
	pL = [correctcoef.(pL, freq) for pL in press]
	ax1 = Axis(fig[1,1], title="Bergh and Tijdeman (1965) figure 11", xlabel="Frequency (Hz)", 
		ylabel="Amplitude")
	ax2 = Axis(fig[1,2], xlabel="Frequency (Hz)", 
		ylabel="Phase")
	
	for i in 1:3
		p = pL[i]
		bA = pbergh11[i][1]
		bϕ = pbergh11[i][2]
		
		A = abs.(p)
		ϕ = phaseangle(p) * 180/π
		lines!(ax1, freq, A, label="R = $(R[i])")
		scatter!(ax1, bA[:,1], bA[:,2])
		lines!(ax2, freq, ϕ)
		scatter!(ax2, bϕ[:,1], bϕ[:,2])
	end
	axislegend(ax1)
	fig
end

# ╔═╡ 4509fed8-9166-46f6-92ff-9fbf1a3ce930
md"""
## Effects of volume with R=0.74 mm and L=1000m
"""

# ╔═╡ 2fe51e05-cb45-4de6-b8fc-f4cd91071354
fig12 = [("bergh65-fig12-R0_75mm-L1000mm-Vv0100mm3-k1.4-s0-ampl.txt", "bergh65-fig12-R0_75mm-L1000mm-Vv0100mm3-k1.4-s0-phase.txt"), ("bergh65-fig12-R0_75mm-L1000mm-Vv0300mm3-k1.4-s0-ampl.txt", "bergh65-fig12-R0_75mm-L1000mm-Vv0300mm3-k1.4-s0-phase.txt"),("bergh65-fig12-R0_75mm-L1000mm-Vv1000mm3-k1.4-s0-ampl.txt", "bergh65-fig12-R0_75mm-L1000mm-Vv1000mm3-k1.4-s0-phase.txt")]	

# ╔═╡ 01dff97b-63ae-47f8-aff8-9791e8a0101f
pbergh12 = loadpoints.(fig12, "points")

# ╔═╡ 82ee815a-a99a-4e90-b5cf-f3639335cc61
press_fig12 = let r=0.75, 
			 L=1000e-3, V=[100,300,1000]*1e-9, 
             σ=0.0, γ=1.4, Pa=101325, Ta=288.15, 
			 Pr=0.707, μ=1.82e-5
	D = r .* 2e-3
	PressureLine.(D, L, V, σ, γ, Ta=Ta, gamma=γ, Pa=Pa, Pr=Pr, mu=μ)
end
	

# ╔═╡ ee02e6a3-81c4-4824-9aab-0c2d9bc245da
let fig=Figure(), press=press_fig12
	pL = [correctcoef.(pL, freq) for pL in press]
	ax1 = Axis(fig[1,1], title="Bergh and Tijdeman (1965) figure 12", xlabel="Frequency (Hz)", 
		ylabel="Amplitude")
	ax2 = Axis(fig[1,2], xlabel="Frequency (Hz)", 
		ylabel="Phase")
	V = [100,300,1000]
	for i in 1:3
		p = pL[i]
		bA = pbergh12[i][1]
		bϕ = pbergh12[i][2]
		
		A = abs.(p)
		ϕ = phaseangle(p) * 180/π
		lines!(ax1, freq, A, label="V = $(V[i])")
		scatter!(ax1, bA[:,1], bA[:,2])
		lines!(ax2, freq, ϕ)
		scatter!(ax2, bϕ[:,1], bϕ[:,2])
	end
	axislegend(ax1)
	fig
end

# ╔═╡ 377438c3-dad2-406a-a624-e688ee39ab7c
md"""
## Influence of polytropic constant 
"""

# ╔═╡ 8ad2f108-89c6-4632-8e17-760e3c9f72b2
fig13 = [("bergh65-fig13-R0_75mm-L1000mm-Vv300mm3-k1.4-s0-ampl.txt", "bergh65-fig13-R0_75mm-L1000mm-Vv300mm3-k1.4-s0-phase.txt"), ("bergh65-fig13-R0_75mm-L1000mm-Vv300mm3-k1.0-s0-ampl.txt", "bergh65-fig13-R0_75mm-L1000mm-Vv300mm3-k1.0-s0-phase.txt")]

# ╔═╡ a0b9a783-7497-47c6-a9e3-21eb2100cbb9
pbergh13 = loadpoints.(fig13, "points")

# ╔═╡ f397fcaa-abd9-4e06-8cc2-11e85f4052b4
press_fig13 = let r=0.75, 
			 L=1000e-3, V=300e-9, 
             σ=0.0, k=[1.4,1.0], γ=1.4, Pa=101325, Ta=288.15, 
			 Pr=0.707, μ=1.82e-5
	D = r .* 2e-3
	[PressureLine(D, L, V, σ, k; Ta=Ta, gamma=γ, Pa=Pa, Pr=Pr, mu=μ) for k in k]
end

# ╔═╡ 18be244d-4e02-4a80-97cc-777c258a1b09
let fig=Figure(), press=press_fig13
	pL = [correctcoef.(pL, freq) for pL in press]
	ax1 = Axis(fig[1,1], title="Bergh and Tijdeman (1965) figure 13", xlabel="Frequency (Hz)", 
		ylabel="Amplitude")
	ax2 = Axis(fig[1,2], xlabel="Frequency (Hz)", 
		ylabel="Phase")
	V = [100,300,1000]
	for i in 1:2
		p = pL[i]
		bA = pbergh13[i][1]
		bϕ = pbergh13[i][2]
		γ = [1.4, 1.0]
		A = abs.(p)
		ϕ = phaseangle(p) * 180/π
		lines!(ax1, freq, A, label="γ = $(γ[i])")
		scatter!(ax1, bA[:,1], bA[:,2])
		lines!(ax2, freq, ϕ)
		scatter!(ax2, bϕ[:,1], bϕ[:,2])
	end
	axislegend(ax1)
	fig
end

# ╔═╡ 1243d194-10e0-4162-916b-3a7c41d93bcc
md"""
## Influence of mean pressure on the dynamic response
"""

# ╔═╡ 87d3db82-29b6-4aef-9805-1a721d29723e
fig14 = [("bergh65-fig14-Pa0_5atm-ampl.txt", "bergh65-fig14-Pa0_5atm-phase.txt"), 
         ("bergh65-fig14-Pa1_0atm-ampl.txt", "bergh65-fig14-Pa1_0atm-phase.txt"),
	     ("bergh65-fig14-Pa2_0atm-ampl.txt", "bergh65-fig14-Pa2_0atm-phase.txt")]

# ╔═╡ 478e20fe-890d-4cc2-ae0f-4eb17ce9d775
pbergh14 = loadpoints.(fig14,"points")

# ╔═╡ a3259d86-dfd1-4863-95a6-8bb713827121
press_fig14 = let r=0.75, 
			 L=1000e-3, V=300e-9, 
             σ=0.0, k=1.4, γ=1.4, atm=[0.5, 1.0, 2.0], Ta=288.15, 
			 Pr=0.707, μ=1.82e-5
	D = r .* 2e-3
	[PressureLine(D, L, V, σ, k; Ta=Ta, gamma=γ, Pa=101325*Pa, Pr=Pr, mu=μ) for Pa in atm]
end

# ╔═╡ 7abe11fe-660e-446b-a2d7-cffac6c279a7
let fig=Figure(), press=press_fig14
	pL = [correctcoef.(pL, freq) for pL in press]
	ax1 = Axis(fig[1,1], title="Bergh and Tijdeman (1965) figure 14", xlabel="Frequency (Hz)", 
		ylabel="Amplitude")
	ax2 = Axis(fig[1,2], xlabel="Frequency (Hz)", 
		ylabel="Phase")
	Pa = [0.5,1.0,2.0]
	for i in 1:3
		p = pL[i]
		bA = pbergh14[i][1]
		bϕ = pbergh14[i][2]
		A = abs.(p)
		ϕ = phaseangle(p) * 180/π
		lines!(ax1, freq, A, label="Pₐ = $(Pa[i]) atm")
		scatter!(ax1, bA[:,1], bA[:,2])
		lines!(ax2, freq, ϕ)
		scatter!(ax2, bϕ[:,1], bϕ[:,2])
	end
	axislegend(ax1)
	fig
end

# ╔═╡ ddcf5078-9f71-4461-a3cd-0e725798b76c
md"""
## Influence of R₂ on the dynamic response of a double pressure measuring system
"""


# ╔═╡ b49ab635-2b7d-48ad-adbe-4dc4fc36cc08
fig16 = [("points/bergh65-fig16-R0_50mm-ampl.txt", "points/bergh65-fig16-R0_50mm-phase.txt"), ("points/bergh65-fig16-R0_75mm-ampl.txt", "points/bergh65-fig16-R0_75mm-phase.txt"), ("points/bergh65-fig16-R1_25mm-ampl.txt", "points/bergh65-fig16-R1_25mm-phase.txt")]

# ╔═╡ c4544166-470f-4f57-96db-92124f060dce
pbergh16 = loadpoints.(fig16,".")

# ╔═╡ 11d06d01-73c7-4cbf-ab5e-52ad0ca33b77
press_fig16 = let R1=0.75, R2=[0.5, 0.75, 1.25],
			 L=[0.25, 0.75], V=[0.0,300e-9], 
             σ=0.0, k=1.4, γ=1.4, Pa=101325.0, Ta=288.15, 
			 Pr=0.707, μ=1.82e-5
	D = [[R1 * 2e-3, r*2e-3] for r in R2]
	[PressureLine(d, L, V, σ, k; Ta=Ta, gamma=γ, Pa=Pa, Pr=Pr, mu=μ) for d in D]
end

# ╔═╡ 27421c02-6c13-4914-8a3a-8660075c19ab
let fig=Figure(), press=press_fig16
	pL = [correctcoef.(pL, freq) for pL in press]
	ax1 = Axis(fig[1,1], title="Bergh and Tijdeman (1965) figure 16", xlabel="Frequency (Hz)", 
		ylabel="Amplitude")
	ax2 = Axis(fig[1,2], xlabel="Frequency (Hz)", 
		ylabel="Phase")
	R2 = [0.5,0.75,1.25]
	for i in 1:3
		p = pL[i]
		bA = pbergh16[i][1]
		bϕ = pbergh16[i][2]
		A = abs.(p)
		ϕ = phaseangle(p) * 180/π
		lines!(ax1, freq, A, label="R₂ = $(R2[i]) mm")
		scatter!(ax1, bA[:,1], bA[:,2])
		lines!(ax2, freq, ϕ)
		scatter!(ax2, bϕ[:,1], bϕ[:,2])
	end
	axislegend(ax1)
	fig
end

# ╔═╡ d79b4d2c-56ff-4139-b0c4-62243b31c273
md"""
## Influence of R₂ on the dynamic response of a double pressure measuring system. Different lengths
"""

# ╔═╡ 4a4f8b3b-0774-43c6-9988-21acadd20cc9
fig18 = [("points/bergh65-fig18-R0_50mm-ampl.txt", "points/bergh65-fig18-R0_50mm-phase.txt"), ("points/bergh65-fig18-R0_75mm-ampl.txt", "points/bergh65-fig18-R0_75mm-phase.txt"), ("points/bergh65-fig18-R1_25mm-ampl.txt", "points/bergh65-fig18-R1_25mm-phase.txt")]

# ╔═╡ c594519b-7e22-4a31-93f9-afb7d16002c0
pbergh18 = loadpoints.(fig18,".")

# ╔═╡ 0f9ef798-2be2-4caa-9dbe-d12be4481610
press_fig18 = let R1=0.75, R2=[0.5, 0.75, 1.25],
			 L=[0.75, 0.25], V=[0.0,300e-9], 
             σ=0.0, k=1.4, γ=1.4, Pa=101325.0, Ta=288.15, 
			 Pr=0.707, μ=1.82e-5
	D = [[R1 * 2e-3, r*2e-3] for r in R2]
	[PressureLine(d, L, V, σ, k; Ta=Ta, gamma=γ, Pa=Pa, Pr=Pr, mu=μ) for d in D]
end

# ╔═╡ a5445dee-a141-4e86-ac88-e423d4baf2a5
let fig=Figure(), press=press_fig18
	pL = [correctcoef.(pL, freq) for pL in press]
	ax1 = Axis(fig[1,1], title="Bergh and Tijdeman (1965) figure 18", xlabel="Frequency (Hz)", 
		ylabel="Amplitude")
	ax2 = Axis(fig[1,2], xlabel="Frequency (Hz)", 
		ylabel="Phase")
	R2 = [0.5,0.75,1.25]
	for i in 1:3
		p = pL[i]
		bA = pbergh18[i][1]
		bϕ = pbergh18[i][2]
		A = abs.(p)
		ϕ = phaseangle(p) * 180/π
		lines!(ax1, freq, A, label="R₂ = $(R2[i]) mm")
		scatter!(ax1, bA[:,1], bA[:,2])
		lines!(ax2, freq, ϕ)
		scatter!(ax2, bϕ[:,1], bϕ[:,2])
	end
	axislegend(ax1)
	fig
end

# ╔═╡ b03ae350-da7c-41ae-b9ab-d1b93715d072
md"""
## Several configurations - Figure 19
"""

# ╔═╡ e9f4e4dd-f727-4558-8757-7824eecd6b5f
# Probably a mistake on the article inverted the phases of II and III
fig19 = [("points/bergh65-fig20-I-ampl.txt", "points/bergh65-fig20-I-phase.txt"), 
("points/bergh65-fig19-II-ampl.txt", "points/bergh65-fig19-III-phase.txt"), #
("points/bergh65-fig19-III-ampl.txt", "points/bergh65-fig19-II-phase.txt"), 
("points/bergh65-fig20-IV-ampl.txt", "points/bergh65-fig20-IV-phase.txt")]

# ╔═╡ 3b869e92-a63c-4c33-8567-d48488176523
pbergh19 = loadpoints.(fig19,".")

# ╔═╡ 55286adf-c25a-435b-8177-2e937187f836
press_fig19 = let ra=0.75, rb=0.5
	La = 0.75
	Lb = 0.25
	Da = ra*2e-3
	Db = rb*2e-3
	V = [0.0, 300e-9]
	k = [1.4,1.4]
	σ = [0.0, 0.0]
	
	[PressureLine(Da, 1.0, 300e-9, 0.0,1.4),
	 PressureLine([Da,Db], [La,Lb], V, σ, k),
	 PressureLine([Db,Da], [Lb,La], V, σ, k),
	 PressureLine(Db, 1.0, 300e-9, 0.0,1.4)]
	
end

# ╔═╡ 94ce520b-dc35-4918-a98c-ff4171c31c83
let fig=Figure(), press=press_fig19
	pL = [correctcoef.(pL, freq) for pL in press]
	ax1 = Axis(fig[1,1], title="Bergh and Tijdeman (1965) figure 19", xlabel="Frequency (Hz)", 
		ylabel="Amplitude")
	ax2 = Axis(fig[1,2], xlabel="Frequency (Hz)", 
		ylabel="Phase")
	config = ["I", "II", "III", "IV"]
	for i in 1:4
		p = pL[i]
		bA = pbergh19[i][1]
		bϕ = pbergh19[i][2]
		A = abs.(p)
		ϕ = phaseangle(p) * 180/π
		lines!(ax1, freq, A, label="$(config[i])")
		scatter!(ax1, bA[:,1], bA[:,2])
		lines!(ax2, freq, ϕ)
		scatter!(ax2, bϕ[:,1], bϕ[:,2])
	end
	axislegend(ax1)
	fig
end

# ╔═╡ cf233fa0-a9a1-4042-bcda-94ca80aa131d
md"""
## Several configurations - Figure 20
"""

# ╔═╡ 712b0917-6bf2-4bb9-943f-15624c4dcb9e
# Probably a mistake on the figure. Inverted the phases of II and III

fig20 = [("points/bergh65-fig20-I-ampl.txt", "points/bergh65-fig20-I-phase.txt"), 
("points/bergh65-fig20-II-ampl.txt", "points/bergh65-fig20-III-phase.txt"), 
("points/bergh65-fig20-III-ampl.txt", "points/bergh65-fig20-II-phase.txt"), 
("points/bergh65-fig20-IV-ampl.txt", "points/bergh65-fig20-IV-phase.txt")]

# ╔═╡ 8fe831b2-6a26-4d85-a52e-49fa41cedc5f
pbergh20 = loadpoints.(fig20,".")

# ╔═╡ 01175395-2ef5-4ed8-be00-b35824f2b1bd
press_fig20 = let ra=0.75, rb=0.5
	Da = ra*2e-3
	Db = rb*2e-3
	L = [0.5,0.5]
	V = [0.0, 300e-9]
	k = [1.4,1.4]
	σ = [0.0, 0.0]
	
	[PressureLine([Da,Da], L, V, σ, k),
	 PressureLine([Da,Db], L, V, σ, k),
	 PressureLine([Db,Da], L, V, σ, k),
	 PressureLine([Db,Db], L, V, σ, k)]
	
end

# ╔═╡ 65b1ab55-af35-4a62-a7d5-840801dec4e9
let fig=Figure(), press=press_fig20
	pL = [correctcoef.(pL, freq) for pL in press]
	ax1 = Axis(fig[1,1], title="Bergh and Tijdeman (1965) figure 20", xlabel="Frequency (Hz)", 
		ylabel="Amplitude")
	ax2 = Axis(fig[1,2], xlabel="Frequency (Hz)", 
		ylabel="Phase")
	config = ["I", "II", "III", "IV"]
	for i in 1:4
		p = pL[i]
		bA = pbergh20[i][1]
		bϕ = pbergh20[i][2]
		A = abs.(p)
		ϕ = phaseangle(p) * 180/π
		lines!(ax1, freq, A, label="$(config[i])")
		scatter!(ax1, bA[:,1], bA[:,2])
		lines!(ax2, freq, ϕ)
		scatter!(ax2, bϕ[:,1], bϕ[:,2])
	end
	axislegend(ax1)
	fig
end

# ╔═╡ Cell order:
# ╠═95d2c466-5910-11ed-2815-0d5b5cb4134b
# ╠═5a31aa00-f652-40b7-989b-0238c8067b87
# ╠═a1b1204a-e892-419d-80bf-901e086d3330
# ╠═49a1a599-4f49-4b39-a5cc-3a8c1f35de53
# ╠═98f0215d-0427-4bfb-a1e3-aadfe08cdccf
# ╠═ca214731-2c87-4dcd-98d1-3c1885673f63
# ╠═8fe1cd47-6717-4391-856a-58e430da800a
# ╠═a636902f-87cc-4432-9678-1aaba7efdc28
# ╠═3062ca4c-c3f6-4a92-8162-c6b40c634a2d
# ╠═4ce7dc80-b1d7-4e30-b434-5388836bdc4b
# ╠═895014f8-a512-4028-be26-c2283f13ec5f
# ╟─3a733a31-90ce-40e3-a497-a5b7aa5add84
# ╠═2a17be96-a5ac-4457-b631-35a61c0b583b
# ╠═e3e0ea67-cb2f-4fdd-a111-fdcafdfe878a
# ╠═ae71115a-567d-4300-a88c-5ea3f3004cf3
# ╠═60098f06-055b-413d-9ff3-919497727050
# ╟─445099a5-e96c-41ac-85b7-b29645222fd6
# ╠═2ea662c3-5eb3-4d29-8172-98df45721b47
# ╠═493f38d9-e876-4436-9ded-4f3c4507a498
# ╠═e2a208d8-c5c3-4ab7-bc83-061a06dc8b43
# ╠═40eebd93-5097-4f6d-b984-d738e5240fd4
# ╟─e881f434-9895-400e-bb49-29b689e37df3
# ╠═4509fed8-9166-46f6-92ff-9fbf1a3ce930
# ╠═2fe51e05-cb45-4de6-b8fc-f4cd91071354
# ╠═01dff97b-63ae-47f8-aff8-9791e8a0101f
# ╠═82ee815a-a99a-4e90-b5cf-f3639335cc61
# ╟─ee02e6a3-81c4-4824-9aab-0c2d9bc245da
# ╠═377438c3-dad2-406a-a624-e688ee39ab7c
# ╠═8ad2f108-89c6-4632-8e17-760e3c9f72b2
# ╠═a0b9a783-7497-47c6-a9e3-21eb2100cbb9
# ╠═f397fcaa-abd9-4e06-8cc2-11e85f4052b4
# ╟─18be244d-4e02-4a80-97cc-777c258a1b09
# ╠═1243d194-10e0-4162-916b-3a7c41d93bcc
# ╠═87d3db82-29b6-4aef-9805-1a721d29723e
# ╠═478e20fe-890d-4cc2-ae0f-4eb17ce9d775
# ╠═a3259d86-dfd1-4863-95a6-8bb713827121
# ╟─7abe11fe-660e-446b-a2d7-cffac6c279a7
# ╠═ddcf5078-9f71-4461-a3cd-0e725798b76c
# ╠═b49ab635-2b7d-48ad-adbe-4dc4fc36cc08
# ╠═c4544166-470f-4f57-96db-92124f060dce
# ╠═11d06d01-73c7-4cbf-ab5e-52ad0ca33b77
# ╟─27421c02-6c13-4914-8a3a-8660075c19ab
# ╠═d79b4d2c-56ff-4139-b0c4-62243b31c273
# ╠═4a4f8b3b-0774-43c6-9988-21acadd20cc9
# ╠═c594519b-7e22-4a31-93f9-afb7d16002c0
# ╠═0f9ef798-2be2-4caa-9dbe-d12be4481610
# ╟─a5445dee-a141-4e86-ac88-e423d4baf2a5
# ╠═b03ae350-da7c-41ae-b9ab-d1b93715d072
# ╠═e9f4e4dd-f727-4558-8757-7824eecd6b5f
# ╠═3b869e92-a63c-4c33-8567-d48488176523
# ╠═55286adf-c25a-435b-8177-2e937187f836
# ╟─94ce520b-dc35-4918-a98c-ff4171c31c83
# ╠═cf233fa0-a9a1-4042-bcda-94ca80aa131d
# ╠═712b0917-6bf2-4bb9-943f-15624c4dcb9e
# ╠═8fe831b2-6a26-4d85-a52e-49fa41cedc5f
# ╠═01175395-2ef5-4ed8-be00-b35824f2b1bd
# ╟─65b1ab55-af35-4a62-a7d5-840801dec4e9
