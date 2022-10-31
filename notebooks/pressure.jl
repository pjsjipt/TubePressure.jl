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
fL500 = [("bergh65-fig9-R0_50mm-L500mm-Vv300mm3-k1.4-s0-ampl.txt", "bergh65-fig9-R0_50mm-L500mm-Vv300mm3-k1.4-s0-phase.txt"), ("bergh65-fig9-R0_75mm-L500mm-Vv300mm3-k1.4-s0-ampl.txt", "bergh65-fig9-R0_75mm-L500mm-Vv300mm3-k1.4-s0-phase.txt"), ("bergh65-fig9-R1_25mm-L500mm-Vv300mm3-k1.4-s0-ampl.txt", "bergh65-fig9-R1_25mm-L500mm-Vv300mm3-k1.4-s0-phase.txt")]

# ╔═╡ a636902f-87cc-4432-9678-1aaba7efdc28
pberghL500 = [(readdlm(f[1], ',', Float64), readdlm(f[2], ',', Float64)) for f in fL500]

# ╔═╡ 3062ca4c-c3f6-4a92-8162-c6b40c634a2d
freq = 1.0:0.2:300.0

# ╔═╡ 4ce7dc80-b1d7-4e30-b434-5388836bdc4b
R = [0.5, 0.75, 1.25]

# ╔═╡ 895014f8-a512-4028-be26-c2283f13ec5f
pressL500 = let r=R, 
			 L=500e-3, Vv=300e-9, 
             σ=0.0, γ=1.4, Pa=101325, Ta=288.15, 
			 Pr=0.707, μ=1.82e-5
	D = r .* 2e-3
	PressureLine.(D, L, Vv; σ=σ, Ta=Ta, γ=γ, Pa=Pa, Pr=Pr, mu=μ)
end
	
	

# ╔═╡ c4f3db1f-58f8-4e69-8ce7-7fba65eb6ca3
pL500 = [pressurecorr.([pL], freq) for pL in pressL500]

# ╔═╡ 3a733a31-90ce-40e3-a497-a5b7aa5add84
let fig=Figure(), pL=pL500
	ax1 = Axis(fig[1,1], title="L=500 mm, V=300 mm³", xlabel="Frequency (Hz)", 
		ylabel="Amplitude")
	ax2 = Axis(fig[1,2], xlabel="Frequency (Hz)", 
		ylabel="Phase")
	
	for i in 1:3
		p = pL[i]
		bA = pberghL500[i][1]
		bϕ = pberghL500[i][2]
		
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
fL1000 = [("bergh65-fig9-R0_50mm-L1000mm-Vv300mm3-k1.4-s0-ampl.txt", "bergh65-fig9-R0_50mm-L1000mm-Vv300mm3-k1.4-s0-phase.txt"), ("bergh65-fig9-R0_75mm-L1000mm-Vv300mm3-k1.4-s0-ampl.txt", "bergh65-fig9-R0_75mm-L1000mm-Vv300mm3-k1.4-s0-phase.txt"), ("bergh65-fig9-R1_25mm-L1000mm-Vv300mm3-k1.4-s0-ampl.txt", "bergh65-fig9-R1_25mm-L1000mm-Vv300mm3-k1.4-s0-phase.txt")]

# ╔═╡ ae71115a-567d-4300-a88c-5ea3f3004cf3
pberghL1000 = [(readdlm(f[1], ',', Float64), readdlm(f[2], ',', Float64)) for f in fL1000]

# ╔═╡ 60098f06-055b-413d-9ff3-919497727050
pressL1000 = let r=R, 
			 L=1000e-3, Vv=300e-9, 
             σ=0.0, γ=1.4, Pa=101325, Ta=288.15, 
			 Pr=0.707, μ=1.82e-5
	D = r .* 2e-3
	PressureLine.(D, L, Vv; σ=σ, Ta=Ta, γ=γ, Pa=Pa, Pr=Pr, mu=μ)
end
	

# ╔═╡ 65aade7d-d24e-4161-bc59-1151c6b1d6d8
pL1000 = [pressurecorr.([pL], freq) for pL in pressL1000]

# ╔═╡ 445099a5-e96c-41ac-85b7-b29645222fd6
let fig=Figure(), pL=pL1000
	ax1 = Axis(fig[1,1], title="L=1000 mm, V=300 mm³", xlabel="Frequency (Hz)", 
		ylabel="Amplitude")
	ax2 = Axis(fig[1,2], xlabel="Frequency (Hz)", 
		ylabel="Phase")
	
	for i in 1:3
		p = pL[i]
		bA = pberghL1000[i][1]
		bϕ = pberghL1000[i][2]
		
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

# ╔═╡ fa3e7448-a35a-438d-80cc-c17f6de0fd69
pbergh11 = [(readdlm(f[1], ',', Float64), readdlm(f[2], ',', Float64)) for f in fig11]

# ╔═╡ 40eebd93-5097-4f6d-b984-d738e5240fd4
press_fig11 = let r=R, 
			 L=3000e-3, Vv=300e-9, 
             σ=0.0, γ=1.4, Pa=101325, Ta=288.15, 
			 Pr=0.707, μ=1.82e-5
	D = r .* 2e-3
	PressureLine.(D, L, Vv; σ=σ, Ta=Ta, γ=γ, Pa=Pa, Pr=Pr, mu=μ)
end
	

# ╔═╡ df142f9c-91ee-4c8a-a503-22835dc4d7cd
pL_fig11 = [pressurecorr.([pL], freq) for pL in press_fig11]

# ╔═╡ e881f434-9895-400e-bb49-29b689e37df3
let fig=Figure(), pL=pL_fig11
	ax1 = Axis(fig[1,1], title="L=3000 mm, V=300 mm³", xlabel="Frequency (Hz)", 
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
pbergh12 = [(readdlm(f[1], ',', Float64), readdlm(f[2], ',', Float64)) for f in fig12]

# ╔═╡ a9145f09-5a29-4382-98d1-b0e61b33596e
pbergh12[1][2]

# ╔═╡ 82ee815a-a99a-4e90-b5cf-f3639335cc61
press_fig12 = let r=0.75, 
			 L=1000e-3, V=[100,300,1000]*1e-9, 
             σ=0.0, γ=1.4, Pa=101325, Ta=288.15, 
			 Pr=0.707, μ=1.82e-5
	D = r .* 2e-3
	PressureLine.(D, L, V; σ=σ, Ta=Ta, γ=γ, Pa=Pa, Pr=Pr, mu=μ)
end
	

# ╔═╡ dedd369e-0809-4569-ad38-323c337399d4
pL_fig12 = [pressurecorr.([pL], freq) for pL in press_fig12]

# ╔═╡ ee02e6a3-81c4-4824-9aab-0c2d9bc245da
let fig=Figure(), pL=pL_fig12
	ax1 = Axis(fig[1,1], title="L=1000 mm, R=0.75 mm²", xlabel="Frequency (Hz)", 
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
pbergh13 = [(readdlm(f[1], ',', Float64), readdlm(f[2], ',', Float64)) for f in fig13]

# ╔═╡ f397fcaa-abd9-4e06-8cc2-11e85f4052b4
press_fig13 = let r=0.75, 
			 L=1000e-3, V=300e-9, 
             σ=0.0, k=[1.4,1.0], γ=1.4, Pa=101325, Ta=288.15, 
			 Pr=0.707, μ=1.82e-5
	D = r .* 2e-3
	[PressureLine(D, L, V; σ=σ, Ta=Ta, γ=γ, k=k, Pa=Pa, Pr=Pr, mu=μ) for k in k]
end

# ╔═╡ 7d76656c-9656-4bb8-975f-1625c915ec4a
pL_fig13 = [pressurecorr.([pL], freq) for pL in press_fig13]

# ╔═╡ 18be244d-4e02-4a80-97cc-777c258a1b09
let fig=Figure(), pL=pL_fig13
	ax1 = Axis(fig[1,1], title="L=1000 mm", xlabel="Frequency (Hz)", 
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
pbergh14 = [(readdlm(f[1], ',', Float64), readdlm(f[2], ',', Float64)) for f in fig14]

# ╔═╡ a3259d86-dfd1-4863-95a6-8bb713827121
press_fig14 = let r=0.75, 
			 L=1000e-3, V=300e-9, 
             σ=0.0, k=1.4, γ=1.4, atm=[0.5, 1.0, 2.0], Ta=288.15, 
			 Pr=0.707, μ=1.82e-5
	D = r .* 2e-3
	[PressureLine(D, L, V; σ=σ, Ta=Ta, γ=γ, k=k, Pa=101325*Pa, Pr=Pr, mu=μ) for Pa in atm]
end

# ╔═╡ d867ef8e-0f63-4fcc-bd6d-9a5d92411e01
pL_fig14 = [pressurecorr.([pL], freq) for pL in press_fig14]

# ╔═╡ 7abe11fe-660e-446b-a2d7-cffac6c279a7
let fig=Figure(), pL=pL_fig14
	ax1 = Axis(fig[1,1], title="L=1000 mm", xlabel="Frequency (Hz)", 
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


# ╔═╡ Cell order:
# ╠═95d2c466-5910-11ed-2815-0d5b5cb4134b
# ╠═5a31aa00-f652-40b7-989b-0238c8067b87
# ╠═a1b1204a-e892-419d-80bf-901e086d3330
# ╠═49a1a599-4f49-4b39-a5cc-3a8c1f35de53
# ╠═98f0215d-0427-4bfb-a1e3-aadfe08cdccf
# ╠═ca214731-2c87-4dcd-98d1-3c1885673f63
# ╠═a636902f-87cc-4432-9678-1aaba7efdc28
# ╠═3062ca4c-c3f6-4a92-8162-c6b40c634a2d
# ╠═4ce7dc80-b1d7-4e30-b434-5388836bdc4b
# ╠═895014f8-a512-4028-be26-c2283f13ec5f
# ╠═c4f3db1f-58f8-4e69-8ce7-7fba65eb6ca3
# ╠═3a733a31-90ce-40e3-a497-a5b7aa5add84
# ╠═2a17be96-a5ac-4457-b631-35a61c0b583b
# ╠═e3e0ea67-cb2f-4fdd-a111-fdcafdfe878a
# ╠═ae71115a-567d-4300-a88c-5ea3f3004cf3
# ╠═60098f06-055b-413d-9ff3-919497727050
# ╠═65aade7d-d24e-4161-bc59-1151c6b1d6d8
# ╠═445099a5-e96c-41ac-85b7-b29645222fd6
# ╠═2ea662c3-5eb3-4d29-8172-98df45721b47
# ╠═493f38d9-e876-4436-9ded-4f3c4507a498
# ╠═fa3e7448-a35a-438d-80cc-c17f6de0fd69
# ╠═40eebd93-5097-4f6d-b984-d738e5240fd4
# ╠═df142f9c-91ee-4c8a-a503-22835dc4d7cd
# ╠═e881f434-9895-400e-bb49-29b689e37df3
# ╠═4509fed8-9166-46f6-92ff-9fbf1a3ce930
# ╠═2fe51e05-cb45-4de6-b8fc-f4cd91071354
# ╠═01dff97b-63ae-47f8-aff8-9791e8a0101f
# ╠═a9145f09-5a29-4382-98d1-b0e61b33596e
# ╠═82ee815a-a99a-4e90-b5cf-f3639335cc61
# ╠═dedd369e-0809-4569-ad38-323c337399d4
# ╠═ee02e6a3-81c4-4824-9aab-0c2d9bc245da
# ╠═377438c3-dad2-406a-a624-e688ee39ab7c
# ╠═8ad2f108-89c6-4632-8e17-760e3c9f72b2
# ╠═a0b9a783-7497-47c6-a9e3-21eb2100cbb9
# ╠═f397fcaa-abd9-4e06-8cc2-11e85f4052b4
# ╠═7d76656c-9656-4bb8-975f-1625c915ec4a
# ╠═18be244d-4e02-4a80-97cc-777c258a1b09
# ╠═1243d194-10e0-4162-916b-3a7c41d93bcc
# ╠═87d3db82-29b6-4aef-9805-1a721d29723e
# ╠═478e20fe-890d-4cc2-ae0f-4eb17ce9d775
# ╠═a3259d86-dfd1-4863-95a6-8bb713827121
# ╠═d867ef8e-0f63-4fcc-bd6d-9a5d92411e01
# ╠═7abe11fe-660e-446b-a2d7-cffac6c279a7
# ╠═ddcf5078-9f71-4461-a3cd-0e725798b76c
# ╠═b49ab635-2b7d-48ad-adbe-4dc4fc36cc08
