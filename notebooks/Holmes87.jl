### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 5c98f2bc-9076-4e10-8b84-cbdaf34e65ca
begin
	import Pkg
	Pkg.activate(".")
end

# ╔═╡ 3db66fc0-6a59-473c-9059-c00251dcacaf
begin
	using CairoMakie
	using DelimitedFiles
	using TubePressure
end

# ╔═╡ fffad682-611f-11ed-3fdf-93926ecdc1a5
md"""
# Reproducing the data from Holmes

In this notebook, we try to reproduce the graphs in the paper

 * Holmes, J. D. and Lewis, R. E., "Optimization of dynamic-pressure-measurement systems. I. Single point measurements, Journal of Wind Engineering and Industrial Aerodynamics, 25, p249-273, 1987.

"""

# ╔═╡ 6c3a2e0a-abe9-4d90-9cd1-b3e0b8373cef
function loadpoints(flst, path="points", sep=',')
	f = joinpath.(path, flst)
	return readdlm.(f, sep, Float64)
end

# ╔═╡ 5f935b56-1f1d-4926-84fd-864f12a3ede6
fig4 = let fampl = "holmes87-fig4-ampl.txt", fphase="holmes87-fig4-phase.txt"
	loadpoints((fampl,fphase), "points")
end

# ╔═╡ f568e0b7-678b-40a4-913d-514a7cbc476e
tube4 = let D=1.5e-3, L=0.43, V=107e-9, 
             σ=0.0,γ=1.4, Pa=101325, Ta=288.15, 
			 Pr=0.707, μ=1.82e-5
	PressureLine.(D*0.98, L, V, σ, γ; Ta=Ta, gamma=γ, Pa=Pa, Pr=Pr, mu=μ)
end

# ╔═╡ c9618ebb-e4b4-47f1-b36e-8826dc3c4870
freq = 0.2:0.2:800.0

# ╔═╡ 88ac849a-f5d7-4d5b-98ab-c469d6dcd58c
let fig=Figure(), press=tube4
	r = press.(freq)
	ax1 = Axis(fig[1,1], title="Holmes and Lewis (1987) figure 4", xlabel="Frequency (Hz)", ylabel="Amplitude")
	ax2 = Axis(fig[1,2], xlabel="Frequency (Hz)", ylabel="Phase")
	bA = fig4[1]
	bϕ = fig4[2]
	A = abs.(r)
	ϕ = phaseangle(r) * 180/π
	lines!(ax1, freq, A)
	scatter!(ax1, bA[:,1], bA[:,2])
	
	lines!(ax2, freq, ϕ)
	scatter!(ax2, bϕ[:,1], bϕ[:,2])
	
	fig
end


# ╔═╡ Cell order:
# ╟─fffad682-611f-11ed-3fdf-93926ecdc1a5
# ╠═5c98f2bc-9076-4e10-8b84-cbdaf34e65ca
# ╠═3db66fc0-6a59-473c-9059-c00251dcacaf
# ╠═6c3a2e0a-abe9-4d90-9cd1-b3e0b8373cef
# ╠═5f935b56-1f1d-4926-84fd-864f12a3ede6
# ╠═f568e0b7-678b-40a4-913d-514a7cbc476e
# ╠═c9618ebb-e4b4-47f1-b36e-8826dc3c4870
# ╠═88ac849a-f5d7-4d5b-98ab-c469d6dcd58c
