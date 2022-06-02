### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ 58a00d34-c172-415c-8a21-8d4e29f390ad
# ╠═╡ show_logs = false
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
		Pkg.add("CairoMakie")
		Pkg.add("CSV")
		Pkg.add("DataFrames")
		Pkg.add("Statistics")
		Pkg.add("FFTW")
		Pkg.add("RollingFunctions")
	end
	
	using PlutoUI
	using CairoMakie
	using CSV
	using DataFrames
	using Statistics
	using FFTW
	using RollingFunctions
end

# ╔═╡ 81588a9b-7066-429a-8061-ba7d4e39706e
TableOfContents()

# ╔═╡ 32ebb54a-f2d4-4cb8-9f5f-85dc274f24d9
md"""
## Load CSVs
"""

# ╔═╡ 6b4f72b8-ef6c-495a-9d6c-0d86506f4d38
begin
	root_path = "/Users/daleblack/Google Drive/UCI/8. Spring 2022/Phys 106W/final/data/"
	paths = []
	for i in range(0, 2, 21)
		if i == 0.0
			path = string(root_path, 0)
		else
			path = string(root_path, i)
		end
		push!(paths, path)
	end
end

# ╔═╡ cf68550c-e1bc-4129-a819-305dd354cabd
begin
	dfs = []
	for path in paths
		df = CSV.read(path, DataFrame; header=false)
		push!(dfs, df)
	end
end

# ╔═╡ 4eeaa080-0096-47ff-b36e-df20e8332f33
let
	f = Figure()
	ax = Axis(f[1, 1])
	df = dfs[1]

	scatter!(range(0, 1, 10000), df[!, :Column1])
	f
end

# ╔═╡ e5f0dbdc-2536-48f7-8b3c-b9f45df41f4c
md"""
## Clean Data
"""

# ╔═╡ 266695f6-f4a7-43a5-a3e9-ec21b5ed1edc
begin
	clean_arrays = []
	for df in dfs
		clean = rollmean(vec(Array(df)), 100)
		m = minimum(clean)
		clean_new = (clean .- m)
		clean_new = clean_new .- 0.02
		push!(clean_arrays, clean_new)
	end
	ts_clean = collect(range(1, length(clean_arrays[1]))) ./ 10000
end

# ╔═╡ b36464bb-adf5-4d0f-8ee4-a83f878eae42
let
	f = Figure()
	ax = Axis(f[1, 1])

	scatter!(ts_clean, clean_arrays[1])
	
	ax.xticks = collect(range(0, 1, 9))
	ax.xlabel = "time (s)"
	ax.ylabel = "height (m)"
	f
end

# ╔═╡ 485d539a-7a6f-45aa-96be-7d1bd2721eb0
md"""
## Crop Data
"""

# ╔═╡ cd4a7d70-4acc-46a6-ad6b-0b0c088a4d73
begin
	diff_arrays = []
	for arr in clean_arrays
		diff_arr = rollmean(diff(arr), 100)
		push!(diff_arrays, diff_arr)
	end
	ts_diff = collect(range(1, length(diff_arrays[1]))) ./ 10000
end

# ╔═╡ 2cd632af-2579-49b4-9dd9-82e2fa0d8eb8
begin
	cropped_arrs = []
	ts_cropped = []
	for i in 1:length(clean_arrays)
		start_ind = findfirst(x->x>0.0000125, diff_arrays[i])
	end_ind = 10000 - findfirst(x->x>0.0000125, reverse(diff_arrays[i]))
		cropped_arr = clean_arrays[i][start_ind:end_ind]
		push!(cropped_arrs, cropped_arr)
		t_cropped = collect(range(1, length(cropped_arr))) ./ 10000
		push!(ts_cropped, t_cropped)
	end
end

# ╔═╡ 760f8f84-2a3b-4a95-8c3f-f9497eefa4e5
let
	f = Figure()
	ax = Axis(f[1, 1])

	scatter!(ts_cropped[4], cropped_arrs[4])
	
	ax.xticks = collect(range(0, 1, 9))
	ax.xlabel = "time (s)"
	ax.ylabel = "height (m)"
	f
end

# ╔═╡ 36aa1a1b-039a-42db-9da7-654b460b5366
md"""
## 1. FFT
"""

# ╔═╡ 6892fe16-663b-43f4-8d78-9e31bade8f67
begin
	ts_cropped_fft = []
	for t in ts_cropped
		maxi = maximum(t) / 2
		t_cropped_fft = t .- maxi
		push!(ts_cropped_fft, t_cropped_fft)
	end
end

# ╔═╡ f35c1c98-042d-413b-8333-199f0500c62a
let
	f = Figure()
	ax = Axis(f[1, 1])
	func = imag

	scatter!(rollmean(ts_cropped_fft[1], 100), rollmean(func.(fft(cropped_arrs[1])), 100))
	scatter!(rollmean(ts_cropped_fft[2], 100), rollmean(func.(fft(cropped_arrs[2])), 100))
	scatter!(rollmean(ts_cropped_fft[3], 100), rollmean(func.(fft(cropped_arrs[3])), 100))
	scatter!(rollmean(ts_cropped_fft[4], 100), rollmean(func.(fft(cropped_arrs[4])), 100))
	scatter!(rollmean(ts_cropped_fft[5], 100), rollmean(func.(fft(cropped_arrs[5])), 100))
	scatter!(rollmean(ts_cropped_fft[6], 100), rollmean(func.(fft(cropped_arrs[6])), 100))
	scatter!(rollmean(ts_cropped_fft[7], 100), rollmean(func.(fft(cropped_arrs[7])), 100))
	scatter!(rollmean(ts_cropped_fft[8], 100), rollmean(func.(fft(cropped_arrs[8])), 100))
	scatter!(rollmean(ts_cropped_fft[9], 100), rollmean(func.(fft(cropped_arrs[9])), 100))
	scatter!(rollmean(ts_cropped_fft[10], 100), rollmean(func.(fft(cropped_arrs[10])), 100))
	scatter!(rollmean(ts_cropped_fft[11], 100), rollmean(func.(fft(cropped_arrs[11])), 100))
	scatter!(rollmean(ts_cropped_fft[12], 100), rollmean(func.(fft(cropped_arrs[12])), 100))
	scatter!(rollmean(ts_cropped_fft[13], 100), rollmean(func.(fft(cropped_arrs[13])), 100))
	scatter!(rollmean(ts_cropped_fft[14], 100), rollmean(func.(fft(cropped_arrs[14])), 100))
	scatter!(rollmean(ts_cropped_fft[15], 100), rollmean(func.(fft(cropped_arrs[15])), 100))
	scatter!(rollmean(ts_cropped_fft[16], 100), rollmean(func.(fft(cropped_arrs[16])), 100))
	scatter!(rollmean(ts_cropped_fft[17], 100), rollmean(func.(fft(cropped_arrs[17])), 100))
	scatter!(rollmean(ts_cropped_fft[18], 100), rollmean(func.(fft(cropped_arrs[18])), 100))
	scatter!(rollmean(ts_cropped_fft[19], 100), rollmean(func.(fft(cropped_arrs[19])), 100))
	scatter!(rollmean(ts_cropped_fft[20], 100), rollmean(func.(fft(cropped_arrs[20])), 100))
	
	ax.xlabel = "f_k"
	ax.ylabel = "FFT"
	f
end

# ╔═╡ 12e0b6d2-a562-4b34-8015-f531ff0a192e
begin
	freqs = zeros(length(ts_cropped_fft))
	for i in 1:length(ts_cropped_fft)
		freqs[i] = ts_cropped_fft[i][end]
	end
end

# ╔═╡ dd41150c-ec92-43b6-a7db-b18aabbb8c2c
let
	f = Figure()
	ax = Axis(f[1, 1])

	scatter!(range(0, 2, 12), rollmean(freqs, 10))
	ax.xlabel = "Coil current (A)"
	ax.ylabel = "Key FFT Frequency (Hz)"
	
	f
end

# ╔═╡ Cell order:
# ╠═58a00d34-c172-415c-8a21-8d4e29f390ad
# ╠═81588a9b-7066-429a-8061-ba7d4e39706e
# ╟─32ebb54a-f2d4-4cb8-9f5f-85dc274f24d9
# ╠═6b4f72b8-ef6c-495a-9d6c-0d86506f4d38
# ╠═cf68550c-e1bc-4129-a819-305dd354cabd
# ╟─4eeaa080-0096-47ff-b36e-df20e8332f33
# ╟─e5f0dbdc-2536-48f7-8b3c-b9f45df41f4c
# ╠═266695f6-f4a7-43a5-a3e9-ec21b5ed1edc
# ╟─b36464bb-adf5-4d0f-8ee4-a83f878eae42
# ╟─485d539a-7a6f-45aa-96be-7d1bd2721eb0
# ╠═cd4a7d70-4acc-46a6-ad6b-0b0c088a4d73
# ╠═2cd632af-2579-49b4-9dd9-82e2fa0d8eb8
# ╠═760f8f84-2a3b-4a95-8c3f-f9497eefa4e5
# ╟─36aa1a1b-039a-42db-9da7-654b460b5366
# ╠═6892fe16-663b-43f4-8d78-9e31bade8f67
# ╟─f35c1c98-042d-413b-8333-199f0500c62a
# ╠═12e0b6d2-a562-4b34-8015-f531ff0a192e
# ╟─dd41150c-ec92-43b6-a7db-b18aabbb8c2c
