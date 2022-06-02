### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ 2cf5ea1c-e12f-11ec-2775-5b259f01aa62
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

# ╔═╡ 06bf6df5-d9f2-4dc0-ac27-e43f96a52675
TableOfContents()

# ╔═╡ 2d6fb5e4-2fb9-4cc4-89b8-124c73b2e62d
md"""
## Load CSVs
"""

# ╔═╡ 56bf3dc7-017b-4086-9cb7-04c18dfccade
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

# ╔═╡ baa710f3-1a9d-4359-b1eb-1e09368f06d1
begin
	dfs = []
	for path in paths
		df = CSV.read(path, DataFrame; header=false)
		push!(dfs, df)
	end
end

# ╔═╡ 2d04c0f1-e909-44cf-85a4-2ee3a6fcabe8
let
	f = Figure()
	ax = Axis(f[1, 1])
	df = dfs[1]

	scatter!(range(0, 1, 10000), df[!, :Column1])
	ax.xlabel = "time (s)"
	ax.ylabel = "height (m)"
	ax.title = "Raw Data"
	f
end

# ╔═╡ 04b77cdd-b3ba-4d14-9504-9bd2d3b22254
md"""
## Clean Data
"""

# ╔═╡ 58de5faf-35af-493e-a382-cba769743f9f
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

# ╔═╡ b9c0cfcb-b9be-4efe-b33c-33f91b61817e
let
	f = Figure()
	ax = Axis(f[1, 1])

	scatter!(ts_clean, clean_arrays[1])
	
	ax.xticks = collect(range(0, 1, 9))
	ax.xlabel = "time (s)"
	ax.ylabel = "height (m)"
	ax.title = "Cleaned Data"
	f
end

# ╔═╡ d1b85c93-d42f-4f21-8c66-977e8100a4a6
md"""
## Crop Data
"""

# ╔═╡ b7eac438-1134-42cb-a7d4-ab238d7e6469
begin
	diff_arrays = []
	for arr in clean_arrays
		diff_arr = rollmean(diff(arr), 100)
		push!(diff_arrays, diff_arr)
	end
	ts_diff = collect(range(1, length(diff_arrays[1]))) ./ 10000
end

# ╔═╡ ba810115-9d62-4fb7-84bc-e27e12cf16a7
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

# ╔═╡ c8902e5b-aba9-4e0b-b928-775efa2f7247
let
	f = Figure()
	ax = Axis(f[1, 1])

	scatter!(ts_cropped[4], cropped_arrs[4])
	
	# ax.xticks = collect(range(0, 1, 9))
	ax.xlabel = "time (s)"
	ax.ylabel = "height (m)"
	ax.title = "Cleaned and Cropped Data"
	f
end

# ╔═╡ 0e662e9f-c402-4e65-8c28-48c75c005052
md"""
## 1. Rest Position Shift
"""

# ╔═╡ a118263f-679b-43dd-a365-96dced810a17
begin
	rest_means = zeros(size(dfs))
	for i in 1:length(dfs)
		rest_means[i] = mean(clean_arrays[i][end-500:end-200])
	end
end

# ╔═╡ 48ff6a41-2ae4-458f-8f85-7b35ed722d84
collect(range(0, 2, 21))

# ╔═╡ a9524661-d945-440a-999c-6ca1a91a085c
let
	f = Figure()
	ax = Axis(f[1, 1])

	scatter!(range(0, 2, 21), rest_means; markersize=10)
	ax.xlabel = "excitation current (A)"
	ax.ylabel = "equilibrium positions (m)"
	f
end

# ╔═╡ c28b9f81-4016-43bb-bc24-72249c5de216
md"""
## 2. Double Well Potential
"""

# ╔═╡ 506dfd2b-3ba3-48e3-b4ec-170db3d84ebf
let
	f = Figure()
	ax = Axis(f[1, 1])

	scatter!(cropped_arrs[1])
	
	ax.xticks = collect(range(0, 10000, 9))
	f
end

# ╔═╡ edb68374-c731-4caf-9fca-178c3090e153
begin
	v01 = zeros(length(cropped_arrs[1]))
	for i in 1:length(v01)
		if i == 1
			v01[i] = 0
		elseif i == length(v01)
			v01[i] = 0
		else
			v01[i] = cropped_arrs[1][i] - cropped_arrs[1][i-1]
		end
	end
end

# ╔═╡ 5d591599-8868-4ee3-80a4-3507e21f4e23
v01_squared = v01.^2

# ╔═╡ 67c325dc-dbd3-4e5d-8504-b49b70339a9e
length(v01_squared), length(v01)

# ╔═╡ 330690b9-bb61-45d5-a655-4a07139c21b3
scatter(rollmean(cropped_arrs[1], 100), rollmean(v01_squared, 100); markersize=4)

# ╔═╡ 85aca39d-7e03-453f-a8ee-adda4d5df54b
let
	f = Figure()
	ax = Axis(f[1, 1])

	scatter!(rollmean(cropped_arrs[1], 100), rollmean(v01_squared, 100); markersize=4)
	
	# ax.xticks = collect(range(0, 1, 9))
	ax.xlabel = "height (m)"
	ax.ylabel = "Velocity^2 (m/s^2)"
	f
end

# ╔═╡ Cell order:
# ╠═2cf5ea1c-e12f-11ec-2775-5b259f01aa62
# ╠═06bf6df5-d9f2-4dc0-ac27-e43f96a52675
# ╟─2d6fb5e4-2fb9-4cc4-89b8-124c73b2e62d
# ╠═56bf3dc7-017b-4086-9cb7-04c18dfccade
# ╠═baa710f3-1a9d-4359-b1eb-1e09368f06d1
# ╟─2d04c0f1-e909-44cf-85a4-2ee3a6fcabe8
# ╟─04b77cdd-b3ba-4d14-9504-9bd2d3b22254
# ╠═58de5faf-35af-493e-a382-cba769743f9f
# ╟─b9c0cfcb-b9be-4efe-b33c-33f91b61817e
# ╟─d1b85c93-d42f-4f21-8c66-977e8100a4a6
# ╠═b7eac438-1134-42cb-a7d4-ab238d7e6469
# ╠═ba810115-9d62-4fb7-84bc-e27e12cf16a7
# ╟─c8902e5b-aba9-4e0b-b928-775efa2f7247
# ╟─0e662e9f-c402-4e65-8c28-48c75c005052
# ╠═a118263f-679b-43dd-a365-96dced810a17
# ╠═48ff6a41-2ae4-458f-8f85-7b35ed722d84
# ╟─a9524661-d945-440a-999c-6ca1a91a085c
# ╟─c28b9f81-4016-43bb-bc24-72249c5de216
# ╠═506dfd2b-3ba3-48e3-b4ec-170db3d84ebf
# ╠═edb68374-c731-4caf-9fca-178c3090e153
# ╠═5d591599-8868-4ee3-80a4-3507e21f4e23
# ╠═67c325dc-dbd3-4e5d-8504-b49b70339a9e
# ╠═330690b9-bb61-45d5-a655-4a07139c21b3
# ╟─85aca39d-7e03-453f-a8ee-adda4d5df54b
