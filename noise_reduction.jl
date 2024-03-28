### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ d655634e-8ee5-48d5-93e3-07f6b0e85567
begin
	using Random
	using Statistics
	using Distributions
	using DSP
	using FFTW
	using PlutoPlotly
	using PlutoUI
end

# ╔═╡ 2ad61290-ec52-11ee-1ed4-bf2c121bd57f
md"""
# Noise Reduction: Hands-On Notebook

Hongyang Zhou, 2024/03/27
"""

# ╔═╡ 23d012e7-3f1f-4205-892e-0a55c284a51f
md"""
## Generating example signal

Imagine we have a signal with 101 sample points. It consists of four parts:

- a sinusoidal wave with period 5
- a sinusoidal wave with period 20
- a sinusoidal wave with period 40
- Gaussian noise with mean 0 and standard deviation 1
"""

# ╔═╡ 729c180b-8f5b-49c5-b991-df081d05ccc7
begin
	"`fs` is the sampling frequency."
	function generate_fake_signal()
		rng = MersenneTwister(1234)
		t = 0:100
		noise = 0.5 .* randn(rng, length(t))

		T1 = 5.0     # period [s]
		f1 = 1 / T1  # frequency [Hz]

		y1 = @. sin(2π*f1*t)

		T2 = 20.0    # period [s]
		f2 = 1 / T2  # frequency [Hz]
		y2 = @. sin(2π*f2*t)		

		T3 = 40.0    # period [s]
		f3 = 1 / T3  # frequency [Hz]
		y3 = @. sin(2π*f3*t)

		y_truth = @. y1 + y2 .+ y3
		y_noise = @. y_truth .+ noise

		t, y_truth, y_noise
	end

	function moving_average(g::AbstractVector{<:AbstractFloat}, n::Int)
   		nbackward = div(n, 2)
   		nforward = isodd(n) ? nbackward : nbackward - 1
   		len = length(g)
   		g_avg = similar(g)
   		@inbounds for i in 1:len
      		lo = max(1, i - nbackward)
      		hi = min(len, i + nforward)
      		g_avg[i] = mean(@view g[lo:hi])
   		end

   		g_avg
	end

	function neighbor_average(g, n=1; iterations=1)
		len = length(g)
   		g_avg = copy(g)
		g_tmp = similar(g)
		for iter in 1:iterations
   			@inbounds for i in 1:len
      			lo = max(1, i - n)
      			hi = min(len, i + n)
      			g_tmp[i] = mean(@view g_avg[lo:hi])
   			end
			g_avg = g_tmp
		end

   		g_avg
	end


	fs = 1 # sampling frequency [Hz]
	t, y_truth, y = generate_fake_signal()
	N = length(t)
end

# ╔═╡ 7de4a2d8-e88d-4ed8-856c-7be58322aa6b
plot(
    scatter(x=t, y=y, mode="lines"),
    Layout(yaxis_title="Outputs", xaxis_title="Samples", title="Signal")
)

# ╔═╡ f54f7a94-a68b-468d-b92e-fe47e69285f6
md"""
## Classic Filters

- High-Pass
- Low-Pass
- Band-Pass
- Moving box average
- Neighbor average
"""

# ╔═╡ 0ab781d9-e250-47f0-948f-e52d4207e385
let
	Tcut = 10.0 # [s]
	f0 = 1 / Tcut  # cutoff frequencies [Hz]
	responsetype = Highpass(f0; fs)
	designmethod = Butterworth(4)
	hpfilt = digitalfilter(responsetype, designmethod)

	y_highpass = filtfilt(hpfilt, y)

	plot([
		scatter(x=t, y=y, name="raw", line = attr(dash="dot")),
		scatter(x=t, y=y_truth, name="truth", line = attr(dash="dashdot")),
    	scatter(x=t, y=y_highpass, mode="lines", name="High-Pass")],
		Layout(yaxis_title="Outputs", xaxis_title="Samples", title="High-Pass Filtered Signal")
	)
end

# ╔═╡ fd7d0d0e-1943-4901-a303-a5774c4bcebc
let
	Tcut = 30
	f0 = 1 / Tcut
	responsetype = Lowpass(f0; fs)
	designmethod = Butterworth(4)
	lpfilt = digitalfilter(responsetype, designmethod)

	y_lowpass = filtfilt(lpfilt, y)

	plot([
		scatter(x=t, y=y, name="raw", line = attr(dash="dot")),
		scatter(x=t, y=y_truth, name="truth", line = attr(dash="dashdot")),
    	scatter(x=t, y=y_lowpass, mode="lines", name="Low-Pass")], 
		Layout(yaxis_title="Outputs", xaxis_title="Samples", title="Low-Pass Filtered Signal")
	)
end

# ╔═╡ 8749dc10-c27b-46a6-b214-687a99c101a7
let
	fl, fh = 1/30, 1/10  # cutoff frequencies [Hz]
	responsetype = Bandpass(fl, fh; fs)
	designmethod = Butterworth(4)
	filter = digitalfilter(responsetype, designmethod)

	y_bandpass = filtfilt(filter, y)

	plot([
		scatter(x=t, y=y, name="raw", line = attr(dash="dot")),
		scatter(x=t, y=y_truth, name="truth", line = attr(dash="dashdot")),
    	scatter(x=t, y=y_bandpass, mode="lines", name="Band-Pass")],
		Layout(yaxis_title="Outputs", xaxis_title="Samples", title="Band-Pass Filtered Signal")
	)
end

# ╔═╡ 97a1f4c0-570e-4035-a028-ec76f61d316f
let
	ȳ3 = moving_average(y, 3)
	ȳ5 = moving_average(y, 5)
	ȳ10 = moving_average(y, 10)

	plot([
		scatter(x=t, y=y, name="raw", line=attr(color="black", dash="dot")),
		scatter(x=t, y=y_truth, name="truth", line = attr(dash="dashdot")),
    	scatter(x=t, y=ȳ3, mode="lines", name="nbox = 3"),
		scatter(x=t, y=ȳ5, mode="lines", name="nbox = 5"),
		scatter(x=t, y=ȳ10, mode="lines", name="nbox = 10")],
		Layout(yaxis_title="Outputs", xaxis_title="Samples", title="Moving Box Averaged Signal")
	)
end

# ╔═╡ e588c06c-6dcc-48d1-bb68-221c19575ca8
let
	ȳ1 = neighbor_average(y, 1; iterations=1)
	ȳ10 = neighbor_average(y, 1; iterations=10)
	ȳ20 = neighbor_average(y, 1; iterations=20)

	plot([
		scatter(x=t, y=y, name="raw", line=attr(color="black", dash="dot")),
		scatter(x=t, y=y_truth, name="truth", line = attr(dash="dashdot")),
    	scatter(x=t, y=ȳ1, mode="lines", name="nIter = 1"),
		scatter(x=t, y=ȳ10, mode="lines", name="nIter = 10"),
		scatter(x=t, y=ȳ20, mode="lines", name="nIter = 20")],
		Layout(yaxis_title="Outputs", xaxis_title="Samples", title="Neighbor-Averaged Signal")
	)
end

# ╔═╡ 47adebc5-306e-4d4a-a2b2-0c9dbd164914
md"""
## Modern Filters

- Kalman Filter
"""

# ╔═╡ f7f6d21b-e830-49e3-8605-2da8ba082cb7
md"""
## Kernel Density Estimation

Imagine you observe 6 particles at a location with velocities [9,12,10,55,62,50]. The question is: **what is the true distribution of particles?** Kernel Density Estimation (KDE) is a method specifically tackling this task.

```math
\hat{f}_h(x) = \frac{1}{n} \sum_{i=1}^n K_h(x - x_i)
```
where n is the number of samples, ``K_h`` is the kernel centered at ``x_i`` with width h. The kernel can take any function that satifies three conditions:

```math
\begin{aligned}
\int K(x) \mathrm{d}x &= 1 \\
K(x) &= K(-x) \\
\int x^2 K(x)\mathrm{d}x &< \infty
\end{aligned}
```

The common choices are:

- Gaussian
- Uniform
- Triangle

As the number of samples n increases, the shape of the kernel function matters less; instead, the width h is the key parameter to find.

Finding the optimal width is equivalent to the well-known *bias-variance trade-off problem* in statistics.

- Larger h => smoother density estimation, larger **bias** error
- Smaller h => spikier density estimation, larger **variance** error
"""

# ╔═╡ 481b9a8c-c0ff-4e50-b0b3-46065558071a
let
	weights = [9, 12, 10, 55, 62, 50]
	n = length(weights)
	x = range(-10, 80, length=500)
	widths = [1, 5, 20]
	f̂ = [Float64[] for _ in eachindex(widths)]
	for iw in eachindex(widths)
		dists = [Normal(weights[i], widths[iw]) for i in eachindex(weights)]
		kernels = [pdf.(dist, x) for dist in dists]
		f̂[iw] = mean(kernels, dims=1)[1]
	end

	plot([
    	scatter(;x, y=f̂[i], mode="lines", name="h = $(widths[i])") for i in eachindex(widths)],
    	Layout(yaxis_title="f", xaxis_title="x", title="Distribution Estimates")
	)
end

# ╔═╡ 80da6802-0585-4fe8-9162-a6989bc48796
md"""
As you can see from above, the width h is a free parameter which exhibits a strong influence on the resulting estimate. How to pick a suitable width for our purpose?

[Wiki: Bandwidth Selection](https://en.wikipedia.org/wiki/Kernel_density_estimation#Bandwidth_selection)

### Error definition

In statistics, we try to minimize the *Mean Integrated Squared Error (MISE)*:

```math
\mathbb{E}\left[ \int_x (\hat{f}_h(x) - f(x))^2\mathrm{d}x \right]
```

The expectation is taken on random variables ``\{X_i | i=1,...,n\}`` of observed distribution ``\hat{f}(x)``. MISE is essentially the sum of squared bias and variance at each point x,

```math
\mathrm{MISE} = \mathrm{bias}^2 + \mathrm{variance}
```

- bias: deviation of the estimated density from the true value, i.e. systematic error
- variance: fluctuation of the density estimation, i.e. noise

### Optimization problem

```math
\frac{\partial}{\partial h}\mathrm{MISE}(h) = 0
```

However, we don't know the true distribution f(x). How can we proceed? A variety of automatic, data-based methods have been developed to select the bandwidth:

- Rules-of-thumb selector: If Gaussian basis functions are used to approximate univariate data, and the underlying density being estimated is Gaussian, the optimal choice for h is given by Silverman

```math
h = \left( \frac{4\hat{\sigma}^5}{3n} \right)^{1/5} \approx 1.06 \hat{\sigma}n^{-1/5}
```
where $\hat{\sigma}$ is the std of samples.

- Plug-in selector: an iterative process which requires ``f^{\prime\prime}``

- Visual inspection: brute force
"""

# ╔═╡ 9ddb45a6-9318-4e33-8f69-f11e5ba9a158
md"""
### Cross-validation (CV) selector

```math
\begin{aligned}
\mathrm{MISE} &= \mathbb{E}\left[ \int_x (\hat{f}_h(x) - f(x))^2\mathrm{d}x \right] \\
&= \mathbb{E}\int \hat{f}(x)^2\mathrm{d}x - 2\mathbb{E}\int\hat{f}(x)f(x)\mathrm{d}x + \int f(x)^2\mathrm{d}x
\end{aligned}
```

- 1st term: direct compute
- 3rd term: constant, independent of h
- 2nd term: involves unknown f, ``\int \hat{f}f\mathrm{d}x =\mathbb{E}[\hat{f}(x)]``
  - ``\sum_{i=1}^{n}\hat{f}(X_i)/n``  =>  ``h = 0``, not good!
  - Leave-one-out strategy
  - Jackknife resampling

```math
\hat{\mathbb{E}}_X\left[ \hat{f}(X) \right] = \frac{1}{n}\sum_{i=1}^n \hat{f}_{-i}(X_i)
```
where
```math
\hat{f}_{-i}(X_i) = \frac{1}{(n-1)h}\sum_{j\neq i, j=1}^n K\left( \frac{X_i - X_j}{h} \right)
```
is the *leave-one-out* estimator.

The CV function of kernel width h is calculated by

```math
\begin{aligned}
\mathrm{CV}(h) &= \int \hat{f}(x)^2\mathrm{d}x - 2\mathbb{E}\int\hat{f}(x)f(x)\mathrm{d}x \\
&= \frac{1}{n^2h^2}\sum_{i=1}^n \sum_{j=1}^n \int K\left( \frac{X_i - x}{h} \right)\, K\left( \frac{X_j - x}{h} \right)\mathrm{d}x \\
& -\frac{2}{n}\sum_{i=1}^n \left[ \frac{1}{(n-1)h}\sum_{j\neq i, j=1}^n K\left( \frac{X_i - X_j}{h} \right) \right] \\
&= \frac{1}{n^2h^2}\sum_{i=1}^n \sum_{j=1}^n \bar{K}\left( \frac{X_i - X_j}{h} \right) \\
& - \frac{2}{n(n-1)h}\sum_{i=1}^n\sum_{j\neq i, j=1}^n K\left( \frac{X_i - X_j}{h} \right)
\end{aligned}
```
where ``\bar{K}(x) = \int K(x-t)K(t)\mathrm{d}t`` is a convolution.

The CV function is then used to estimate ``\hat{h}`` that minimize it.
"""

# ╔═╡ c8b068ef-0c62-4ea7-bbea-1c15defcbecf
md"""

## PIC Noise Reduction Techniques

- High-order shape functions
- Kernel width optimization
- Adaptive width adjustment

### High-order shape functions

FLEKS use a uniform box shape function. High-order (>2) shape functions can mitigate the noise effect.

### Kernel width optimization

```math
\hat{f}_h(x) = \frac{1}{nh} \sum_{i=1}^n K(\frac{x - x_i}{h})
```

- Larger width h => lower and flatter shape function
- Smaller width h => higher and narrower shape function

### Adaptive width adjustment

If the density has sharp peaks or narrow valleys, it is more likely that the density estimator that minimizes the mean and integrated error will *cut the peaks and fill the valleys*.

If the bandwidth is not held fixed, but is varied depending upon the location of either the estimate (balloon estimator) or the samples (pointwise estimator), this produces a particularly powerful method termed adaptive or variable bandwidth kernel density estimation.

- Priori density distribution of the CV-optimal width
- Shrinks the widths where the density is high.
- Enlarges the widths where the density is low.

Algorithm:

1. Find a pilot estimate ``\tilde{f}(x)`` by minimizing the CV function.
2. Define local wdith factor ``\lambda_i`` by

```math
\lambda_i = \left( \frac{f(X_i)}{g} \right)^{-\alpha}
```
where ``g`` is the geometric mean of ``\tilde{f}``, i.e.

```math
\log(g) = \frac{1}{n}\sum \log\left( \tilde{f}(X_i) \right)
```
and ``\alpha`` is a sensitive factor, which is usually set to 0.5. Larger ``\alpha`` corresponds to more sensitive dependence of the corrected density on the prior density, rather than the samples.

### When to apply the adjustment

- Density estimation is not cheap.
- When the distribution is close to uniform, the optimal width -> ``\infty``, so we shouldn't even perform the estimation!

Criterion:
- Anderson-Darling test, to tell whether or not samples obey uniform distribution.
"""

# ╔═╡ a111afa2-0039-4551-a26d-0f4f8e539de4
md"""
$(Resource("https://raw.githubusercontent.com/henry2004y/pluto_playground/master/figures/density_estimate_KDEvsPIC.png"))
"""

# ╔═╡ c94516da-72bb-42d6-b7fe-9fdaab1fc3c9
md"""
$(Resource("https://raw.githubusercontent.com/henry2004y/pluto_playground/master/figures/noise_reduction_shapefunction_kernelwidth.png"))
"""

# ╔═╡ 3099f6b6-dc46-46b6-b7a5-63907dfa32a4
md"""
## Using Diffusion to Detect PIC Noise

Let us think of a simple scenario. We are trying to represent a Maxwellian distribution in a PIC cell by sampling 10 points. If we keep increasing the number of samples, we can approximate the Maxwellian distribution. The goal is to learn the mapping from few particles to many particles to better approximate the physical distributions.

State-of-the-art!
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
DSP = "~0.7.9"
Distributions = "~0.25.107"
FFTW = "~1.8.0"
PlutoPlotly = "~0.4.6"
PlutoUI = "~0.7.58"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.2"
manifest_format = "2.0"
project_hash = "04e31178c49ef47686bcb61a4d6ce4f165dcdba6"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

    [deps.AbstractFFTs.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "0f748c81756f2e5e6854298f11ad8b2dfae6911a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BaseDirs]]
git-tree-sha1 = "3e93fcd95fe8db4704e98dbda14453a0bfc6f6c3"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.2.3"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "c955881e3c981181362ae4088b35995446298b80"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.14.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.0+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "260fd2400ed2dab602a7c15cf10c1933c59930a2"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.5"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "f7f4319567fe769debfcf7f8c03d8da1dd4e2fb0"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.9"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "0f4b5d62a88d8f59003e43c25a8a90de9eb76317"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.18"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "7c302d7a5fec5214eb8a5a4c466dcf7a51fcf169"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.107"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "5b93957f6dcd33fc343044af3d48c215be2562f1"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.9.3"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5fdf2fe6724d8caabf43b557b84ce53f3b7e2f6b"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.0.2+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "18144f3e9cbe9b15b070288eef858f71b291ce37"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.27"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "72dc3cf284559eb8f53aa593fe62cb33f83ed0c0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.0.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "BaseDirs", "Colors", "Dates", "Downloads", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "Pkg", "PlotlyBase", "Reexport", "TOML"]
git-tree-sha1 = "1ae939782a5ce9a004484eab5416411c7190d3ce"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.4.6"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "71a22244e352aa8c5f0f2adde4150f62368a3f2e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.58"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "a9c7a523d5ed375be3983db190f6a5874ae9286d"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.6"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9b23c31e76e333e6fb4c1595ae6afa74966a729e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.4"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "1d77abd07f617c4868c33d4f5b9e1dbb2643c9cf"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.2"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "cef0472124fab0695b58ca35a77c6fb942fdab8a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.1"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─2ad61290-ec52-11ee-1ed4-bf2c121bd57f
# ╠═d655634e-8ee5-48d5-93e3-07f6b0e85567
# ╟─23d012e7-3f1f-4205-892e-0a55c284a51f
# ╠═729c180b-8f5b-49c5-b991-df081d05ccc7
# ╠═7de4a2d8-e88d-4ed8-856c-7be58322aa6b
# ╠═f54f7a94-a68b-468d-b92e-fe47e69285f6
# ╠═0ab781d9-e250-47f0-948f-e52d4207e385
# ╠═fd7d0d0e-1943-4901-a303-a5774c4bcebc
# ╠═8749dc10-c27b-46a6-b214-687a99c101a7
# ╠═97a1f4c0-570e-4035-a028-ec76f61d316f
# ╠═e588c06c-6dcc-48d1-bb68-221c19575ca8
# ╟─47adebc5-306e-4d4a-a2b2-0c9dbd164914
# ╟─f7f6d21b-e830-49e3-8605-2da8ba082cb7
# ╠═481b9a8c-c0ff-4e50-b0b3-46065558071a
# ╟─80da6802-0585-4fe8-9162-a6989bc48796
# ╟─9ddb45a6-9318-4e33-8f69-f11e5ba9a158
# ╟─c8b068ef-0c62-4ea7-bbea-1c15defcbecf
# ╟─a111afa2-0039-4551-a26d-0f4f8e539de4
# ╟─c94516da-72bb-42d6-b7fe-9fdaab1fc3c9
# ╟─3099f6b6-dc46-46b6-b7a5-63907dfa32a4
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
