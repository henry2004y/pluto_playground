### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ e9c03142-0848-4317-907c-c0702d6da74a
using PlutoUI

# ╔═╡ a5f670e0-66b1-11ee-2e93-bb1e0b2857e6
md"
# Test Particle Model for Ganymede

Ganymede is believed to have both an ionosphere and an exosphere.

The ionospheric model is based on a Monte Carlo approach, where test particles, i.e., macro particles representing a certain physical number of ions, are created and followed in the presence of electric and magnetic fields. The ionosphere is created from ionization of the neutral exosphere.

## Carnielli's Model

* Spherical grid, 100x90x180 in $r, \theta, \phi$, $r \in [1.0, 6.7]$

* number of particles: $7.1\times 10^7$

Inputs:

1. E and B, Cartesian mesh $\pm 8$ RG
2. Neutral species from an exosphere model, i.e. number density, bulk velocity and temperature of different species $O_2, H_20, H_2, H, O$, and $OH$.

The test particles are created at random positions in every cell of the exospheric grid where they get assigned a weight and an initial velocity.
The weight, $W_j$, equals the number of physical particles represented by the test particle, and is assigned according to

```math
W_j = \mathrm{d}t \times V_{\mathrm{cell},\mathrm{exo}}\times\sum_n n_n \times \nu_{i,n}
```
where $\mathrm{d}t$ is the time step, $V_{\mathrm{cell},\mathrm{exo}}$ is the volume of the cell of the exospheric grid, $n_n$ is the number density of the neutral species n in the cell where the macro-particle is produced and $\nu_{i,n}$ is the ionization frequency of the ion species i generated from the ionization of the neutral n.

The initial velocity of the test particle, $\mathbf{v}_j(t=0)$, is assigned in relation to the average velocity $\mathbf{v}_n$ of the neutral species in the cell where the particle is produced, according to

```math
\mathbf{v}_j(t=0) = \frac{\sum_n n_n \times \nu_{i,n} \times \left( \mathbf{v}_n \pm R\sqrt{2k_B\, T_n/m_j} \right)}{\sum_n n_n \times \nu_{i,n}}
```
where the sum is over all neutral species n whose ionization can lead to the ion species i, $T_n$ is the temperature of the neutral species n and R is a uniformly distributed random number between 0 and 1 that generates the velocity dispersion for the ionospheric species.

It was found that the effect of initial velocity is negligible, since the ions are quickly accelerated by the EM fields with speeds significantly higher compared to that of the neutral species. The effect of gravity is also negligible.

Tracing criteria:

1. The distance travelled by the test particle must not exceed the size of the surrounding ionospheric grid cells.
2. The distance travelled by the test particle must not exceed the resolution of the field grid.
3. The time step is smaller than $\frac{1}{20}T$ of the instant gyroperiod.
4. The test particles are followed until either they impact the moon's surface or they cross the outer boundaries of the simulation grid.

At the end of the simulation, the number density $n_i$, bulk velocity $\mathbf{u}_i$, and temperature $T_i$ of the ion species are given as

```math
\begin{aligned}
n_i &= \frac{\sum_j W_j}{V_\mathrm{iono}\, N_\mathrm{stat}} \\
\mathbf{u}_i &= \frac{\sum_j \mathbf{v}_j\, W_j}{\sum_j W_j} \\
T_i &= \frac{m_i}{3k_B}\left( \frac{\sum_j \mathbf{v}_j^2 W_j}{\sum_j W_j} - \mathbf{u}_j^2 \right)
\end{aligned}
```
where $V_\mathrm{iono}$ is the volume of the ionospheric grid cell, $N_\mathrm{stat}$ is a statistical parameter representing the number of test particles injected initially per exospheric grid cell.

### Parameters

#### Sources for the ionization processes

* Photo-ionization

Ganymede's sunlit exosphere is constantly photo-ionized by solar EUV radiation. The photo-ionization frequency of the ion species i, $\nu_i^\gamma$, is calculated as follows:

```math
\nu_i^\gamma = \sum_n \int_0^{\lambda_\mathrm{th}}\mathrm{d}\lambda I^\infty(\lambda)\times\sigma_{i,n}^\mathrm{ion}(\lambda)
```
where $\lambda$ is the radiation wavelength, $\lambda_\mathrm{th}$ is the threshold wavelength for ionization, $\sigma_{i,n}^\mathrm{ion}(\lambda)$ is the ionization cross-section of the neutral $n$ producing the ion species, $I^\infty(\lambda)$ is the unattenuated solar radiation spectral flux at Ganymede's orbit and the sum runs over all neutral species which produce i.

The solar flux may be obtain from the [TIMED/SEE data base](http://lasp.colorado.edu/see/) at 1 AU and extrapolated to Jupiter's orbital distance at 5.46 AU.

* Electron impact

The Jovian plasma is able to partially penetrate inside Ganymede's magnetosphere, and the energetic electrons (>tens of eV) are able to ionize the neutral exosphere.
It is natural to expect an asymmetry in the energy distribution between the open and closed magnetic field lines: in the region of closed magnetic field lines mainly low energy electrons would be present from photo-ionization of the neutral atmosphere, and not energetic electrons from the Jovian plasma sheet which are not able to penetrate.
Ionization from Jovian electrons is assumed to take place only within the region of open magnetic field lines. The energy distribution of the Jovian electrons is set to be spatially constant within this simulation volume. The ionization frequency, νie, is calculated as follows:

```math
\nu_i^e = \sum_n \int_{E_\mathrm{th}}^{E_\mathrm{max}}\mathrm{d}E\, I(E)\sigma_{e,n}^\mathrm{ion}(E)
```
where $E$ is the electron energy, $E_\mathrm{th}$ is the energy threshold for electron-impact ionization, $E_\mathrm{max}$ is the highest electron energy available from data, $I(E)$ is the electron differential flux, and $\sigma_{e,n}^\mathrm{ion}(E)$ is the energy-dependent ionization cross-section for the electron impact on the neutral species n producing the ion i. The sum runs over all neutral species n whose electron-impact ionization can lead to i.

Secondary electron contribution to the ionization was neglected.

The ionization frequencies calculated for electron-impact are given in the following table. For all neutral-ion pair considered, the electron-impact ionization dominates over photo-ionization in the region of open magnetic field lines, where both ionizing sources are active (electron impact is not included in the region of closed magnetic field lines).

| Ionization by solar photons (electrons) | $\nu^{h\nu}\,[10^{-8}\,\mathrm{s}^{-1}]$ | $\nu^e\,[\mathrm{s}^{-1}]$ |
|-----------------------------------------|------------------------------------------------|----------------------------|
| $H + h\nu(e^-) \rightarrow H^+ + e^-(+e^-)$ | 0.24                                           | 2.41                       |
| $H_2 + h\nu(e^-) \rightarrow H_2^+ + e^-(+e^-)$ | 0.23                                           | 3.02                       |
| $H_2 + h\nu(e^-) \rightarrow H + H^+ + e^-(+e^-)$ | 0.01                                           | 0.24                       |
| $H_2O + h\nu(e^-) \rightarrow H_2O^+ + e^-(+e^-)$ | 1.13                                           | 4.25                       |
| $H_2O + h\nu(e^-) \rightarrow H + OH^+ + e^-(+e^-)$ | 0.23                                           | 1.29                       |
| $H_2O + h\nu(e^-) \rightarrow O^+ + H_2 + e^-(+e^-)$ | 0.02                                           | 0.19                       |
| $H_2O + h\nu(e^-) \rightarrow H^+ + OH + e^-(+e^-)$ | 0.11                                           | 1.03                       |
| $O + h\nu(e^-) \rightarrow O^+ + e^-(+e^-)$ | 0.87                                           | 4.90                       |
| $O_2 + h\nu(e^-) \rightarrow O_2^+ + e^-(+e^-)$ | 1.75                                           | 9.05                       |
| $O_2 + h\nu(e^-) \rightarrow O + O^+ + e^-(+e^-)$ | 0.44                                           | 0.90                       |
| $OH + h\nu(e^-) \rightarrow OH^+ + e^-(+e^-)$ | 1.75                                           | 5.76                       |
"

# ╔═╡ 8f0f14a3-7af7-4fa5-9e38-3ee37130b15f
md"""
#### Configuration for the exosphere

Three species are considered: $H_2$, $O_2$, and $H_2O$.

$(LocalResource("./figures/exosphere_ganymede.jpg"))

#### Configuration for the magnetosphere

The lifetime of ionospheric species is typically a few minutes, in which the EM field does not change significantly.
"""

# ╔═╡ 230453df-7f81-4827-82ef-6ea98db3197f
md"""
## Our Model

We use Hall MHD G8 spherical output as a starting point.

### Procedures to obtain surface flux

1. Obtain magnetic field, density, velocity and current from the unstructured grid.

Requires a conversion from BATSRUS unstructured *.out format to VTK *.vtu format.

2. Calculate the electric field including convective and Hall terms.
```math
\mathbf{E} = -\mathbf{U}\times\mathbf{B} + \frac{\mathbf{J}\times\mathbf{B}}{en}
```

3. Interpolate the E and B fields onto a Cartesian grid.

```python
from paraview.simple import *

# create a new 'XML Unstructured Grid Reader'
outvtu = XMLUnstructuredGridReader(registrationName='out.vtu', FileName=['out.vtu'])

# Properties modified on outvtu
outvtu.PointArrayStatus = ['B_x [nT]', 'B_y [nT]', 'B_z [nT]', 'J_x [`mA/m^2]', 'J_y [`mA/m^2]', 'J_z [`mA/m^2]', 'Rho [amu/cm^3]', 'U_x [km/s]', 'U_y [km/s]', 'U_z [km/s]']
outvtu.TimeArray = 'None'

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=outvtu)
calculator1.Function = ''

# Properties modified on calculator1
calculator1.ResultArrayName = 'B'
calculator1.Function = '"B_x [nT]"*iHat+"B_y [nT]"*jHat+"B_z [nT]"*kHat'

# create a new 'Calculator'
calculator2 = Calculator(registrationName='Calculator2', Input=calculator1)
calculator2.Function = ''

# Properties modified on calculator2
calculator2.ResultArrayName = 'U'
calculator2.Function = '"U_x [km/s]"*iHat+"U_y [km/s]"*jHat+"U_z [km/s]"*kHat'

# create a new 'Calculator'
calculator3 = Calculator(registrationName='Calculator3', Input=calculator2)
calculator3.Function = ''

# Properties modified on calculator3
calculator3.ResultArrayName = 'J'
calculator3.Function = '"J_x [`mA/m^2]"*iHat+"J_y [`mA/m^2]"*jHat+"J_z [`mA/m^2]"*kHat'

# create a new 'Calculator'
calculator4 = Calculator(registrationName='Calculator4', Input=calculator3)
calculator4.Function = ''

# Properties modified on calculator4
calculator4.ResultArrayName = 'E'
calculator4.Function = 'cross(B, U)*1e-6 + cross(J, B) / ("Rho [amu/cm^3]"*1.6)*1e1'

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(registrationName='ResampleToImage1', Input=calculator4)
resampleToImage1.UseInputBounds = 0
resampleToImage1.SamplingBounds = [-3.0, 3.0, -3.0, 3.0, -3.0, 3.0]
resampleToImage1.SamplingDimensions = [100, 100, 100]

# save data
SaveData('test.vti', proxy=resampleToImage1, ChooseArraysToWrite=1,
    PointDataArrays=['B', 'E'])
```

4. Read the field info into Julia.

```julia
using ReadVTK

vtk = VTKFile("test.vti")
point_data = get_point_data(vtk)
B = get_data(point_data["B"]) .* 1e-9
E = get_data(point_data["E"])
datasize = (3, 100, 100, 100)
B = reshape(B, datasize)
E = reshape(E, datasize)
```

5. Tracing with a given seed.

```julia
using TestParticle
using OrdinaryDiffEq
using Meshes

function trace(x, y, z, E, B; trajectories::Int=1)
   ## Initialize particles
   x0 = [0.0, 0.0, 0.0]
   u0 = [0.0, 0.0, 0.0]
   stateinit = [x0..., u0...]

   param = prepare(x, y, z, E, B, species=Proton)
   tspan = (0.0, 10.0)
   
   prob = ODEProblem(trace!, stateinit, tspan, param)
   sols = Vector{ODESolution}(undef, trajectories)
   particles = sample_particles(trajectories)

   for i in 1:trajectories
      u0 = [p[i] for p in particles]
      prob = remake(prob; u0)

      sol = solve(prob, Tsit5(); isoutofdomain, verbose=false)
      sols[i] = sol
   end
   
   return sols
end   

const RG = 2634e3 # [m]

origin = Tuple(ReadVTK.get_origin(vtk))
Δ = Tuple(get_spacing(vtk))

grid = Meshes.CartesianGrid(datasize[2:4], origin, Δ)

extent = extrema(grid)
x = range(extent[1].coords[1], extent[2].coords[1], length=datasize[2]) .* RG
y = range(extent[1].coords[2], extent[2].coords[2], length=datasize[3]) .* RG
z = range(extent[1].coords[3], extent[2].coords[3], length=datasize[4]) .* RG

## Solve for the trajectories

@time sols = trace(x, y, z, E, B; trajectories=10)
```

Steps 2-3 are performed in ParaView:



Questions:

* How to seed the initial test particles?
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.52"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "f5c06f335ceddc089c816627725c7f55bb23b077"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "91bd53c39b9cbfb5ef4b015e8b582d344532bd0a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e47cd150dbe0443c3a3651bc5b9cbd5576ab75b7"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.52"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "b7a5e99f24892b6824a954199a45e9ffcc1c70f0"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═e9c03142-0848-4317-907c-c0702d6da74a
# ╟─a5f670e0-66b1-11ee-2e93-bb1e0b2857e6
# ╟─8f0f14a3-7af7-4fa5-9e38-3ee37130b15f
# ╟─230453df-7f81-4827-82ef-6ea98db3197f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
