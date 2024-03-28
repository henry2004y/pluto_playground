### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ e9c03142-0848-4317-907c-c0702d6da74a
using PlutoUI

# ╔═╡ 4e71cdfb-1584-4c49-9b9a-b7cb86e4da78
TableOfContents()

# ╔═╡ a5f670e0-66b1-11ee-2e93-bb1e0b2857e6
md"
# Test Particle Model for Ganymede

Ganymede is believed to have both an ionosphere and an exosphere.

The ionospheric model is based on a Monte Carlo approach, where test particles, i.e., macro particles representing a certain physical number of ions, are created and followed in the presence of electric and magnetic fields. The ionosphere is created from ionization of the neutral exosphere.

"

# ╔═╡ 07edf9ea-8b85-40ef-a13c-b6907da2f546
md"
## Carnielli 2019

This is an ionosphere model of Ganymede assuming a known exosphere.

* Spherical grid, $100\times 90\times 180$ in $r, \theta, \phi$, $r \in [1.0, 6.7]$

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

# ╔═╡ fb6f9c35-c3b9-455f-9a5e-a73d4d2a3976
md"""
## Carnielli 2020

This is an exosphere model built on top of ion sputtering. The three ion sources considered are: thermal Jovian ions, energetic Jovian ions, and ionospheric ions.


* Sputtering of Ganymede’s ionospheric ions provides a significant contribution to the moon’s exosphere ($\sim 10%$)

* Ionospheric $\mathrm{O}_2^+$ is the major contributor for surface sputtering

* Ionospheric sputtering occurs mainly in the leading hemisphere at low latitudes

* Ionospheric sputtering could explain the asymmetry in the $\mathrm{H}_2\mathrm{O}$ density between Ganymede’s leading and trailing hemispheres.

The sputtering rate is calculated by multiplying the impact rate by Y, the energy-dependent sputtering yield. The sputtering yield comes from [Fama+ 2008] for impact energies less than 100 keV and [Johnson+ 2009] for impact energies larger than 100 keV:

```math
Y =
    \begin{cases}
      ... & E < 100\,\mathrm{keV}\\
      ... & E \ge 100\,\mathrm{keV}
    \end{cases} 
```

θ is the angle of incidence of the impacting ion ($0^\circ$ normal to the surface) and T is the surface temperature. The incidence angle is assumed to be $45^\circ$.

Open boundaries are applied.

### Thermal Jovian ions

$\mathrm{O}^+,\, n = 1.75\,\mathrm{cm}^{-3},\,T=360\,\mathrm{eV},\, U=140\,\mathrm{km/s}$

### Energetic Jovian ions

$\mathrm{H}^+,\, \mathrm{O}^{n+},\, \mathrm{and}\quad \mathrm{S}^{n+}.$

Analytical fits from [Mauk+ 2004] in the energy range from 20 keV to 30 MeV, with charge state of $q=2$ for atomic oxygen and $q=3$ for atomic sulphur. Divided into 10 equally spaced energy bins (in logarithmic scale).

Test particles were initialized with a random velocity direction, a random energy within the energy bin and a weight determined by the mean intensity.
For energetic ions the corotation flow speed – approximately 140 km/s in Ganymede’s frame of reference – contributes negligibly (a few keV) to the ion’s kinetic energy, hence the injection flux is uniform across all boundaries of the simulation box.

### Ionospheric ions

3D distribution of the neutral exosphere obtained from the model of [Leblanc+ 2017], and ionized it via photo- and electron-impact ionization. The ion species simulated are $\mathrm{O}_2^+, \mathrm{O}^+, \mathrm{H}_2\mathrm{O}^+, \mathrm{H}_2^+, \mathrm{H}^+,$ and $\mathrm{OH}^+$.

"""

# ╔═╡ a43f12de-ac3b-44ea-8b5e-fe719f0784da
md"""
## Fatemi 2016

[Fatemi+ 2016](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016GL068363) performed a hybrid PIC simulation and then a test particle for energetic ions using the EM field output. One important thing to notice is that in [Jia+ 2008](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2007JA012748) the total thermal pressure (3.8 nPa for G8, 1.9 nPa for other Galileo flybys) includes both thermal and non-thermal components. In the hybrid model of [Fatemi+ 2016], energetic ion component is not included because it will significantly reduce the time steps. However, in terms of magnetic field comparisons, this is not important at all.

To compensate for this lack of energetic ions in the hybrid model, they performed an additional test particle simulation for $\mathrm{H}^+,\, \mathrm{O}^{2+},\, \mathrm{and}\quad \mathrm{S}^{3+}$ ranging from energies of 1 keV to $10^4$ keV, again in logarithmic bins. Whereas the thermal Jovian plasma contain the bulk of the flux at Ganymede's orbit, the energetic particle spectrum (i.e., $E\ge 10$ keV) is thought to be the main driver for particle surface weathering, i.e., radiolysis and neutral sputtering. The distribution is again from [Mauk+ 2004] from Galileo Energetic Particle Detector (EPD) data for the G8 flyby.

Furthermore, in order to calculate the high-resolution spatial distribution of precipitating ions onto Ganymede's surface, we have used the background electric and magnetic fields from the hybrid model to trace ions from the upstream boundary through the model until each ion either struck Ganymede or left the simulation. Each ion trajectory is advanced in time at urn:x-wiley:grl:media:grl54330:grl54330-math-0003, where Ωg is the instantaneous ion gyrofrequency. If an ion strikes the surface of Ganymede, we record the local flux in a two-dimensional grid covering Ganymede's longitude and latitude. We have individually simulated H+, O++, and S+++ ions ranging from energies of 1 keV to 104 keV in 100.25 sized logarithmic bins. While the thermal Jovian plasma consists of singly charged protons and oxygen [Kivelson et al., 2004], the energetic component of Jupiter's magnetospheric particles is dominated mainly by a combination of H+, O++, and S+++ [Collier and Hamilton, 1995; Keppler and Krupp, 1996; Mauk et al., 2004]. Whereas the thermal Jovian plasma contain the bulk of the flux at Ganymede's orbit [Kivelson et al., 2004], the energetic particle spectrum (i.e., urn:x-wiley:grl:media:grl54330:grl54330-math-0004 keV) is thought to be the main driver for particle surface weathering, i.e., radiolysis and neutral sputtering [Johnson, 1990; Shi et al., 1995; Baragiola et al., 2013]. Thus, we focus on the energetic particle population for the particle-tracing model. To absolutely normalize the flux of incident energetic particles, we use the distributions for each species as determined by Mauk et al. [2004] from Galileo Energetic Particle Detector (EPD) data for the G8 flyby specifically. For each specific pair of ion composition and energy, we have simulated >108 particles to ensure sufficient statistics.

## Poppe 2018

They applied a backwards Liouville tracing technique to trace ion trajectories from a given endpoint in time and space (whether on the surface of Ganymede or in near-Ganymede space) backwards through the EM fields. For each selected endpoint, they traced ion trajectories with discrete values of energy and angle backwards in time until the trajectory either strikes the surface of Ganymede or exits the simulation domain. Particles whose back-traced trajectory intersects Ganymede are assigned zero distribution function, while trajectories that exit the simulation domain are assigned a distribution function value based on their velocity at the simulation boundary assuming that the boundaries represent the undisturbed Jovian plasma and energetic particle distributions.

The backwards Liouville technique provides the advantage of high spatial and velocity space resolution of plasma distribution functions, since only particles that hit the surface are counted.

## Potini 2021

This is the Swedish team.

They took the EM field from [Fatemi+ 2016] as input, first did test particle ions and then sputtering estimation based on empirical formulas. My preliminary thermal $O^+$ result is very similar to theirs.
"""

# ╔═╡ 9425d497-7b5f-4317-8895-f0cdbc19901c
md"""
## Plainaki 2015 & 2020

This is an Italian team.

They mentioned a very good point about test particle tracing. The larger the energy is, the more important is the finite-Larmor-radius effect. For ions at $\sim \mathrm{MeV}$, the ion gyroradius is at Ganymede's radius scale, which is also what I observe by looking at single ion trajectories.

### Test particle seeds

They placed a $1\, \mathrm{R}_\mathrm{G}$ thick planar source surface perpendicular to the moon’s orbit, located between $X = −3\,\mathrm{R}_\mathrm{G}$ and $X = −4\, $\mathrm{R}_\mathrm{G}$ upstream of Ganymede. The source surface is subdivided in $0.2\times 0.2\,\mathrm{R}_\mathrm{G}$ cells and 1000 test particles are launched with a defined initial energy but random initial direction from each cell, simulating a total of about $10^7$ ions in each run.

hyzhou comment: the division is completely unnecessary. Random direction is also a tricky part, but I assume they did it right, i.e. sampling from a unit sphere or equivalently even sampling for $(\theta, \phi)$.

### Ion sputtering and radiolysis

Once the magnetospheric ions precipitate on the surface of Ganymede, they can cause sputtering, ionization and excitation of water–ice molecules. Following electronic excitations and ionizations, water–ice molecules can get dissociated; chemical reactions among the water-dissociation products result in the formation of new molecules (e.g. O2, H2, OH and minor species) that are finally ejected from the surface into the moon’s exosphere. Laboratory measurements of ice irradiation experiments have shown that H2O molecules dominate the total release yield at lower temperatures (<120 K) and O2 and H2 at higher (>120 K) temperatures (Johnson, 2001). Nevertheless, any H2 formed in ice diffuses and escapes much more efficiently than O2 at the relevant temperatures in the outer Solar System; moreover, H2 escapes from the icy moons because of its low mass and the relatively weak gravitational fields. Therefore, the irradiation of Ganymede’s surface can preferentially populate the magnetosphere with hydrogen, as is the case at Europa, leaving behind an oxygen-rich satellite surface.

### Exosphere

They developed a collisionless MC model that uses the simulated magnetospheric ion fluxes impacting the surface of the moon as inputs. The surface is considered to be composed of pure ice.

The modeling technique simulating $\mathrm{O}_2$ radiolysis and sputtering is based on the Europa Global model of Exospheric Outgoing Neutrals (EGEON).
The H2O sputtering simulation technique is based on the exospheric MC model of [Mura+ 2009] firstly applied to Mercury.
Those models provided the 3D spatial density distribution of the main species released from the bodies’ surfaces.

"""

# ╔═╡ 230453df-7f81-4827-82ef-6ea98db3197f
md"""
## Our Model

We use Hall MHD G8 spherical output as a starting point.

### Seeding

An important question is how to set seeds.

### Procedures to obtain surface flux

1. Obtain magnetic field, density, velocity and current from the unstructured grid.

Requires a conversion from BATSRUS unstructured *.out format to VTK *.vtu format.

2. Interpolate MHD quantities onto a Cartesian grid and calculate the derived electric field including the convective term and Hall term.
```math
\mathbf{E} = -\mathbf{U}\times\mathbf{B} + \frac{\mathbf{J}\times\mathbf{B}}{en}
```
If it is an ideal MHD output, we should only consider the convective term. In Xianzhe's ideal MHD model, there is an additional artificial resistivity. I am not sure if it is also counted.

```python
from paraview.simple import *

# create a new 'XML Unstructured Grid Reader'
outvtu = XMLUnstructuredGridReader(registrationName='out.vtu', FileName=['out.vtu'])

# Properties modified on outvtu
outvtu.PointArrayStatus = ['B_x [nT]', 'B_y [nT]', 'B_z [nT]', 'J_x [`mA/m^2]', 'J_y [`mA/m^2]', 'J_z [`mA/m^2]', 'Rho [amu/cm^3]', 'U_x [km/s]', 'U_y [km/s]', 'U_z [km/s]']
outvtu.TimeArray = 'None'

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=outvtu)
calculator1.ResultArrayName = 'B'
calculator1.Function = '"B_x [nT]"*iHat+"B_y [nT]"*jHat+"B_z [nT]"*kHat'

# create a new 'Calculator'
calculator2 = Calculator(registrationName='Calculator2', Input=calculator1)
calculator2.ResultArrayName = 'U'
calculator2.Function = '"U_x [km/s]"*iHat+"U_y [km/s]"*jHat+"U_z [km/s]"*kHat'

# create a new 'Calculator'
calculator3 = Calculator(registrationName='Calculator3', Input=calculator2)
calculator3.ResultArrayName = 'J'
calculator3.Function = '"J_x [`mA/m^2]"*iHat+"J_y [`mA/m^2]"*jHat+"J_z [`mA/m^2]"*kHat'

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(registrationName='ResampleToImage1', Input=calculator3)
resampleToImage1.UseInputBounds = 0
resampleToImage1.SamplingBounds = [-4.0, 4.0, -4.0, 4.0, -4.0, 4.0]
resampleToImage1.SamplingDimensions = [160, 160, 160]

# create a new 'Calculator', assume average mass = 14 amu
calculator4 = Calculator(registrationName='Calculator4', Input=resampleToImage1)
calculator4.ResultArrayName = 'Econv'
calculator4.Function = 'cross(B, U)*1e-6'

# create a new 'Calculator', assume average mass = 14 amu
calculator5 = Calculator(registrationName='Calculator5', Input=calculator4)
calculator5.ResultArrayName = 'E'
calculator5.Function = 'cross(B, U)*1e-6 + cross(J, B) * 14 / ("Rho [amu/cm^3]"*1.6)*1e-2'

# save data
SaveData('test.vti', proxy=calculator5, ChooseArraysToWrite=1,
    PointDataArrays=['B', 'E', 'U', 'n', 'Uth'])
```

Note the ordering of interpolations. It would be better to first convert the basic quantities from spherical coordinates to Cartesian coordinates, then calculate the derived electric field. If we first calculate the electric field and then interpolate to Cartesian coordinates, then around the inner boundary region there might be some artifacts.

4. Read the field info into Julia.

```julia
using ReadVTK

const RG = 2634e3 # [m]

function getEM(file::String)
   vtk = VTKFile(file)
   point_data = get_point_data(vtk)
   B = get_data(point_data["B"]) .* 1e-9
   #E = get_data(point_data["E"])
   E = get_data(point_data["Econv"])
   datasize = (3, 160, 160, 160)
   B = reshape(B, datasize)
   E = reshape(E, datasize)

   origin = Tuple(ReadVTK.get_origin(vtk))
   Δ = Tuple(get_spacing(vtk))
   
   grid = Meshes.CartesianGrid(datasize[2:4], origin, Δ)
   
   extent = extrema(grid)
   x = range(extent[1].coords[1], extent[2].coords[1], length=datasize[2]) .* RG
   y = range(extent[1].coords[2], extent[2].coords[2], length=datasize[3]) .* RG
   z = range(extent[1].coords[3], extent[2].coords[3], length=datasize[4]) .* RG

   x, y, z, E, B
end
```

5. Tracing with a given seed.

```julia
using TestParticle
using OrdinaryDiffEq
using StaticArrays
using Meshes
using Vlasiator: qᵢ, mᵢ, kB
using Random
using ProgressMeter
using ChunkSplitters

"Initial state perturbation for EnsembleProblem."
function prob_func(prob, i, repeat)
   #r = rand(3)
   #x = ((r[1] .* dx) .- 3.9) .* RG
   #y = dy .* (r[2] .- 0.5) .* RG
   #z = dz .* (r[3] .- 0.5) .* RG
   # Sample from Maxwellian
   #v = TestParticle.sample(vdf, 1)
   #prob = remake(prob; u0=SVector{6}([x, y, z, v...]))

   # How about energetic particles?
   r = randn(6)
   rmag1 = hypot(r[1], r[2], r[3])
   energy = 1e5
   rseed = 3.5RG
   x = r[1]/rmag1*rseed
   y = r[2]/rmag1*rseed
   z = r[3]/rmag1*rseed
   vmag = @. √(2*energy*qᵢ/(mass*mᵢ))
   # Sample on a sphere
   rmag2 = hypot(r[4], r[5], r[6])
   vx = r[4]/rmag2*vmag
   vy = r[5]/rmag2*vmag
   vz = r[6]/rmag2*vmag

   prob = remake(prob; u0=SVector{6}([x, y, z, vx, vy, vz]))

   prob
end

function output_func(sol, i)
   return last(sol), false
   #return sol, false
end

function trace_trajectory_fast(x, y, z, E, B; trajectories::Int=1, mass::Float64=1.0,
   tmax::Float64=500.0)
   ## Initialize particles
   x0 = [0.0, 0.0, 0.0] # initial position, [m]
   u0 = [0.0, 0.0, 0.0] # initial velocity, [m/s]
   stateinit = SVector{6}([x0..., u0...])

   param = prepare(x, y, z, E, B, species=User, q=TestParticle.qᵢ, m=mass*TestParticle.mᵢ)
   tspan = (0.0, tmax)
   
   prob = ODEProblem(trace, stateinit, tspan, param)
   ensemble_prob = EnsembleProblem(prob; prob_func, output_func, safetycopy=false)

   # p = (q, m, E, B)
   dtFE(u, p, t) = p[2] / (2π * abs(p[1]) * hypot(p[4](u, t)...))
   cb = StepsizeLimiter(dtFE; safety_factor=1 // 15, max_step=true)

   sols = solve(ensemble_prob, Vern6(), EnsembleThreads(); callback=cb, dt=0.1,
      trajectories, isoutofdomain, verbose=false);
   #sols = solve(ensemble_prob, Tsit5(), EnsembleThreads();
   #   trajectories, isoutofdomain, verbose=false);

   sols
end


function find_hits(sols; nchunks=Threads.nthreads())
   xu_t = [Vector{SVector{6, Float64}}(undef, 0) for _ in 1:nchunks]
   nhits, nouts = zeros(Int, nchunks), zeros(Int, nchunks)

   Threads.@threads for (irange, ichunk) in chunks(sols, nchunks)
      for i in irange
         rfinal = hypot(sols[i][1:3]...) / RG
         if rfinal < 1.08
            nhits[ichunk] += 1
            push!(xu_t[ichunk], sols[i])
         elseif rfinal > 3
            nouts[ichunk] += 1
         end
      end
   end

   xu = vcat(xu_t...)
   nhit = sum(nhits)
   nout = sum(nouts)
   finish_rate = (nhit + nout) / length(sols)
   if finish_rate ≥ 0.99
      @info "Finished rate: $finish_rate"
   else
      @error "Finished rate: $finish_rate too low! Consider increase tmax!"
   end
   @info "Number of particles hit the exosphere: $(length(xu))"

   xu
end

#######################################

const U = 140e3 # [m/s]
# assume a constant thermal speed upstream (CHECK!)
# Uth = √(2 p/ρ), P = n kB T 
# G8: √(2.0 * 3.6e-9 / (60e6 * mᵢ))
# Carnielli has very cold O+, only at 360 eV
# If we assume H+, 28e6 * 1.38e-23 * 360 * 11606 = 1.6 nPa makes sense
# If we assume O+, 1.75e6 * 1.38e-23 * 360 * 11606 = 0.1 nPa does not make sense
# Let's skip this point for now, and assume a super cold upstream species!!!
# G2: √(2.0 * 1.8e-9 / (60e6 * mᵢ))
#const Uth = 2.7e5
#const Uth = 0.0 # test
const Uth = 44629.0 # Carnielli's value
# particle mass [amu]
const mass = 1.0
# number density [m^-3]
const ni = 28e6
# number of seeds
trajectories = 1_000_00
# seeding box widths [RG]
const dx = 0.2
const dy = 4.0
const dz = 4.0
# tracing method
method = :accurate
# tracing max time [s]
tmax = 1000.0
# EM field file
file = "data/test_ideal_G2.vti"

## Data processing

# Create a Maxwellian distribution
const vdf = Maxwellian([U, 0.0, 0.0], Uth)

x, y, z, E, B = getEM(file)

## Solve for the trajectories

Random.seed!(1234)

#@time xu, tinterval = trace_trajectory(x, y, z, E, B;
#   trajectories, mass, tmax, method, dx, dy, dz)

@time sols = trace_trajectory_fast(x, y, z, E, B; trajectories, mass, tmax)

xu = find_hits(sols)

#tinterval = dx * RG / U
tinterval = 1.0

# particle weight
#wi = (ni/mass/(trajectories/(dx*dy*dz*RG^3)))
diffflux = 1e4 # [cm² s sr keV]^-1
energy = 1e2 # [keV]
wi = 4π*diffflux*1e4*energy*4π*(3.5RG)^2*tinterval / trajectories

# If all are H+, then one test particle represents (56e6 / (n/V)) physical particles;
# If all are O+, then one test particle represents (3.5e6 / (n/V)) physical particles.
@info "One test particle = $(round(wi, digits=4)) physical particles"
@info "Number of threads: $(Threads.nthreads())"
@info "Tracing method: $method"
@info "Particle mass: $mass [amu]"
@info "Number of particles: $trajectories"
@info "Sample time interval: $(round(tinterval, digits=2)) [s]"


using JLD2
@time jldsave("particle_$(Int(mass)).jld2"; xu, mass, RG, mᵢ, trajectories, tinterval, wi)
```

Step 2 is performed in ParaView with Python.
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

julia_version = "1.10.0-rc1"
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
version = "1.0.5+1"

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
version = "2.28.2+1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+2"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

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
deps = ["SHA"]
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
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
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
# ╠═e9c03142-0848-4317-907c-c0702d6da74a
# ╟─4e71cdfb-1584-4c49-9b9a-b7cb86e4da78
# ╟─a5f670e0-66b1-11ee-2e93-bb1e0b2857e6
# ╟─07edf9ea-8b85-40ef-a13c-b6907da2f546
# ╟─8f0f14a3-7af7-4fa5-9e38-3ee37130b15f
# ╟─fb6f9c35-c3b9-455f-9a5e-a73d4d2a3976
# ╟─a43f12de-ac3b-44ea-8b5e-fe719f0784da
# ╟─9425d497-7b5f-4317-8895-f0cdbc19901c
# ╠═230453df-7f81-4827-82ef-6ea98db3197f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
