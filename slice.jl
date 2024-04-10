### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ 524af530-61f5-11ee-1216-c527b94397de
md"
# SLICE2D

## Nomenclature

* **ECV**: Eulerian control volume. $\mathrm{ECV}_{i,j}$ denotes the Eulerian Control Volume surrounding the Eulerian grid point $(x_i,y_j)$ whose area is $\mathrm{EA}_{i,j} = (x_{i+1/2}-x_{i-1/2})\times(y_{j+1/2} - y_{j-1/2})$.
* **IECV**: Intermediate Eulerian control volume, associated with surrounding intersection points and aligned along the first cascade direction. $\mathrm{IECV}_{i,j}$ denotes the CV associated with the intersection point $(x_{i,j}^\prime, y_j)$ whose area is $\mathrm{IEA}_{i,j} = (x_{i+1/2}^\prime - x_{i-1/2}^\prime)\times(y_{j+1/2} - y_{j-1/2})$.
* **ILCV**: Intermediate Lagrangian control volume, associated with surrounding intersection points and aligned along the second cascade direction.
* Coordinates with *lower cases* (e.g. $x_{i,j}$) represent locations in the Eulerian mesh.
* Coordinates with a *prime* (e.g. $x_{i+1/2,j}^\prime$) represent locations for the IECV.
* Coordinates with a *tilde* (e.g. $\tilde{x}_{i,j-1/2}$) represent locations for the ILCV.
"

# ╔═╡ 66a18629-719e-4b0e-bc8e-21e5c1604404
md"
## Description

The arrival points at time $t^{n+1}$ collectively define an Eulerian $(x,y)$ mesh which subdivides the entire computational domain into a finite number of contiguous and non-overlapping regular 2D ECVs; similarly, the departure points at time $t^{n}$ collectively define a Lagrangian $(X,Y)$ mesh which subdivides the same domain into an equal number of contiguous and non-overlapping irregular 2D LCVs.

Consider a 2D domain $\Omega=[x_\mathrm{min}, x_\mathrm{max}]\times [y_\mathrm{min}, y_\mathrm{max}]$ subdivided into $N_x\times N_y$ ECVs with centers $(x_i, y_j), i=1,...,N_x, j=1,...,N_y$, and whose boundaries are the straight lines that join the four corner points $(x_{i\pm 1/2}, y_{j\pm 1/2}),i=1,...,N_x, j=1,...,N_y$. The curvilinear Lagrangian coordinate system, which corresponds to the transformed Eulerian $(x,y)$ system by the action of the velocity field over the interval $[t^n, t^{n+1}]$, is denoted by $(X, Y)$. A $y_j$-line of the Eulerian grid is the horizontal line defined by the set of points ${(x_{i-1/2}, y_j), i=1,...,N_x+1}$ and the Lagrangian $Y_i$ is the curve defied by the set of points ${(\tilde{x}_{i, j-1/2}, \tilde{y}_{i,j-1/2}), j=1,...,N_y+1}$. These latter points are the center points of the bottom edges of $\mathrm{LCV}_{i,j}$ along the Lagrangian $X_{j-1/2}$, i.e. $(\tilde{x}_{i, j-1/2}, \tilde{y}_{i,j-1/2}) = [(x_{i-1/2}^d, y_{j-1/2}^d) + (x_{i+1/2}^d, y_{j-1/2}^d)]/2$, where $(x^d, y^d)$ denotes the departure point of $(x, y)$.

The conservative semi-lagrangian scheme in 2D is a combination of the cascade approach and the SLICE-1D algorithm. It comprises a multiple sweep of SLICE-1D along one of the two Eulerian directions followed by a similar multiple sweep along one of the two Lagrangian ones. The order of sweeps does not matter: it can be a first cascade along Eulerian $x$ (holding $y$ constant), followed by one along the Lagrangian $Y$ (holding $X$ constant), or a first cascade along Eulerian $y$, then a second one along the Lagrangian $X$.

Consistent with the adopted CV approach, the problem considered here is the evaluation of the time evolution of the density averaged over an Eulerian grid box, i.e. given $\bar{\rho}(x,y,t^n)$, compute $\bar{\rho}(x,y,t^{n+1})$. The scheme advectively transports the mass to be conserved, then computes the grid-box-averaged density at the new time step $t^n$ from the mass that arrives at the ECV at time $t^{n+1}$. This is equivalent to computing the elementary masses of LCVs, then translating them to their corresponding ECVs.
"

# ╔═╡ 6abad483-131d-4458-b7fd-d4fa3f00f226
md"
### First stage

#### Computation of intersection points

The Lagrangian grid $(X_j, Y_i)$ is the linkage of Lagrangian points by straight lines in Eulerian space. Thus $Y_i$ is linear (in $x-y$ coordinates) in the region of an intersection. The intersection points $(x_{i-1/2, j}^\prime, y_j)$ between $Y_{i-1/2}$ and $y_j$ are computed using linear interpolation, and similarly for the intersection points $(x_{i,j}^\prime, y_j)$.

#### Transfer of mass to intersection control volumes (IECVs)

Since $\mathrm{ECV}_{i,j}$ and $\mathrm{IECV}_{i,j}$ share the same boundaries in the y-direction, $\mathrm{IECV}_{i,j}$ has an inherently Eulerian flavour and therefore the name of Intersection Eulerian CV seems appropriate. It is worth noting that for each column $Y_i$ two intersection points are computed:

* $(x_{i-1/2,j}^\prime, y_j)$ which is the intersection between $Y_{i-1/2}$ and $y_j$ and this gives directly the left boundary of $\mathrm{IECV}_{i,j}$
* $(x_{i,j}^\prime, y_j)$ which is the intersection between $Y_i$ and $y_j$ and this gives directly the center of mass for $\mathrm{IECV}_{i,j}$.

The first stage of the scheme is to redistribute conservatively the given masses
from $\mathrm{ECV}_{i,j}$ to $\mathrm{IECV}_{i,j}$. This is done using SLICE1D.

Given the masses
```math
\mathrm{ME}_{i,j}^n = \bar{\rho}_{i,j}^n \Delta x_i \Delta y_j
```
of each $\mathrm{ECV}_{i,j}$, 
"

# ╔═╡ b8ac1dec-aa06-443d-ac88-c43cb3769195
md"
### Second Stage

The second stage consists of redistributing the masses $\mathrm{MI}_{i,j}^n$ amongst the $\mathrm{LCV}_{i,j}$ in a conservative way. Although in physical space $(x,y)$ the $\mathrm{ILCV}_{i,j}$ and $\mathrm{LCV}_{i,j}$ are complex distorted areas, they are nevertheless regular in the Lagrangian space $(X, Y)$.

First the masses associated with the ILCVs need to be evaluated. This is done by
making the approximation that each ILCV actually has the same mass as its corresponding IECV. This is a reasonable approximation provided that the velocity field is smooth enough to give a sufficiently smooth distorted Lagrangian mesh. It is consistent with the assumption that the cascade approach is only accurate provided the velocity field has a sufficient degree of smoothness. The second stage of the remapping consists of 1D integrations over $\mathrm{LCV}_{i,j}$ of local pseudo-density functions $\tilde{\rho}(s)$ (mass per unit length along the curved distance of the Lagrangian $Y_i$). Similarly as the first stage, the $\tilde{\rho}(s)$ are defined such that their integrals over $\mathrm{ILCV}_{i,j}$ has the associated mass $\mathrm{MI}_{i,j}^n$ computed in the first stage.  This requires specification
of how $\mathrm{ILCV}_{i,j}$ and $\mathrm{LCV}_{i,j}$ overlap.

...

Having defined the topology, a pseudo-density $\tilde{\rho}(s)$ (or mass per unit length along the second cascade direction $s$, which defines the distance along the Lagrangian $Y_i$) can be fitted to the $\mathrm{ILCV}_{i,j}$, where $\mathrm{MI}_{i,j}^n$ is known.

Finally $\bar{\rho}_{i,j}^{n+1}$ is evaluated as
```math
\bar{\rho}_{i,j}^{n+1} \equiv \frac{M_{i,j}^{n+1}}{\Delta x_i \Delta y_j} = \frac{(M_{i,j}^d)^n}{\Delta x_i \Delta y_j} = \frac{ML_{i,j}^n}{\Delta x_i \Delta y_j}
```
"

# ╔═╡ 8eff7d2a-b948-4201-bde8-0e6ee388d6f7
md"
### Algorithm

1. Define the corner points of $\mathrm{ECV}_{i,j}$ to be $(x_{i\pm 1/2}, y_{j\pm 1/2}), i=1,...,N_x,\, j=1,...,N_y$.

2. Given $\bar{\rho}_{i,j}^n\equiv \bar{\rho}(x_i,y_j, t^n)$, compute the mass $ME_{i,j}^n = \bar{\rho}_{i,j}^n EA_{i,j}$ of $\mathrm{ECV}_{i,j}$.

3. Represent the density for each $\mathrm{ECV}_{i,j}$.

4. Locate the corners of $\mathrm{LCV}_{i,j}$, which are the departure points for the corners of $\mathrm{ECV}_{i,j}$, i.e. $(x_{i\pm 1/2}^d, y_{j\pm 1/2}^d)$.

5. Locate intersection points $(x_{i,j}^\prime, y_j)$.

6. Locate the $\mathrm{IECV}_{i,j}$ boundaries $[x_{i-1/2,j}^\prime, x_{i+1/2,j}^\prime]$.

7. Compute the mass $\mathrm{MI}_{i,j}^n$ of each $\mathrm{IECV}_{i,j}$ using SLICE1D.

8. Set the mass of each $\mathrm{ILCV}_{i,j}$ equal to $\mathrm{MI}_{i,j}^n$.

9. Compute the distances $s_{i,j}$ that define the boundaries of $\mathrm{LCV}_{i,j}$ along $Y_i$.

10. Compute the distances $s_{i,j}^\prime$ that define the boundaries of $\mathrm{ILCV}_{i,j}$ along $Y_i$.

11. Construct the representation for each $\mathrm{ILCV}_{i,j}$.

12. Compute the mass $\mathrm{ML}_{i,j}^n$ of each $\mathrm{LCV}_{i,j}$ using SLICE1D.

13. Transport $\mathrm{ML}_{i,j}^n$ to its corresponding $\mathrm{ECV}_{i,j}$ and compute $\bar{\rho}_{i,j}^{n+1}=ML_{i,j}^n / EA_{i,j}$.
"

# ╔═╡ 74d33f87-6fb7-4567-9b9e-ccff956425fe
md"
## Computational Efficiency

For a problem in $d$ dimensions, the cost of SLICE with a kth-order piecewise polynomial representation scales as $\mathcal{O}(kd)$, whereas it is $\mathcal{O}(k^d)$ for a standard SL scheme. SLICE has an additional overhead associated with the computation of intersections and piecewise polynomial coefficients.
"

# ╔═╡ Cell order:
# ╟─524af530-61f5-11ee-1216-c527b94397de
# ╟─66a18629-719e-4b0e-bc8e-21e5c1604404
# ╟─6abad483-131d-4458-b7fd-d4fa3f00f226
# ╟─b8ac1dec-aa06-443d-ac88-c43cb3769195
# ╟─8eff7d2a-b948-4201-bde8-0e6ee388d6f7
# ╟─74d33f87-6fb7-4567-9b9e-ccff956425fe
