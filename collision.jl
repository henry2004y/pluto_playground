### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ bfa63490-6151-11ee-2f43-45656d4df512
md"
# Collision and the Raise of π

## Problem description

Imagine there are two blocks with mass M and m in an ideal system: all collisions are elastic, no resistive force. The M block has an initial velocity of u, while the m block has an initial velocity of 0. There is a wall on one side. What is the total number of collisions between the two blocks?

## Solution

Before and after the nth collision, we have the convervation of momentum and energy:
```math
\begin{aligned}
M u_n + m v_n &= M u_{n+1} + m v_{n+1} \\
\frac{1}{2}Mu_n^2 + \frac{1}{2}mv_n^2 &= \frac{1}{2}Mu_{n+1}^2 + \frac{1}{2}mv_{n+1}^2
\end{aligned}
```

The solution at time step n+1 can be decided from time step n:
```math
\begin{pmatrix}
u_{n+1} \\ v_{n+1}
\end{pmatrix} = 
\begin{pmatrix}
\frac{M-m}{M+m} & - \frac{2m}{M+m} \\
\frac{2M}{M+m} & \frac{M-m}{M+m}
\end{pmatrix}
\begin{pmatrix}
u_n \\ v_n
\end{pmatrix}
```

Let $\alpha=m/M$, then
```math
\begin{pmatrix}
u_{n+1} \\ v_{n+1}
\end{pmatrix} = (1+\alpha)^{-n}
\begin{pmatrix}
1-\alpha & -2\alpha \\
2        & 1-\alpha
\end{pmatrix}^n
\begin{pmatrix}
u \\ 0
\end{pmatrix}
```

Let matrix $A = \begin{pmatrix}
1-\alpha & -2\alpha \\
2        & 1-\alpha
\end{pmatrix}$. Our goal is to compute $A^n$.

_Eigenvalue Decomposition_
```math
\mathrm{det}(A - \lambda I) = 0 \rightarrow \lambda_\pm = (1\pm i\sqrt{\alpha})^2
```

For $\lambda_+$, can take right eigenvector $x_+ = \begin{pmatrix} i\sqrt{\alpha} \\ 1 \end{pmatrix}$.

For $\lambda_-$, can take right eigenvector $x_- = \begin{pmatrix} -i\sqrt{\alpha} \\ 1 \end{pmatrix}$.

If we choose $P = \begin{pmatrix} x_+ & x_- \end{pmatrix}$ and $\Lambda = \begin{pmatrix} \lambda_+ & \\ & \lambda_- \end{pmatrix}$, then the eigenvalue decomposition gives
```math
A = P\Lambda P^{-1}
```
where
```math
P^{-1} = \frac{1}{\mathrm{det}P}\begin{pmatrix} 1 & i\sqrt{\alpha} \\ -1 & i\sqrt{\alpha} \end{pmatrix} = \frac{1}{2i\sqrt{\alpha}} \begin{pmatrix} 1 & i\sqrt{\alpha} \\ -1 & i\sqrt{\alpha} \end{pmatrix}
```

Therefore $A^n = (P\Lambda P^{-1})^n = P \Lambda^n P^{-1}$, and finally we have
```math
\begin{pmatrix}
u_{n+1} \\ v_{n+1}
\end{pmatrix} = (1+\alpha)^{-n}\, A^n \begin{pmatrix} u \\ 0 \end{pmatrix} =
\frac{(1+\alpha)^{-n}}{2}u \begin{pmatrix} (1+i\sqrt{\alpha})^{2n} + (1-i\sqrt{\alpha})^{2n} \\ \frac{1}{i\sqrt{\alpha}}\left[ (1+i\sqrt{\alpha})^{2n} - (1-i\sqrt{\alpha})^{2n} \right] \end{pmatrix}
```
from which one can read off
```math
\begin{aligned}
u_n &= (1+\alpha)^{-n} \left[ (1+i\sqrt{\alpha})^{2n} + (1-i\sqrt{\alpha})^{2n} \right] u \\
v_n &= (1+\alpha)^{-n}\frac{1}{i\sqrt{\alpha}} \left[ (1+i\sqrt{\alpha})^{2n} - (1-i\sqrt{\alpha})^{2n} \right] u
\end{aligned}
```

Although the expressions for $u_n$ and $v_n$ contain a bunch of $i\sqrt{\alpha}$, it can be shown that the imaginary parts in both expressions vanish, so $u_n$ and $v_n$ are always real as they should be.

After the last collision, we would expect $u_n < 0$ and $|u_n| > v_n$. This gives the criterion
```math
-u_n > v_n \rightarrow - \left[ (1+i\sqrt{\alpha})^{2n} + (1-i\sqrt{\alpha})^{2n} \right] > \frac{1}{i\sqrt{\alpha}} \left[ (1+i\sqrt{\alpha})^{2n} - (1-i\sqrt{\alpha})^{2n} \right]
```

Assume that $M \gg m$, so $\alpha = m/M \approx 0$. We can use the approximation
```math
1\pm i\sqrt{\alpha} \approx e^{\pm i\sqrt{\alpha}}
```
while inserting into the inequity above gives
```math
\begin{aligned}
- \left[ e^{i2n\sqrt{\alpha}} + e^{-i2n\sqrt{\alpha}} \right] &> \frac{1}{i\sqrt{\alpha}} \left[ e^{i2n\sqrt{\alpha}} - e^{-i2n\sqrt{\alpha}} \right] \\
-\cos(2n\sqrt{\alpha}) &> \frac{1}{\sqrt{\alpha}}\sin(2n\sqrt{\alpha})
\end{aligned}
```

Note that $\cos(2n\sqrt{\alpha}) < 0$, when we divide both sides of the inequility by $\cos(2n\sqrt{\alpha})$, we should flip the sign and obtain
```math
\tan(2n\sqrt{\alpha}) > -\sqrt{\alpha}
```

Again let's notice that we shall restrict ourselves to the region where $\cos(2n\sqrt{\alpha}) < 0$, i.e. $\frac{\pi}{2}<2n\sqrt{\alpha}<\frac{3}{2}\pi$.

Since $\sqrt{\alpha} \approx 0$, so $2n\sqrt{\alpha} = 2n\sqrt{\frac{m}{M}}\approx \pi$, we get the number of collisions between blocks
```math
n \approx \frac{1}{2}\sqrt{\frac{M}{m}}\pi
```

If we also count the collisions between the m block and the wall,
```math
N = 2n \approx \sqrt{\frac{M}{m}}\pi
```

"

# ╔═╡ f738fb5d-a31a-484d-95ff-d2fc394affb8
begin
	m = 1.0
	M = 100.0
	α = M / m
	u = 1.0
	v = 0.0
end

# ╔═╡ d37eea2a-5337-4ba9-98b8-e0f4bed3ca58
function trial(α, u, v, M, m)
	n = 0
	while -u < v
		u_next = (M-m)/(M+m)*u - 2m/(M+m)*v
		v_next = 2M/(M+m)*u + (M-m)/(M+m)*v
		u, v = u_next, v_next
		n += 1
	end

	n
end

# ╔═╡ 4d6b99e4-d374-4787-ac08-1a65a15f5619
begin
	n = trial(α, u, v, M, m)
	println("number of collisions: ", n)
	println("mass ratio: ", α)
	println("2n/sqrt(α) = ", 2n/√(α))
end

# ╔═╡ dd888f9c-d997-4cec-978d-05979c5f8af5
md"
## Results from numerical simulation

M/m     | n | N | difference to $\pi$
:------- | :----: | :--------: | :----:
1        | 2  | 4       | 27%
100      | 16 | 3.2     | 1.9%
10000    | 157 | 3.14   | 0.051% 
1000000  | 1571 | 3.142 | 0.013%

"

# ╔═╡ de1ef8dd-7f98-4063-8f5c-89746e96e1cd
md"""
## Thinking in the phase space

This is not the end of the story. Given the fact that $\pi$ arises from the solution, we should immediately think of any periodic motion in a circle. In this specific case,
it is the phase space consists with x-axis being $x=\sqrt{M}u$, and y-axis being $y=\sqrt{m}v$. Then the conservation of energy means that a state of the system consists of two blocks must be sitting on a circle, and the conservation of momentum means that the state can only jump along a line with a given slope that depends on the mass ratio!

See [the video of 3Blue1Brown](https://youtu.be/jsYwFizhncE?si=r27PS1C9-51PPA6W) for a detailed explanation.
"""

# ╔═╡ Cell order:
# ╟─bfa63490-6151-11ee-2f43-45656d4df512
# ╠═f738fb5d-a31a-484d-95ff-d2fc394affb8
# ╠═d37eea2a-5337-4ba9-98b8-e0f4bed3ca58
# ╠═4d6b99e4-d374-4787-ac08-1a65a15f5619
# ╟─dd888f9c-d997-4cec-978d-05979c5f8af5
# ╟─de1ef8dd-7f98-4063-8f5c-89746e96e1cd
