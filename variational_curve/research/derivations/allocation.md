# 固定位数和熵约束下的码点分配

固定码本规模、熵约束和固定体积胞元是三种约束。逐项变分得到密度指数和局部曲率修正。

[研究索引](../README.md) · [符号和引用所用书目](../../ref.bib)

## Variational Codeword Allocation

<a id="eq-2-1"></a>

We next translate local distortion into an allocation of reproduction points. Let $\lambda(x)>0$ denote the number of codewords per unit Riemannian volume. A local cell then has volume approximately $\lambda(x)^{-1}$. For a specified cell family or local task weight, write $a(x)>0$ for its normalized distortion coefficient. The leading high-rate functional is

<a id="eq-2-2"></a>

$$
\mathcal D_0[\lambda]
=\int_\mathcal{M}p(x)a(x)\lambda(x)^{-2/d}\,\mathrm{d}V_{g}(x).
$$

For isotropic cells under geodesic squared error, $a(x)$ is the Euclidean cell-shape coefficient. The manifold quantization theorem gives the corresponding optimal leading constant [`iacobelli2016quantization`]; the allocation calculation below is the continuous high-rate design used to place a finite codebook [`gersho1979asymptotically`, `gray1998quantization`].

### Fixed-cardinality and entropy-constrained designs

**Theorem (Optimal high-rate codeword densities).**

<a id="eq-2-3"></a>

Let $p$ and $a$ be bounded and strictly positive.

1.  Under the fixed-cardinality constraint $\int_\mathcal{M}\lambda\,\mathrm{d}V_{g}=N$, the unique minimizer of [2.2](#eq-2-2) is

<a id="eq-2-4"></a> <a id="eq-2-5"></a>

$$
\begin{aligned}
\lambda_{\rm F}(x)
&=\frac{N}{Z_a}\bigl(p(x)a(x)\bigr)^{d/(d+2)},\\
Z_a&=\int_\mathcal{M}(pa)^{d/(d+2)}\,\mathrm{d}V_{g},\\
\mathcal D_0[\lambda_{\rm F}]
&=Z_a^{(d+2)/d}N^{-2/d}.
\end{aligned}
$$

2.  In the high-rate entropy model

<a id="eq-2-6"></a>

$$
H(\hat S)=h_g(p)+\mathbb{E}_p\ln\lambda(S)+o(1),
$$

the minimizer at entropy budget $R$ nats is

<a id="eq-2-7"></a> <a id="eq-2-8"></a>

$$
\begin{aligned}
\lambda_{\rm H}(x)
&=\exp\!\left(R-h_g(p)-\tfrac d2\mathbb{E}_p\ln a\right)a(x)^{d/2},
\\
\mathcal D_0[\lambda_{\rm H}]
&=\exp\!\left[-\tfrac2d(R-h_g(p))+\mathbb{E}_p\ln a\right].
\end{aligned}
$$

Here the $o(1)$ term in [2.6](#eq-2-6) is omitted in defining the high-rate optimization problem.

**Proof.** The fixed-cardinality functional is strictly convex in $\lambda$. Its Euler--Lagrange equation is $(2/d)pa\lambda^{-2/d-1}=\eta$; normalization gives [2.4](#eq-2-4) and [2.5](#eq-2-5). For entropy coding, a cell near $x$ has probability $p(x)/\lambda(x)$, which yields [2.6](#eq-2-6). Set $u=\ln\lambda$. The objective $\mathbb{E}_p[a\exp(-2u/d)]$ is strictly convex in $u$, and $\mathbb{E}_pu\le R-h_g(p)$ is affine. The active-constraint equation is $a\exp(-2u/d)=\text{constant}$. Solving for the constant gives [2.7](#eq-2-7) and [2.8](#eq-2-8).

For constant $a$, fixed-cardinality coding allocates more codewords where $p$ is large, with exponent $d/(d+2)$. Entropy coding gives equal local cell distortion and constant codeword density. These are distinct operational designs: $\lambda$ is a codeword count density, and its normalization under entropy coding is set by the entropy budget.

### Curvature in a local-cell allocation model

The next result isolates the geometric contribution at finite cell volume. It also fixes whether curvature is compared at a common radius or at a common communication resource.

**Proposition (Geodesic-ball distortion at fixed volume).**

<a id="eq-2-9"></a>

Let $V$ be the intrinsic volume of a small geodesic ball centered at $x$, with a source uniform inside that ball. Its mean squared geodesic distance from $x$ satisfies

<a id="eq-2-10"></a>

$$
D_{\rm ball}(x,V)=a_dV^{2/d}+b_d\mathrm{Scal}(x)V^{4/d}
+o(V^{4/d}),
$$

where

<a id="eq-2-11"></a>

$$
a_d=\frac d{d+2}\omega_d^{-2/d},
\qquad
b_d=\frac{\omega_d^{-4/d}}{3(d+2)(d+4)}.
$$

**Proof.** Integrating the normal-coordinate volume element gives 

$$
D_{\rm ball}(x,r)=\frac{dr^2}{d+2}
\left(1-\frac{2\mathrm{Scal}(x)r^2}{3d(d+2)(d+4)}+O(r^4)\right).
$$

Equation [1.10](geometry_and_rate.md#eq-1-10) implies 

$$
r^2=(V/\omega_d)^{2/d}
\left(1+\frac{\mathrm{Scal}(x)}{3d(d+2)}
(V/\omega_d)^{2/d}+O(V^{4/d})\right).
$$

Substitution yields [2.10](#eq-2-10).

For a uniform source and the local geodesic-ball cell model, [2.10](#eq-2-10) gives the second-order allocation functional, with $q=2/d$,

<a id="eq-2-12"></a>

$$
\mathcal D_{\rm ball}[\lambda]
=\int_\mathcal{M}p
\left[a_d\lambda^{-q}
+b_d\mathrm{Scal}\,\lambda^{-2q}\right]\mathrm{d}V_{g}.
$$

The cells are modeled locally by balls of the specified volume; a particular tessellation can use its own measured moments in the same optimization. Spatial source-density derivatives and cell shape must be included when calibrating that model for a nonuniform source.

At sufficiently high density, the first term in [2.12](#eq-2-12) dominates its second variation. Let $\lambda_0$ denote [2.4](#eq-2-4) with $a=a_d$ and define

$$
h_\lambda(x)=\frac{b_d}{a_d}\mathrm{Scal}(x)\lambda_0(x)^{-q},
\qquad
\bar h_\lambda=\frac1N\int_\mathcal{M}\lambda_0h_\lambda\,\mathrm{d}V_{g}.
$$

The local stationary solution with $\int\lambda=N$ has expansion

<a id="eq-2-13"></a>

$$
\frac{\lambda_{\rm C}(x)}{\lambda_0(x)}
=1+\frac{2}{1+q}
\bigl(h_\lambda(x)-\bar h_\lambda\bigr)+O(N^{-2q}).
$$

Indeed, differentiating [2.12](#eq-2-12) and writing $\lambda=\lambda_0(1+\varepsilon)$ gives $(1+q)\varepsilon=2h_\lambda-c$ to first order. The normalization condition determines $c=2\bar h_\lambda$.

For the uniform source, $\lambda_0=N/\mathrm{Vol}_g(\mathcal{M})$, so [2.13](#eq-2-13) allocates additional codewords where scalar curvature exceeds its volume average. Equation [2.10](#eq-2-10) explains the sign: positive curvature increases the radius needed to enclose a fixed intrinsic volume. The resulting density can initialize a manifold Lloyd quantizer, followed by nearest-codeword assignment and Fréchet-mean updates. Exact minimization in either step does not increase the measured quantization distortion. This provides an operational use of the local curvature calculation while retaining the classical leading high-rate allocation.

## Extrinsic Curvature and Local Cell Moments

<a id="eq-2-14"></a>

### The quartic chordal correction

Fix $x\in\mathcal{M}$ and a tangent vector $u$, and let $\gamma(s)=\operatorname{Exp}_x(su)$, viewed in the ambient Euclidean space. Because $\gamma$ is a geodesic of constant speed, 

$$
\begin{gathered}
\gamma'(0)=u,\qquad \gamma''(0)=\mathrm{I\!I}_x(u,u),\\
\langle u,\gamma'''(0)\rangle=-\|\mathrm{I\!I}_x(u,u)\|^2.
\end{gathered}
$$

The last equality follows by differentiating $\langle\gamma',\gamma''\rangle=0$. Taylor expansion therefore gives

<a id="eq-2-15"></a>

$$
\begin{aligned}
&\|\operatorname{Exp}_xu-x\|^2=|u|^2+\tfrac14\|\mathrm{I\!I}_x(u,u)\|^2\\
&\qquad+\tfrac13\langle u,\gamma'''(0)\rangle+O(|u|^5)\\
&\quad=|u|^2-\tfrac1{12}\|\mathrm{I\!I}_x(u,u)\|^2+O(|u|^5).
\end{aligned}
$$

For $U\sim\mathcal N(0,(2t)^{-1}I_d)$, the fourth-moment identity is

<a id="eq-2-16"></a>

$$
\mathbb{E}\|\mathrm{I\!I}(U,U)\|^2
=\frac{\|\mathbf H\|^2+2\|\mathrm{I\!I}\|_{\mathrm{HS}}^2}{4t^2}.
$$

Consequently the quartic change of exponent contributes $(\|\mathbf H\|^2+2\|\mathrm{I\!I}\|_{\mathrm{HS}}^2)/(48t)$ to the normalized kernel. Adding the volume correction yields

$$
\begin{aligned}
A_E&=-\frac{\mathrm{Scal}}{12}
+\frac{\|\mathbf H\|^2+2\|\mathrm{I\!I}\|_{\mathrm{HS}}^2}{48}\\
&=\frac{2\|\mathrm{I\!I}\|_{\mathrm{HS}}^2-\|\mathbf H\|^2}{16},
\end{aligned}
$$

which establishes [1.19](geometry_and_rate.md#eq-1-19) in every codimension. All odd terms integrate to zero in the local expansion. Smoothness and compactness give the uniform $O(t^{-2})$ remainder in [1.21](geometry_and_rate.md#eq-1-21).

### Radius and volume in a local cell model

<a id="eq-2-17"></a>

Let $C\subset T_x\mathcal{M}$ be a centrally symmetric reference cell with Euclidean volume $V_0$. Assume its second and radial fourth moments are isotropic, and write 

$$
m_2=\frac1{V_0}\int_C|u|^2du,\qquad
m_4=\frac1{V_0}\int_C|u|^4du.
$$

For $T_r=\operatorname{Exp}_x(rC)$, denote its intrinsic volume by $V_r$ and its mean squared radius under the uniform measure by 

$$
M_r=\frac1{V_r}\int_{T_r}d_{g}(x,y)^2\,\mathrm{d}V_{g}(y).
$$

Their local expansions are

<a id="eq-2-18"></a> <a id="eq-2-19"></a>

$$
\begin{aligned}
V_r&=V_0r^d\left[1-\frac{\mathrm{Scal}(x)m_2}{6d}r^2+O(r^4)\right],
\\
M_r&=m_2r^2-\frac{\mathrm{Scal}(x)}{6d}(m_4-m_2^2)r^4+O(r^6).
\end{aligned}
$$

These formulas follow by integrating the volume expansion against $1$ and $|u|^2$, respectively. In particular, $\int_C|u|^2\mathrm{Ric}(u,u)du=V_0m_4\mathrm{Scal}/d$.

To express the cell moment at a prescribed intrinsic volume $V_r$, put $w=(V_r/V_0)^{2/d}$. Inverting [2.18](#eq-2-18) gives $r^2=w+\mathrm{Scal}(x)m_2w^2/(3d^2)+O(w^3)$, and hence

<a id="eq-2-20"></a>

$$
M_r=m_2w+\frac{\mathrm{Scal}(x)}{6d}
\left[\frac{d+2}{d}m_2^2-m_4\right]w^2+O(w^3).
$$

For the unit-ball reference cell, $V_0=\omega_d$, $m_2=d/(d+2)$, and $m_4=d/(d+4)$, so this becomes

<a id="eq-2-21"></a>

$$
\frac d{d+2}\left(\frac{V_r}{\omega_d}\right)^{2/d}
+\frac{\mathrm{Scal}(x)}{3(d+2)(d+4)}
\left(\frac{V_r}{\omega_d}\right)^{4/d}
+O(V_r^{6/d}).
$$

The fixed-volume expression is the appropriate local input when codeword density specifies a cell volume. The coefficient also identifies the reference-cell moments required by a chosen allocation model.
