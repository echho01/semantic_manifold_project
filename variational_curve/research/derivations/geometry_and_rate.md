# 参考几何、局部积分核与率失真展开

从源的正则性假设推到核的零、二、四阶矩，再用匹配上下界得到曲率系数。下文的附加推导分别给出内在损失与弦长损失的适用条件。

[研究索引](../README.md) · [符号和引用所用书目](../../ref.bib)

## System Model

<a id="eq-1-1"></a>

### Knowledge-assisted semantic representation

A fixed semantic encoder maps an observation to a representation $S\in\mathcal{M}$. The shared knowledge base specifies the support $\mathcal{M}$ and a metric for comparing representations. We study an independent and identically distributed sequence of these representations, with the encoder and the knowledge base held fixed. This makes the source distribution, the reconstruction criterion, and the communication resources explicit.

**Assumption (Semantic source).**

<a id="eq-1-2"></a>

The following conditions hold for the source-coding results.

1.  $\mathcal{M}^d\subset\mathbb{R}^{D_a}$ is a compact, connected, smooth embedded manifold without boundary. Its reach $\tau>0$ is the largest radius within which the nearest-point projection is unique [`federer1959curvature`].

2.  The metric $g$ is induced by the embedding. All geometric derivatives used below are bounded; in particular, sectional curvature is bounded and $\|\mathrm{I\!I}_x(u,u)\|\le\tau^{-1}$ for unit tangent vectors $u$.

3.  The source has a smooth density $p$ relative to $\mathrm{d}V_{g}$, with $0<p_{\min}\le p\le p_{\max}<\infty$.

Here $d$ is the intrinsic dimension and $D_a$ is the ambient dimension. The differential entropy relative to the Riemannian volume is

<a id="eq-1-3"></a>

$$
h_g(p)=-\int_\mathcal{M}p(x)\ln p(x)\,\mathrm{d}V_{g}(x).
$$

For an orthonormal tangent frame $\{e_i\}_{i=1}^d$, let $\mathbf H=\sum_i\mathrm{I\!I}(e_i,e_i)$ be the trace of the second fundamental form. We use the trace convention, without division by $d$. The Gauss equation gives $\mathrm{Scal}=\|\mathbf H\|^2-\|\mathrm{I\!I}\|_{\mathrm{HS}}^2$. Expectations of geometric quantities are taken under $p$.

### Source coding and reconstruction criteria

Reconstructions belong to $\mathcal{M}$. We consider the intrinsic loss $\delta_g(x,y)=d_{g}(x,y)^2$ and the ambient loss $\delta_E(x,y)=\|x-y\|^2$. Their rate--distortion functions are

<a id="eq-1-4"></a>

$$
R_\delta(D)=
\inf_{P_{\hat S|S}:\,\mathbb{E}\delta(S,\hat S)\le D}I(S;\hat S).
$$

All logarithms and mutual informations are in nats unless a subscript $2$ is used; $R_{\delta,2}=R_\delta/\ln2$ is in bits per source symbol. The high-resolution limit is $D\downarrow0$ with the source and metric fixed. The reproduction constraint $\hat S\in\mathcal{M}$ is part of [1.4](#eq-1-4) for both criteria.

The information rate in [1.4](#eq-1-4) describes asymptotically long source codes. For an implementable single-vector quantizer with $N$ reproduction points, the index requires $\lceil\log_2N\rceil$ bits. Section [2.1](allocation.md#eq-2-1) gives the corresponding high-rate allocation problem and its solution. For a memoryless physical channel of capacity $C$ bits per use and $b$ channel uses per source symbol, the usual separation condition is $R_{\delta,2}(D)<bC$ [`cover2006elements`].

### Projection reconstruction after a noisy channel

For the channel analysis we use the explicit analog reconstruction model

<a id="eq-1-5"></a>

$$
\begin{gathered}
\tilde S=S+N,\qquad N\sim\mathcal N(0,vI_{D_a}),\\
\sigma^2:=D_a v=\mathbb{E}\|N\|^2.
\end{gathered}
$$

The noise is independent of $S$. The receiver projects $\tilde S$ onto the support supplied by its knowledge base. Projection is unique inside the reach tube. A measurable nearest-point selection is used elsewhere; compactness bounds its intrinsic squared loss in the connected-source model. The small-noise limit is $v\downarrow0$ at fixed $D_a$. Thus the dimension-ratio gain below compares reconstruction errors at the same per-coordinate noise variance.

Section [3.1](projection_and_task.md#eq-3-1) also considers several separated components as a model of distinct semantic concepts. Component errors are then evaluated separately, with a finite task penalty. In Section [3.21](projection_and_task.md#eq-3-21), a task map specifies the information that the receiver needs, and a relative metric bound quantifies differences between the transmitter and receiver knowledge bases.

## Geometric Preliminaries

<a id="eq-1-6"></a>

Let $\operatorname{Exp}_x:T_x\mathcal{M}\to\mathcal{M}$ be the exponential map. In a normal neighborhood of $x$, a tangent vector $u$ represents the point $\operatorname{Exp}_xu$, and $d_{g}(x,\operatorname{Exp}_xu)=\|u\|$. Compactness permits a common neighborhood radius and uniform remainder estimates.

**Lemma (Normal-coordinate volume).**

<a id="eq-1-7"></a>

In geodesic normal coordinates,

<a id="eq-1-8"></a> <a id="eq-1-9"></a>

$$
\begin{aligned}
g_{ij}(u)&=\delta_{ij}-\tfrac13 R_{ikjl}u^ku^l+O(\|u\|^3),
\\
\mathrm{d}V_{g}(\operatorname{Exp}_xu)&=
\left(1-\tfrac16\mathrm{Ric}_x(u,u)+O(\|u\|^3)\right)\,\mathrm du.
\end{aligned}
$$

Consequently,

<a id="eq-1-10"></a>

$$
\mathrm{Vol}_g(B_r(x))
=\omega_d r^d
\left(1-\frac{\mathrm{Scal}(x)}{6(d+2)}r^2+O(r^4)\right),
$$

where $\omega_d$ is the Euclidean unit-ball volume.

**Proof.** The normal-coordinate metric and volume expansions are standard [`gray2004tubes`]. Rotational symmetry gives 

$$
\int_{\|u\|\le r}\mathrm{Ric}_x(u,u)\,\mathrm du
=\frac{\mathrm{Scal}(x)\omega_dr^{d+2}}{d+2}.
$$

Integrating [1.9](#eq-1-9) yields [1.10](#eq-1-10); odd Taylor terms integrate to zero on the centered ball.

**Lemma (Local ambient distortion).**

<a id="eq-1-11"></a>

Uniformly for $u$ near zero,

<a id="eq-1-12"></a>

$$
\|\operatorname{Exp}_xu-x\|^2
=\|u\|^2-\frac1{12}\|\mathrm{I\!I}_x(u,u)\|^2+O(\|u\|^5).
$$

Hence intrinsic and ambient squared distances have the same quadratic term, and their first difference is quartic.

**Proof.** For the unit-speed geodesic $\gamma$ with initial velocity $e$, $\gamma''(0)=\mathrm{I\!I}_x(e,e)$ and $\langle\gamma'(0),\gamma'''(0)\rangle
=-\|\gamma''(0)\|^2$. Expanding $\|\gamma(r)-\gamma(0)\|^2$ gives the fourth-order coefficient $1/4-1/3=-1/12$. Setting $u=re$ proves [1.12](#eq-1-12).

These two expansions produce the intrinsic and ambient kernel coefficients in Section [1.16](#eq-1-16). The normal-coordinate volume contributes scalar curvature; the quartic term in [1.12](#eq-1-12) adds the extrinsic curvature contribution.

**Lemma (Ambient Gaussian norm tail).**

<a id="eq-1-13"></a>

For $N\sim\mathcal N(0,vI_{D_a})$ and $\sigma^2=D_av$, define

<a id="eq-1-14"></a>

$$
T_{D_a}(r,\sigma)
=\Pr\{\chi^2_{D_a}\ge D_ar^2/\sigma^2\}.
$$

For $r>\sigma$,

<a id="eq-1-15"></a>

$$
T_{D_a}(r,\sigma)\le
\exp\!\left[-\frac{D_a}{2}
\psi\!\left(\frac{r^2}{\sigma^2}\right)\right],
\quad \psi(u)=u-1-\ln u.
$$

For $r\le\sigma$, the upper bound $1$ applies.

**Proof.** $\|N\|^2/v$ has a chi-square distribution with $D_a$ degrees of freedom. Its moment-generating function is $(1-2t)^{-D_a/2}$. Optimizing the Chernoff bound over $0<t<1/2$ yields [1.15](#eq-1-15).

## Curvature-Corrected Semantic Rate--Distortion

<a id="eq-1-16"></a>

### Source and reproduction geometry

Both the source $S$ and the reproduction $\hat S$ take values in $\mathcal{M}$. For the criteria $\delta_g(x,y)=d_{g}(x,y)^2$ and $\delta_E(x,y)=\|x-y\|^2$, let $R_g$ and $R_E$ denote the rate--distortion functions in [1.4](#eq-1-4). We use the Riemannian source entropy in [1.3](#eq-1-3).

For an orthonormal tangent frame $e_1,\ldots,e_d$, write

<a id="eq-1-17"></a>

$$
\mathbf H=\sum_{i=1}^d\mathrm{I\!I}(e_i,e_i),\qquad
\|\mathrm{I\!I}\|_{\mathrm{HS}}^2=
\sum_{i,j=1}^d\|\mathrm{I\!I}(e_i,e_j)\|^2.
$$

Thus $\mathbf H$ is the trace of the second fundamental form, without division by $d$. The Euclidean Gauss equation gives $\mathrm{Scal}=\|\mathbf H\|^2-\|\mathrm{I\!I}\|_{\mathrm{HS}}^2$. These invariants determine how the two distortion criteria depart from their common quadratic approximation.

### A common local kernel

For $t>0$ and a smooth function $f$, let

<a id="eq-1-18"></a>

$$
(K_{\delta,t}f)(x)=\int_\mathcal{M}e^{-t\delta(x,y)}f(y)\,\mathrm{d}V_{g}(y).
$$

Here $t$ controls the reconstruction scale, which is of order $t^{-1/2}$. Introduce the two scalar fields

<a id="eq-1-19"></a>

$$
A_g=-\frac{\mathrm{Scal}}{12},\qquad
A_E=\frac{2\|\mathrm{I\!I}\|_{\mathrm{HS}}^2-\|\mathbf H\|^2}{16}.
$$

**Lemma (Weighted local kernel expansion).**

<a id="eq-1-20"></a>

Under A1--A3, for either distortion criterion,

<a id="eq-1-21"></a>

$$
K_{\delta,t}f=
\left(\frac{\pi}{t}\right)^{d/2}
\left[f+\frac{1}{t}\left(\frac14\Delta_g f+A_\delta f\right)
+O(t^{-2})\right],
$$

uniformly on $\mathcal{M}$. We use $\Delta_g=\operatorname{div}_g\nabla_g$. Consequently, for $Z_{\delta,t}=K_{\delta,t}1$,

<a id="eq-1-22"></a>

$$
Z_{\delta,t}(x)=
\left(\frac{\pi}{t}\right)^{d/2}
\left[1+\frac{A_\delta(x)}{t}+O(t^{-2})\right],
$$

and the corresponding local mean distortion is

<a id="eq-1-23"></a>

$$
-\partial_t\ln Z_{\delta,t}(x)
=\frac{d}{2t}+\frac{A_\delta(x)}{t^2}+O(t^{-3}).
$$

**Proof.** In normal coordinates, the volume correction is $-\mathrm{Ric}(v,v)/6$. For geodesic distortion the exponent is exactly $-t|v|^2$. For ambient distortion the additional quartic term is $t\|\mathrm{I\!I}(v,v)\|^2/12$. Gaussian integration yields [1.21](#eq-1-21); the complete calculation is in Appendices [1.36](#eq-1-36) and [2.14](allocation.md#eq-2-14).

**Theorem (Curvature-corrected rate--distortion).**

<a id="eq-1-24"></a>

Under A1--A3, let $\overline A_\delta=\int_\mathcal{M}A_\delta(x)p(x)\,\mathrm{d}V_{g}(x)$. Then

<a id="eq-1-25"></a>

$$
R_\delta(D)=h_g(p)-\frac d2\ln\frac{2\pi eD}{d}
-\frac{2\overline A_\delta}{d}D+O(D^2),
\qquad D\downarrow0.
$$

In particular,

<a id="eq-1-26"></a> <a id="eq-1-27"></a>

$$
\begin{aligned}
R_g(D)&=h_g(p)-\frac d2\ln\frac{2\pi eD}{d}
+\frac{\overline{\mathrm{Scal}}_p}{6d}D+O(D^2),
\\
R_E(D)&=h_g(p)-\frac d2\ln\frac{2\pi eD}{d}\\
&\quad+\frac{\mathbb{E}_p[\|\mathbf H\|^2-2\|\mathrm{I\!I}\|_{\mathrm{HS}}^2]}{8d}D
+O(D^2),
\end{aligned}
$$

where $\overline{\mathrm{Scal}}_p=\mathbb{E}_p\mathrm{Scal}(S)$. The geodesic formula holds on any compact smooth Riemannian manifold without boundary with a smooth positive source density.

**Proof.** The kernel expansion supplies matching bounds. For the converse, choose a positive function $a_t$ satisfying $\int p(x)a_t(x)e^{-t\delta(x,y)}\mathrm{d}V_{g}(x)\le1$ for every $y$. Conditional relative entropy then gives $I(S;\hat S)\ge\mathbb{E}\ln a_t(S)-tD$. The choice 

$$
a_t(x)=\frac{(t/\pi)^{d/2}}{p(x)}
\exp\left(-\frac{A_\delta(x)}{t}-\frac{C}{t^2}\right)
$$

is feasible for a sufficiently large constant $C$. For achievability, use the source-weighted channel

<a id="eq-1-28"></a>

$$
q_t(y|x)=\frac{p(y)e^{-t\delta(x,y)}}{(K_{\delta,t}p)(x)}.
$$

Its distortion is $D_t=d/(2t)+\overline A_\delta/t^2+O(t^{-3})$. The density term in this expansion integrates to zero because $\int_\mathcal{M}\Delta_g p\,\mathrm{d}V_{g}=0$. Eliminating $t$ gives [1.25](#eq-1-25) from both bounds. Appendix [1.36](#eq-1-36) gives the feasibility argument and the mutual information calculation.

The leading logarithmic cost is determined by the intrinsic dimension. The next term distinguishes intrinsic curvature from bending in the feature space. The source weights these corrections through its probability density; its derivatives cancel from the averaged coefficient.

**Corollary (Constant-curvature sources).**

<a id="eq-1-29"></a>

For a compact space form with constant sectional curvature $\kappa$, the geodesic coefficient is $(d-1)\kappa/6$. If the source is uniform,

$$
R_g(D)=\ln\mathrm{Vol}_g(\mathcal{M})-\frac d2\ln\frac{2\pi eD}{d}
+\frac{(d-1)\kappa}{6}D+O(D^2).
$$

For the standard sphere $S^d(\rho)$, the ambient coefficient is $(d-2)/(8\rho^2)$.

**Proof.** A space form has $\mathrm{Scal}=d(d-1)\kappa$. On $S^d(\rho)$, $\|\mathrm{I\!I}\|_{\mathrm{HS}}^2=d/\rho^2$ and $\|\mathbf H\|^2=d^2/\rho^2$. Substitute these identities into Theorem [1.24](#eq-1-24).

### The rate cost of changing the distortion criterion

**Theorem (Intrinsic--ambient rate gap).**

<a id="eq-1-30"></a>

Under A1--A3,

<a id="eq-1-31"></a>

$$
R_g(D)-R_E(D)=
\frac{\mathbb{E}_p[\|\mathbf H\|^2+2\|\mathrm{I\!I}\|_{\mathrm{HS}}^2]}{24d}D
+O(D^2).
$$

Thus the two criteria have the same logarithmic and constant terms, and their difference is $O(D/\tau^2)$ for fixed intrinsic dimension.

**Proof.** Subtract [1.27](#eq-1-27) from [1.26](#eq-1-26) and use $\mathrm{Scal}=\|\mathbf H\|^2-\|\mathrm{I\!I}\|_{\mathrm{HS}}^2$. The reach bound on the second fundamental form controls both terms in the numerator. The nonnegative coefficient agrees with the pointwise ordering $\delta_E\le\delta_g$.

**Proposition (Surface coefficient).**

<a id="eq-1-32"></a>

For a compact smooth surface $\mathcal{M}^2\subset\mathbb{R}^3$ with principal curvatures $\kappa_1,\kappa_2$,

<a id="eq-1-33"></a>

$$
R_E(D)=h_g(p)-\ln(\pi eD)
-\frac1{16}\mathbb{E}_p[(\kappa_1-\kappa_2)^2]D+O(D^2).
$$

The correction vanishes on a sphere and is negative for a surface whose principal curvatures differ on a set of positive source measure.

**Proof.** For a surface, $\|\mathbf H\|^2-2\|\mathrm{I\!I}\|_{\mathrm{HS}}^2
=-(\kappa_1-\kappa_2)^2$. The result follows from [1.27](#eq-1-27).

### The fixed-codebook limit

A quantizer with $N$ codewords has the distinct operational objective 

$$
D_N=\inf_{|\mathcal C|\le N}
\int_\mathcal{M}p(x)\min_{c\in\mathcal C}d_{g}(x,c)^2\,\mathrm{d}V_{g}(x).
$$

Its leading asymptotics determine the codeword allocation in Section [2.1](allocation.md#eq-2-1).

**Proposition (Manifold quantization asymptotics).**

<a id="eq-1-34"></a>

Under A1--A3,

<a id="eq-1-35"></a>

$$
\lim_{N\to\infty}N^{2/d}D_N
=Q_d\left(\int_\mathcal{M}p^{d/(d+2)}\,\mathrm{d}V_{g}\right)^{(d+2)/d},
$$

where $Q_d$ is the Euclidean squared-error quantization constant.

**Proof.** Apply the Riemannian quantization theorem of Iacobelli [`iacobelli2016quantization`]. Compactness makes its moment and growth conditions automatic.

The fixed-codebook law has exponent $2/d$, with a constant determined by source density and cell geometry. Equation [1.25](#eq-1-25) optimizes mutual information over all test channels. These operational definitions provide the respective bases for source coding and codebook design [`gray1998quantization`].

## Matching Bounds for the Rate--Distortion Expansion

<a id="eq-1-36"></a>

### Uniform Laplace expansion

Fix a normal-coordinate radius below the injectivity radius of the compact manifold. The integral outside that neighborhood is exponentially small in $t$: geodesic distortion is bounded away from zero there, and the same holds for chordal distortion by compactness and injectivity of the embedding. Inside the neighborhood, write $y=\operatorname{Exp}_xu$. The normal-coordinate volume expansion [`gray2004tubes`] and the Taylor expansion of $f$ give

$$
\begin{aligned}
\mathrm{d}V_{g}(\operatorname{Exp}_xu)&=
\left(1-\tfrac16\mathrm{Ric}_x(u,u)+O(|u|^3)\right)du,\\
f(\operatorname{Exp}_xu)&=f(x)+\mathrm df_x(u)
+\tfrac12\mathrm{Hess}_g f_x(u,u)+O(|u|^3).
\end{aligned}
$$

The terms of odd total degree integrate to zero. Smoothness permits expansion through the fourth degree with a uniform remainder. The identities 

$$
\frac{\int u_i u_j e^{-t|u|^2}\,du}
{\int e^{-t|u|^2}\,du}=\frac{\delta_{ij}}{2t}
$$

give $\Delta_g f/(4t)$ and $-f\mathrm{Scal}/(12t)$. The ambient quartic term is calculated in Appendix [2.14](allocation.md#eq-2-14). Thus, with $c_t=(\pi/t)^{d/2}$ and $L_\delta f=\Delta_g f/4+A_\delta f$,

<a id="eq-1-37"></a>

$$
K_{\delta,t}f=c_t[f+t^{-1}L_\delta f+O(t^{-2})].
$$

The expansion holds uniformly for functions in bounded smooth families, which also allows its use below with $e^{-A_\delta/t}$. Applying the same expansion with an extra factor $\delta(x,y)$ gives

<a id="eq-1-38"></a>

$$
\begin{aligned}
&\int_\mathcal{M}\delta(x,y)e^{-t\delta(x,y)}f(y)\,\mathrm{d}V_{g}(y)\\
&\quad=c_t\left[\frac{d}{2t}f(x)
+\frac{d+2}{2t^2}L_\delta f(x)+O(t^{-3})\right].
\end{aligned}
$$

This establishes Lemma [1.20](#eq-1-20) and the moment estimate used in the achievable channel.

### A feasible dual function

<a id="eq-1-39"></a>

Let $a_t$ be positive and suppose

<a id="eq-1-40"></a>

$$
b_t(y):=\int_\mathcal{M}p(x)a_t(x)e^{-t\delta(x,y)}\,\mathrm{d}V_{g}(x)\le1
\quad\text{for all }y.
$$

For any reconstruction law of finite mutual information, apply nonnegativity of relative entropy to $P_{S|\hat S=y}$ and the probability density $p(x)a_t(x)e^{-t\delta(x,y)}/b_t(y)$. After averaging over $\hat S$, this gives

<a id="eq-1-41"></a>

$$
\begin{aligned}
I(S;\hat S)
&\ge\mathbb{E}\ln a_t(S)-t\mathbb{E}\delta(S,\hat S)
-\mathbb{E}\ln b_t(\hat S)\\
&\ge\mathbb{E}\ln a_t(S)-tD.
\end{aligned}
$$

The bound also holds when the mutual information is infinite.

Set 

$$
a_t(x)=\frac{c_t^{-1}}{p(x)}
\exp\left[-\frac{A_\delta(x)}{t}-\frac{C}{t^2}\right].
$$

The symmetry of the distortion and [1.37](#eq-1-37) imply 

$$
c_t^{-1}K_{\delta,t}(e^{-A_\delta/t})(y)
=1+O(t^{-2})
$$

uniformly in $y$: the $A_\delta/t$ term from the kernel cancels the one from the function. Choose $C$ to dominate the uniform upper remainder. Then $b_t(y)\le e^{-C/t^2}(1+C/t^2)\le1$ for all sufficiently large $t$. The resulting lower bound is

<a id="eq-1-42"></a>

$$
R_\delta(D)\ge h_g(p)+\frac d2\ln\frac t\pi
-\frac{\overline A_\delta}{t}-tD-\frac C{t^2}.
$$

Taking $t=d/(2D)$ gives the right-hand side of [1.25](#eq-1-25) with an $O(D^2)$ remainder from below.

### An achievable source-weighted channel

Use $q_t$ from [1.28](#eq-1-28). Since $p$ is strictly positive on a compact manifold, division in [1.37](#eq-1-37) is uniform. Equations [1.37](#eq-1-37) and [1.38](#eq-1-38) give

<a id="eq-1-43"></a>

$$
\mathbb{E}[\delta(S,\hat S)\mid S=x]
=\frac d{2t}+\frac{1}{t^2}
\left[\frac{\Delta_g p(x)}{4p(x)}+A_\delta(x)\right]+O(t^{-3}).
$$

Since $\mathcal{M}$ has no boundary, integration by parts gives $\int_\mathcal{M}\Delta_g p\,\mathrm{d}V_{g}=0$, and hence

<a id="eq-1-44"></a>

$$
D_t=\frac d{2t}+\frac{\overline A_\delta}{t^2}+O(t^{-3}).
$$

Let $r_t$ denote the reproduction density. The definition of $q_t$ gives the exact identity

<a id="eq-1-45"></a>

$$
\begin{aligned}
I_t+tD_t
&=-\int_\mathcal{M}p\ln(K_{\delta,t}p)\,\mathrm{d}V_{g}
-D_{\mathrm{KL}}(r_t\|p)\\
&\le h_g(p)+\frac d2\ln\frac t\pi
-\frac{\overline A_\delta}{t}+O(t^{-2}).
\end{aligned}
$$

Here the source-derivative term again integrates to zero. The distortion $D_t$ is continuous and strictly decreasing for large $t$; its derivative is the negative average conditional variance of $\delta(S,\hat S)$. Thus it can be chosen equal to any sufficiently small $D$, with

<a id="eq-1-46"></a>

$$
t(D)=\frac d{2D}+\frac{2\overline A_\delta}{d}+O(D).
$$

Substitution into [1.45](#eq-1-45) gives 

$$
I_{t(D)}\le h_g(p)-\frac d2\ln\frac{2\pi eD}{d}
-\frac{2\overline A_\delta}{d}D+O(D^2).
$$

Together with [1.42](#eq-1-42), this proves Theorem [1.24](#eq-1-24).
