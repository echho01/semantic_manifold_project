# 投影信道、任务约化与失配

从法管内投影的扰动展开算到四阶系数，再给出任务充分性和度量失配界。

[研究索引](../README.md) · [符号和引用所用书目](../../ref.bib)

## Projection Decoding and Semantic Reliability

<a id="eq-3-1"></a>

The channel model is [1.5](geometry_and_rate.md#eq-1-5): $v$ is the per-coordinate variance and $\sigma^2=D_av$ is the total noise energy. We first derive the local reconstruction error, then consider transmission between separated semantic components.

### Projection distortion and dimension gain

**Theorem (Projection-noise expansion).**

<a id="eq-3-2"></a>

Under A1--A3, let $\hat S$ be a measurable nearest-point projection of $S+N$. As $v\downarrow0$,

<a id="eq-3-3"></a> <a id="eq-3-4"></a>

$$
\begin{aligned}
\mathbb{E}d_{g}(S,\hat S)^2
&=dv+B_pv^2+o(v^2),\\
B_p&=\mathbb{E}_p\left[\|\mathrm{I\!I}\|_{\mathrm{HS}}^2-\tfrac23\mathrm{Scal}\right].
\end{aligned}
$$

Consequently, the reconstruction-error ratio is

<a id="eq-3-5"></a>

$$
G(v):=\frac{\mathbb{E}\|N\|^2}{\mathbb{E}d_{g}(S,\hat S)^2}
=\frac{D_a}{d}
\left(1-\frac{B_p}{d}v+o(v)\right).
$$

**Proof.** The differential of nearest-point projection at $x\in\mathcal{M}$ is the orthogonal projection onto $T_x\mathcal{M}$, yielding $dv$. The cubic expansion of the projection map and the quartic expansion of intrinsic squared distance give $\|\mathrm{I\!I}_x\|_{\mathrm{HS}}^2-\frac23\mathrm{Scal}(x)$ after Gaussian averaging. Appendix [3.14](#eq-3-14) gives the calculation. For any fixed $r<\tau$, the event $\|N\|\ge r$ contributes an exponentially small term because the intrinsic loss is bounded on compact connected $\mathcal{M}$. Averaging over the source gives [3.3](#eq-3-3).

The leading ratio $D_a/d$ compares the raw and projected reconstruction errors in the specified analog model. The normal noise components are removed to first order. The correction $B_p$ describes the interaction between tangent noise and the embedding curvature.

**Corollary (Round spheres).**

<a id="eq-3-6"></a>

For $\mathcal{M}=S^d(\rho)$, whose reach is $\rho$,

<a id="eq-3-7"></a>

$$
\mathbb{E}d_{g}(S,\hat S)^2
=dv+\frac{d(5-2d)}{3\rho^2}v^2+o(v^2).
$$

In particular, the fourth-order coefficient is positive for $d=1,2$ and negative for $d\ge3$.

**Proof.** On the sphere, $\|\mathrm{I\!I}\|_{\mathrm{HS}}^2=d/\rho^2$ and $\mathrm{Scal}=d(d-1)/\rho^2$. Substitution into [3.4](#eq-3-4) proves the result.

In decibels, [3.5](#eq-3-5) becomes 

$$
10\log_{10}(D_a/d)-\frac{10B_pv}{d\ln10}+o(v).
$$

The sphere family therefore supplies a direct test of both the leading dimension ratio and the sign of the curvature term.

### Separated concepts and component errors

For this subsection, let $\mathcal{M}=\bigcup_{i=1}^{N_c}\mathcal{M}_i$ be a finite union of smooth compact $d$-dimensional components without boundary. Each component satisfies the geometric regularity in A1--A2, and the union has positive reach $\tau$. The pairwise Euclidean separation is at least $\Delta_0>0$. The decoder chooses a nearest point in the union. Write $\mathcal E_{\rm c}$ for a decoded component different from the source component and $\mathcal E_{\rm t}$ for leaving the reach tube. Component error carries a fixed finite penalty $d_\times$.

**Theorem (Component and projection-reliability bounds).**

<a id="eq-3-8"></a>

Under [1.5](geometry_and_rate.md#eq-1-5),

<a id="eq-3-9"></a> <a id="eq-3-10"></a>

$$
\begin{aligned}
\Pr(\mathcal E_{\rm c})
&\le T_{D_a}(\Delta_0/2,\sigma),\\
\Pr(\mathcal E_{\rm c}\cup\mathcal E_{\rm t})
&\le T_{D_a}(m,\sigma),\qquad
m=\min\{\tau,\Delta_0/2\}.
\end{aligned}
$$

Here $T_{D_a}$ is the chi-square tail in [1.14](geometry_and_rate.md#eq-1-14), and [1.15](geometry_and_rate.md#eq-1-15) gives an exponential upper bound when its radius exceeds $\sigma$.

**Proof.** Since $S$ is a candidate reconstruction, $\|\hat S-\tilde S\|\le\|S-\tilde S\|=\|N\|$. Thus $\|\hat S-S\|\le2\|N\|$. A component error requires $\|\hat S-S\|\ge\Delta_0$, and hence $\|N\|\ge\Delta_0/2$. Also $\operatorname{dist}(\tilde S,\mathcal{M})\le\|N\|$, so leaving the reach tube requires $\|N\|\ge\tau$. These event inclusions prove the bounds.

The half-separation is the decoding margin between components. It gives an explicit noise budget for a prescribed error probability: solve $T_{D_a}(\Delta_0/2,\sigma)\le\epsilon_{\rm c}$. The geometry of a particular decision boundary can yield a tighter expression. For two concentric circles, component selection depends only on the received radius; Section [numerical evidence](../EXPERIMENTS.md) evaluates its exact noncentral chi-square probability.

**Corollary (Rayleigh-fading bound).**

<a id="eq-3-11"></a>

Let $U=|h|^2\sim\mathrm{Exp}(1)$ be independent of $N$, and let the equalized perturbation be $N/\sqrt U$. For any radius $r>0$,

<a id="eq-3-12"></a>

$$
\Pr\{\|N\|/\sqrt U\ge r\}
=1-\left(1+\frac{2\sigma^2}{D_ar^2}\right)^{-D_a/2}.
$$

Using $r=\Delta_0/2$ bounds component error; using $r=m$ bounds component error or departure from the reach tube. At high SNR, the right-hand side is $\sigma^2/r^2+O(\sigma^4/r^4)$.

**Proof.** Condition on $N$. The event has probability $1-\exp(-\|N\|^2/r^2)$. The Laplace transform of the chi-square distribution gives [3.12](#eq-3-12). Expanding at $\sigma^2=0$ gives the last claim.

### Combining quantization and reconstruction error

Let $S_{\rm q}$ be a quantized representative in the source component, and transmit $S_{\rm q}$ through [1.5](geometry_and_rate.md#eq-1-5). Define $D_{\rm q}=\mathbb{E}d_{g}(S,S_{\rm q})^2$ and 

$$
D_{\rm ch}=\mathbb{E}\!\left[
d_{g}(S_{\rm q},\hat S)^2\mathbf1_{\mathcal E_{\rm c}^{\,c}}\right],
$$

using the within-component distance on the event of correct decoding. For the loss equal to squared within-component distance on success and $d_\times$ on a component error, we obtain

<a id="eq-3-13"></a>

$$
D_{\rm total}\le
\left(\sqrt{D_{\rm q}}+\sqrt{D_{\rm ch}}\right)^2
+d_\times\Pr(\mathcal E_{\rm c}).
$$

The triangle inequality on each component followed by Minkowski's inequality proves [3.13](#eq-3-13). At high SNR, $D_{\rm ch}=dv+O(v^2)$.

Equations [2.5](allocation.md#eq-2-5), [3.3](#eq-3-3), and [3.9](#eq-3-9) give a joint resource calculation. The codebook size sets $D_{\rm q}$; the channel variance sets $D_{\rm ch}$; and component separation sets the allowable probability of a semantic decision error. For a target distortion, these equations allow the codebook and noise budgets to be selected together.

## Fourth-Order Projection-Noise Expansion

<a id="eq-3-14"></a>

Fix $x\in\mathcal{M}$ and decompose the ambient noise as $N=(\xi,\eta)$ in $T_x\mathcal{M}\oplus N_x\mathcal{M}$, with independent components $\xi\sim\mathcal N(0,vI_d)$ and $\eta\sim\mathcal N(0,vI_{D_a-d})$. All quantities in this appendix are evaluated at $x$. Choose a fixed sufficiently small radius $r<\tau$. On $\{\|N\|<r\}$ the nearest-point projection is smooth. The complementary event contributes an exponentially small amount to any bounded reconstruction loss, uniformly in $x$. Local-coordinate expectations below are taken on this event; replacing truncated Gaussian polynomial moments by their full moments changes only an exponentially small term. Write $\hat S_x$ for the projection receiver with the bounded extension specified in Section [3.1](#eq-3-1).

### Expansion of the nearest point

Use tangent coordinates to write the manifold locally as a graph $F(u)=(u,\phi(u))$, where

<a id="eq-3-15"></a>

$$
\phi_\alpha(u)=\tfrac12u^TB_\alpha u
+\tfrac16C_\alpha(u,u,u)+O(|u|^4).
$$

Each $B_\alpha$ is the symmetric matrix of the second fundamental form in normal direction $\alpha$; $C_\alpha$ is a symmetric cubic tensor. Define 

$$
B_\eta=\sum_\alpha\eta_\alpha B_\alpha,\qquad
C_\eta=\sum_\alpha\eta_\alpha C_\alpha.
$$

For a cubic tensor, $C_\eta(\xi,\xi)$ denotes the vector obtained by leaving its first index uncontracted.

Let $u$ be the tangent coordinate of $\pi_{\mathcal{M}}(x+N)$. The first-order condition for nearest-point projection is 

$$
u-\xi+D\phi(u)^T[\phi(u)-\eta]=0.
$$

Solving by homogeneous degree in $(\xi,\eta)$ gives

<a id="eq-3-16"></a>

$$
\begin{aligned}
u={}&\xi+B_\eta\xi+B_\eta^2\xi
+\tfrac12C_\eta(\xi,\xi)\\
&-\tfrac12\sum_\alpha B_\alpha\xi
(\xi^TB_\alpha\xi)+O(\|N\|^4).
\end{aligned}
$$

The quadratic and cubic terms both contribute to the fourth-order mean squared error.

### Tangent-coordinate mean squared error

Write $B=\|\mathrm{I\!I}\|_{\mathrm{HS}}^2$ and $H_2=\|\mathbf H\|^2$. Gaussian moment identities yield

$$
\begin{aligned}
\mathbb{E}|B_\eta\xi|^2&=v^2B,\\
2\mathbb{E}\langle\xi,B_\eta^2\xi\rangle&=2v^2B,\\
\sum_\alpha\mathbb{E}(\xi^TB_\alpha\xi)^2&=v^2(H_2+2B).
\end{aligned}
$$

The term involving $C_\eta$ has zero expectation by independence and centering of $\eta$. The total cubic term also has zero expectation. It follows from [3.16](#eq-3-16) that

<a id="eq-3-17"></a>

$$
\mathbb{E}|u|^2=dv+(B-H_2)v^2+O(v^3).
$$

To justify the stated remainder, expand the smooth squared error through degree five. Every odd-degree polynomial has zero expectation on the symmetric event $\{\|N\|<r\}$, and the remainder is bounded by a constant times $\|N\|^6$. Gaussian tails account for the complement.

### Conversion to geodesic error

Let $w$ be the geodesic normal coordinate of the same point $F(u)$. The tangent part of the geodesic Taylor expansion used in [2.15](allocation.md#eq-2-15) is 

$$
u=w-\tfrac16\sum_\alpha B_\alpha w(w^TB_\alpha w)+O(|w|^4).
$$

Consequently,

<a id="eq-3-18"></a>

$$
d_{g}(x,F(u))^2=|u|^2+\tfrac13\|\mathrm{I\!I}(u,u)\|^2+O(|u|^5).
$$

The fourth-order term can be evaluated at $u=\xi$, giving 

$$
\mathbb{E}\|\mathrm{I\!I}(\xi,\xi)\|^2=v^2(H_2+2B).
$$

Combining this identity with [3.17](#eq-3-17) and applying the same parity argument to the remainder yields

<a id="eq-3-19"></a>

$$
\begin{aligned}
\mathbb{E}d_{g}(x,\hat S_x)^2
&=dv+\left(\tfrac53B-\tfrac23H_2\right)v^2+O(v^3)\\
&=dv+\left(\|\mathrm{I\!I}\|_{\mathrm{HS}}^2
-\tfrac23\mathrm{Scal}\right)v^2+O(v^3).
\end{aligned}
$$

This proves the expansion in Section [3.1](#eq-3-1). Compactness makes the remainder uniform, so the formula can also be averaged over any distribution of $x$.

For $S^d(\rho)$, $B=d/\rho^2$ and $\mathrm{Scal}=d(d-1)/\rho^2$, hence

<a id="eq-3-20"></a>

$$
\mathbb{E}d_{g}(S,\hat S)^2=dv+
\frac{d(5-2d)}{3\rho^2}v^2+O(v^3),
\qquad v=\frac{\sigma^2}{D_a}.
$$

The coefficient depends on the intrinsic dimension, while adding orthogonal ambient coordinates changes only the relation between per-coordinate variance $v$ and total variance $\sigma^2$.

## Task-Oriented Compression and Metric Mismatch

<a id="eq-3-21"></a>

### Exact reduction to the task source

Let $f:\mathcal{M}^d\to(\mathcal N^{d_f},h)$ be a smooth surjective submersion, where $\mathcal N$ is compact and without boundary. The receiver needs $Y_f=f(S)$. The task distortion is 

$$
\delta_f(s,\hat s)=d_h(f(s),f(\hat s))^2.
$$

The pushforward source $Y_f$ has density $p_f$ relative to the volume of $h$. Assume $p_f$ is smooth and strictly positive.

**Theorem (Task-source rate--distortion).**

<a id="eq-3-22"></a>

For every distortion level,

<a id="eq-3-23"></a>

$$
R_{\rm task}(D)=R^{\mathcal N}_{p_f}(D).
$$

For $d_f\ge1$, its high-resolution expansion is

<a id="eq-3-24"></a>

$$
\begin{aligned}
R_{\rm task}(D)
={}&h_h(p_f)-\frac{d_f}{2}
\ln\frac{2\pi eD}{d_f}\\
&+\frac{\mathbb{E}_{p_f}[\mathrm{Scal}_h]}{6d_f}D+o(D).
\end{aligned}
$$

Consequently,

<a id="eq-3-25"></a>

$$
\begin{gathered}
R_g(D)-R_{\rm task}(D)
=\frac{d-d_f}{2}\ln\frac1D+C_f\\
+\frac D6\left(\frac{\mathbb{E}_p[\mathrm{Scal}_g]}d
-\frac{\mathbb{E}_{p_f}[\mathrm{Scal}_h]}{d_f}\right)+o(D),
\end{gathered}
$$

where 

$$
C_f=h_g(p)-h_h(p_f)
-\frac d2\ln\frac{2\pi e}{d}
+\frac{d_f}{2}\ln\frac{2\pi e}{d_f}.
$$

**Proof.** For any reconstruction $\hat S$, data processing gives $I(S;\hat S)\ge I(f(S);f(\hat S))$, proving the lower bound in [3.23](#eq-3-23). Conversely, use a test channel from $Y_f$ to $\hat Y_f$ and reconstruct a fixed measurable representative of the fiber $f^{-1}(\hat Y_f)$. This preserves the task distortion and introduces no information about the original fiber coordinate. Thus $I(S;\hat S)=I(Y_f;\hat Y_f)$. The intrinsic expansion of Theorem [1.24](geometry_and_rate.md#eq-1-24), applied on $\mathcal N$, gives [3.24](#eq-3-24); subtraction gives [3.25](#eq-3-25). Appendix [3.30](#eq-3-30) records the measurable construction and the comparison argument.

The fibers identify representations that are equivalent for the task. The pullback tensor $f^*h$ has those fiber directions as its kernel, and compression can discard them. Equation [3.25](#eq-3-25) expresses the saving through the task dimension, task-source entropy, and task-space curvature. For example, on a product $\mathcal{M}=\mathcal N\times\mathcal F$ with task $f(y,z)=y$, only the $\mathcal N$ coordinate needs to be described. Task-oriented information-bottleneck communication uses the same principle of retaining task-relevant information [`shao2022taskoriented`]; [3.23](#eq-3-23) specifies the distortion-based source-coding limit.

### Relative metric mismatch

Suppose the source support and coordinates are shared but the reconstruction criterion used to design the codec is $g_R$, while performance is evaluated using $g_T$. Assume

<a id="eq-3-26"></a>

$$
(1-\varepsilon_g)g_T\preceq g_R\preceq(1+\varepsilon_g)g_T,
\qquad 0\le\varepsilon_g<1.
$$

This is a relative comparison of quadratic forms on each tangent space. It is invariant to a change of local coordinates.

**Theorem (Finite-rate metric comparison).**

<a id="eq-3-27"></a>

Under [3.26](#eq-3-26),

<a id="eq-3-28"></a>

$$
R_R((1+\varepsilon_g)D)\le R_T(D)
\le R_R((1-\varepsilon_g)D).
$$

For a fixed codebook $\mathcal C$, let $Q_R$ and $Q_T$ assign a source point to a nearest codeword under $g_R$ and $g_T$, respectively. Then

<a id="eq-3-29"></a>

$$
D_T(Q_T)\le D_T(Q_R)
\le\frac{1+\varepsilon_g}{1-\varepsilon_g}D_T(Q_T).
$$

The same factor bounds the true distortion of a codec optimized under $g_R$ relative to the optimal true distortion at the same codebook cardinality or information rate.

**Proof.** Integrating the quadratic-form comparison along a path and minimizing over paths gives $(1-\varepsilon_g)d_T^2\le d_R^2\le(1+\varepsilon_g)d_T^2$. This inclusion of feasible test channels proves [3.28](#eq-3-28). For each $x$, 

$$
d_T(x,Q_Rx)^2
\le\frac{d_R(x,Q_Rx)^2}{1-\varepsilon_g}
\le\frac{1+\varepsilon_g}{1-\varepsilon_g}d_T(x,Q_Tx)^2.
$$

Integration proves [3.29](#eq-3-29). Comparing the $g_R$-optimal codec with a $g_T$-optimal codec gives the last statement; arbitrarily near-optimal codecs cover cases where the infimum is not attained.

The excess distortion is bounded by $2\varepsilon_gD_T^*/(1-\varepsilon_g)$ at the specified resource budget. Thus improving the codebook or increasing the information rate reduces both the matched distortion and the excess caused by this metric perturbation. At high resolution, a sufficient additional rate for the factor in [3.29](#eq-3-29) is 

$$
\frac d2\ln\frac{1+\varepsilon_g}{1-\varepsilon_g}+o(1)
\quad\text{nats per source symbol}.
$$

This provides a quantitative tolerance for updating a learned metric while keeping the representation support and codebook interface fixed.

## Task Reduction and Metric Comparison

<a id="eq-3-30"></a>

### Exact task rate--distortion reduction

<a id="eq-3-31"></a>

Let $f:\mathcal{M}\to\mathcal N$ be the smooth surjective submersion of Theorem [3.22](#eq-3-22), and put $Z=f(S)$ and $\hat Z=f(\hat S)$. For every reconstruction channel, 

$$
I(S;\hat S)\ge I(Z;\hat Z),\qquad
\mathbb{E}d_{\mathcal N}(Z,\hat Z)^2
=\mathbb{E}d_{\mathcal N}(f(S),f(\hat S))^2.
$$

Taking the infimum establishes $R_{\mathrm{task}}(D)\ge R_{f_\#p}^{\mathcal N}(D)$.

For the reverse bound, a smooth submersion admits a smooth local section around every point of $\mathcal N$. Choose a finite cover by such neighborhoods and select the first applicable local section. This defines a Borel measurable map $s:\mathcal N\to\mathcal{M}$ satisfying $f(s(z))=z$. For any admissible task-space test channel $P_{\hat Z|Z}$, generate $\hat Z$ from $Z=f(S)$ and reconstruct $\hat S=s(\hat Z)$. The Markov chain $S\longrightarrow Z\longrightarrow\hat Z$ and the deterministic relations between $\hat Z$ and $s(\hat Z)$ give 

$$
I(S;\hat S)=I(S;\hat Z)=I(Z;\hat Z).
$$

The distortion is preserved. Taking test channels approaching the task-space infimum proves the exact identity

<a id="eq-3-32"></a>

$$
R_{\mathrm{task}}(D)=R_{f_\#p}^{\mathcal N}(D)
\quad\text{for every }D.
$$

The reconstruction rule supplies a representative of each fiber without communicating a fiber coordinate.

To apply Theorem [1.24](geometry_and_rate.md#eq-1-24) on $\mathcal N$, its source density with respect to the task-space volume is given by the coarea formula [`federer1959curvature`]:

<a id="eq-3-33"></a>

$$
p_Z(z)=\int_{f^{-1}(z)}\frac{p(x)}{J_f(x)}\,
\mathrm dV_{f^{-1}(z)}(x),
$$

where $J_f(x)=\sqrt{\det(\mathrm df_x\mathrm df_x^*)}$ and the adjoint uses the source and task metrics. A proper smooth submersion and a smooth positive $p$ give a smooth positive $p_Z$. Compactness makes its upper and lower bounds uniform. Applying [1.26](geometry_and_rate.md#eq-1-26) with dimension $d_f=\dim\mathcal N$ proves the task expansion.

In particular, subtraction of the full-source and task expansions gives

<a id="eq-3-34"></a>

$$
\begin{aligned}
R_g(D)-R_{\mathrm{task}}(D)
={}&\frac{d-d_f}{2}\ln\frac1D+C_f\\
&+\left[\frac{\overline{\mathrm{Scal}}_p}{6d}
-\frac{\mathbb{E}_{p_Z}\mathrm{Scal}_{\mathcal N}}{6d_f}\right]D+O(D^2),
\end{aligned}
$$

where

$$
\begin{aligned}
C_f={}&h_g(p)-h_{\mathcal N}(p_Z)
-\frac d2\ln\frac{2\pi e}{d}
+\frac{d_f}{2}\ln\frac{2\pi e}{d_f}.
\end{aligned}
$$

### A uniform metric comparison

<a id="eq-3-35"></a>

Suppose $0\le\varepsilon<1$ and $(1-\varepsilon)g_T\preceq g_R\preceq(1+\varepsilon)g_T$. Integrating the corresponding norm inequalities along every curve and taking infima over curves gives

<a id="eq-3-36"></a>

$$
(1-\varepsilon)d_{g_T}(x,y)^2
\le d_{g_R}(x,y)^2
\le(1+\varepsilon)d_{g_T}(x,y)^2.
$$

For a fixed codebook $\mathcal C$, let $Q_T(x)$ and $Q_R(x)$ be nearest-codeword assignments under the respective metrics. Then

$$
\begin{aligned}
d_{g_T}(x,Q_R(x))^2
&\le\frac{d_{g_R}(x,Q_R(x))^2}{1-\varepsilon}\\
&\le\frac{d_{g_R}(x,Q_T(x))^2}{1-\varepsilon}\\
&\le\frac{1+\varepsilon}{1-\varepsilon}
d_{g_T}(x,Q_T(x))^2.
\end{aligned}
$$

Averaging proves the codeword-selection bound. The same factor applies to an entire design optimized for $g_R$. Indeed, let $\mathcal A_R$ denote a common set of admissible reconstruction laws at rate $R$, and let $D_T^*(R)$ be the minimum true distortion over this set. Compare a $g_R$-optimal law to a $g_T$-optimal law using [3.36](#eq-3-36) twice. The resulting true distortion satisfies

<a id="eq-3-37"></a>

$$
D_T^*(R)\le D_{\mathrm{mis}}(R)
\le\frac{1+\varepsilon}{1-\varepsilon}D_T^*(R).
$$

Infimizing sequences give the same conclusion if a minimizer is not attained. This proves Theorem [3.27](#eq-3-27).

For the mutual-information rate constraint, write $\chi=(1+\varepsilon)/(1-\varepsilon)$. A sufficient rate for true distortion $D$ is $R_g(D/\chi)$. The high-resolution rate penalty therefore satisfies

<a id="eq-3-38"></a>

$$
\Delta R(D)\le R_g(D/\chi)-R_g(D)
=\frac d2\ln\chi+O(D)
$$

in nats for fixed $\varepsilon$. Division by $\ln2$ gives the corresponding expression in bits.
