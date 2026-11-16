# 覆盖、可识别性与有限字典码率

来自 topology 分支的证明。先约定滤过的半径尺度、覆盖和 reach，再计算概念码率、描述符编解码与双尺度持久像。

[研究索引](../README.md) · [符号和引用所用书目](../../ref.bib)

## System Model and Communication Resources

<a id="eq-1-1"></a>

### Representations and structural distortion

A feature block is an indexed cloud $X=(x_i)_{i=1}^n\subset\mathbb{R}^{D_a}$, allowing repeated points. Noise results hold for every deterministic $X$; manifold inference additionally requires:

**Assumption (Geometric source).**

<a id="eq-1-2"></a>

**(A1)** The support $\mathcal{M}\subset\mathbb{R}^{D_a}$ is a compact smooth $d$-dimensional embedded manifold without boundary, with positive reach $\tau$. **(A2)** A certified lower bound $r\le\tau$ is available when a reconstruction guarantee is requested. **(A3)** When source probabilities are used, $S$ has a positive $C^1$ density $p$ with respect to induced volume, bounded above and below. Its connected-component index is denoted by $C(S)$.

A knowledge base (KB) specifies a representation and its filtration. For a fixed homological degree $k$, structural distortion is

<a id="eq-1-3"></a>

$$
d_{\rm topo}(X,\widehat X)
=d_{\mathrm{B}}\big(\mathrm{Dgm}_k(X),\mathrm{Dgm}_k(\widehat X)\big).
$$

For several degrees we use their maximum. Both sides use the same filtration (Section [1.5](#eq-1-5)). Component-label error and whole-block partition recovery are separate metrics: a correct component count can still group the wrong items.

### Dictionary and descriptor messages

A shared dictionary $P=\{p_1,\ldots,p_N\}$ gives a nearest-atom quantizer $Q_P(S)$. Each atom index occupies $R=\lceil\log_2N\rceil$ bits with a fixed-length code. The initial distribution of the dictionary is a setup cost; it is not included in this per-symbol rate. A KB update instead has length $B$ bits for a whole diagram, with source class and side information specified in Section [3.1](synchronization_and_reeb.md#eq-3-1).

### Calibrated feature channel

For a received feature block we analyze

<a id="eq-1-4"></a>

$$
Y_i=x_i+\xi_i,\qquad
\xi_i\overset{\rm iid}{\sim}
\mathcal N(0,\sigma^2 I_{D_a}/D_a).
$$

Thus $\sigma^2=\mathbb{E}\|\xi_i\|^2$ is the total noise energy per feature vector, and its estimate is available to a debiasing receiver. One physical realization sends $\alpha x_i$ through a real AWGN link with per-coordinate noise variance $\nu^2$, then divides by $\alpha$: $\sigma^2=D_a\nu^2/\alpha^2$. A block consumes $nD_a$ real channel uses; its mean transmit power per use is $P=\alpha^2\sum_i\|x_i\|^2/(nD_a)$, so $P/\nu^2=\sum_i\|x_i\|^2/(n\sigma^2)$. Dictionary and descriptor bits are source-message lengths; digital transport requires a channel code. Correctly decoded digital messages have no additive feature noise.

## Persistence and Three Geometric Lemmas

<a id="eq-1-5"></a>

### Filtrations, diagrams, distances

For finite $X\subset\mathbb{R}^{D_a}$, the Čech complex $\mathrm{C}_r(X)$ is the nerve of radius-$r$ balls and is homotopy equivalent to their union $X^r$. The Rips complex $\mathrm{R}_r(X)$ contains simplices with all edges at most $2r$. Both filtrations use the *radius convention* and have the same $H_0$; Jung's theorem gives

<a id="eq-1-6"></a>

$$
\mathrm{C}_r(X)\subseteq\mathrm{R}_r(X)\subseteq\mathrm{C}_{\sqrt2 r}(X)
$$

in any dimension [`edelsbrunner2010computational`]. For a compact manifold we use the offset filtration $\mathcal{M}^r$ directly. We work over $\mathbb F_2$. Finite-cloud modules are pointwise finite-dimensional and interval-decomposable. Compact-support offset and Rips diagrams are understood in the q-tame sense [`chazal2014persistence`]; zero-length and ephemeral classes are omitted. Intrinsic homology below the reach is finite-dimensional.

The module $(H_k(F_r))_{r\ge0}$ gives a diagram $\mathrm{Dgm}_k$ of bars $[b,d)$ with persistence $d-b$. Its map from scale $a$ to $b$ has rank equal to the number of bars born by $a$ and dying after $b$. Bottleneck distance minimizes the largest $L^\infty$ matching cost; an unmatched bar costs $(d-b)/2$. It equals interleaving distance $d_{\mathrm{I}}$ [`chazal2016structure`]. For functions $f,g$ on a common space,

<a id="eq-1-7"></a>

$$
d_{\mathrm{B}}\big(\mathrm{Dgm}_k(f),\mathrm{Dgm}_k(g)\big)\le\|f-g\|_\infty
$$

(stability [`cohensteiner2007stability`]). For point clouds with a common index set, [1.7](#eq-1-7) applied to the Rips filtration values on the full simplex gives $d_{\mathrm{B}}(\mathrm{Dgm}(\mathrm{R}X),\mathrm{Dgm}(\mathrm{R}Y))\le\max_i\|x_i-y_i\|$.

### Intrinsic features persist to the reach

**Lemma (Intrinsic persistence).**

<a id="eq-1-8"></a>

Let $\mathcal{M}$ be compact with reach $\tau$. In the offset filtration $(\mathcal{M}^r)_{r\ge0}$, every class of $H_*(\mathcal{M})$ is born at $0$ and survives at least to $\tau$, and every class born after $0$ is born at $r\ge\tau$. If $\mathcal{M}$ has components at pairwise distance at least $\Delta$, every finite $H_0$ bar dies at radius at least $\Delta/2$. If there are at least two components and their minimum distance is $\Delta_{\min}$, the first such death is $\Delta_{\min}/2$ and $\tau\le\Delta_{\min}/2$.

**Proof.** Below the reach, nearest-point projection gives compatible deformation retractions [`federer1959curvature`, `niyogi2008finding`], so all offset inclusions induce isomorphisms. Component offsets first meet at half their minimum separation. The minimizing pair's midpoint has two nearest points, giving the reach bound.

### Homology inference from a Rips complex

**Lemma (Rips inference).**

<a id="eq-1-9"></a>

Let $X\subset\mathbb{R}^{D_a}$ be finite with Hausdorff distance $d_H(X,\mathcal{M})\le\varepsilon$ to a compact manifold $\mathcal{M}$ of reach $\tau$, and let $\eta>\varepsilon$ with $(3+2\sqrt2)\eta<\tau$. Put $c_R=2+\sqrt2$. Then for every $k$, 

$$
\operatorname{rank}\Big(H_k\big(\mathrm{R}_\eta(X)\big)\to H_k\big(\mathrm{R}_{c_R\eta}(X)\big)\Big)=\beta_k(\mathcal{M}).
$$

If $\tilde X$ is any filtration on the same vertex set with $d_{\mathrm{B}}(\mathrm{Dgm}_k\tilde X,\mathrm{Dgm}_k\mathrm{R}X)\le\zeta$ for all $k$, the same holds for $\tilde X$ whenever $\varepsilon+2\zeta<\eta$ and $(3+2\sqrt2)\eta<\tau$.

**Proof.** Appendix [1.18](#eq-1-18) first bounds the persistent image of the Cech filtration using inclusions between offsets of $X$ and $\mathcal{M}$, then applies [1.6](#eq-1-6) and stability to obtain both statements.

The receiver can set $\eta=r/6$ from a reach lower bound $r$, since $3+2\sqrt2<6$. Persistent-image estimators of this type appear in [`chazal2008towards`].

### Covering numbers

$N_{\mathrm{cov}}(\mathcal{M},s)$ is the minimal number of points of $\mathcal{M}$ whose $s$-balls cover $\mathcal{M}$, and $N_{\mathrm{pack}}(\mathcal{M},s)$ the maximal number of points of $\mathcal{M}$ at pairwise distance greater than $s$; a maximal packing is a cover, so $N_{\mathrm{cov}}\le N_{\mathrm{pack}}$.

**Lemma (Covering numbers of a manifold with reach).**

<a id="eq-1-10"></a>

For $s\le\tau/2$, with $\omega_d$ the volume of the unit $d$-ball,

$$
\begin{aligned}
\frac{\mathrm{Vol}_g(\mathcal{M})}{\omega_d s^d}\,(1.1)^{-d}
&\le N_{\mathrm{cov}}(\mathcal{M},s)\\
&\le N_{\mathrm{pack}}(\mathcal{M},s)\le
\frac{\mathrm{Vol}_g(\mathcal{M})}{\omega_d (s/2)^d\cos^d\theta},
\end{aligned}
$$

with $\theta=\arcsin\big(s/(4\tau)\big)$.

**Proof.** Upper bound: the balls of radius $s/2$ around a packing are disjoint, and $\mathrm{Vol}_g(\mathcal{M}\cap B(x,s/2))\ge\cos^d\theta\,\omega_d(s/2)^d$ by [`niyogi2008finding`]. Lower bound: points of $\mathcal{M}$ within ambient distance $s\le\tau/2$ of $x$ are within geodesic distance $\rho=2\tau\arcsin(s/(2\tau))\le1.0108\,s$ [`boissonnat2019reach`]; by the Gauss equation the sectional curvature is at least $-2/\tau^2$, so Bishop's comparison [`petersen2016riemannian`] bounds the geodesic ball volume by $\omega_d\rho^d(\sinh(\sqrt2\rho/\tau)/(\sqrt2\rho/\tau))^{d-1}
\le\omega_d s^d(1.1)^d$.

In particular, $\log_2N_{\mathrm{cov}}(\mathcal{M},s)=\log_2\frac{\mathrm{Vol}_g(\mathcal{M})}{\omega_ds^d}+O(d)$.

## Topological Distortion and Dictionary Rate

<a id="eq-1-11"></a>

### A stable distortion for semantic structure

We use the following established maximality property to select a structural distortion. Here universal function stability means $\delta(\mathrm{Dgm}_k(f),\mathrm{Dgm}_k(g))\le\|f-g\|_\infty$ over all common spaces and functions with pointwise finite-dimensional modules, in a positive homological degree.

**Theorem (Bottleneck universality [`lesnick2015theory`]).**

<a id="eq-1-12"></a>

On diagrams of pointwise finite-dimensional one-parameter modules over a prime field, every universally function-stable extended pseudometric is bounded above by $d_{\mathrm{B}}$.

**Proof.** Lesnick's geometric lifting realizes the modules by functions whose sup-norm distance approaches their interleaving distance. Apply universal stability and $d_{\mathrm{I}}=d_{\mathrm{B}}$.

This motivates [1.3](#eq-1-3). Bottleneck distance is invariant under module isomorphism and is direct-sum consistent: unions of matchings bound a sum's distance by the maximum summand distance. For finite diagrams of sizes $n,m$, the unnormalized $p$-Wasserstein distance obeys $d_{\mathrm{B}}\le W_p\le(n+m)^{1/p}d_{\mathrm{B}}$, introducing a cardinality-dependent stability factor.

### The concept rate and a finite diagram code

**Proposition (Concept rate).**

<a id="eq-1-13"></a>

Let $C\in\{1,\ldots,N_c\}$ be the source concept, with $N_c\ge2$. For i.i.d. source blocks, a fixed-length code using at most $TR$ bits for $T$ symbols and attaining average concept error $P_c$ satisfies 

$$
R\ge H(C)-h_2(P_c)-P_c\log_2(N_c-1).
$$

For equiprobable concepts $H(C)=\log_2N_c$; for $N_c=1$ the rate is zero.

**Proof.** Data processing, the entropy chain rule, and symbolwise Fano give $TR\ge I(C^T;\widehat C^T)$ and $H(C^T\mid\widehat C^T)\le\sum_t[h_2(P_{c,t})+
P_{c,t}\log_2(N_c-1)]$. Concavity of $h_2$ yields the stated average bound [`cover2006elements`]. Lossless block source coding approaches $H(C)$ bits per symbol.

**Proposition (A finite diagram code).**

<a id="eq-1-14"></a>

Suppose the source diagram has at most $K_k$ finite bars in each $k\in\mathcal H$, with endpoints in $[0,L]$, $L>0$, where $K_k$ and $L$ are known to both terminals. For $\varepsilon>0$, put $J=\lceil L/(2\varepsilon)\rceil$ and $s=\lceil\log_2\max(1,J)\rceil$. There is a fixed-length code with error at most $\varepsilon$ in every degree and length

<a id="eq-1-15"></a>

$$
B_{\mathrm{list}}
=\sum_{k\in\mathcal H}\left[
\lceil\log_2(K_k+1)\rceil+2K_k s\right].
$$

Essential bars are agreed separately; their infinite endpoints are never rounded to finite values.

**Proof.** Quantize endpoints to the midpoints of $\max(1,J)$ equal bins. The error is at most $L/(2\max(1,J))\le\varepsilon$, and monotonicity preserves endpoint order. Remove bars whose indices coincide; their diagonal cost is at most $\varepsilon$. Send the remaining count and pad its endpoint list to $K_k$ slots. Matching retained bars by identity proves the claim.

For $K=\sum_k K_k$ and fixed degree set, the high-resolution length is $2K\log_2(L/\varepsilon)+O(K+\sum_k\log(K_k+1))$. Section [3.1](synchronization_and_reeb.md#eq-3-1) sharpens this upper bound when the receiver has a separated reference diagram.

### Dictionaries with a coverage certificate

Let $P=\{c_1,\ldots,c_N\}\subset\mathcal{M}$ be a dictionary and encode each symbol by its nearest atom. Its uniform reconstruction error is 

$$
e(P)=\sup_{x\in\mathcal{M}}\min_{c\in P}\|x-c\|.
$$

The receiver knows $e(P)\le\varepsilon$ and $r\le\tau(\mathcal{M})$. An index costs $\lceil\log_2N\rceil$ bits; dictionary distribution is a separate setup operation. Persistence can be computed at setup or after observing the activated atoms.

**Theorem (Certified dictionary rate).**

<a id="eq-1-16"></a>

Fix $0<\varepsilon<r/6$ and set $\eta=r/6$.

1.  Every dictionary with $e(P)\le\varepsilon$ satisfies 

$$
\beta_k(\mathcal{M})=\operatorname{rank}\Big(
H_k(\mathrm{R}_\eta(P))\to
H_k(\mathrm{R}_{c_R\eta}(P))\Big).
$$

2.  The minimum fixed-length index rate among dictionaries satisfying this uniform-error requirement is exactly 

$$
R_\varepsilon^*=\left\lceil\log_2N_{\mathrm{cov}}(\mathcal{M},\varepsilon)\right\rceil.
$$

Writing $A_\varepsilon=\log_2(\mathrm{Vol}_g(\mathcal{M})/(\omega_d\varepsilon^d))$ and $\theta_\varepsilon=\arcsin(\varepsilon/(4\tau))$, it obeys

$$
\begin{aligned}
A_\varepsilon-d\log_2(1.1)&\le R_\varepsilon^*,\\
R_\varepsilon^*&\le A_\varepsilon+d
-d\log_2\cos\theta_\varepsilon+1.
\end{aligned}
$$

Consequently, at a fixed $\varepsilon=cr$ with $0<c<1/6$, 

$$
R_\varepsilon^*=
\log_2\frac{\mathrm{Vol}_g(\mathcal{M})}{\omega_d r^d}+O(d)+O(1).
$$

**Proof.** Part (a) follows from $d_H(P,\mathcal{M})=e(P)$, Lemma [1.9](#eq-1-9), and $(3+2\sqrt2)r/6<r\le\tau$. For (b), the error requirement is exactly the covering condition. A minimum cover attains the necessary cardinality $N_{\mathrm{cov}}(\mathcal{M},\varepsilon)$; Lemma [1.10](#eq-1-10) and integer rounding give the bounds.

Part (b) concerns uniform reconstruction. A transmitter knowing the invariant can encode it directly at a lower rate. Atom indices alone need a coverage promise to identify the source topology.

**Proposition (Identifiability requires a coverage promise).**

<a id="eq-1-17"></a>

Let $d\ge1$, $D_a\ge d+1$, and let $P\subset\mathcal{M}$ be any finite set of distinct atoms on a compact smooth $d$-manifold of reach at least $r$. Among manifolds of reach at least $r$ containing $P$, there are two with different $\beta_0$ and strictly positive smooth source densities whose nearest-atom index sequences have exactly the same law. Under equal prior probabilities, every receiver has error at least $1/2$, even with infinitely many indices. More generally, for any $k\ge1$ there are $2^k$ such alternatives for which the largest possible probability of recovering $\beta_0$ is at most $2^{-k}\binom{k}{\lfloor k/2\rfloor}$.

**Proof.** Appendix [1.18](#eq-1-18) constructs separated spheres carrying extra components and transfers each sphere's atom probabilities to smooth densities inside the original Voronoi cells. This leaves the entire index law unchanged while changing the number of components.

**Remark (Finite observation time).** If the receiver sees only activated atoms and each atom has probability at least $q_{\min}>0$, after $T$ independent symbols the probability that any atom is absent is at most $N e^{-Tq_{\min}}$. Thus $T\ge q_{\min}^{-1}\log(N/\delta_s)$ exposes the full dictionary with probability at least $1-\delta_s$, a separate session-length cost.

## Proofs for Sections [1.5](#eq-1-5) and [1.11](#eq-1-11)

<a id="eq-1-18"></a>

### Lemma [1.9](#eq-1-9)

*Cech step.* If $d_H(X,\mathcal{M})\le\varepsilon$, then for $a\ge\varepsilon$, $\mathcal{M}^{a-\varepsilon}\subseteq X^a\subseteq\mathcal{M}^{a+\varepsilon}$. Choose $a'>\varepsilon$, $b'\ge a'+2\varepsilon$, and $b'+\varepsilon<\tau$. The inclusion chain

$$
\begin{aligned}
\mathcal{M}^{a'-\varepsilon}\subseteq X^{a'}
&\subseteq\mathcal{M}^{a'+\varepsilon}
\subseteq X^{a'+2\varepsilon}\\
&\subseteq X^{b'}\subseteq\mathcal{M}^{b'+\varepsilon}
\end{aligned}
$$

contains an inclusion between manifold offsets below the reach. Its homology map is an isomorphism, so $\operatorname{rank}(H_k(X^{a'})\to H_k(X^{b'}))\ge\beta_k(\mathcal{M})$. The same map factors through $H_k(\mathcal{M}^{a'+\varepsilon})$, giving the opposite inequality. The natural nerve equivalences transfer the rank statement to Cech complexes.

*Rips step.* The two inclusion chains

$$
\begin{aligned}
\mathrm{C}_a\subseteq\mathrm{R}_a\subseteq\mathrm{R}_b
&\subseteq\mathrm{C}_{\sqrt2b},\\
\mathrm{R}_a\subseteq\mathrm{C}_{\sqrt2a}\subseteq\mathrm{C}_b
&\subseteq\mathrm{R}_b
\end{aligned}
$$

give lower and upper bounds on $\operatorname{rank}(\mathrm{R}_a\to\mathrm{R}_b)$. The first is $\beta_k(\mathcal{M})$ when $a>\varepsilon$, $\sqrt2b\ge a+2\varepsilon$, and $\sqrt2b+\varepsilon<\tau$; the second is $\beta_k(\mathcal{M})$ when $b\ge\sqrt2a+2\varepsilon$ and $b+\varepsilon<\tau$. Set $a=\eta>\varepsilon$ and $b=(2+\sqrt2)\eta$. All these inequalities follow from $(3+2\sqrt2)\eta<\tau$.

*Perturbed version.* Write $\rho_X(a,b)$ for the persistent rank of $\mathrm{R}X$. If $d_{\mathrm{B}}\le\zeta$ and $b-a>2\zeta$, matching bars gives 

$$
\rho_X(a-\zeta,b+\zeta)
\le\rho_{\tilde X}(a,b)
\le\rho_X(a+\zeta,b-\zeta).
$$

A bar counted on either relevant side has length greater than $2\zeta$, so it cannot disappear into the diagonal. Each match changes birth and death by at most $\zeta$, giving the stated inequalities with the convention "birth at most $a$, death greater than $b$." For $a=\eta$, $b=(2+\sqrt2)\eta$, the Rips-step conditions for both outer pairs follow from

$$
\begin{aligned}
\eta-\zeta&>\varepsilon,\\
2\eta&>(1+\sqrt2)\zeta+2\varepsilon,\\
(2+2\sqrt2)\eta+\sqrt2\zeta+\varepsilon&<\tau.
\end{aligned}
$$

Indeed, $\varepsilon+2\zeta<\eta$ implies the first two and $\sqrt2\zeta+\varepsilon<\eta$; the last then follows from $(3+2\sqrt2)\eta<\tau$. Both outer ranks equal $\beta_k(\mathcal{M})$, proving the perturbed statement. $\blacksquare$

### Proposition [1.17](#eq-1-17)

Place $k$ round $d$-spheres $S_1,\ldots,S_k$ of radius $r$ in the ambient space, mutually separated and separated from $\mathcal{M}$ by more than $2r$. This is possible because $\mathcal{M}$ is compact and $D_a\ge d+1$. For every subset $A\subseteq\{1,\ldots,k\}$, define 

$$
\mathcal{M}_A=\mathcal{M}\cup\bigcup_{i\in A}S_i.
$$

A disjoint union of sets of reach at least $r$ and pairwise distance at least $2r$ has reach at least $r$: a point at distance less than $r$ cannot have nearest points in two components, and its nearest point within one component is unique. Thus all $\mathcal{M}_A$ satisfy the geometric promise, contain the same dictionary $P$, and have $\beta_0(\mathcal{M}_A)=\beta_0(\mathcal{M})+|A|$.

Fix a deterministic rule for nearest-atom ties, giving Voronoi cells $V_1,\ldots,V_N$. Since the atoms are distinct, $V_j\cap\mathcal{M}$ contains a neighborhood of $c_j$. Choose a nonnegative smooth bump $h_j$ supported strictly inside this neighborhood with unit integral over $\mathcal{M}$. Let $u_i$ be the uniform density on $S_i$, and set 

$$
\pi_{ij}=\int_{S_i\cap V_j}u_i\,dV,
\qquad h_i^{\circ}=\sum_{j=1}^{N}\pi_{ij}h_j.
$$

The density $h_i^{\circ}$ on $\mathcal{M}$ has exactly the same atom probabilities as $u_i$ on $S_i$. Some $\pi_{ij}$ may vanish, so these bumps need not be strictly positive.

Let $p_0$ be a strictly positive smooth probability density on $\mathcal{M}$. Choose weights $w_i>0$ with $w=\sum_iw_i<1$, and on $\mathcal{M}_A$ define

$$
\begin{aligned}
p_A|_{\mathcal{M}}&=(1-w)p_0+\sum_{i\notin A}w_i h_i^{\circ},\\
p_A|_{S_i}&=w_i u_i,\qquad i\in A.
\end{aligned}
$$

These densities integrate to one and are smooth and strictly positive on each compact component. Their density bounds may depend on the chosen alternative, as allowed by the source model. The probability of atom $j$ under every alternative is exactly 

$$
(1-w)\int_{\mathcal{M}\cap V_j}p_0\,dV
+\sum_{i=1}^{k}w_i\pi_{ij},
$$

which is independent of $A$. All i.i.d. atom-index sequences therefore have the same law. With $k=1$, the two possible values of $\beta_0$ are indistinguishable and the equal-prior error is at least $1/2$. For uniform $A$, the component-count increment has the $\mathrm{Binomial}(k,1/2)$ distribution and is independent of every observation. Its largest point probability is $2^{-k}\binom{k}{\lfloor k/2\rfloor}$, which bounds any receiver's success probability.$\blacksquare$
