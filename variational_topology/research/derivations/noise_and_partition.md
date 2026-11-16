# 点云噪声下的滤过重参数化与分区恢复

推导距离集中、去偏与接触阈值，并说明什么条件足以恢复一个有身份的分区。

[研究索引](../README.md) · [符号和引用所用书目](../../ref.bib)

## Channel Noise: Reparametrization and Receiver Reliability

<a id="eq-2-1"></a>

Condition on the clean indexed cloud in [1.4](coverage_and_dictionary.md#eq-1-4), allowing repeated points, and write $y_i=Y_i=x_i+N_i$. Throughout, $D_k(Z)=\mathrm{Dgm}_k(\mathrm{R}Z)$ and $a=\log(2n^2/\delta)/D_a$, with $n\ge2$, $0<\delta<1$, and natural logarithms. Following the high-dimensional persistence analyses in [`damrich2024persistent`, `hiraoka2026curse`], we retain the clean distances in a shifted reference to obtain finite-sample receiver guarantees.

### Uniform stability and the shifted reference filtration

**Theorem (Worst-case stability).**

<a id="eq-2-2"></a>

With probability at least $1-\delta$, simultaneously in every degree $k$,

$$
\begin{aligned}
d_{\mathrm{B}}\big(D_k(Y),D_k(X)\big)
&\le\max_i\|N_i\|\\
&\le\sigma\Big(1+\sqrt{2\log(n/\delta)/D_a}\Big).
\end{aligned}
$$

**Proof.** Each edge radius changes by at most $\max_i\|N_i\|$, so filtration stability applies. Gaussian concentration of the norm, with $\mathbb{E}\|N_i\|\le\sigma$, and a union bound give the second inequality.

This displacement bound remains of order $\sigma$. Define the sharper noise-dependent reference

<a id="eq-2-3"></a>

$$
\begin{gathered}
\varphi(r)=\sqrt{r^2+\sigma^2/2},\\
F_r^\sigma=\left\{\omega:\max_{i,j\in\omega,\,i\ne j}
\frac{\sqrt{\|x_i-x_j\|^2+2\sigma^2}}2\le r\right\},
\end{gathered}
$$

with all vertices present at zero. Let $D_k^\sigma=\mathrm{Dgm}_k(F^\sigma)$. For $k\ge1$, it transforms each bar $(b,d)$ to $(\varphi(b),\varphi(d))$. In $H_0$, births stay at zero and deaths become $\varphi(d)$. This requires all $n-1$ merge radii, including the zero radii from duplicate samples that an ordinary diagram omits.

**Theorem (Reparametrization with a finite-sample error).**

<a id="eq-2-4"></a>

With probability at least $1-\delta$, simultaneously for every degree $k$,

$$
\begin{aligned}
d_{\mathrm{B}}\big(D_k(Y),D_k^\sigma\big)&\le\zeta_{\mathrm{raw}},\\
\zeta_{\mathrm{raw}}&=\sigma\big((2+\sqrt2)\sqrt a+\sqrt2\,a\big).
\end{aligned}
$$

On the same event, for every pair, with $l_{ij}=\|x_i-x_j\|$,

<a id="eq-2-5"></a>

$$
\left|\|y_i-y_j\|^2-l_{ij}^2-2\sigma^2\right|
\le\alpha l_{ij}+\beta,
$$

where $\alpha=4\sigma\sqrt a$ and $\beta=4\sigma^2(\sqrt a+a)$.

**Proof.** Appendix [2.16](#eq-2-16) combines Gaussian and chi-square concentration [`laurent2000adaptive`], bounds the edge-radius error, and applies filtration stability.

The reference compresses a positive-dimensional bar's lifetime to $(d^2-b^2)/(\varphi(d)+\varphi(b))$ and raises every finite $H_0$ death to at least $\sigma/\sqrt2$. A fixed threshold can therefore fail even when the error relative to $D_k^\sigma$ is small.

### Debiasing and its precision requirement

The receiver removes the known squared-distance offset using

<a id="eq-2-6"></a>

$$
\hat d_{ij}=\sqrt{(\|y_i-y_j\|^2-2\sigma^2)_+},\qquad \hat d_{ii}=0.
$$

Although this dissimilarity need not satisfy the triangle inequality, its clique filtration obeys simplexwise stability.

**Corollary (Debiased persistence).**

<a id="eq-2-7"></a>

On the event of Theorem [2.4](#eq-2-4),

$$
\begin{aligned}
d_{\mathrm{B}}\big(D_k(\hat d),D_k(X)\big)&\le\zeta_{\mathrm{deb}},\\
\zeta_{\mathrm{deb}}&=\frac{\alpha+\sqrt{\alpha^2+4\beta}}4.
\end{aligned}
$$

As $a\to0$, $\zeta_{\mathrm{deb}}=\sigma a^{1/4}(1+O(a^{1/4}))$. If a variance estimate is used and $|\hat\sigma^2-\sigma^2|\le v$, the same result holds with $\beta$ replaced by $\beta+2v$.

**Proof.** The positive-part map is $1$-Lipschitz, so $|\hat d_{ij}^2-l_{ij}^2|\le\alpha l_{ij}+\beta$. Set $r_*=(\alpha+\sqrt{\alpha^2+4\beta})/2$. For $\hat d\ge l$, the error is at most $\max(\alpha/2,\sqrt\beta)\le r_*$. For $\hat d<l$, it is at most $\min(l,\alpha+\beta/l)\le r_*$, with $l=0$ treated separately. Halving the distance error gives the radius error $\zeta_{\mathrm{deb}}$. The estimated offset introduces at most $2v$ of additional squared-distance error before clipping, which proves the last statement.

The fourth-root rate comes from square-root inversion near zero clean distance; positive clean separations have a smaller inversion error. Calibration uncertainty $2v$ consumes the same budget as random fluctuations.

**Corollary (Certified dictionary rate over the analog channel).**

<a id="eq-2-8"></a>

Let the transmitted cloud contain every atom of the dictionary in Theorem [1.16](coverage_and_dictionary.md#eq-1-16). With $\eta=r/6$, the debiased rank receiver reports $\beta_*(\mathcal{M})$ with probability at least $1-\delta$ if $\varepsilon+2\zeta_{\mathrm{deb}}<\eta$. Any $0<\varepsilon<\eta-2\zeta_{\mathrm{deb}}$ is therefore admissible, with source index length at most $\lceil\log_2N_{\mathrm{pack}}(\mathcal{M},\varepsilon)\rceil$.

**Proof.** Apply Corollary [2.7](#eq-2-7) and the perturbed version of Lemma [1.9](coverage_and_dictionary.md#eq-1-9). Repetitions do not change clean coverage; $n$ counts all independently perturbed transmitted occurrences.

This radius contraction adds $d\log_2[\eta/(\eta-2\zeta_{\mathrm{deb}})]$ to the radius-dependent term of the volume rate bound, with its constants and rounding retained. For $2\zeta_{\mathrm{deb}}\ge\eta$ the certificate gives no positive covering budget.

### Fragmentation, cross-contact, and the inferred partition

An $H_0$ receiver counts threshold-graph components: one plus the number of finite deaths exceeding $\theta$. Let the clean merge radii be at most $s$ for sampling bars and at least $\gamma$ for concept bars, with $s<\theta<\gamma$.

**Theorem (Threshold reliability).**

<a id="eq-2-9"></a>

With probability at least $1-\delta$:

1.  If $\sigma/\sqrt2-\zeta_{\mathrm{raw}}>\theta$, the raw receiver reports $n$ singletons. It reports the clean concept count if $\varphi(s)+\zeta_{\mathrm{raw}}<\theta<\varphi(\gamma)-\zeta_{\mathrm{raw}}$.

2.  The debiased receiver reports the clean concept count if $\zeta_{\mathrm{deb}}<\min(\theta-s,\gamma-\theta)$.

**Proof.** The edge-radius error bounds each ordered $H_0$ merge radius by the same amount, including zero radii (Appendix [2.16](#eq-2-16)). The first condition in (i) puts every raw edge above $\theta$.

At fixed $n,\theta$ and $D_a\to\infty$, the raw fragmentation scale approaches $\sigma=\sqrt2\theta$ when the clean sampling radii are zero. Debiased thresholding at $\theta\ge0$ selects exactly the raw edges at $\varphi(\theta)$, correcting threshold calibration while retaining pairwise fluctuations.

To test the partition itself, consider $m$ samples at each of two centers $x_A,x_B$ separated by $\Delta>0$, with $n=2m$. Let $r_{\rm in}$ be the larger within-concept connectivity radius (the maximum minimum-spanning-tree edge radius), and $r_{\rm out}=\frac12\min_{i\in A,j\in B}\hat d_{ij}$. The true partition is recovered exactly when

<a id="eq-2-10"></a>

$$
r_{\rm in}\le\theta<r_{\rm out}.
$$

A cross edge violates this partition even if two components remain; zero-radius collapse additionally requires every vertex to be connected.

**Proposition (Finite-sample cross-contact bounds).**

<a id="eq-2-11"></a>

For $\sigma>0$ and $0\le\theta<\Delta/2$, define

<a id="eq-2-12"></a>

$$
p_\theta=F_{D_a,\lambda}\!\left(D_a(1+2\theta^2/\sigma^2)\right),
$$

where $F_{D,\lambda}$ is the noncentral chi-square CDF with $D$ degrees of freedom and $\lambda=D_a\Delta^2/(2\sigma^2)$. Then

<a id="eq-2-13"></a>

$$
1-(1-p_\theta)^m
\le\Pr(r_{\rm out}\le\theta)
\le\min(1,m^2p_\theta).
$$

Suppose $m\to\infty$, $\log m=o(D_a)$, and $\theta/\Delta<1/2$ is fixed. Writing 

$$
S_\theta=\sqrt{\Delta^2-4\theta^2}\,(D_a/\log m)^{1/4},
$$

we have $\Pr(r_{\rm out}\le\theta)\to0$ if $\sigma=cS_\theta$ with $c<2^{-5/4}$, and this probability tends to one if $c>1/2$.

**Proof.** For one cross pair, $D_a\|y_i-y_j\|^2/(2\sigma^2)$ has the distribution in [2.12](#eq-2-12). A union bound over all $m^2$ pairs gives the upper bound. A matching of $m$ cross pairs uses disjoint noise vectors, so those events are independent and gives the lower bound. The asymptotics follow from the moderate-deviation estimate [2.17](#eq-2-17) proved in Appendix [2.16](#eq-2-16).

The zero-contact scale is $\Delta(D_a/\log m)^{1/4}$. Assuming all $m^2$ pairs independent would select $2^{-5/4}$ as a transition constant; the proven bounds retain their dependence. At $\theta=\Delta/4$, $\zeta_{\mathrm{deb}}<\Delta/4$ ensures the true partition, while Proposition [2.11](#eq-2-11) supplies a cross-contact obstruction. Together these establish the fourth-root noise scale for this receiver.

### The value of a known knowledge base

**Corollary (Matched reliability comparison).**

<a id="eq-2-14"></a>

The nearest-center receiver knowing $x_A,x_B$ has per-point error $q(\sigma)=Q(\Delta\sqrt{D_a}/(2\sigma))$ and block error $1-(1-q(\sigma))^n$. For $0<\delta<1/2$, its largest noise amplitude for block error at most $\delta$ is

<a id="eq-2-15"></a>

$$
\sigma_{\rm KB}(\delta)=
\frac{\Delta\sqrt{D_a}}
{2Q^{-1}(1-(1-\delta)^{1/n})}.
$$

For fixed $0<\delta<1/2$ in the joint limit of Proposition [2.11](#eq-2-11), this scales as $\Delta\sqrt{D_a/\log n}$. Relative to the threshold graph receiver, the ratio of the noise scales at block reliability is therefore of order $(D_a/\log n)^{1/4}$.

**Proof.** An error requires a noise displacement at least $\Delta/2$ along $x_B-x_A$. Independence gives the block expression; solving it and using $Q^{-1}(u)\sim\sqrt{2\log(1/u)}$ proves the remaining claims.

At fixed *per-point* reliability the known-KB scale is instead $\Delta\sqrt{D_a}$, giving ratio $(D_a\log n)^{1/4}$. These comparisons concern the specified threshold graph; spectral persistence [`damrich2024persistent`] and distance-to-measure constructions [`chazal2011geometric`] can use additional global information.

## Channel Concentration and Contact Asymptotics

<a id="eq-2-16"></a>

### Theorem [2.4](#eq-2-4)

Put $v=x_i-x_j$, $l=\|v\|$, and $G=N_i-N_j\sim\mathcal N(0,2\sigma^2I_{D_a}/D_a)$. Expanding $\|v+G\|^2$, the cross-term bound is 

$$
\Pr\{|2\langle v,G\rangle|>4l\sigma\sqrt{t/D_a}\}\le2e^{-t},
$$

with the zero vector handled identically. Since $D_a\|G\|^2/(2\sigma^2)\sim\chi^2_{D_a}$, the two tails in [`laurent2000adaptive`] give 

$$
\Pr\!\left\{\left|\|G\|^2-2\sigma^2\right|>
4\sigma^2\left(\sqrt{t/D_a}+t/D_a\right)\right\}\le2e^{-t}.
$$

A union bound over all pairs with $t=\log(2n^2/\delta)$ gives [2.5](#eq-2-5), with failure probability $2n(n-1)e^{-t}<\delta$. For $\sigma=0$ the bounds are deterministic.

Writing $L=\|y_i-y_j\|$ and $L_\varphi=\sqrt{l^2+2\sigma^2}$,

$$
\begin{aligned}
|L-L_\varphi|
&\le\frac{4l\sigma\sqrt a+4\sigma^2(\sqrt a+a)}{L_\varphi}\\
&\le4\sigma\sqrt a+2\sqrt2\,\sigma(\sqrt a+a)=2\zeta_{\mathrm{raw}}.
\end{aligned}
$$

Halving distances gives edge radii; taking maxima gives simplex radii. Thus the filtration functions differ by at most $\zeta_{\mathrm{raw}}$, proving the bottleneck bound.

The strictly increasing reference transformation preserves nonvertex simplex order and persistence pairings. Kruskal's $n-1$ edges give the transformed $H_0$ deaths, including originally zero radii. In higher degrees, zero-length bars have both endpoints transformed together and remain on the diagonal.

### Ordered merge radii and partition recovery

For a symmetric nonnegative dissimilarity matrix $D$ with zero diagonal, let $G_D(r)$ contain edges with $D_{ij}/2\le r$. Maximum edge-radius error $\zeta$ gives $G_D(r)\subseteq G_{D'}(r+\zeta)$ and the reverse inclusion with $D,D'$ exchanged. Component-count jumps therefore shift by at most $\zeta$, proving stability of all $n-1$ ordered merge radii, including zeros, and hence Theorem [2.9](#eq-2-9).

The graphs induced by $A,B$ are connected exactly for $r\ge r_{\rm in}$, and there is no cross edge exactly for $r<r_{\rm out}$, proving [2.10](#eq-2-10).

### Moderate deviations for Proposition [2.11](#eq-2-11)

Write $D=D_a$, $v=\Delta^2-4\theta^2>0$, and $\sigma=c\sqrt v(D/\log m)^{1/4}$ with fixed $c>0$. Let $Z\sim\chi^2_D(\lambda)$, $\lambda=D\Delta^2/(2\sigma^2)$, and 

$$
b=\frac{Dv}{2\sigma^2}=\frac{\sqrt{D\log m}}{2c^2}.
$$

The event in [2.12](#eq-2-12) is $Z-(D+\lambda)\le-b$. Because $\log m=o(D)$, $\lambda/D\to0$, $b/D\to0$, and $b/\sqrt D\to\infty$. The centered cumulant generating function is

$$
\begin{aligned}
K(t)&=-\tfrac D2\log(1-2t)+\frac{\lambda t}{1-2t}
-t(D+\lambda)\\
&=(D+2\lambda)t^2+O((D+\lambda)|t|^3)
\end{aligned}
$$

uniformly as $t\to0$. The negative solution of $K'(t_*)=-b$ therefore satisfies $t_*=-\frac{b}{2(D+2\lambda)}(1+o(1))$. Chernoff's bound gives 

$$
\log p_\theta\le K(t_*)+t_*b
=-\frac{b^2}{4D}(1+o(1)).
$$

For completeness the same exponent is a lower bound. Under exponential tilting with density $\exp(t_*[Z-D-\lambda]-K(t_*))$, the random variable $Z$ is a scaled noncentral chi-square with $D$ degrees of freedom, noncentrality $\lambda/(1-2t_*)$, and variance asymptotic to $2D$. Its standardized distribution converges to a standard normal: rotate the mean equally over the $D$ coordinates, express $Z$ as a sum of independent squared normals, and apply the Lyapunov central limit theorem, since their noncentrality per coordinate tends to zero. Consequently the tilted probability of $-b-\sqrt D\le Z-D-\lambda\le-b$ is bounded below by a positive constant. Changing measure on this interval yields 

$$
\log p_\theta\ge K(t_*)+t_*b-|t_*|\sqrt D+O(1).
$$

Here $|t_*|\sqrt D=O(\sqrt{\log m})=o(\log m)$, so

<a id="eq-2-17"></a>

$$
\begin{aligned}
\log p_\theta
&=-\frac{D(\Delta^2-4\theta^2)^2}{16\sigma^4}(1+o(1))\\
&=-\frac{\log m}{16c^4}(1+o(1)).
\end{aligned}
$$

If $c<2^{-5/4}$, then $m^2p_\theta\to0$; if $c>1/2$, then $mp_\theta\to\infty$. Substitution in [2.13](#eq-2-13) completes the proof without assuming independent cross pairs.
