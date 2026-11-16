# 编码曲面的连续变分与有限实现

在固定参考任务上改变编码几何和质量。下文保留光滑曲面上的物质导数、核失真界、有限码本梯度与接收似然推导。

[研究索引](../README.md) · [符号和引用所用书目](../../ref.bib)

## A Fixed Semantic Source and a Learned Surface

### Knowledge, task specification, and local dimension

Let $C$ denote declared knowledge available to both ends. For each fixed $C=c$, the following analysis applies to $(\mathcal M_c,g_c,\mu_c)$; notation suppresses $c$. Across independent source/knowledge units, conditional source coding minimizes $I(U;\widehat U\mid C)$ with an average distortion budget. Knowledge stored only at the receiver requires a different side-information coding theorem [`guo_side`]; this is not inferred from the shared-knowledge case. A local chart $\varphi_c:\mathbb R^m\to\mathbb R^d$ has $\operatorname{rank}J_{\varphi_c}=m\ll d$. This is local tangent rank, not global covariance rank: a curved one-dimensional circle already spans two ambient directions. Data experiments must measure tangent residuals and dimension stability rather than infer a manifold from a two-dimensional visualization. The task family and reconstruction interface are fixed in advance. When only one known decision is required, directly coding that decision or its posterior is an essential baseline; a full-state representation is not automatically necessary.

### Task calibration

Let $U=E_0(X)$ be a fixed reference and $T$ a finite task target with alphabet size $K_T$. A frozen calibrated posterior $q_t(u)$ is positive. The theoretical task sufficiency condition is $I(T;X\mid U)=0$; residual task loss must be measured when the learned reference is not sufficient. On a compact connected boundaryless $m$-manifold $\mathcal M$, fix a smooth immersion $\psi$ and $\eta>0$. Define

<a id="eq-4-1"></a>

$$
\begin{aligned}
\Phi(u)&=(2\sqrt{q_1(u)},\ldots,2\sqrt{q_{K_T}(u)},\sqrt\eta\psi(u)),\\
g&=J_q^\top\operatorname{diag}(q)^{-1}J_q+\eta J_\psi^\top J_\psi.
\end{aligned}
$$

The second term prevents degeneracy in directions not distinguished by the task posterior. We assume a $C^8$ positive metric and a positive $C^8$ density bounded above and below, $d\mu=p_gdV_g$. All logarithms are natural; division by $\log2$ converts nats to bits. The squared semantic distortion is the fixed $d_g^2(u,v)$, not a distance that the codec is free to shrink.

The Hellinger chord is bounded by the length induced by $g$, and Cauchy--Schwarz gives

$$
\operatorname{TV}(q(u),q(v))\le\left\lVert\sqrt{q(u)}-\sqrt{q(v)}\right\rVert\le\tfrac12d_g(u,v).
$$

For decision losses in $[0,1]$, the excess risk of using the Bayes action at $v$ is at most $2\operatorname{TV}$, hence at most $d_g(u,v)$. Therefore expected squared distortion at most $D$ implies mean regret at most $\sqrt D$. If $q_t\ge q_{\min}>0$ and $T-U-\widehat U$, then

<a id="eq-4-2"></a>

$$
I(T;U)-I(T;\widehat U)\le\mathbb ED_{\rm KL}(q(U)\Vert q(\widehat U))\le D/q_{\min}.
$$

For the last step use $D_{\rm KL}(p\Vert q)\le\sum_t(p_t-q_t)^2/q_t$ and $\left\lVert q(u)-q(v)\right\rVert_2\le2\left\lVert\sqrt{q(u)}-\sqrt{q(v)}\right\rVert\le d_g(u,v)$. For the first, the true conditional posterior given $\widehat U$ minimizes conditional cross-entropy. These are conditional guarantees on the declared posterior, not automatic guarantees for any pretrained feature extractor.

### Encoding geometry and anti-collapse constraints

A sufficiently smooth global embedding $F:\mathcal M\to\mathbb R^d$ induces $h=J_F^\top J_F$. Fix constants $0<a\le b$ and constrain

<a id="eq-4-3"></a>

$$
ag\preceq h\preceq bg.
$$

Then $\sqrt a\,d_g\le d_h\le\sqrt b\,d_g$ for intrinsic distances. Local rank or this local metric condition alone does not prove global injectivity; the embedding requirement is separate. For empirical charts, both conditions require validation rather than being inferred from a sampled penalty.

Let $A_F=\int dV_h$ (area only if $m=2$), $J_F^{\rm vol}=\sqrt{\det h/\det g}$, and $p_h=p_g/J_F^{\rm vol}$. The source mass is fixed under all surface variations. Writing $h_g(\mu)=-\int\log p_g\,d\mu$, we obtain

<a id="eq-4-4"></a> <a id="eq-4-5"></a>

$$
\begin{aligned}
h_h(\mu)-h_g(\mu)&=\mathbb E\log J_F^{\rm vol},\\
h_h(\mu)&=\log A_F-D_{\rm KL}(\mu\Vert dV_h/A_F).
\end{aligned}
$$

Thus volume contributes together with density and semantic scale.

**Theorem (Representation invariance).**

<a id="eq-4-6"></a>

For an invertible measurable embedding and the transported distortion $\widetilde d(F(u),F(v))=d_g(u,v)$, the unrestricted rate--distortion function is unchanged. Surface area alone cannot be its operational rate.

**Proof.** Push every test channel $P_{\widehat U|U}$ through $F$ at both endpoints. Mutual information and distortion are unchanged by the measurable bijection. Apply the inverse map for the reverse inequality between the infima. In particular, $F\mapsto cF$ changes $A_F$ to $c^mA_F$ without changing this coding problem. For a two-dimensional surface, its usual Willmore energy is scale invariant, so adding that energy does not repair area-based rate accounting.

## Rate Bounds and Their Channel Interpretation

Define

$$
R_g(D)=\inf_{\mathbb Ed_g^2(U,\widehat U)\le D}I(U;\widehat U).
$$

For $s\ge0$, let

<a id="eq-4-7"></a>

$$
\begin{aligned}
Z_s^h(u)&=\int_\mathcal Me^{-sd_h^2(u,v)}dV_h(v),\\
Q_s^h(dv|u)&=(Z_s^h(u))^{-1}e^{-sd_h^2(u,v)}dV_h(v).
\end{aligned}
$$

Let $p_{h,s}$ be its output density and set

<a id="eq-4-8"></a>

$$
\begin{aligned}
D_h(s)&=-\partial_s\mathbb E\log Z_s^h(U),\\
I_h(s)&=h_h(p_{h,s})-sD_h(s)-\mathbb E\log Z_s^h(U),\\
D_{g|h}(s)&=\iint d_g^2(u,v)Q_s^h(dv|u)d\mu(u).
\end{aligned}
$$

**Theorem (Computable rate sandwich).**

<a id="eq-4-9"></a>

With

$$
L_g(D)=\max\{0,\sup_{s\ge0}[h_g(\mu)-sD-\log\sup_vZ_s^g(v)]\},
$$

we have $L_g(D)\le R_g(D)\le I_h(s)$ whenever $D_{g|h}(s)\le D$. The more conservative condition $D_h(s)\le aD$ is sufficient. Over a known memoryless complex Gaussian channel $Y=HS+N$, $N\sim\mathcal{CN}(0,N_0I)$, with power $P$, let

$$
C(H)=\max_{Q\succeq0,\operatorname{tr}Q\le P}\log\det(I+HQH^*/N_0).
$$

Every asymptotically achievable state throughput $v=n/k$ satisfies $vR_g(D)\le C(H)$. Every strict $vI_h(s)<C(H)$ is achievable with arbitrary positive distortion slack. For positive rates, closure gives $C(H)/I_h(s)\le v^*(D)=C(H)/R_g(D)\le C(H)/L_g(D)$.

**Proof.** Compare the backward conditional law $P_{U|\widehat U=v}$ with $e^{-sd_g^2(u,v)}dV_g(u)/Z_s^g(v)$. Nonnegative KL gives $h_g(U|\widehat U=v)\le s\mathbb E[d_g^2(U,v)|v]+\log Z_s^g(v)$. Average, subtract from $h_g(\mu)$, and optimize in $s$. This is the reference-measure Shannon lower-bound construction [`riegler`].

For the forward channel, its conditional entropy is $sD_h+\mathbb E\log Z_s^h$. Subtracting it from the output entropy proves [4.8](#eq-4-8). The distortion inequality follows pointwise from [4.3](#eq-4-3). Compactness and bounded continuous distortion permit finite-partition approximation and the ordinary random-covering source-coding theorem, yielding achievability with positive slack. A finite quadrature calculation is exact for its discretized source; a continuous-source certificate additionally needs controlled quadrature error.

For i.i.d. states, the entropy chain rule, conditioning, and convexity imply $I(U^n;\widehat U^n)\ge\sum_i I(U_i;\widehat U_i)\ge nR_g(\bar D)$. Data processing and the power constraint give $I(U^n;\widehat U^n)\le I(S^k;Y^k)\le kC(H)$. Conversely concatenate source and channel codes, with rates strictly inside both limits. Bounded distortion makes the contribution of vanishing channel error vanish. This proves the throughput statement; a fading outage channel or a finite block needs its own reliability analysis.

## Dimension, Volume, Density, and Curvature

Define $S_h=\operatorname{Scal}_h$, $\bar S_h=\mathbb ES_h(U)$,

$$
\begin{aligned}
\mathcal J_h&=\int\left\lVert\nabla_h\log p_h\right\rVert_h^2d\mu,\\
\Gamma_F&=\tfrac12\mathbb E\log\frac{\det h}{a^m\det g},\\
B_g(D)=h_g(\mu)-\tfrac m2\log(2\pi eD/m).
\end{aligned}
$$

By [4.3](#eq-4-3), $0\le\Gamma_F\le(m/2)\log(b/a)$.

**Theorem (Expansion of the specified achievable kernel).**

<a id="eq-4-10"></a>

For sufficiently small $D>0$, choose $s_D$ from $D_h(s_D)=aD$. Under the smoothness assumptions,

<a id="eq-4-11"></a> <a id="eq-4-12"></a>

$$
\begin{aligned}
I_h(s_D)&=B_g(D)+\Gamma_F+\frac{aD}{6m}(\bar S_h+3\mathcal J_h)+O(D^2),\\
R_g(D)&\ge B_g(D)+\frac{D}{6m}\min_\mathcal MS_g-O(D^2).
\end{aligned}
$$

For $h=g$, these imply $R_g(D)=B_g(D)+O(D)$. On a homogeneous manifold with uniform source, the exact parametric RDF is $R_g(D_g(s))=\log A_g-\log Z_s^g-sD_g(s)$, and its linear curvature coefficient is $S_g/(6m)$.

**Proof.** Use normal coordinates below the positive injectivity radius. With $c_s=(\pi/s)^{m/2}$ and volume density $1-\operatorname{Ric}_{ij}x^ix^j/6+O(|x|^3)$, uniform Taylor expansion and Gaussian moments yield

<a id="eq-4-13"></a> <a id="eq-4-14"></a>

$$
\begin{aligned}
Z_s^h(u)&=c_s[1-S_h(u)/(12s)+O(s^{-2})],\\
\mathcal A_sf(v)
&=f(v)+\Delta_hf(v)/(4s)\\
&\quad-S_h(v)f(v)/(12s)+O(s^{-2}).
\end{aligned}
$$

Here $\mathcal A_sf(v)=c_s^{-1}\int e^{-sd_h^2(u,v)}f(u)dV_h(u)$. Odd moments vanish. Contributions outside the normal ball are exponentially small; compactness bounds the remainders uniformly. Apply [4.14](#eq-4-14) to $p_h(u)/Z_s^h(u)$ after substituting [4.13](#eq-4-13). The two curvature terms cancel, giving $p_{h,s}=p_h+\Delta_hp_h/(4s)+O(s^{-2})$. Positivity, probability conservation, and integration by parts imply

$$
h_h(p_{h,s})=h_h(\mu)+\mathcal J_h/(4s)+O(s^{-2}).
$$

Expand the distance moment separately, rather than differentiating an unspecified remainder:

$$
D_h(s)=m/(2s)-\bar S_h/(12s^2)+O(s^{-3}).
$$

Its derivative is $-\mathbb E\operatorname{Var}_{Q_s^h}(d_h^2|U)<0$. Hence $s_D=m/(2aD)-\bar S_h/(6m)+O(D)$. Substitution into [4.8](#eq-4-8) cancels the constant shift between $(m/2)\log s_D$ and $-s_DaD$, leaving $h_h(\mu)-(m/2)\log(2\pi eaD/m)+aD(\bar S_h+3\mathcal J_h)/(6m)+O(D^2)$. Equation [4.4](#eq-4-4) proves [4.11](#eq-4-11).

Also $\log\sup_u Z_s^g(u)=(m/2)\log(\pi/s)-\min S_g/(12s)+O(s^{-2})$. Evaluate the converse at $s=m/(2D)$ to obtain [4.12](#eq-4-12). Homogeneity makes $Z_s^g$ independent of its center, and the symmetric kernel preserves uniform volume. Its information equals the converse at $D_g(s)$, proving the exact specialization.

The term $\Gamma_F$ is the leading excess of this construction relative to $B_g$, using the prescribed sufficient distortion conversion. It depends on $a$ and is not a universal penalty for all codes on $F(\mathcal M)$. The coefficient $\mathcal J_h/(2m)$ belongs to this forward kernel; it is not an asserted exact RDF coefficient for a nonuniform source. Curvature here is intrinsic scalar curvature. For $m=2$, $S=2K$; extrinsic mean curvature is a different object. In Euclidean ambient space $S_h=\left\lVert H_{\rm tr}\right\rVert^2-\left\lVert\mathrm{II}\right\rVert_h^2$. A large curvature need not identify task importance or immunity to noise.

The corresponding small-ball calculation is

$$
\frac{\mu(B_h(u,r))}{v_mr^m}=p_h(u)+\frac{r^2}{m+2}\left(\tfrac12\Delta_hp_h-\tfrac16S_hp_h\right)+O(r^4).
$$

It explains the roles of occurrence probability and finite-resolution geometry without treating a ball covering count as an entropy-coded bit rate.

## Variational Formulation and Realizable Codec

### Material variations and exact discrete gradients

For $F_t=F+tV$ and fixed $a,\mu$, $\dot h_{ab}=\langle\partial_aV,\partial_bF\rangle+\langle\partial_aF,\partial_bV\rangle$. Thus

<a id="eq-4-15"></a>

$$
\begin{aligned}
\delta\Gamma_F[V]&=\tfrac12\int\operatorname{tr}(h^{-1}\dot h)d\mu
=\int h^{ab}\langle\partial_aF,\partial_bV\rangle d\mu.
\end{aligned}
$$

Integration by parts shows its $L^2(\mu)$ gradient is $-(p_g\sqrt{\det g})^{-1}\partial_b(p_g\sqrt{\det g}\,h^{ab}\partial_aF)$. For comparison, the unweighted area variation under normal velocity is $\delta A_F=-\int\langle H_{\rm tr},V^\perp\rangle dV_h$. These are not the same functional.

For a finite approximation let fixed masses be $w_i$, positive quadrature weights $v_j(F)$, and $c^h_{ij}=d_h^2(u_i,u_j)$. Define

$$
Q_{ij}=\frac{v_je^{-sc^h_{ij}}}{\sum_l v_le^{-sc^h_{il}}},\qquad \pi_j=\sum_iw_iQ_{ij}.
$$

$$
I_F=\sum_{ij}w_iQ_{ij}\log(Q_{ij}/\pi_j).
$$

With $A_{ij}=\delta\log v_j-s\delta c^h_{ij}-c^h_{ij}\delta s$,

<a id="eq-4-16"></a>

$$
\begin{aligned}
\delta Q_{ij}&=Q_{ij}(A_{ij}-\mathbb E_{Q_i}A_i),\\
\delta I_F&=\sum_iw_i\operatorname{Cov}_{Q_i}(\log(Q_{ij}/\pi_j),A_{ij}).
\end{aligned}
$$

The marginal derivatives cancel because rows remain normalized. With the fixed semantic cost $c^g_{ij}$, the variation of $I_F+\lambda\sum w_iQ_{ij}c^g_{ij}$ replaces the first covariance argument by $\log(Q_{ij}/\pi_j)+\lambda c^g_{ij}$. There is no $\delta c^g$ term. If $s$ is optimized implicitly to satisfy a distortion equality, its variation must be retained; the simpler fixed-$s$ gradient sets $\delta s=0$.

An implementable fourth-order smoother is

$$
\mathcal B_g(F)=\tfrac12\int\left\lVert\Delta_gF\right\rVert^2d\mu,\quad
\delta\mathcal B_g[V]=\int\langle\Delta_gF,\Delta_gV\rangle d\mu.
$$

Its $L^2(\mu)$ gradient is $p_g^{-1}\Delta_g(p_g\Delta_gF)$. This fixed-reference Sobolev penalty is not identical to Willmore energy. With a differentiable coding functional $\mathcal J_0$ containing rate, distortion, and metric constraints, stationarity is exactly

$$
\delta\mathcal J_0[V]+\beta\int\langle\Delta_gF,\Delta_gV\rangle d\mu=0
$$

for every admissible variation. When its Riesz derivative exists this becomes $\operatorname{grad}_{L^2(\mu)}\mathcal J_0+\beta p_g^{-1}\Delta_g(p_g\Delta_gF)=0$. Power and spectral constraints add their normal-cone/KKT terms. It is a stationarity condition, not a uniqueness or global-optimality claim. Boundary charts require specified boundary conditions instead of discarding boundary terms.

### Rate accounting and training

The digital codec shares its reference, atlas, network weights, and entropy model. It transmits state indices $K=Q_\theta(F_\theta(U))$, not an uncharged new manifold for every sample. The operational surrogate is

$$
\mathcal R_{\rm ent}=\mathbb E[-\log\pi_\omega(K)]=H(K)+D_{\rm KL}(P_K\Vert\pi_\omega).
$$

Actual arithmetic-coded lengths include framing, chart labels, and other side information. If a sample-specific surface is transmitted, its parameters must also be coded. A deterministic continuous encoder does not in general have a finite mutual information with its unquantized output.

Train at a fixed reference and power convention using

$$
\begin{aligned}
\mathcal L={}&\mathcal R_{\rm ent}+\lambda_D\mathbb Ed_g^2(U,\widehat U)
+\lambda_T\mathbb E\ell(T,\widehat T)\\
&+\lambda_\Gamma\Gamma_F+\lambda_B\mathcal B_g(F)+\lambda_M\mathcal L_{\rm met}.
\end{aligned}
$$

A sampled eigenvalue hinge implements $\mathcal L_{\rm met}$ for the eigenvalues of $g^{-1/2}hg^{-1/2}$ relative to $[a,b]$. It is not a global certificate. True extrinsic bending can alternatively be differentiated using $\mathrm{II}_{ab}=(I-J_Fh^{-1}J_F^\top)\partial_{ab}F$; it should be named and measured separately from $\mathcal B_g$. Analog fixed-symbol JSCC instead fixes the channel-use budget and omits the entropy-rate term; the two interfaces are not mixed in comparisons.

### Likelihood-matched reception and semantic sensitivity

For a continuous analog interface with known $A_H$ and positive noise covariance $\Sigma$, use the negative log posterior relative to the fixed $dV_g$ measure,

<a id="eq-4-17"></a>

$$
\mathcal E_y(u)=\tfrac12\left\lVert y-A_HF(u)\right\rVert_{\Sigma^{-1}}^2-\log p_g(u).
$$

A minimizer exists by compactness and positivity, but may not be unique. This is a MAP estimator under the stated reference measure; it is not automatically the Bayes estimator for squared semantic loss. Its coordinate gradient is

$$
\nabla_g\mathcal E_y=g^{-1}J_F^\top A_H^\top\Sigma^{-1}(A_HF(u)-y)-\nabla_g\log p_g.
$$

Iterate $u^+=\operatorname{Exp}_u(-\gamma\nabla_g\mathcal E_y)$, or a declared retraction. Under an $L$-smooth pullback on the step neighborhood,

$$
\mathcal E_y(u^+)\le\mathcal E_y(u)-\gamma(1-L\gamma/2)\left\lVert\nabla_g\mathcal E_y\right\rVert^2,
$$

by the second-order descent inequality. For fixed $0<\gamma<2/L$ and a lower-bounded energy, summing proves square-summability of gradients and stationarity of accumulation points. An arbitrary learned unrolling does not inherit this guarantee. CSI and pilot requirements are unchanged by the geometric prior.

In a local quadratic approximation let equalized mode variances be $c_j/p_j$ and task sensitivities be $w_j>0$. Then $D_{\rm loc}\simeq\sum_jw_jc_j/p_j$. Cauchy--Schwarz gives

$$
p_j^*=P\frac{\sqrt{w_jc_j}}{\sum_l\sqrt{w_lc_l}},\quad
D_{\rm loc}^{\min}=\frac{(\sum_j\sqrt{w_jc_j})^2}{P}.
$$

Here $w_j$ comes from the fixed task metric composed with the receiver Jacobian, not from curvature alone. This solves the declared local surrogate, not arbitrary nonlinear fading-channel resource allocation.

## Geometry and Coding-Measure Co-Design

The volume-reference Gibbs kernel is valuable for a geometric expansion, but it need not be an efficient finite-distortion code for a nonuniform source. Introduce a positive coding density $r_\omega$ relative to $dV_h$. In a fixed finite quadrature,

<a id="eq-4-18"></a>

$$
Q_{ij}=\frac{v_j r_{\omega,j}\exp(-s c^h_{ij})}{\sum_l v_l r_{\omega,l}\exp(-s c^h_{il})},\quad
D_g=\sum_{ij}w_iQ_{ij}c^g_{ij}.
$$

The source weights $w_i$ and semantic costs $c^g_{ij}$ remain fixed. Optimize $J=I(Q)+\lambda D_g$ subject to the declared geometry constraints. For a variation define $A_{ij}=\delta\log(v_jr_{\omega,j})-s\delta c^h_{ij}-c^h_{ij}\delta s$. Row normalization and the marginal cancellation used above yield

<a id="eq-4-19"></a>

$$
\delta J=\sum_i w_i\operatorname{Cov}_{Q_i}
\left(\log\frac{Q_{ij}}{\pi_j}+\lambda c^g_{ij},A_{ij}\right).
$$

This is an implementable finite-dimensional variational derivative. It is not a global convergence result. The earlier curvature--density expansion applies to its stated volume-reference kernel $r_\omega=1$, not automatically to [4.18](#eq-4-18).

If $h=g$, $s=\lambda$, and the reproduction weights are unrestricted, the self-consistency equation for $v_jr_{\omega,j}=\pi_j$ recovers the classical Blahut--Arimoto (BA) fixed point. We do not claim that rediscovering this fixed point is novel. Empirical RDF bounds and reproduction-measure optimization already exist [`yang_sandwich`, `yang_wgd`]. The research objective is a controlled low-dimensional parameterization and its approximation gap under fixed task geometry. A geometry-only advantage over a density-only control is an empirical question.

The mismatch issue also requires a conditional, task-relevant statement. Let $p(T|y,c)$ be the true posterior, $q(T|y,c)$ the receiver model, and losses lie in $[0,1]$. Comparing the two Bayes risks by adding and subtracting the model risk gives

$$
\begin{aligned}
\mathcal R_p(a_q)-\mathcal R_p(a_p)
&\le 2\mathbb E_{Y,C}\operatorname{TV}(p(\cdot|Y,C),q(\cdot|Y,C))\\
&\le\sqrt{2\mathbb E_{Y,C}D_{\rm KL}(p(\cdot|Y,C)\Vert q(\cdot|Y,C))}.
\end{aligned}
$$

The last step is Pinsker followed by Jensen. Matching unconditional feature marginals is neither necessary nor sufficient for this task guarantee. This observation constrains our interpretation; it does not supply an unvalidated domain-adaptation method.
