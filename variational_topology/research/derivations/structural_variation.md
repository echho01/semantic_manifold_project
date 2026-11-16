# 几何与关系双失真、结构变分与功率分配

结构类型是可解的特例。主稿使用有坐标身份的连通关系。这里推导联合约束、滤过稳定性、充分功率和可容错的乘积信源。

[研究索引](../README.md) · [符号和引用所用书目](../../ref.bib)

## Structural Types, Continuous States, and Knowledge

Unless explicitly stated otherwise, logarithms and information rates use nats; divide by $\log 2$ for bits. Let the semantic source be $S=(B,U)$ with $B\in\{1,\ldots,K\}$ and $U\in\mathcal M_B$. Thus the source support is a disjoint union $\bigsqcup_b\{b\}\times\mathcal M_b$, not necessarily one connected smooth manifold. Within each stratum fix a task metric $g_b$ and a conditional probability law. Let $C$ be side information available to both source encoder and decoder. We assume i.i.d. source/knowledge tuples across coding units. Decoder-only knowledge requires a Wyner--Ziv formulation and is not covered by the following identity.

Define type error $d_T(S,\widehat S)=\mathbf1\{B\ne\widehat B\}$. For $\widehat B=B$, define $d_G=d_{g_B}^2(U,\widehat U)$; for a wrong type assign a fixed bounded penalty $D_{\rm bad}$ larger than all within-stratum squared diameters. The information rate is

<a id="eq-5-1"></a>

$$
R(D,\epsilon\mid C)=\inf_{\mathbb Ed_G\le D,\,\Pr(B\ne\widehat B)\le\epsilon}I(S;\widehat S\mid C).
$$

The source type can be a task-defined connection class. It is not assumed to equal a Betti number or an unlabelled persistence diagram.

**Theorem (Structural rate and shared knowledge).**

<a id="eq-5-2"></a>

For exact type preservation at the test-channel level,

<a id="eq-5-3"></a>

$$
\begin{aligned}
R(D,0\mid C)=H(B\mid C)+\inf_{\{D_{bc}\}}\sum_{b,c}p(b,c)R_{bc}(D_{bc}),\\
\text{subject to}\quad \sum_{b,c}p(b,c)D_{bc}\le D,
\end{aligned}
$$

where $R_{bc}$ is the within-stratum conditional RDF. Sums are replaced by expectations for a general $C$. For $0\le\epsilon\le1-1/K$,

<a id="eq-5-4"></a>

$$
R(D,\epsilon\mid C)\ge\max\{0,H(B\mid C)-h_2(\epsilon)-\epsilon\log(K-1)\}.
$$

If a structure entropy code uses a positive model $q(b|c)$ instead of the true law, its expected ideal length is

<a id="eq-5-5"></a>

$$
\mathbb E[-\log q(B|C)]=H(B|C)+\mathbb E_CD_{\rm KL}(p(B|C)\Vert q(B|C)).
$$

**Proof.** If $\widehat B=B$ almost surely, the chain rule yields

$$
\begin{aligned}
I(B,U;\widehat B,\widehat U\mid C)
&=I(B;\widehat B,\widehat U\mid C)
+I(U;\widehat B,\widehat U\mid B,C)\\
&=H(B\mid C)+I(U;\widehat U\mid B,C).
\end{aligned}
$$

Each conditional mutual information is at least $R_{bc}(D_{bc})$. Conversely, select conditional test channels approaching the displayed infimum and reproduce $B$ exactly. The resulting test channel attains the same expression up to arbitrary slack. The ordinary conditional source-coding theorem supplies operational block codes: encode the structural sequence given shared $C^n$, then its continuous states conditionally, with vanishing error and positive rate/distortion slack. This is not a zero-error guarantee on a noisy finite block.

For nonzero type error, data processing gives $I(S;\widehat S|C)\ge I(B;\widehat B|C)$. Set $E=\mathbf1\{B\ne\widehat B\}$. The entropy chain rule gives $H(B|\widehat B,C)\le H(E)+\Pr(E=1)\log(K-1)$. Thus $I(B;\widehat B|C)\ge H(B|C)-h_2(P_e)-P_e\log(K-1)$. This penalty is increasing up to $1-1/K$, establishing [5.4](#eq-5-4). Finally add and subtract $-\log p(B|C)$ to prove [5.5](#eq-5-5).

Equation [5.3](#eq-5-3) identifies the structural and within-structure costs only under exact type preservation and shared side information. The structural-bit saving from knowledge is $I(B;C)$; it is not necessarily the total saving of the continuous-plus-discrete codec. Equation [5.5](#eq-5-5) is an exact penalty for that mismatched entropy code, not a universal change of the optimal RDF. Unknown task types and a misspecified support are separate modeling errors. For a known memoryless channel, source--channel separation gives the asymptotic benchmark $vR(D,\epsilon|C)\le C_{\rm ch}$. The practical certificate below concerns a specified representation and is not equated to this RDF.

### Connection to the geometric small-distortion law

Assume finitely many positive-probability pairs $(b,c)$ and the compact, positive-density smooth-source assumptions within each stratum. Let its dimension be $m_b>0$, its volume-relative entropy be $h_{bc}$, and $\bar m=\sum_{b,c}p(b,c)m_b$. The within-stratum geometric law gives

$$
R_{bc}(D_{bc})=h_{bc}-\frac{m_b}{2}\log\frac{2\pi eD_{bc}}{m_b}+O(D_{bc}).
$$

Consequently,

<a id="eq-5-6"></a>

$$
\begin{aligned}R(D,0|C)={}&H(B|C)+\sum_{b,c}p(b,c)h_{bc}\\&-\frac{\bar m}{2}\log\frac{2\pi eD}{\bar m}+O(D).\end{aligned}
$$

To verify this, minimize the leading convex logarithmic expression with a multiplier for $\sum pD_{bc}=D$. Its stationarity condition is $-m_b/(2D_{bc})+\lambda=0$, hence $D_{bc}=m_bD/\bar m$. Every feasible $D_{bc}\le D/p(b,c)$ is small, and the finite collection gives a uniform total $O(D)$ remainder. Sandwiching the true objective between its leading expression plus/minus that uniform remainder proves [5.6](#eq-5-6). The displayed allocation attains the leading asymptotic value; an exact finite-distortion optimizer need not coincide with it.

Thus structure contributes conditional type entropy, whereas continuous precision contributes a probability-weighted intrinsic dimension and entropy. If each conditional source is uniform on a homogeneous manifold, its exact RDF has linear coefficient $S_{bc}/(6m_b)$. Under the corresponding smooth expansion the mixture coefficient is $\mathbb ES_{BC}/(6\bar m)$: inserting the leading allocation yields this term, while the first-order change of the logarithmic objective vanishes under the distortion constraint. This specialization does not assert that a general nonuniform source has that exact coefficient.

## A Task-Registered Filtration

### Episode and fixed reference

A transmitted episode is $W=(U_1,\ldots,U_N)$, with declared grouping and correspondence. States inside it may be dependent. Let $Z=(z_1,\ldots,z_N)$, $z_i=F_\theta(U_i)\in\mathbb R^d$, use a fixed scale, and constrain reference metric distortion. A batch of unrelated images is not silently reinterpreted as the topology of one scene.

For the point-cloud implementation, use the edge-length convention

<a id="eq-5-7"></a>

$$
a_Z(\sigma)=\max_{i,j\in\sigma}\left\lVert z_i-z_j\right\rVert,\qquad K_r(Z)=\{\sigma:a_Z(\sigma)\le r\}.
$$

The full two-skeleton suffices for $H_0,H_1$ over $\mathbb F_2$. Compare finite bars and pair the common essential $H_0$ interval with itself. Let $\mathcal D(Z)$ collect these diagrams and $d_{\rm PH}$ be the maximum of their bottleneck distances, using the sup norm in birth--death coordinates.

A frozen reference episode supplies $\mathcal D^0(W)$. Register task-selected source bars with minimum lifetime $\ell_0>0$, not arbitrary topological noise. For fixed diagram prototypes $\{\mathcal P_c\}$ define

<a id="eq-5-8"></a>

$$
\begin{aligned}
b_0&=\operatorname*{arg\,min}_c d_{\rm PH}(\mathcal D^0,\mathcal P_c),\\
\Delta_0&=\min_{c\ne b_0}d_{\rm PH}(\mathcal D^0,\mathcal P_c)
-d_{\rm PH}(\mathcal D^0,\mathcal P_{b_0}),\\
\tau_0&=\tfrac12\min\{\ell_0,\Delta_0\},\qquad
\epsilon_{\rm enc}=d_{\rm PH}(\mathcal D^0,\mathcal D(Z)).
\end{aligned}
$$

Assume a positive prototype margin when certifying that decision. A structure-only certificate uses $\tau_0=\ell_0/2$. The reference is fixed: redefining it after the encoder has collapsed cannot improve its certificate.

### What this topology does not identify

Persistent homology is not a complete invariant of labeled task relations. Two labeled graphs with identical node coordinates but different edge labels have identical Euclidean Rips filtrations; Betti numbers also fail to identify which endpoints are connected. A relation-sensitive experiment must therefore use a task-justified graph/image filtration or retain a labeled relation head. Our prototype theorem certifies only the specified diagram-based decision.

For a separate edge-score representation $w_e$, if $\max_e|\widehat w_e-w_e|\le\xi$ and $\min_e|w_e-r_0|>\xi$, then every thresholded edge at $r_0$ is identical. This follows immediately because no edge score can cross the threshold; hence the labeled graph and all its connectivity queries are unchanged. Such a bound must be proved for the actual graph head, for example through a verified Lipschitz constant. It cannot be imported from the Rips theorem without that link.

### Variation versus changing topological type

For a smooth scalar potential $f$ on a manifold, sublevel sets provide another filtration. If $\left\lVert f-\widehat f\right\rVert_\infty\le\xi$, the two filtrations are $\xi$-interleaved. For a Morse function, a compact band without critical points is transported by $\nabla f/\left\lVert\nabla f\right\rVert^2$, since $df/dt=1$; regular levels in the band are equivalent. Topological changes require crossing critical values. On a closed oriented surface, Gauss--Bonnet fixes $\int K\,dA=2\pi\chi$ at fixed topology. Its derivative there is zero, so it is not a useful smooth genus-changing loss. We optimize filtration values and code selection rather than asserting a smooth genus-changing flow.

## A Reference-Anchored Wireless Certificate

In real-equivalent coordinates, let

$$
\begin{aligned}
y_i&=A_H(z_i+q_i)+b_H+r_i+n_i,\\
\widehat z_i&=L_Hy_i,\quad \left\lVert q_i\right\rVert\le q,\quad \left\lVert L_Hr_i\right\rVert\le c.
\end{aligned}
$$

Conditional on the channel, $n_i$ are independent zero-mean Gaussians and $e_i=L_Hn_i\sim\mathcal N(0,\Sigma_H)$. The channel/equalizer is common to this episode, and $b_H$ is a common offset. Noncommon CSI and synchronization errors require the declared residual bound. Define

<a id="eq-5-9"></a>

$$
\begin{aligned}
\delta_H(Z)&=\max_{i<j}\left|\left\lVert L_HA_H(z_i-z_j)\right\rVert-\left\lVert z_i-z_j\right\rVert\right|,\\
\eta_0&=\delta_H(Z)+2\left\lVert L_HA_H\right\rVert_{\rm op}q+2c,\\
t_\alpha&=\log\frac{N(N-1)}{2\alpha},\quad 0<\alpha<1,\quad N\ge2,\\
\nu_\alpha&=\sqrt{2\operatorname{tr}\Sigma_H}+2\sqrt{\left\lVert\Sigma_H\right\rVert_{\rm op}t_\alpha},\\
\eta_\alpha&=\eta_0+\nu_\alpha.
\end{aligned}
$$

**Theorem (Reference structure and decision preservation).**

<a id="eq-5-10"></a>

With conditional probability at least $1-\alpha$,

<a id="eq-5-11"></a>

$$
\begin{aligned}
d_{\rm PH}(\mathcal D(Z),\mathcal D(\widehat Z))&\le\eta_\alpha,\\
d_{\rm PH}(\mathcal D^0,\mathcal D(\widehat Z))&\le\epsilon_{\rm enc}+\eta_\alpha.
\end{aligned}
$$

On this event, every registered reference bar of lifetime greater than $2(\epsilon_{\rm enc}+\eta_\alpha)$ has an off-diagonal matching partner. If $\Delta_0>2(\epsilon_{\rm enc}+\eta_\alpha)$, the nearest-prototype decision equals $b_0$. Thus $\epsilon_{\rm enc}+\eta_\alpha<\tau_0$ suffices for both targets.

**Proof.** Common offsets cancel. For every pair, the reverse triangle inequality gives

<a id="eq-5-12"></a>

$$
\left|\left\lVert\widehat z_i-\widehat z_j\right\rVert-\left\lVert z_i-z_j\right\rVert\right|\le\eta_0+\left\lVert e_i-e_j\right\rVert.
$$

For $v\sim\mathcal N(0,C)$, diagonalize $C$ and put $X=\left\lVert v\right\rVert^2-\operatorname{tr}C=\sum_j\lambda_j(g_j^2-1)$. For $0<s<(2\left\lVert C\right\rVert_{\rm op})^{-1}$,

$$
\begin{aligned}
\log\mathbb Ee^{sX}
&=\tfrac12\sum_j[-\log(1-2s\lambda_j)-2s\lambda_j]\\
&\le\frac{s^2\operatorname{tr}C^2}{1-2s\left\lVert C\right\rVert_{\rm op}}.
\end{aligned}
$$

The scalar inequality follows from the logarithm power series. Chernoff's inequality with $s=\sqrt t/(\sqrt{\operatorname{tr}C^2}+2\left\lVert C\right\rVert_{\rm op}\sqrt t)$ yields

$$
\Pr\{\left\lVert v\right\rVert^2>\operatorname{tr}C+2\sqrt{\operatorname{tr}C^2\,t}+2\left\lVert C\right\rVert_{\rm op}t\}\le e^{-t}.
$$

Since $\operatorname{tr}C^2\le\left\lVert C\right\rVert_{\rm op}\operatorname{tr}C$, a looser convenient radius is $\sqrt{\operatorname{tr}C}+\sqrt{2\left\lVert C\right\rVert_{\rm op}t}$. Apply it with $C=2\Sigma_H$ to each $e_i-e_j$ and union bound over $N(N-1)/2$ pairs. Pair events need not be independent. The zero-covariance case is deterministic. This proves that the maximum edge-value discrepancy in [5.12](#eq-5-12) is at most $\eta_\alpha$ with the stated probability.

A simplex value in [5.7](#eq-5-7) is the maximum of its edges, a one-Lipschitz operation in the sup norm. Thus $K_r(Z)\subseteq K_{r+\eta_\alpha}(\widehat Z)$ and conversely. On the common labeled vertex set these maps commute, giving an $\eta_\alpha$-interleaving. Algebraic stability for finite persistence modules gives the bottleneck bound [`stability`, `modules`]. The triangle inequality adds $\epsilon_{\rm enc}$.

The sup-norm distance of $(b,d)$ to the diagonal is $(d-b)/2$, so a bar longer than twice the matching bound cannot disappear into the diagonal. Each prototype distance changes by at most the same matching bound. A competitor-minus-selected gap therefore changes by at most twice that bound, proving decision preservation.

Matching survival does not ensure the identity of a particular spatial cycle without additional feature separation and correspondence assumptions. Homology groups of the noiseless support are another distinct object: if $A_H$ is injective on a compact embedded manifold, its restriction is a homeomorphism onto its image. Conversely, independent full-rank Gaussian perturbations produce a full-dimensional distribution; it is inappropriate to describe each such draw as a smooth deformation of the original support. The finite-sample filtration certificate avoids this ambiguity.

For correlated token noise use $C_{ij}=\operatorname{Cov}(e_i-e_j|H)$ and replace $\nu_\alpha$ by

$$
\max_{i<j}\left[\sqrt{\operatorname{tr}C_{ij}}+\sqrt{2\left\lVert C_{ij}\right\rVert_{\rm op}t_\alpha}\right].
$$

For different deterministic maps per token use their actual maximum noiseless pair-distance deformation. For a random channel and model-validity event $\mathcal V$, a union bound gives

$$
P_{\rm failure}\le\Pr\{\epsilon_{\rm enc}+\eta_\alpha\ge\tau_0\}+\alpha+\Pr(\mathcal V^c).
$$

Add clean task error for a task-accuracy statement, and add control/channel decoding errors when those are outside the model. A sampled covariance estimate without a confidence guarantee does not automatically yield the advertised coverage.

## Minimum Power for a Sufficient Certificate

Consider known nonzero diagonal fading gains, fixed coordinate second moments, and zero-forcing reception. Let $p_j>0$ be actual mode power under the specified normalization. If source-coordinate second moment is $s_j$, include it explicitly in $c_j=\sigma_{n,j}^2s_j/|h_j|^2$; for unit second moments $s_j=1$. Then $[\Sigma_H]_{jj}=c_j/p_j$.

Let $S_p=\sum_jc_j/p_j$, $K_\alpha=1+\sqrt{2t_\alpha}$, and choose strict slack $\zeta>0$. Since $\left\lVert\Sigma_H\right\rVert_{\rm op}\le\operatorname{tr}\Sigma_H=S_p$, a sufficient noise bound is $\nu_\alpha\le K_\alpha\sqrt{2S_p}$. Define

$$
b_{\rm avail}=\tau_0-\epsilon_{\rm enc}-\eta_0-\zeta>0,\qquad T_0=b_{\rm avail}^2/(2K_\alpha^2).
$$

Here $\eta_0,\epsilon_{\rm enc},\tau_0,c_j$ are fixed with respect to the allocation. This holds for the stated fixed quantizer, fixed representation, and perfect-equalization model, not arbitrary joint representation/power changes.

**Theorem (Minimum sufficient power).**

<a id="eq-5-13"></a>

For $c_j>0$, the problem $\min_{p_j>0}\sum_jp_j$ subject to $\sum_jc_j/p_j\le T_0$ has the unique solution

<a id="eq-5-14"></a>

$$
p_j^*=\frac{\sqrt{c_j}\sum_l\sqrt{c_l}}{T_0},\qquad P_{\rm cert}^*=\frac{(\sum_l\sqrt{c_l})^2}{T_0}.
$$

It ensures $\epsilon_{\rm enc}+\eta_\alpha\le\tau_0-\zeta$. With caps $p_j\le\bar p_j$, feasibility is exactly $\sum_jc_j/\bar p_j\le T_0$. A feasible solution is $p_j^*=\min\{\bar p_j,\sqrt{\lambda c_j}\}$ with the reciprocal constraint active, including the all-capped feasibility boundary.

**Proof.** Cauchy--Schwarz gives $(\sum_j\sqrt{c_j})^2\le(\sum_jp_j)(\sum_jc_j/p_j)\le T_0\sum_jp_j$. Equality requires $p_j$ proportional to $\sqrt{c_j}$; saturation gives [5.14](#eq-5-14) and uniqueness. Under caps, the smallest reciprocal sum is attained at all caps, proving feasibility. Convexity and KKT stationarity give $1-\lambda c_j/p_j^2+\beta_j=0$ and complementary slackness, which imply the clipped square-root rule. The boundary solution follows directly or by a limit.

The optimum is exact for this sufficient trace constraint, not the minimum power of all semantic codes. Evaluate the tighter covariance certificate and actual failure rates to quantify conservatism. A nonpositive $b_{\rm avail}$ is an encoding/quantization/margin failure that cannot be repaired by unlimited power in this model. Zero gains and power caps can cause outage; no pilot reduction is inferred.

## Variational Topological Implementation

Use fixed task references and shared prototypes. A digital source interface has entropy cost $\mathcal R_{\rm ent}=\mathbb E[-\log\pi_\omega(K)]$ including type, chart, and required control indices. A fixed-symbol analog interface instead fixes channel uses explicitly. The training objective is

<a id="eq-5-15"></a>

$$
\begin{aligned}
\mathcal L={}&\lambda_G\mathbb Ed_G+\lambda_T\mathbb E\ell(T,\widehat T)
+\lambda_R\mathcal R_{\rm ent}+\lambda_P\mathbb E\sum_jp_j\\
&+\lambda_{\rm PH}\sum_{k=0}^1W_{2,E}^2(\mathcal D_k^0,\mathcal D_k(Z))\\
&+\lambda_C\mathbb E[\epsilon_{\rm enc}+\eta_\alpha-\tau_0+\zeta]_+^2\\
&+\lambda_M\mathcal L_{\rm met}.
\end{aligned}
$$

$W_{2,E}$ uses Euclidean birth--death costs and diagonal augmentation; the certificate uses bottleneck distance with sup-norm costs. These metrics are named separately. If a smooth surrogate replaces $\epsilon_{\rm enc}$ for training, evaluation must recompute the true bottleneck quantity. Fixed task references and scale constraints prevent diagram collapse or arbitrary dilation from manufacturing robustness.

Within a region of fixed persistence and matching pairings,

<a id="eq-5-16"></a>

$$
\begin{aligned}
\nabla_\theta\mathcal L_{\rm PH}=2\sum_j\{(b_j-b_j^0)\nabla_\theta a(\sigma_{b_j})
+(d_j-d_j^0)\nabla_\theta a(\sigma_{d_j})\}.
\end{aligned}
$$

For a uniquely maximal Rips edge $(i,j)$ with nonzero length, $\nabla_\theta a=(z_i-z_j)^\top(\nabla_\theta z_i-\nabla_\theta z_j)/\left\lVert z_i-z_j\right\rVert$. For a diagonal match, its Euclidean squared cost is $(d-b)^2/2$ and its derivatives with respect to $(b,d)$ are $(b-d,d-b)$. At ties and pairing changes the objective is nonsmooth; use a specified subgradient/tie convention, not an assertion of everywhere differentiability. The certificate hinge differentiates its encoder and channel terms while the reference margin remains fixed.

An implementation iteration is: obtain fixed-reference episode targets; encode and quantize; sample a calibrated channel; compute task and rate losses; compute the episode filtration and encoder discrepancy; calculate the physical certificate; optimize the declared loss; then project the scale/power constraints. At deployment the sender may choose a protection level from the margin and CSIT, but any required signaling and model-sharing cost must be counted. Full two-skeleton storage is $O(N^3)$; matrix-reduction PH can be substantially more expensive. Small episode landmarks and periodic PH evaluation control training cost, but introduce a new approximation whose certificate must be validated. A certificate for the landmarks is not automatically a certificate for all tokens.

## An Error-Tolerant Product-Source Benchmark

To quantify an additional structural rate beyond the zero-type-error result, consider a deliberately separable benchmark: $B$ is a uniform binary type, $U$ is an independent uniform circular state, and $C$ reveals $B$ with probability $\rho$ and otherwise reports an erasure. Both ends know $C$. Here attribute distortion is $d_g^2(U,\widehat U)$ even when the type is wrong, unlike the cross-stratum penalty in [5.1](#eq-5-1). This distinction is essential.

Let $R_U(D)$ be the circular RDF. Conditional independence and the mutual-information chain rule imply $I(B,U;\widehat B,\widehat U|C)\ge I(B;\widehat B|C)+I(U;\widehat U|C)$. Convexity of $R_U$ then yields a lower bound $R_U(D)$ on the attribute contribution. Known types need no additional information; on erased units the binary Hamming RDF is $1-h_2(e)$ bits at error $e\le1/2$. Independent optimal test channels attain the sum. Therefore

<a id="eq-5-17"></a>

$$
\begin{aligned}
R_{\rm sep}(D,\epsilon|C)&=R_U(D)\\
&+(1-\rho)\left[1-h_2\left(\min\left\{\frac{\epsilon}{1-\rho},\frac12\right\}\right)\right].
\end{aligned}
$$

The extra term is zero for $\rho=1$, interpreted separately. Here $h_2$ and all rates in [5.17](#eq-5-17) use base-two logarithms. At $\rho=0.5,\epsilon=0.05$, the additional type rate is $0.26550$ bit/episode. This is an exact special-case source-coding statement, not a universal topology penalty or an analog protection-power formula.

An elementary negative control also separates task semantics from distribution matching. A uniform circle rotated by $\pi$ has exactly the same probability law and pairwise metric, hence identical Rips persistence. However, a fixed binary semantic partition into two semicircles has error probability one after that rotation. Thus equal feature marginals and equal persistence diagrams do not by themselves ensure correct task interpretation.
