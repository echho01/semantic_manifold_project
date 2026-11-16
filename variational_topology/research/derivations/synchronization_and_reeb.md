# 知识字典失配、同步与 Reeb 分支信息

先算不同知识库之间的极限误差，再给局部同步的消息数和任务分支的条件信息量。

[研究索引](../README.md) · [符号和引用所用书目](../../ref.bib)

## Knowledge-Base Mismatch and Synchronization

<a id="eq-3-1"></a>

Transmitter and receiver represent the same items $t\in\mathcal T$ by continuous onto maps $\Phi_T:\mathcal T\to\mathcal{M}_T$ and $\Phi_R:\mathcal T\to\mathcal{M}_R$. A dictionary $P_K\subset\mathcal T$ contains $K$ shared item identities. An *interpretation-constrained receiver* reconstructs their geometry as $\Phi_R(P_K)$; it does not learn an alignment map from additional messages. Write $\varepsilon_R(K)=d_H(\Phi_R(P_K),\mathcal{M}_R)$. We first use offset diagrams, denoted by $\mathcal D_{T,k}$ and $\mathcal D_{R,k}$, and define the KB discrepancy $\delta_{{\rm KB},k}=d_{\mathrm{B}}(\mathcal D_{T,k},\mathcal D_{R,k})$.

Define the finite-dictionary distortion by $d_{\rm topo}(K)=d_{\mathrm{B}}(\mathrm{Dgm}_k(\Phi_R(P_K)),\mathcal D_{T,k})$.

### The dense-dictionary limit

**Theorem (Mismatch limit and finite-resolution certificate).**

<a id="eq-3-2"></a>

For every dictionary size $K$ and degree $k$,

<a id="eq-3-3"></a>

$$
\big|d_{\rm topo}(K)-\delta_{{\rm KB},k}\big|
\le\varepsilon_R(K).
$$

Consequently $\varepsilon_R(K)\to0$ implies $d_{\rm topo}(K)\to\delta_{{\rm KB},k}$. For finite reference samples $X_T,X_R$ with coverage errors $\varepsilon_T,\varepsilon_R$, the estimate $\widehat\delta_k=d_{\mathrm{B}}(\mathrm{Dgm}_k(X_T),\mathrm{Dgm}_k(X_R))$ satisfies

<a id="eq-3-4"></a>

$$
|\widehat\delta_k-\delta_{{\rm KB},k}|
\le\varepsilon_T+\varepsilon_R.
$$

Both statements also hold for radius-convention Rips diagrams when Rips filtrations are used consistently, including on the supports.

**Proof.** Offset stability gives $d_{\mathrm{B}}(\mathrm{Dgm}_k(A),\mathrm{Dgm}_k(B))\le d_H(A,B)$. The reverse triangle inequality proves [3.3](#eq-3-3); applying the triangle inequality twice proves [3.4](#eq-3-4). Radius-convention Rips stability has the same Hausdorff constant [`chazal2014persistence`].

The persistent mismatch is an asymptotic obstruction for this interpretation rule: $d_{\rm topo}(K)\ge(\delta_{{\rm KB},k}-\varepsilon_R(K))_+$. It is not a lower bound at every finite $K$. For example, concentric circles of radii $1$ and $2$ have offset $H_1$ diagrams $\{(0,1)\}$ and $\{(0,2)\}$, with $\delta_{\rm KB}=1$. A one-atom receiver dictionary has no $H_1$ bar and hence distortion $1/2$, below the limit. Increasing dictionary resolution can reveal more of the receiver's mismatch. Synchronization changes the interpretation rule and can remove the limit.

### Metric deformation and homological mismatch

**Proposition (Two mechanisms of KB discrepancy).**

<a id="eq-3-5"></a>

**(a)** For corresponding clouds $X_T,X_R$, put $l_{ij}=\|x_i-x_j\|$. If $|\|y_i-y_j\|-l_{ij}|\le e l_{ij}$ for all pairs, where $0\le e<1$, their Rips diagrams obey $d_{\mathrm{B}}\le eL$, with $L=\frac12\operatorname{diam}(X_T)$.

**(b)** For compact positive-reach supports, if $\beta_k(\mathcal{M}_T)\ne\beta_k(\mathcal{M}_R)$, their *offset* diagrams satisfy

<a id="eq-3-6"></a>

$$
\delta_{{\rm KB},k}\ge
\tfrac12\min\{\tau_T,\tau_R\}.
$$

**Proof.** (a) The filtration value of each simplex is half its largest edge length, so it changes by at most $eL$; apply stability. (b) Suppose $\beta_k(\mathcal{M}_T)>\beta_k(\mathcal{M}_R)$. By Lemma [1.8](coverage_and_dictionary.md#eq-1-8), the transmitter has $\beta_k(\mathcal{M}_T)$ bars born at zero and surviving to $\tau_T$. The receiver has only $\beta_k(\mathcal{M}_R)$ such bars; every other bar is born at least at $\tau_R$. An excess transmitter bar must be matched either to the diagonal, at cost at least $\tau_T/2$, or to a later-born receiver bar, at cost at least $\tau_R$. Interchange $T$ and $R$ for the other case.

A small metric deformation, informally an *accent*, has a proportional Rips discrepancy. A homology change, informally a *language* mismatch, has the offset lower bound [3.6](#eq-3-6). These are geometric mechanisms, not a classification of natural languages or learned models. In particular, a small nonzero diagram distance need not indicate different Betti numbers, and a vanishing reach makes the lower bound uninformative.

### An optimal local synchronization code

The receiver holds $K$ finite descriptor points $p_i=(b_i,d_i)$ and knows their identities. The transmitter's descriptor belongs to the local class 

$$
\mathcal K=\big\{\{p_i+o_i\}_{i=1}^K:
o_i\in[-\varepsilon_0,\varepsilon_0]^2\big\}.
$$

This is a source class of diagrams with fixed feature correspondence; it permits independent endpoint updates.

**Theorem (Exact fixed-length synchronization cost).**

<a id="eq-3-7"></a>

Assume $b_i\ge\varepsilon_0$, $d_i-b_i>6\varepsilon_0$, and $\|p_i-p_j\|_\infty>4\varepsilon_0$ for $i\ne j$. For $0<\varepsilon_1<\varepsilon_0$, let $m=\lceil\varepsilon_0/\varepsilon_1\rceil$. The minimum fixed binary message length for reconstructing every $\mathcal D_T\in\mathcal K$ with bottleneck error at most $\varepsilon_1$ is

<a id="eq-3-8"></a>

$$
B^*=\left\lceil\log_2(m^{2K})\right\rceil.
$$

**Proof.** Partition each offset interval into $m$ equal cells and encode the $2K$ cell indices as one base-$m$ integer in $[0,m^{2K}-1]$. Decoding to cell centers gives coordinate error $\varepsilon_0/m\le\varepsilon_1$. For the converse, use $m$ equally spaced values including both endpoints of the interval. Their spacing exceeds $2\varepsilon_1$. The separation assumptions force every nonidentity matching or diagonal deletion to cost more than $2\varepsilon_0$. The resulting $m^{2K}$ diagrams are therefore pairwise more than $2\varepsilon_1$ apart and require distinct messages. Appendix [3.16](#eq-3-16) gives the endpoint and rounding details.

The joint integer code improves on separately coding each coordinate, which uses $2K\lceil\log_2m\rceil$ bits. For $\varepsilon_1=\varepsilon_0/2^s$, each additional refinement level costs exactly two bits per feature. Feature creation, deletion, and uncertain correspondence require a larger source class and additional signaling. Equation [3.8](#eq-3-8) is exact for the stated local class, which need not be realizable by a fixed-size Euclidean Rips cloud.

## Task-Oriented Topological Invariants

<a id="eq-3-9"></a>

Let $f:\mathcal{M}\to\mathbb{R}$ be a Morse function with distinct critical values. Its Reeb quotient $q:\mathcal{M}\to\mathcal{R}_f$ identifies points in the same connected component of a level set. Away from critical values, $q(S)$ is represented by its edge $E$ and score $T=f(S)$ [`dey2013reeb`]. We distinguish the scalar task $T$ from the *branch-aware task* $(E,T)$, which also identifies the region producing that score. For the latter, use

<a id="eq-3-10"></a>

$$
d((e,t),(\widehat e,\widehat t))=
\begin{cases}
(t-\widehat t)^2,&e=\widehat e,\\
d_\times,&e\ne\widehat e,
\end{cases}\qquad d_\times>0.
$$

Rates below are Shannon rates in bits per iid source symbol. Every edge with positive mass has a bounded score interval; assume its conditional score density has finite differential entropy. These conditions allow the integrable density singularities near Morse critical values.

**Theorem (Rate overhead of branch information).**

<a id="eq-3-11"></a>

Under these conditions, as $D\downarrow0$,

<a id="eq-3-12"></a> <a id="eq-3-13"></a> <a id="eq-3-14"></a>

$$
\begin{aligned}
R_{\rm task}(D)
&=H(E)+h(T\mid E)\\
&\quad-\tfrac12\log_2(2\pi eD)+o(1),\\
R_f(D)&=h(T)-\tfrac12\log_2(2\pi eD)+o(1),\\
R_{\rm task}(D)-R_f(D)&\longrightarrow H(E\mid T).
\end{aligned}
$$

If $b(t)$ is the number of components of the regular level $f^{-1}(t)$, then $H(E\mid T)\le\mathbb{E}\log_2 b(T)\le\log_2 b_{\max}$.

**Proof.** Transmit the finite edge index losslessly and apply asymptotically optimal scalar coding conditionally on each edge. The scalar Shannon lower bound is asymptotically tight for bounded-support densities of finite entropy [`linder1994asymptotic`, `koch2016shannon`]. For a converse, distortion $D$ implies edge-error probability at most $D/d_\times$. Fano's inequality bounds the remaining edge uncertainty, and a conditional Gaussian entropy bound on correctly decoded edges bounds the score uncertainty. The error events contribute $O(D|\log D|)$, leaving the constant in [3.12](#eq-3-12) unchanged; Appendix [3.17](#eq-3-17) proves this explicitly. Finally, $H(E)+h(T\mid E)-h(T)=H(E\mid T)$.

This overhead is required by the branch-aware output specification. A scalar-score task retains the rate [3.13](#eq-3-13) even when its level sets are disconnected. For entropy-coded scalar quantizers the exact finite-resolution identity is 

$$
H(E,Q_\Delta(T))-H(Q_\Delta(T))=H(E\mid Q_\Delta(T)).
$$

For refining partitions this tends to $H(E\mid T)$. It verifies the branch overhead while leaving the scalar-quantizer gap to the Shannon rate explicit.

**Proposition (Task sufficiency and homological screening).**

<a id="eq-3-15"></a>

For any measurable distortion that depends on $S$ only through $q(S)$, the source and its pushforward $q_\#p$ have identical task rate--distortion functions. For the Reeb quotient of a smooth Morse function, $q_*:H_1(\mathcal{M})\to H_1(\mathcal{R}_f)$ is surjective; hence $\beta_1(\mathcal{R}_f)\le\beta_1(\mathcal{M})$ and $H_k(\mathcal{R}_f)=0$ for $k\ge2$.

**Proof.** Data processing gives $I(S;\widehat U)\ge I(q(S);\widehat U)$. Conversely, any conditional reconstruction law given $q(S)$ defines a law given $S$ with the same distortion and mutual information, proving equality of the rate--distortion functions. The quotient induces a surjection of fundamental groups on each connected component for $\dim\mathcal{M}\ge2$ [`kaluba2015reeb`]; abelianization and direct sums give the $H_1$ statement. In dimension one, every level component is a point and the quotient is a homeomorphism. A graph has no higher-dimensional homology.

Thus a task code can operate on the quotient distribution without reconstructing the homology killed by $q_*$. A persistence diagram is a stable summary of this structure, not a complete encoding of Reeb graph branches or their probabilities. If the receiver approximates $f$ on the same domain with $\|\widehat f-f\|_\infty\le u$, sublevel diagrams satisfy $d_{\mathrm{B}}(\mathrm{Dgm}(\widehat f),\mathrm{Dgm}(f))\le u$. This alone does not identify the correct branch.

Near the finite set $\mathcal C$ of critical values, the probability that a score perturbation of size $u$ can cross a critical level is bounded by $\Pr[\operatorname{dist}(T,\mathcal C)\le u]$. This is $o(1)$ for an absolutely continuous score and at most $2|\mathcal C|Lu$ if its density is bounded by $L$ there. The torus example has a saddle singularity, so a universal linear critical-level error law is not assumed.

## Synchronization and Task-Rate Proofs

### Exact synchronization message count

<a id="eq-3-16"></a>

Let $m=\lceil\varepsilon_0/\varepsilon_1\rceil\ge2$. Each of the $2K$ offset coordinates lies in $[-\varepsilon_0,\varepsilon_0]$. Divide this interval into $m$ equal cells, assigning a boundary to the cell on its right, except the last endpoint. The cell centers remain in the interval, and their maximum rounding error is $\varepsilon_0/m\le\varepsilon_1$. For cell indices $a_0,\ldots,a_{2K-1}\in\{0,\ldots,m-1\}$, send $J=\sum_{\ell=0}^{2K-1}a_\ell m^\ell$ as a $\lceil\log_2(m^{2K})\rceil$-bit integer. Repeated division by $m$ recovers every index. Because $b_i\ge\varepsilon_0$ and $d_i-b_i>6\varepsilon_0$, all encoded and decoded points have valid nonnegative births and positive persistence.

For the converse take $G=\{-\varepsilon_0+2j\varepsilon_0/(m-1):0\le j<m\}$. Its spacing is strictly greater than $2\varepsilon_1$ because $m-1<\varepsilon_0/\varepsilon_1$. Consider all $m^{2K}$ diagrams with offsets in $G^2$. Any matching of different diagrams that pairs each reference feature with itself has cost greater than $2\varepsilon_1$. A cross-feature pair has cost at least $\|p_i-p_j\|_\infty-2\varepsilon_0>2\varepsilon_0$; a diagonal deletion costs at least $(d_i-b_i-2\varepsilon_0)/2>2\varepsilon_0$. Every matching thus costs more than $2\varepsilon_1$. Two of these diagrams cannot share a reconstruction of error at most $\varepsilon_1$, by the triangle inequality. A $B$-bit message offers at most $2^B$ reconstructions, proving the matching lower bound in [3.8](#eq-3-8).

### A converse with a finite branch-error penalty

<a id="eq-3-17"></a>

Put $U=(E,T)$ and $\widehat U=(\widehat E,\widehat T)$, and let $p_\times=\Pr(E\ne\widehat E)\le D/d_\times$. There are $J$ positive-mass edges, and all score intervals have length at most $A>0$. With entropy measured relative to counting measure on edges and Lebesgue measure on scores, 

$$
I(U;\widehat U)=H(E)+h(T\mid E)
-H(E\mid\widehat U)-h(T\mid E,\widehat U).
$$

Fano's inequality gives $H(E\mid\widehat U)\le h_2(p_\times)+p_\times\log_2 J$. On correctly decoded edges, the conditional mean-square error, averaged over the conditioning variables, is at most $D/(1-p_\times)$. The maximum-entropy bound and Jensen's inequality therefore give

$$
\begin{aligned}
h(T\mid E,\widehat U)
&\le\frac{1-p_\times}{2}
\log_2\frac{2\pi eD}{1-p_\times}
+p_\times\log_2 A .
\end{aligned}
$$

On error events we used only the bounded support of $T$. Uniformly for $0\le p_\times\le D/d_\times$, these two bounds yield 

$$
I(U;\widehat U)\ge H(E)+h(T\mid E)
-\tfrac12\log_2(2\pi eD)-o(1).
$$

Indeed, $h_2(p_\times)=O(D|\log D|)$ and $p_\times|\log D|=o(1)$.

For achievability, reproduce $E$ exactly and, on each edge, use its scalar rate--distortion test channel with distortion at most $D$. Projection of a scalar reconstruction onto that edge's closed score interval cannot increase squared error. The finite-entropy, bounded-support scalar asymptotic [`linder1994asymptotic`] then gives 

$$
I(U;\widehat U)\le H(E)+h(T\mid E)
-\tfrac12\log_2(2\pi eD)+o(1).
$$

There are finitely many edges, so their remainder terms average to $o(1)$. The same scalar result proves [3.13](#eq-3-13). Given a regular score $t$, the edge has at most $b(t)$ possible values, which proves the entropy bound. For a simple Morse function with $c$ critical points, every vertex has degree at most three and there are at most $3c/2$ edges; this also bounds $b_{\max}$.
