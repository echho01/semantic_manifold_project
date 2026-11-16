# 结构裕量、描述符精度与传输资源

从特征误差传到持久图裕量，分别推导区间计数可靠性、量化打包和联合精度分配。保留构造性的编码及高斯界的证明。

[研究索引](../README.md) · [符号和引用所用书目](../../ref.bib)

## System Model and Structural Tasks

<a id="eq-4-1"></a>

### Feature Blocks and Persistence Descriptors

A source block is an indexed finite feature set $X=(x_1,\ldots,x_n)$ with $x_i\in\mathbb{R}^p$. The index identifies feature correspondence before and after reconstruction. Each coordinate lies in a specified interval of length at most $2L$. The metric and its physical or normalized units are shared by transmitter and receiver. Euclidean distance is used for the explicit source quantizer and experiments. A learned metric can be used after its distance-error envelope is calibrated; metric learning is a separate representation stage [`bamberger2026metric`].

For $r\geq0$, the Vietoris--Rips complex $\operatorname{VR}_r(X)$ contains a simplex precisely when all vertex pairs have distance at most $r$. Thus, $r$ denotes the maximum *edge length*. Homology is computed over a fixed field, here $\mathbb F_2$. For degree $k$, let

<a id="eq-4-2"></a>

$$
D_k(X)=\operatorname{Dgm}_k\bigl(\operatorname{VR}_\bullet(X)\bigr)
=\{(u_j,v_j):u_j<v_j\}
$$

denote its finite persistence intervals, represented by their birth and death coordinates. The lifetime is $\ell_j=v_j-u_j$. The essential $H_0$ interval of a nonempty finite Rips filtration is recorded separately; it does not enter the finite endpoint code. All comparisons use the same degree, filtration, and scale units.

The bottleneck distance between two finite diagrams is

<a id="eq-4-3"></a>

$$
d_{\mathrm{B}}(D,D')=\inf_\pi\sup_{z\in D\cup\Delta}\|z-\pi(z)\|_\infty ,
$$

where the diagonal $\Delta=\{(t,t):t\in\mathbb{R}\}$ is included with unlimited multiplicity. Matching $(u,v)$ to the diagonal costs $(v-u)/2$. We use the maximum of [4.3](#eq-4-3) across selected degrees when several degrees are transmitted. This descriptor retains the chosen multiscale information, including its scale; homology and persistence computations follow [`edelsbrunner2010computational`, `chazal2016structure`].

### Source Formats and Physical Channel

In the *feature format*, a uniform scalar quantizer uses $b$ bits per coordinate, with $2^b$ midpoint levels over the known range. Its source length and worst-case feature error satisfy

<a id="eq-4-4"></a>

$$
\begin{gathered}
B_X=npb,\qquad C_X=\sqrt p L,\\
q_b=\max_i\|x_i-Q_b(x_i)\|_2\leq C_X2^{-b}.
\end{gathered}
$$

The range and block dimensions belong to the agreed source format. For unequal coordinate ranges, $C_X$ is one half of the Euclidean norm of the range-length vector. In the *descriptor format*, the encoder sends a quantized list or enumerative code for [4.2](#eq-4-2). Section [4.23](#eq-4-23) develops its source length. Descriptor transmission supports the specified structural queries; feature transmission also supports later coordinate-based processing.

The physical link has $N$ real AWGN channel uses,

<a id="eq-4-5"></a>

$$
Y_t=U_t+W_t,\quad W_t\sim\mathcal N(0,\sigma_w^2),
\qquad \frac1N\sum_{t=1}^N\mathbb{E}U_t^2\leq P .
$$

Source length $B$ is measured in bits per block, while $N$ counts channel uses; the channel noises are independent. For digital delivery, the decoded source message equals the encoded message outside a block-error event of probability $\delta_{\mathrm{ch}}$. On successful decoding, the feature error is bounded by [4.4](#eq-4-4).

We also analyze a linear feature-transmission realization to evaluate structural robustness. The transmitter sends $\alpha Q_b(x_i)$ componentwise over [4.5](#eq-4-5); inverse scaling gives

<a id="eq-4-6"></a>

$$
\widehat x_i=Q_b(x_i)+Z_i,\qquad
Z_i\sim\mathcal N(0,\sigma^2 I_p),\quad
\sigma=\sigma_w/\alpha .
$$

Here $N=np$ and $\alpha$ satisfies the power constraint. The source alphabet has at most $2^{B_X}$ elements; the receiver in [4.6](#eq-4-6) directly reconstructs the features. The finite Gaussian perturbation bound below applies to this realization. More generally, the receiver can supply an envelope

<a id="eq-4-7"></a>

$$
\Pr\!\left\{\max_i\|\widehat x_i-Q_b(x_i)\|_2>a\right\}
\leq\delta_{\mathrm{ch}},
$$

which includes digital delivery with $a=0$.

### Receiver Tasks and Knowledge

A feature-count query reports

<a id="eq-4-8"></a>

$$
c_\theta(D)=\#\{(u,v)\in D:v-u>\theta\},\qquad\theta>0 .
$$

The parameter $\theta$ selects the scale of features relevant to the task. In degree zero, adding the essential interval gives the number of groups retained by this persistence threshold. The task can equivalently identify a separated group of long-lived intervals and permit any threshold in its separation gap.

A knowledge-based structural decision compares the payload descriptor $D_X$ with a reference $K_T$ available at the transmitter:

<a id="eq-4-9"></a>

$$
A_T=\mathbf{1}\{d_{\mathrm{B}}(D_X,K_T)>\tau\},\qquad \tau>0 .
$$

The output flags structural departure from the reference beyond the task tolerance $\tau$. The receiver uses reconstructed descriptors $\widehat D_X$ and $\widehat K$. Its local cache $K_R$ is updated when necessary. The reference can serve $H$ payload blocks, so an update of $B_K$ bits has amortized source cost $B_K/H$. Thresholds and degree identifiers belong to the shared task schema. The communicated object is the specified finite-block descriptor.

## Task Margins and Reliable Feature Precision

<a id="eq-4-10"></a>

### From Feature Errors to Diagram Errors

The following standard consequence of persistence stability fixes the constant required by our edge-length convention.

**Lemma (Rips stability [`chazal2014geometric`]).**

<a id="eq-4-11"></a>

For indexed feature sets $X$ and $\widehat X$, if $\max_i\|x_i-\widehat x_i\|_2\leq e$, then, in every degree,

<a id="eq-4-12"></a>

$$
d_{\mathrm{B}}\bigl(D_k(X),D_k(\widehat X)\bigr)\leq2e .
$$

**Proof.** The reverse triangle inequality bounds every pairwise-distance change by $2e$. The identity on vertex indices maps $\operatorname{VR}_r(X)$ into $\operatorname{VR}_{r+2e}(\widehat X)$ and conversely. These maps commute with filtration inclusions and give a $2e$-interleaving. Stability for finite persistence modules yields [4.12](#eq-4-12).

For a reference diagram and threshold, define the lifetime margin

<a id="eq-4-13"></a>

$$
m(D,\theta)=\min\bigl(\{\theta\}\cup
\{|v-u-\theta|:(u,v)\in D\}\bigr).
$$

The term $\theta$ also controls intervals matched to the diagonal. The margin is positive whenever no reference interval lies exactly on the decision threshold.

**Theorem (Feature-count reliability).**

<a id="eq-4-14"></a>

If $d_{\mathrm{B}}(D,\widehat D)<m(D,\theta)/2$, then $c_\theta(D)=c_\theta(\widehat D)$. Consequently, a feature source code and channel satisfying

<a id="eq-4-15"></a>

$$
C_X2^{-b}+a<\frac{m(D_k(X),\theta)}4
$$

preserve the count with probability at least $1-\delta_{\mathrm{ch}}$ under [4.7](#eq-4-7).

**Proof.** Choose a matching of cost $\eta<m(D,\theta)/2$. For matched off-diagonal intervals, the lifetime changes by at most $2\eta$, which is smaller than the distance to $\theta$. Their positions relative to the threshold therefore agree. Every interval matched to the diagonal has lifetime at most $2\eta<\theta$, so it contributes to neither count. This gives a bijection between counted intervals. On the event in [4.7](#eq-4-7), total feature displacement is at most $C_X2^{-b}+a$. Lemma [4.11](#eq-4-11) completes the proof.

Theorem [4.14](#eq-4-14) controls a discrete task output through a continuous descriptor error. Feature survival and feature classification both enter the margin. An interval can survive while crossing the selected threshold; its distance from $\theta$ determines whether the receiver reports it.

### Finite-Block Gaussian Reliability

For [4.6](#eq-4-6), define

<a id="eq-4-16"></a>

$$
a_\delta=\sigma\left(\sqrt p+\sqrt{2\ln(n/\delta)}\right).
$$

**Corollary (Gaussian feature transmission).**

<a id="eq-4-17"></a>

Under [4.6](#eq-4-6), with probability at least $1-\delta$,

<a id="eq-4-18"></a>

$$
d_{\mathrm{B}}(D_k(X),D_k(\widehat X))
\leq2\bigl(C_X2^{-b}+a_\delta\bigr).
$$

If $a_\delta<m/4$, where $m=m(D_k(X),\theta)$, the smallest integer precision satisfying this bound is

<a id="eq-4-19"></a>

$$
b_{\min}=\max\!\left\{1,
\left\lfloor \log_2\!\frac{C_X}{m/4-a_\delta}\right\rfloor+1\right\}.
$$

For a digital link, the same expression applies with $a_\delta=0$ and failure probability $\delta_{\mathrm{ch}}$.

**Proof.** For a standard Gaussian $G\in\mathbb{R}^p$, $\Pr\{\|G\|_2>\sqrt p+\sqrt{2t}\}\leq e^{-t}$; Appendix [4.47](#eq-4-47) provides a derivation. Apply this bound to each feature with $t=\ln(n/\delta)$ and take a union bound. The result follows from Theorem [4.14](#eq-4-14); the floor plus one enforces the strict inequality in [4.15](#eq-4-15).

Equation [4.16](#eq-4-16) exhibits the block-size cost of a simultaneous structural guarantee. Feature dimension enters through $\sqrt p$, while the number of features enters through a logarithmic maximum-error term. Increasing source precision reduces quantization error. The noise term is controlled through link gain, power, or coding design.

### Threshold Selection and Precision Reduction

Suppose the task separates a nonempty set of retained intervals, whose minimum lifetime is $v$, from discarded intervals, whose maximum lifetime is $u$, with $0\leq u<v$. Set $u=0$ when the discarded set is empty. For any $\theta\in(u,v)$, the reference count is the same and

<a id="eq-4-20"></a>

$$
m(D,\theta)=\min\{\theta-u,v-\theta\}.
$$

**Proposition (Margin-maximizing task threshold).**

<a id="eq-4-21"></a>

Among thresholds selecting this fixed group of intervals,

<a id="eq-4-22"></a>

$$
\theta^*=\frac{u+v}{2},\qquad m^*=\frac{v-u}{2}
$$

maximize the sufficient feature-error budget in [4.15](#eq-4-15), which becomes $(v-u)/8$.

**Proof.** The two affine terms in [4.20](#eq-4-20) sum to $v-u$. Their minimum is maximized when they are equal.

For several degrees, use a threshold for each and the minimum of the resulting margins. The encoder computes the interval lifetimes once and selects precision by [4.19](#eq-4-19). Sorting costs $O(M\log M)$ for $M$ intervals; persistence computation precedes this allocation step. Threshold adaptation is permitted within the task's accepted gap and preserves the same requested reference features over a larger perturbation range.

## Descriptor Rate and Knowledge Synchronization

<a id="eq-4-23"></a>

### Finite Description Bounds

Let $\mathcal{D}_{M,T}$ be the class of diagrams containing at most $M$ finite intervals with $0\leq u<v\leq T$. The interval-count bound $M$ and scale range $T$ are part of the source model. Define $B^*(\epsilon;M,T)$ as the smallest fixed message length for which an encoder and decoder achieve bottleneck error at most $\epsilon$ for every diagram in this class.

**Theorem (Descriptor covering and packing).**

<a id="eq-4-24"></a>

Let $J=\left\lceil T/(2\epsilon)\right\rceil$ and $N_\epsilon=J(J+1)/2$. Then

<a id="eq-4-25"></a>

$$
B^*(\epsilon;M,T)\leq
\left\lceil \log_2\binom{N_\epsilon+M}{M}\right\rceil.
$$

For $0<\epsilon\leq T/64$, let $K_\epsilon=\left\lfloor T/(16\epsilon)\right\rfloor^2$. Then

<a id="eq-4-26"></a>

$$
B^*(\epsilon;M,T)\geq
\left\lceil \log_2\binom{K_\epsilon+M}{M}\right\rceil.
$$

For fixed $M$ and $\epsilon/T\to0$,

<a id="eq-4-27"></a>

$$
B^*(\epsilon;M,T)
=2M\log_2(T/\epsilon)-\log_2(M!)+O(M).
$$

The upper bound uses an unordered multiset of quantized intervals; the $M!$ term in [4.27](#eq-4-27) expresses the leading effect of order invariance. Multiplicities remain meaningful and are encoded. The packing family consists of realizable degree-one diagrams of finite filtered complexes, built from disjoint triangles whose boundaries and faces enter at prescribed times. The converse therefore concerns this descriptor source class. For Euclidean Rips diagrams of a fixed point count, the upper bound still applies; additional geometric constraints can reduce the source class. Appendix [4.49](#eq-4-49) provides the code, packing, and realization.

For implementation, an endpoint list avoids enumerative ranking. Use $2^s$ midpoint levels over $[0,T]$ for each endpoint and sort the intervals lexicographically. The message contains the rounded endpoint pairs; the receiver discards pairs with coincident endpoints. If $m\leq M$ intervals are present,

<a id="eq-4-28"></a> <a id="eq-4-29"></a>

$$
\begin{aligned}
B_K(s)&=h_M+2ms, &h_M&=\left\lceil \log_2(M+1)\right\rceil,
\\
d_{\mathrm{B}}(K_T,Q_s(K_T))&\leq C_K2^{-s}, &C_K&=T/2.
\end{aligned}
$$

The count field uses $h_M$ bits. A collapsed interval has lifetime at most $2C_K2^{-s}$ and diagonal cost at most $C_K2^{-s}$. Endpoint lists provide a direct finite-length update message. Scale, degree, version, and precision fields have fixed lengths determined by the shared protocol and can be added to $h_M$.

### AWGN Channel Uses

The source bounds translate to a physical transmission cost. Real AWGN capacity is

<a id="eq-4-30"></a>

$$
C_{\mathrm{AWGN}}=\frac12\log_2(1+P/\sigma_w^2)
\quad\text{bits per channel use}.
$$

**Corollary (Channel-use lower bound).**

<a id="eq-4-31"></a>

Suppose the receiver reconstructs every diagram in $\mathcal{D}_{M,T}$ to bottleneck error at most $\epsilon$, where $0<\epsilon\leq T/64$, with probability at least $1-\delta$, where $0<\delta<1/2$, using $N$ channel uses of [4.5](#eq-4-5). With $K_\epsilon$ from Theorem [4.24](#eq-4-24),

<a id="eq-4-32"></a>

$$
NC_{\mathrm{AWGN}}\geq(1-\delta)
\log_2\binom{K_\epsilon+M}{M}-h_2(\delta),
$$

where $h_2$ is binary entropy.

**Proof.** Choose the packing diagrams equiprobably. Any reconstruction within $\epsilon$ identifies the transmitted member uniquely. Fano's inequality bounds the required mutual information from below, while the AWGN power constraint bounds it above by $NC_{\mathrm{AWGN}}$ [`cover2006elements`].

For an explicit finite-block construction, transmit a $B$-bit source message by binary phase-shift keying with energy $P$ per use and repeat each bit $r$ times. Averaging gives bit error $Q(\sqrt{rP/\sigma_w^2})$. A Gaussian tail bound and union bound give

<a id="eq-4-33"></a>

$$
P_{\mathrm{block}}\leq\frac B2
\exp\!\left(-\frac{rP}{2\sigma_w^2}\right).
$$

Thus $N=rB$ uses suffice with

<a id="eq-4-34"></a>

$$
r=\max\!\left\{1,\left\lceil \frac{2\sigma_w^2}{P}
\ln\frac{B}{2\delta}\right\rceil\right\}.
$$

Channel codes can replace this explicit construction. Source precision requirements remain unchanged when the required block-error probability is met. In the long-block regime, separated coding approaches the source length divided by [4.30](#eq-4-30) [`cover2006elements`].

### Knowledge Updates and Structural Decisions

Knowledge synchronization targets the descriptor used by the task. Both nodes retain the last acknowledged reference version, so the transmitter knows the cache $K_R$ used by the receiver. It evaluates $e_R=d_{\mathrm{B}}(K_T,K_R)$. A cache within its assigned error budget can be reused with a version indicator. Otherwise, [4.28](#eq-4-28) updates it to $\widehat K=Q_s(K_T)$. Successful delivery gives [4.29](#eq-4-29), independent of the old discrepancy.

For decision [4.9](#eq-4-9), define its source margin

<a id="eq-4-35"></a>

$$
\mu=\left|d_{\mathrm{B}}(D_X,K_T)-\tau\right|.
$$

**Theorem (Joint payload and knowledge reliability).**

<a id="eq-4-36"></a>

If $d_{\mathrm{B}}(D_X,\widehat D_X)\leq\eta_X$ and $d_{\mathrm{B}}(K_T,\widehat K)\leq\eta_K$, then

<a id="eq-4-37"></a>

$$
\left|d_{\mathrm{B}}(\widehat D_X,\widehat K)
-d_{\mathrm{B}}(D_X,K_T)\right|\leq\eta_X+\eta_K.
$$

The receiver decision equals $A_T$ whenever $\eta_X+\eta_K<\mu$. If the two bounds fail with probabilities at most $\delta_X$ and $\delta_K$, the decision error probability is at most $\delta_X+\delta_K$.

**Proof.** Apply the triangle inequality in both directions in [4.37](#eq-4-37). A perturbation smaller than $\mu$ cannot move the decision score across $\tau$. A union bound gives the probability statement without requiring independence.

For $H$ blocks sharing one update, the probability that any decision fails is at most $\delta_K+\sum_{h=1}^H\delta_{X,h}$ when each block satisfies its margin condition. An acknowledged update and its version identifier determine which payloads use that reference.

The decision uses the sum of payload and knowledge errors. Increasing payload precision is most useful after the reference precision is compatible with the task margin. When the reference changes, its endpoint update creates a new operating point with a chosen $\eta_K$. Reuse of a correct cache reduces subsequent update traffic.

## Joint Precision Allocation Protocol

<a id="eq-4-38"></a>

We minimize source cost for the feature format and endpoint knowledge code. Let $a$ satisfy the payload link envelope, and choose positive engineering slack $\zeta$. For a block with margin $\mu$, define

<a id="eq-4-39"></a>

$$
E=\mu-2a-\zeta>0.
$$

The finite-precision design is

<a id="eq-4-40"></a>

$$
\begin{aligned}
\min_{b,s\in\mathbb N_{\geq1}}\quad&
npb+\frac{h_M+2ms}{H}
\\
\text{subject to}\quad&
2C_X2^{-b}+C_K2^{-s}\leq E.
\end{aligned}
$$

The slack enforces the strict decision condition in Theorem [4.36](#eq-4-36). For blocks sharing an update, use a common minimum margin or one constraint per block with its own $b$. The number $H$ is the planned reuse count; version agreement is maintained by the protocol.

**Theorem (Optimal continuous allocation).**

<a id="eq-4-41"></a>

Write $a_X=np$ and $a_K=2m/H$. In the interior precision regime with $m>0$, the continuous relaxation of [4.40](#eq-4-40) assigns source-induced diagram errors

<a id="eq-4-42"></a>

$$
e_X^*=\frac{a_X}{a_X+a_K}E,\qquad
e_K^*=\frac{a_K}{a_X+a_K}E,
$$

and precisions

<a id="eq-4-43"></a>

$$
b^*=\log_2\frac{2C_X}{e_X^*},\qquad
s^*=\log_2\frac{C_K}{e_K^*}.
$$

If $E\geq C_X+C_K/2$, the optimal precisions are $(1,1)$. Otherwise, if one unconstrained precision is below one, set it to one and allocate the remaining error budget to the other source. Rounding both interior precisions upward is feasible and adds less than $np+2m/H$ bits per block to the relaxed objective.

**Proof.** For $e_X=2C_X2^{-b}$ and $e_K=C_K2^{-s}$, the variable cost is

<a id="eq-4-44"></a>

$$
a_X\log_2\frac{2C_X}{e_X}
+a_K\log_2\frac{C_K}{e_K}.
$$

It is strictly convex for positive errors and decreases in each error, so the constraint is active. Lagrange stationarity gives $a_X/e_X=a_K/e_K$, yielding [4.42](#eq-4-42). Constraints $b,s\geq1$ bound the allowed errors above; convexity gives the boundary rule. Rounding upward decreases quantization error and increases each precision by less than one, proving the overhead bound.

The larger error share is assigned to the object whose precision is more expensive. A knowledge descriptor with few intervals can be encoded accurately at small source cost. Increasing its reuse count further favors a precise update, because the additional endpoint bits are shared across blocks.

An exact integer solution is inexpensive. Starting with a feasible rounded pair of cost $J_0$, every better solution has $b\leq\left\lfloor J_0/(np)\right\rfloor$. For each such $b$ satisfying $E-2C_X2^{-b}>0$, choose

<a id="eq-4-45"></a>

$$
s(b)=\max\!\left\{1,
\left\lceil \log_2\frac{C_K}{E-2C_X2^{-b}}\right\rceil\right\}
$$

and evaluate the objective. This bounded search returns the exact minimum for the stated formats. A nonempty source and an empty knowledge diagram require only the count field for the latter; the optimization then reduces to payload precision alone.

<figure data-latex-placement="t">
<p><a id="eq-4-46"></a></p>
<figcaption>The protocol selects source precision from the structural task margin and amortizes knowledge precision over its reuse period.</figcaption>
</figure>

If the task includes both counts and a reference decision, impose [4.15](#eq-4-15) and [4.40](#eq-4-40) together. The search starts at the precision required by the count task. For descriptor payloads, replace $npb$ and $2C_X2^{-b}$ by the descriptor length and error from [4.28](#eq-4-28)--[4.29](#eq-4-29); the same allocation argument applies. Persistence can be computed using Rips algorithms such as Ripser [`bauer2021ripser`]; coding and allocation operate on the resulting finite intervals.

## Gaussian Maximum-Error Bound

<a id="eq-4-47"></a>

For $V=\|G\|_2^2$ with $G\sim\mathcal N(0,I_p)$, and $0<\lambda<1/2$, the moment-generating function gives

<a id="eq-4-48"></a>

$$
\begin{aligned}
\ln\mathbb{E}e^{\lambda(V-p)}
&=-\frac p2\ln(1-2\lambda)-p\lambda\\
&\leq\frac{p\lambda^2}{1-2\lambda}.
\end{aligned}
$$

Expand $-\ln(1-z)-z$ and bound its terms by $z^2/[2(1-z)]$ for $0<z<1$ to obtain the inequality. Chernoff's inequality with $\lambda=\sqrt{t/p}/(1+2\sqrt{t/p})$ gives

$$
\Pr\{V-p\geq2\sqrt{pt}+2t\}\leq e^{-t}.
$$

Since $(\sqrt p+\sqrt{2t})^2\geq p+2\sqrt{pt}+2t$, the norm tail used in Corollary [4.17](#eq-4-17) follows. The union bound over features requires their marginal norm tails, without requiring independence between features.

## Proof of the Descriptor Coding Bounds

<a id="eq-4-49"></a>

### Cover and Encoder

Let the endpoint grid contain $0,2\epsilon,4\epsilon,\ldots,T$, with its last spacing shortened if necessary. It has $J+1$ distinct values. Round each endpoint to a nearest grid value, using a monotone tie rule. Each coordinate moves by at most $\epsilon$. If rounded endpoints coincide, the original lifetime is at most $2\epsilon$, so matching to the diagonal costs at most $\epsilon$. Every other interval moves by at most $\epsilon$ in the $\ell_\infty$ norm. Thus reconstruction error is at most $\epsilon$.

There are $N_\epsilon=\binom{J+1}{2}$ pairs of distinct grid endpoints. A diagram is a multiset of at most $M$ such pairs. Introduce one additional category for unused interval slots. There are $\binom{N_\epsilon+M}{M}$ nonnegative integer count vectors summing to $M$. An agreed lexicographic ranking yields the fixed-length code [4.25](#eq-4-25); the decoder un-ranks the vector. Multiple intervals with identical endpoints are represented by their count.

### Packing and Realization

Let $q=\left\lfloor T/(16\epsilon)\right\rfloor$, and choose the $q^2$ pairs

<a id="eq-4-50"></a>

$$
\mathcal G=\{(4\epsilon i,\,3T/4+4\epsilon j):
0\leq i,j<q\}.
$$

They lie in $[0,T]^2$ above the diagonal. Different points are separated by at least $4\epsilon$ in $\ell_\infty$, and each has diagonal distance at least $T/4>2\epsilon$. Form all multisets of at most $M$ points from $\mathcal G$. If two diagrams differ, some excess point must match to a different grid point or to the diagonal, at cost exceeding $2\epsilon$. An $\epsilon$-accurate code cannot assign the diagrams the same message. Counting multisets gives [4.26](#eq-4-26).

Every packing diagram is realized in $H_1$ by a finite filtered simplicial complex. For each $(u,v)$, create a separate triangle. Put its vertices into the complex at time zero, its three edges at $u$, and its face at $v$. The boundary creates one degree-one class at $u$, and the face kills it at $v$. Disjoint union yields the selected multiset, including repeated intervals. A single vertex realizes the empty diagram.

Finally, for fixed $M$ and $z\to\infty$, $\log_2\binom{z+M}{M}=M\log_2z-\log_2(M!)+o(1)$. Both $N_\epsilon$ and $K_\epsilon$ are constant multiples of $(T/\epsilon)^2$ up to rounding, which gives [4.27](#eq-4-27).
