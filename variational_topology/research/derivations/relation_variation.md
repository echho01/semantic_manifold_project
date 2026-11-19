# 关系泛函的一阶变分：临界载荷、像素影响与变分功率分配

把"带位置身份的连通询问"写成表示上的泛函，求它的精确一阶变分，再由变分推出资源分配。与持久图泛函的对比放在第 3 节。

[研究索引](../README.md) · [连通风险与 MST](connectivity.md) · [对称轨道率失真与衰落逆定理](orbit_rdf_and_fading.md) · [实现](../experiments/relation_variation.py) · [测试](../tests/test_variation.py) · [数值](../results/relation_variation.json)

## 0 Setting

$G=(\mathcal V,\mathcal E)$ is a finite connected graph with $n$ vertices, $f:\mathcal V\to\mathbb R$ a score field, $m_f(a,b)=\min_{\pi:a\leadsto b}\max_{v\in\pi}f(v)$ the merge level and $B_f(a,b,t)=\mathbf 1\{m_f(a,b)\le t\}$. The query law $\nu$ draws the endpoints and, independently, a threshold with continuous CDF $F$. The **relation functional** is

$$
\mathcal R[\hat f;f]=d_\nu(f,\hat f)=\Pr_\nu[B_f\ne B_{\hat f}]=\mathbb E_{(a,b)}\big|F(m_f(a,b))-F(m_{\hat f}(a,b))\big| .
$$

## 1 Critical loads at a generic field

Let $f$ be injective with smallest gap $\gamma_f$. Activate vertices in increasing order of $f$. When $v$ enters, the active neighbors of $v$ lie in components of sizes $n_{v,1},n_{v,2},\dots$ (possibly none), which merge with $v$ into a component of size $N_v=1+\sum_jn_{v,j}$.

**Lemma.** For every pair $(a,b)$ there is exactly one vertex $v^\star(a,b)$ with $m_f(a,b)=f(v^\star(a,b))$, and it is the vertex whose activation first puts $a$ and $b$ in one component.

*Proof.* $B_f(a,b,t)$ switches from 0 to 1 at $t=m_f(a,b)$, which is a value of $f$ and hence the value of a single vertex. Connectivity of $G$ makes $m_f$ finite. $\square$

**Theorem.** Let $c_f(v)=\Pr_{(a,b)}[v^\star(a,b)=v]$.

1. $\sum_vc_f(v)=1$. For uniform ordered pairs with replacement, $c_f(v)=(N_v^2-\sum_jn_{v,j}^2)/n^2$.
2. If $\|\delta\|_\infty<\gamma_f/2$ then $\mathcal R[f+\delta;f]=\sum_vc_f(v)\,|F(f(v)+\delta(v))-F(f(v))|$.
3. For injective $\hat f$ and a uniform threshold, at points where $m_{\hat f}\ne m_f$ almost surely, $\partial\mathcal R/\partial\hat f(v)=\mathbb E[\operatorname{sgn}(m_{\hat f}-m_f)\mathbf 1\{v^\star_{\hat f}=v\}]$, and $\sum_v|\partial\mathcal R/\partial\hat f(v)|\le1$.

*Proof.* (1) The events $\{v^\star=v\}$ partition the pairs. After $v$ enters, the ordered pairs inside its component number $N_v^2$; those already inside one earlier component number $\sum_jn_{v,j}^2$. (2) For $u\ne w$ with $f(u)<f(w)$, $f(w)+\delta(w)-f(u)-\delta(u)>\gamma_f-2\|\delta\|_\infty>0$, so $f+\delta$ is injective with the same order. The activation sequence, the component history and $v^\star$ are unchanged, hence $m_{f+\delta}(a,b)=f(v^\star)+\delta(v^\star)$. Insert in the definition of $\mathcal R$ and group pairs by $v^\star$. (3) By the same argument $m_{\hat f}$ is locally the coordinate function $\hat f(v^\star_{\hat f})$; differentiate under the expectation. The bound follows from $\sum_v\Pr[v^\star_{\hat f}=v]=1$. $\square$

**Remarks.**

- Every vertex has positive load: the pair $(v,v)$ gives $c_f(v)\ge1/n^2$. A vertex that only extends one component has load $(2n_{v,1}+1)/n^2$. A vertex joining two components of sizes $n_1,n_2$ has load $(2n_1n_2+2n_1+2n_2+1)/n^2$.
- The functional is not differentiable at $\hat f=f$: it has a kink, and its one-sided variation is the weighted $\ell^1$ norm $\sum_vc_f(v)|\delta(v)|$. The local geometry of the relation pseudometric is of Finsler type. It is not Riemannian, so a quadratic (Fisher-type) expansion does not apply to this task.
- One union–find sweep in sorted order computes all loads in $O(n\,\alpha(n))$ after sorting.
- Beyond $\gamma_f/2$ the order can change and the formula can fail; the test suite exhibits failures at $3\gamma_f$.

**数值核对。** 100 个随机 $6\times6$ 单射场：载荷公式与逐对枚举的最大差为 0，$\sum_vc_f=1$ 的误差 $2\times10^{-16}$，局部律与直接计算的 $\mathbb E|m_{f+\delta}-m_f|$ 最大差 $2\times10^{-17}$。

## 2 General variation and pixel influence

**Directional derivative with ties.** $m_f=\min_\pi\max_{v\in\pi}f(v)$ is a finite min of finite maxima of coordinates. The directional derivative of a finite maximum is the maximum of the directional derivatives of its active terms, and dually for a minimum. Hence

$$
m_f'(a,b;\delta)=\min_{\pi\in\Pi^\star(a,b)}\ \max_{v\in\pi:\,f(v)=m_f(a,b)}\delta(v),
$$

where $\Pi^\star(a,b)$ are the paths attaining the minimax value. Choosing one active branch, as the MST implementation does, gives an element of this directional structure. (Checked by path enumeration on $3\times3$ fields with tied values.)

**Binary images.** For $x\in\{0,1\}^{\mathcal V}$ let $\mathrm{Inf}_x(v)=d_\nu(x,x\oplus e_v)$. Let the pixels flip independently with probabilities $\varepsilon_v$, $S=\sum_v\varepsilon_v$, $\Lambda=\sum_v\varepsilon_v\mathrm{Inf}_x(v)$.

**Proposition.** $-S\Lambda\le\mathbb E\,d_\nu(x,x\oplus Z)-\Lambda\le S^2/2$.

*Proof.* $\Pr[\text{only }v\text{ flips}]=\varepsilon_v\prod_{u\ne v}(1-\varepsilon_u)\in[\varepsilon_v(1-S),\varepsilon_v]$. Lower bound: keep only the single-flip events. Upper bound: $\Pr[\ge2\text{ flips}]\le\sum_{u<v}\varepsilon_u\varepsilon_v\le S^2/2$ and $d_\nu\le1$. $\square$

The constant $3/2$ of the first plan is not needed; $1/2$ suffices for the upper side and the lower side is multiplicative. The bound is a worst case: it is informative only when $S\ll\Lambda/S$, which fails for $32\times32$ images at practical error rates. The measured law is much better (below).

Under uniform ordered pairs with replacement and a threshold in $(0,1)$:

- erasing a foreground pixel $v$ of a component of size $s$ that splits into parts of sizes $s_1,\dots,s_k$ ($\sum s_i=s-1$) has influence $(s^2-\sum_is_i^2)/n^2$, largest at cut pixels;
- adding a background pixel that touches components of sizes $n_1,\dots,n_k$ has influence $((1+\sum n_i)^2-\sum n_i^2)/n^2$; an isolated one has $1/n^2$.

## 3 Contrast with diagram functionals

The sublevel $H_0$ diagram of a generic $f$ has one birth per local minimum and one death per component absorbed at a merge vertex (elder rule). A functional of the diagram therefore has a first variation supported on local minima and merge vertices, with weights that depend on the diagram point and not on the sizes $n_{v,j}$ of the merged components. Three consequences:

1. it is constant on automorphism orbits of the readout graph (Proposition on automorphisms in the paper), where the relation functional is positive;
2. it gives zero weight to vertices that extend a component, although each such vertex decides $2n_{v,1}+1$ ordered pairs;
3. a merge of two pixels and a merge of two large strokes receive weights unrelated to the $2n_1n_2$ pairs they connect.

Batchwise latent persistence adds a fourth limit: it describes the sampled batch, and identifying the homology of an underlying manifold needs sampling and reach conditions that training batches do not provide.

## 4 Variational allocation

**Proposition.** For $w_v>0$ and $\varphi$ convex, decreasing, differentiable, the minimizer of $\sum_vw_v\varphi(p_v)$ subject to $\sum_vp_v=P$, $p_v>0$, satisfies $w_v\varphi'(p_v)=-\mu$ (convex separable program; the Lagrange condition is sufficient).

(a) *Analog, generic field.* $\delta(v)\sim\mathcal N(0,N_0/p_v)$, $\mathbb E|\delta(v)|=\sqrt{2N_0/(\pi p_v)}$. In the regime of Section 1 the risk is $\sqrt{2N_0/\pi}\sum_vc_f(v)p_v^{-1/2}$. Stationarity: $c_f(v)p_v^{-3/2}=\text{const}$, so $p_v=P\,c_f(v)^{2/3}/\sum_uc_f(u)^{2/3}$ with value $\sqrt{2N_0/(\pi P)}\,(\sum_vc_f(v)^{2/3})^{3/2}$. Uniform power gives $\sqrt{2N_0n/(\pi P)}$. The ratio is $\big(\frac1n\sum_v(n\,c_f(v))^{2/3}\big)^{3/2}\le1$ by Jensen's inequality for $t\mapsto t^{2/3}$, with equality only for uniform loads. The Gaussian perturbation is unbounded, so the formula is the linearized risk; it is exact only on the event $\|\delta\|_\infty<\gamma_f/2$.

(b) *Variance certificate.* $\varphi(p)=1/p$ gives $p_v\propto\sqrt{w_v}$. With $w_j=\sigma_j^2/|h_j|^2$ this is the allocation of the two-mode experiment: gains $1$ and $0.22$ give shares $(0.1803,0.8197)$.

(c) *BPSK, binary image.* $\varepsilon_v=Q(\sqrt{2p_v/N_0})$ and first-order risk $\sum_v\mathrm{Inf}_x(v)Q(\sqrt{2p_v/N_0})$. Since $\frac{d}{dp}Q(\sqrt{2p/N_0})=-e^{-p/N_0}/(2\sqrt{\pi pN_0})$ is negative and increasing, the objective is convex. Stationarity reads $\mathrm{Inf}_x(v)e^{-p_v/N_0}/(2\sqrt{\pi p_vN_0})=\mu$. With $u=2p_v/N_0$: $u\,e^{u}=\mathrm{Inf}_x(v)^2/(2\pi\mu^2N_0^2)$, i.e. $p_v=\tfrac{N_0}2W_0\big(\mathrm{Inf}_x(v)^2/(2\pi\mu^2N_0^2)\big)$. The multiplier is found by bisection on the power constraint. A sign detector does not need to know the allocation.

The protection exponents are $2/3$ (analog, $\ell^1$ risk), $1/2$ (variance certificate) and logarithmic (digital), since $W_0(y)\sim\ln y$.

**规则 (a) 的适用范围（2026-11-18 实测）。** (a) 只在 $\|\delta\|_\infty<\gamma_f/2$ 内是精确线性化的最优解。在 40 个 $16\times16$ 灰度字形场上用精确全点对风险和公共噪声测得：噪声标准差 $10^{-6}$ 时分配/均匀风险比为 0.76（线性化预测 0.73）；$10^{-3}$、$10^{-2}$、0.03、0.1 时为 1.84、1.08、1.06、1.12，即比均匀功率更差。主稿命题因此不再陈述模拟规则，只保留方差型与 BPSK 两条，并把这一结果作为对照写入实验节。[脚本](../experiments/allocation_controls.py) · [数值](../results/allocation_controls.json)

**数值核对（合成）。** $6\times6$ 随机单射场上 (a) 的风险比均值 0.771，范围 0.684–0.883。测试验证了 (a) 的闭式不被 200 组随机可行分配超过，(c) 的解满足等斜率条件并不被随机扰动分配超过。

## 5 Omniglot 上的实测（2026-11-16）

协议与镜像实验相同：20 个 evaluation 字母表各 16 张，种子 20261115，$32\times32$ 最近邻缩放，阈值 128，四邻接，全部像素为合法端点。管线核对：ID 哈希一致，镜像均值 0.012096 与 0.778214，逐图差为 0。

| 量 | 结果 |
|---|---|
| 总影响 $\sum_v\mathrm{Inf}_x(v)$ | 均值 0.08265（95% 区间 0.07385–0.09145），中位数 0.05887 |
| 前景份额 | 前景占像素 8.43%，承担影响的 70.83% |
| 集中度 | 影响最大的 10 个像素（1%）占 30.94%，前 51 个（5%）占 65.59% |
| 一阶律，$\varepsilon=10^{-4}$（每图 1000 次） | 实测/预测 0.955，配对差区间含零 |
| 一阶律，$\varepsilon=10^{-3}$（每图 300 次） | 0.992，配对差区间含零 |
| 一阶律，$\varepsilon=10^{-2}$（每图 100 次） | 0.951，配对差区间不含零（每图约 10 个翻转，出现交互） |

未编码 BPSK，每图 100 次公共噪声，平均每像素能量相同：

| $E_b/N_0$ | 均匀功率的关系风险 | 变分分配 | 比值 | 像素误码率：均匀 → 分配 | 图内功率跨度 |
|---|---|---|---|---|---|
| 4 dB | $9.685\times10^{-4}$ | $4.18\times10^{-5}$ | 0.043 | 0.01249 → 0.02913 | 7.1 dB |
| 6 dB | $1.971\times10^{-4}$ | $7.08\times10^{-6}$ | 0.036 | 0.00240 → 0.00548 | 5.2 dB |
| 8 dB | $1.59\times10^{-5}$ | $7.2\times10^{-7}$ | 0.045 | 0.00019 → 0.00043 | 3.7 dB |

这是对读出的受控实验：发端用自己的图像计算影响图并分配功率，链路未编码，接收端逐像素符号判决。它不构成对神经编解码器的排名，也没有计入把分配告知接收端的任何开销（符号判决不需要）。8 dB 处分配方案的风险估计区间为 $[0.30,1.14]\times10^{-6}$，相对误差较大。
