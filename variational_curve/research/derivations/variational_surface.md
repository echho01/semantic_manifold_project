> **2026-11-18 更正。** 推论 1 的最优曲面就是与 $g_\eta$ 等距的曲面，零空间方向的 $\sqrt\eta$ 尺度来自评价度量自身的权重，任意 $\eta>0$ 下曲面仍是 $m$ 维。下文中“涌现维数”“收缩到任务秩”的说法已从主稿删除，推论改述为忽略各向异性的码率代价。定理 2 需要流形无边界；由标量曲率给出的 reach 上界只对标量曲率在某处为正的曲面成立。

# 编码曲面的变分泛函：等距亏量、维数、面积与 reach

把编码曲面当作率失真问题的变量，推出两个泛函：数字编码器的等距亏量 $\Lambda[h]$ 和模拟映射的面积–reach–功率界。本文件给出主稿定理 1、命题 1、推论 1、定理 2 及曲率说明的完整证明。有限字母表上的对应结论见 [encoder_identifiability.md](encoder_identifiability.md) 与 [interaction_rank.md](interaction_rank.md)。

[研究索引](../README.md) · [计算脚本](../experiments/surface_functional.py) · [数值](../results/surface_functional.json) · [测试](../tests/test_variational.py)

## 0 Setting

$(\mathcal M,g)$ is a compact smooth manifold of dimension $m$ with semantic metric $g$; the source has a continuous positive density $p$ relative to $V_g$, and $h_g(p)=-\int p\ln p\,dV_g$. A cost $\delta$ is locally Riemannian with metric $g$ when $|\delta-d_g^2|\le c_0d_g^3$ near the diagonal and $\delta\ge c_1>0$ elsewhere. The reference law is $R_g(D)=h_g(p)-\tfrac m2\ln\frac{2\pi eD}m+o(1)$ [`linder1999high`].

A coding surface is a smooth embedding $F:\mathcal M\to\mathbb R^k$ with induced metric $h=J_F^\top J_F$, volume element $dV_h=\varsigma\,dV_g$, $\varsigma=\sqrt{\det(g^{-1}h)}$, and area $V_h=\int\varsigma\,dV_g$. The entropy of the source relative to the surface area is $h_h(p)=h_g(p)+\mathbb E\ln\varsigma=h_g(p)+\tfrac12\mathbb E\ln\det(g^{-1}h)$.

**Digital coder.** $W(da\mid u)=\mu(a)e^{-s\|F(u)-F(a)\|^2}dV_h(a)/Z_s(u)$ with one global $s>0$ and $\mu$ in a bounded subset of $C^1(\mathcal M)$ with $\mu\ge\mu_0>0$. $R_F(D)$ is the least $I(U;A)$ in the family subject to $\mathbb E\delta(U,A)\le D$.

**Analog mapping.** $Y=F(U)+N$, $N\sim\mathcal N(0,\sigma^2I_k)$; the receiver returns the nearest point of $F(\mathcal M)$.

## 1 Isometry defect (Theorem 1)

**Statement.** As $D\downarrow0$, $R_F(D)=h_g(p)-\tfrac m2\ln\frac{2\pi eD}m+\Lambda[h]+o(1)$ with

$$
\Lambda[h]=\frac m2\ln\Big(\frac1m\mathbb E\operatorname{tr}(h^{-1}g)\Big)+\frac12\mathbb E\ln\det(g^{-1}h)\ \ge0 .
$$

$\Lambda=0$ iff $h=c\,g$ $p$-a.e. for a constant $c$; and $a\,g\preceq h\preceq b\,g$ implies $\Lambda\le\tfrac m2\ln(b/a)$.

**Proof.** *Local form.* For an embedding, $\|F(u)-F(a)\|^2=d_h(u,a)^2+O(d_h^4)$ and the chord is bounded below away from the diagonal, so the coding cost is locally Riemannian with metric $h$. In $h$-normal coordinates $v$ about $u$, $dV_h=(1+O(|v|^2))dv$ and $\mu(a)=\mu(u)(1+O(|v|))$ uniformly over the family. Laplace's method gives, uniformly in $u$ and in $\mu$,

$$
Z_s(u)=\mu(u)\Big(\frac\pi s\Big)^{m/2}(1+O(s^{-1/2})),\qquad A\mid U=u\ \approx\ \mathcal N\big(0,\tfrac1{2s}I_m\big)\ \text{in }v .
$$

*Information.* The conditional density of $A$ relative to $V_h$ is $\mu(a)e^{-s\|\cdot\|^2}/Z_s(u)$, hence

$$
h_h(A\mid U=u)=-\mathbb E[\ln\mu(A)\mid u]+s\,\mathbb E[\|F(u)-F(A)\|^2\mid u]+\ln Z_s(u)=\frac m2\ln\frac{\pi e}s+o(1).
$$

The terms in $\ln\mu$ cancel to leading order. The output density relative to $V_h$ is an approximate identity applied to the source density and converges uniformly, so $h_h(A)\to h_h(p)$ and

$$
I(U;A)=h_h(p)-\frac m2\ln\frac{\pi e}s+o(1).
$$

*Distortion.* In the same coordinates $\delta(u,\exp_uv)=v^\top\tilde g(u)v+O(|v|^3)$, where $\tilde g$ is the matrix of $g$ in an $h$-orthonormal frame, whose trace is $\operatorname{tr}(h^{-1}g)$. Therefore $\mathbb E\delta=\frac1{2s}\mathbb E\operatorname{tr}(h^{-1}g)(1+o(1))$; contributions away from the diagonal are exponentially small.

*Elimination.* Both quantities are monotone in $s$ for large $s$ and their leading terms do not depend on $\mu$. With $s=\mathbb E\operatorname{tr}(h^{-1}g)/(2D)\,(1+o(1))$,

$$
R_F(D)=h_h(p)-\frac m2\ln\frac{2\pi eD}{\mathbb E\operatorname{tr}(h^{-1}g)}+o(1)=h_g(p)-\frac m2\ln\frac{2\pi eD}m+\Lambda[h]+o(1).
$$

*Sign.* Let $\zeta_1,\dots,\zeta_m$ be the eigenvalues of $h^{-1}g$. Jensen for the logarithm and the inequality of arithmetic and geometric means give

$$
\ln\Big(\frac1m\mathbb E\sum_i\zeta_i\Big)\ \ge\ \mathbb E\ln\Big(\frac1m\sum_i\zeta_i\Big)\ \ge\ \frac1m\mathbb E\sum_i\ln\zeta_i=-\frac1m\mathbb E\ln\det(g^{-1}h).
$$

Equality in the second step requires equal eigenvalues at each point, and in the first a constant trace, so $h^{-1}g=c^{-1}I$ $p$-a.e. If $\zeta_i\in[1/b,1/a]$ then $\Lambda\le\tfrac m2[\ln(1/a)-\ln(1/b)]$. $\square$

**Remarks.** (i) $\Lambda[\kappa h]=\Lambda[h]$ for a constant $\kappa$: the digital coder is indifferent to the area of the surface. (ii) The second fundamental form of $F$ enters only through the $O(d_h^4)$ difference between chord and geodesic, i.e. at order $D$ in the rate. (iii) The leading term does not depend on $\mu$: a smooth codeword density cannot repair a metric mismatch. An anisotropic or position-dependent scale could, but that is a different coder, equivalent to changing the surface.

**数值核对（720 态圆，$g=d\theta^2$，$h=(1+0.5\cos\theta)^2d\theta^2$）。** 闭式 $\Lambda=\tfrac12\ln(1-b^2)^{-3/2}+\ln\tfrac{1+\sqrt{1-b^2}}2$，$b=0.5$ 时 0.21125 bit。与 Blahut 最优在相同 $D_g$ 处的码率差：

| $D_g$ | 1 | 0.3 | 0.1 | 0.03 | 0.01 | 0.005 | $\Lambda$ |
|---|---|---|---|---|---|---|---|
| 均匀源，均匀质量 | 0.1334 | 0.1862 | 0.2037 | 0.2095 | 0.2108 | 0.2110 | 0.2112 |
| 均匀源，优化质量 | 0.1558 | 0.1978 | 0.2073 | 0.2101 | 0.2108 | 0.2110 | 0.2112 |
| 倾斜源，优化质量 | 0.0560 | 0.1240 | 0.1355 | 0.1391 | 0.1401 | 0.1404 | 0.1406 |
| 倾斜源，均匀质量 | 0.3698 | 0.3061 | 0.2310 | 0.1772 | 0.1544 | 0.1478 | 0.1406 |
| 倾斜源，匹配曲面 $h=g$、均匀质量 | 0.1946 | 0.0842 | 0.0319 | 0.0101 | 0.0034 | 0.0017 | 0 |

最后一行是纯密度效应：低分辨率时显著，高分辨率时消失。前四行是几何效应：高分辨率时趋于 $\Lambda$，优化质量也去不掉。Blahut 参考的对偶间隙最大 0.0019 bit（出现在低斜率处）。

## 2 First variation (Proposition 1)

With $T=\mathbb E\operatorname{tr}(h^{-1}g)$, $\delta(h^{-1})=-h^{-1}\delta h\,h^{-1}$ and $\delta\ln\det h=\operatorname{tr}(h^{-1}\delta h)$:

$$
\delta\Lambda=-\frac m{2T}\mathbb E\operatorname{tr}(h^{-1}\delta h\,h^{-1}g)+\frac12\mathbb E\operatorname{tr}(h^{-1}\delta h)=\frac12\mathbb E\operatorname{tr}\Big[h^{-1}\delta h\Big(I-\frac mTh^{-1}g\Big)\Big].
$$

For a variation $\delta F$ of the surface, $\delta h=J_F^\top\delta J_F+\delta J_F^\top J_F$. As a functional of the metric, $\Lambda$ is stationary for all symmetric $\delta h$ iff $h^{-1}g=(T/m)I$ $p$-a.e., and every such metric attains the global minimum zero by Theorem 1. For surface variations the set of admissible $\delta h$ is smaller (a tangential field $\xi$ and a normal field $\nu$ give $\delta h_{ab}=\nabla_a\xi_b+\nabla_b\xi_a-2\langle\nu,\mathrm{I\!I}_{ab}\rangle$), so a surface may be stationary without being a minimizer; the statement in the letter is therefore made for the metric. The gradient direction is explicit: stretch where $h^{-1}g$ exceeds its mean, contract where it is below. The finite covariance formula $\delta L=\sum_ip_i\operatorname{Cov}_{W_i}(A_{ij},S_{ij})$ is the discrete form used in training.

有限差分核对：随机对称正定场上公式与中心差分相差 $<10^{-8}$。

## 3 Emergent dimension (Corollary 1)

Let $g=g_\eta$ have, in a frame adapted to the task, $r$ sensitive directions on which $h^{-1}g$ has eigenvalues $\zeta_1,\dots,\zeta_r$ and $m-r$ null directions of $g$-weight $\eta$ that the surface scales by $\varepsilon$, so their eigenvalue is $\nu=\eta/\varepsilon^2$. With $\bar\zeta=\tfrac1r\mathbb E\sum_{i\le r}\zeta_i$,

$$
\Lambda(\nu)=\frac m2\ln\frac{r\bar\zeta+(m-r)\nu}m-\frac12\mathbb E\sum_{i\le r}\ln\zeta_i-\frac{m-r}2\ln\nu,\qquad
\Lambda'(\nu)=0\iff\nu=\bar\zeta .
$$

The second derivative is positive there. Hence $\varepsilon^2=\eta/\bar\zeta$, and $\Lambda(\eta)-\Lambda(\bar\zeta)=\frac{m-r}2\ln\frac{\bar\zeta}\eta+\frac m2\ln\frac{r\bar\zeta+(m-r)\eta}{m\bar\zeta}$. As $\eta\downarrow0$, $\bar\zeta$ tends to its task value, $\varepsilon\propto\sqrt\eta\to0$ and the penalty of keeping unit scale is $\frac{m-r}2\ln(1/\eta)+O(1)$.

**数值核对（平环面，后验只依赖 $\theta$，$m=2$，$r=1$，$\eta=0.01$）。** 原始图卡中敏感方向特征值 $\tfrac12g_F(\theta)+\eta$，均值 $\bar\zeta=0.11$，最优 $\varepsilon^2=0.0909$。

| 曲面 | $\Lambda$（bit） | $D=0.03$ | $0.003$ | $0.0005$ nat |
|---|---|---|---|---|
| 原始 $h=I$ | 0.8689 | 0.8075 | 0.8643 | 0.8681 |
| 仅调维数 $h=\mathrm{diag}(1,\varepsilon^2)$ | 0.0136 | 0.0238 | 0.0144 | 0.0137 |
| 完全匹配 $h=g_\eta$ | 0 | 0.0075 | 0.0001 | 0.0000 |

维数一项解释 0.8689 − 0.0136 = 0.8553 bit，与上式一致。

## 4 Area, reach and power (Theorem 2)

**Weak-noise distortion.** Inside the reach the differential of the nearest-point projection at a point of the surface is the orthogonal projector onto the tangent space. The error in a chart is $J_F^+N$ with covariance $\sigma^2h^{-1}$, so $D_{\rm w}[F]=\sigma^2\mathbb E\operatorname{tr}(h^{-1}g)$ to leading order. (椭圆上的 Monte Carlo 核对：与公式之比在 8% 以内，4000 次抽样。)

**Area bound.** With $\varsigma=\sqrt{\det(g^{-1}h)}$, the inequality of arithmetic and geometric means gives $\operatorname{tr}(h^{-1}g)\ge m\det(h^{-1}g)^{1/m}=m\varsigma^{-2/m}$. Put $\alpha=m/(m+2)$ and write $p^\alpha=(p\varsigma^{-2/m})^\alpha\varsigma^{2/(m+2)}$. Hölder with exponents $(m+2)/m$ and $(m+2)/2$:

$$
\int p^\alpha dV_g\le\Big(\int p\varsigma^{-2/m}dV_g\Big)^{\frac m{m+2}}\Big(\int\varsigma\,dV_g\Big)^{\frac2{m+2}},
\qquad\text{so}\qquad
\mathbb E\,\varsigma^{-2/m}\ \ge\ \|p\|_\alpha\,V_h^{-2/m}.
$$

Hence $D_{\rm w}\ge m\sigma^2\|p\|_{m/(m+2)}V_h^{-2/m}$. Equality requires $h=c\,g$ pointwise (conformal) and $\varsigma\propto p^{m/(m+2)}$, i.e. $c\propto p^{2/(m+2)}$. Since the Rényi entropy of order $\alpha<1$ dominates the Shannon entropy, $\ln\|p\|_\alpha=\tfrac2mH_\alpha(p)\ge\tfrac2mh_g(p)$, with equality for the uniform density.

**Packing bound.** Let $F(\mathcal M)\subset B(0,\rho)$ have reach at least $\tau$. For a unit normal $n$ and $t<\tau$ the normal exponential map $(x,tn)\mapsto x+tn$ is injective on the tube of radius $t$ and has Jacobian $\det(I-tA_n)$, where the shape operator satisfies $\|A_n\|\le1/\tau$ [`federer1959curvature`, `gray2004tubes`]. Thus the Jacobian is at least $(1-t/\tau)^m$ and

$$
\operatorname{Vol}_k(T_{\tau/2})\ \ge\ V_h\,\omega_{k-m}(\tau/2)^{k-m}\,2^{-m}=V_h\,\omega_{k-m}\tau^{k-m}2^{-k}.
$$

The tube lies in $B(0,\rho+\tau/2)$, which gives $V_h\le2^k\omega_k(\rho+\tau/2)^k/(\omega_{k-m}\tau^{k-m})$. 核对：平面圆、$\mathbb R^3$ 中球面、$\mathbb R^4$ 中 Clifford 环面和两个椭圆（reach $=b^2/a$）都满足该界。

**Scaling.** The weak-noise formula needs the noise inside the reach, so $\tau\ge c\,\sigma$ for a constant $c$ that depends on $k$ and the tolerated threshold probability. With $\tau=c\sigma\ll\rho$ the two bounds give $D_{\rm w}\gtrsim\mathrm{const}\cdot\sigma^2\,(\sigma^{k-m}/\rho^k)^{2/m}=\mathrm{const}\cdot(\sigma^2/\rho^2)^{k/m}$, the exponent of the Shannon limit for $k$ real channel uses per $m$-dimensional source.

**Contrast with the digital coder.** The digital functional is scale-invariant and is minimized by an isometry whatever the density. The analog functional rewards area and is minimized by a conformal surface whose area density is $p^{m/(m+2)}$; for a uniform source the two optimal shapes coincide.

## 5 Which curvature a surface controls

For a submanifold of Euclidean space the Gauss equation gives, in an $h$-orthonormal frame,

$$
\mathrm{Scal}_h=\sum_{i\ne j}\big[\langle\mathrm{I\!I}_{ii},\mathrm{I\!I}_{jj}\rangle-\|\mathrm{I\!I}_{ij}\|^2\big]\le\sum_{i\ne j}\|\mathrm{I\!I}_{ii}\|\,\|\mathrm{I\!I}_{jj}\|\le\frac{m(m-1)}{\tau^2},
$$

because reach at least $\tau$ implies $\|\mathrm{I\!I}(e,e)\|\le1/\tau$ for unit tangent $e$. So $\tau\le\sqrt{m(m-1)/\max\mathrm{Scal}_h}$ wherever the scalar curvature is positive; the round sphere attains equality. An isometric surface has $\mathrm{Scal}_h=\mathrm{Scal}_g/c$, so the positive curvature of the semantic metric caps the noise margin of the surface at a given scale. The intrinsic term $\mathbb E[\mathrm{Scal}_g]D/(6m)$ of the semantic rate depends on $g$ only.

## 6 与有限模型和已有实验的关系

- 有限 Gibbs 编码器中，代价差可被质量吸收当且仅当它可加可分；两个在对角线上一阶为零的光滑代价之差若可分则恒为零。所以双中心化残差非零是 $h\ne c\,g$ 的有限判据。
- 96 点圆周（主稿表 II）：固定曲面就是评价度量（$\Lambda=0$，$\epsilon=0$），工作点约 1.2 bit 属于低分辨率区，所以密度解释了 0.08179 bit 的节省而几何只有 0.00155 bit。
- SEMASIA 上的固定曲面 $z=(2\sqrt q,\sqrt\eta u)$ 对评价代价等距；有界形变 $F=z+\alpha B\tanh(Az)$ 满足 $(1-\alpha)^2g\preceq h\preceq(1+\alpha)^2g$，亏量至多 $m\ln\frac{1+\alpha}{1-\alpha}$。在定理假设成立的范围内，形变在高分辨率下不可能有益；可能的收益只在低码率或受限质量模型下出现。

## 7 未解决的问题

- 定理 1 是 $D\downarrow0$ 的结论；低分辨率下密度与几何的联合最优没有闭式。
- 定理 2 的两个界分别紧（面积界由共形 companding 曲面取等，packing 界由密堆管取等），二者同时取等的曲面是否存在未讨论。
- 曲面变分 $\delta F$ 下的驻点方程（含 $\mathrm{I\!I}$ 的法向项）未求解。
- 真实特征上的 $\Lambda$ 需要 $h^{-1}g$ 的局部估计，尚未测量。
