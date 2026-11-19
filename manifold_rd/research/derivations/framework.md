# 流形建模、语义失真与三层语义率失真

主稿第 II–IV 节的完整推导：语义度量与校准、信道引起的失真、广义流形率失真（含曲率项）、知识库条件率失真、任务导向率失真和联合判据。数值核对见 [framework_checks.py](../framework_checks.py) 与 [测试](../tests/test_theory.py)。接收机字母表、损失预算与信道极限见 [structure_and_limits.md](structure_and_limits.md)。

[研究索引](../README.md) · [对 Word 与开题 PPT 的逐条对照](../notes/06_framework_critique.md) · [计算结果](../results/framework_checks.json)

## 1 Semantic metric, calibration and the distortion a channel induces

**Model.** $U=E_0(X,C)\in\mathcal M_c\subset\mathbb R^d$, $C=(K,t)$. Assumption 1: $\mathcal M_c$ compact, connected, smooth, boundaryless, of dimension $m_c$; $\mu_c$ has a continuous positive density $p_c$ relative to the volume of a Riemannian metric $g_c$.

**Semantic distortion.** $d_g(u,\hat u)^2$. A cost $\delta$ is locally Riemannian with semantic metric $g$ if $|\delta-d_g^2|\le c_0d_g^3$ for $d_g\le\epsilon_0$ and $\delta\ge c_1$ otherwise. Instances: generalized $d_{g_0}^2$; knowledge $d_{g_k}^2$; task $\delta_T=D_{\mathrm{KL}}(q(u)\Vert q(\hat u))$ with local form $\tfrac12v^\top G_qv$, $G_q=J_q^\top\operatorname{diag}(1/q)J_q$ of rank $r$; joint $\delta_\eta=\delta_T+\eta d_{g_k}^2$ with semantic metric $g_\eta=\tfrac12G_q+\eta g_k$.

*Why $\delta_\eta$ is locally Riemannian.* With $\iota(q)=2\sqrt q$, $D_{\mathrm{KL}}(q\Vert q')=\tfrac12\|\iota(q)-\iota(q')\|^2+O(\|\cdot\|^3)$ and $\|\iota q(\exp_uv)-\iota q(u)\|^2=v^\top G_q(u)v+O(|v|^3)$, because the Fisher metric is the pullback of the Euclidean one under $\iota$. A smooth metric satisfies $d_{g_\eta}(u,\exp_uv)^2=v^\top g_\eta(u)v+O(|v|^3)$. Away from the diagonal $\delta_\eta\ge\eta d_{g_k}^2$.

**Calibration.** (i) $a\,g\preceq g'\preceq b\,g$ gives $a\,d_g^2\le d_{g'}^2\le b\,d_g^2$, since the inequalities hold for the length of every curve. (ii) $D_{\mathrm{KL}}(q\Vert q')\le\chi^2(q\Vert q')\le\|q-q'\|^2/q_{\min}\le\ell_q^2d_g^2/q_{\min}$.

**Channel-induced distortion.** Let $v=\exp_z^{-1}(\hat z)$ with conditional mean $b_H$ and covariance $\Sigma_H$.

1. Inside the injectivity radius $d_g(z,\exp_zv)^2=v^\top Gv$ exactly, hence $\mathbb E[d_g^2\mid z,H]=b_H^\top Gb_H+\operatorname{tr}(G\Sigma_H)$. This is the formula on slide 19 with equality in normal coordinates. In another chart the identity holds up to $O(\mathbb E|v|^3)$.
2. For $y=Hf(z)+n$, $n\sim\mathcal N(0,\Sigma_n)$, the Fisher information is $G_{\rm obs}=J_f^\top H^\top\Sigma_n^{-1}HJ_f$. A locally unbiased receiver has $\Sigma_H\succeq G_{\rm obs}^{-1}$ to leading order, so the distortion is at least $\operatorname{tr}(GG_{\rm obs}^{-1})=\sum_i1/\lambda_i$, $\lambda_i$ the eigenvalues of $G^{-1}G_{\rm obs}$. A zero eigenvalue removes local identifiability of one direction; the Bayesian receiver then keeps the prior spread along it. The source RDF is unaffected.
3. Ambient noise $N\sim\mathcal N(0,\sigma^2I_d)$ followed by nearest-point projection: the differential of the projection at a point of $\mathcal M$ is the tangent projector, so $\mathbb Ed_g^2=m\sigma^2+O(\sigma^4)$. On a circle of radius one the next term is $\sigma^4$: with $x=n_2/(1+n_1)$ and $\arctan^2x=x^2-\tfrac23x^4+\dots$, $\mathbb E\theta^2=\sigma^2+3\sigma^4-2\sigma^4$.

**数值核对。** 圆映射 $\rho(\cos\theta,\sin\theta)$，$G_{\rm obs}=\rho^2$。$\rho=4,8,16$ 时精确角度均方误差 0.06743、0.015880、0.0039217，与 $1/\rho^2$ 之比 1.079、1.016、1.004，二阶式 $v+v^2$ 给出 0.06641、0.015869、0.0039215。$\rho=2$ 时比值 1.47，局部分析失效。任务加权损失 $\tfrac12\mathbb E[g]/\rho^2$ 在 $\rho=16$ 为 $3.906\times10^{-4}$ nat，精确值 $3.918\times10^{-4}$。

## 2 Generalized manifold RDF

Let $Z_s(a)=\int e^{-s\delta(u,a)}dV_g(u)$ and $\tilde Z_s(u)=\int e^{-s\delta(u,a)}dV_g(a)$.

**(i) Brackets at every $D$.** Relative entropy of $P_{U|A=a}$ against $e^{-s\delta(\cdot,a)}/Z_s(a)$ gives $h_g(U\mid A=a)\le s\mathbb E[\delta\mid A=a]+\ln Z_s(a)$ and $R(D)\ge h_g(p)-sD-\ln\sup_aZ_s(a)$. The kernel $e^{-s\delta(u,a)}dV_g(a)/\tilde Z_s(u)$ attains $D_s=-\mathbb E\partial_s\ln\tilde Z_s(U)$ with information $h_g(A_s)-\mathbb E\ln\tilde Z_s(U)-sD_s$.

**(ii) Leading law.** For a locally Riemannian cost, $Z_s,\tilde Z_s=(\pi/s)^{m/2}(1+o(1))$ uniformly, $sD_s\to m/2$, and the output density converges uniformly to $p$. Both bounds give $R(D)=h_g(p)-\tfrac m2\ln(2\pi eD/m)+o(1)$.

**(iii) Curvature term.** For $\delta=d_g^2$ and $p\in C^4$. Kernel expansion: with $T_tf(x)=\int e^{-t\delta(x,y)}f(y)dV_g(y)$,

$$
T_tf=\Big(\frac\pi t\Big)^{m/2}\Big[f+\frac1t\Big(\frac{\Delta_gf}4-\frac{\mathrm{Scal}}{12}f\Big)+O(t^{-2})\Big].
$$

Normal coordinates give $dV=(1-\tfrac16\mathrm{Ric}(v,v)+O(|v|^3))dv$ and $f(\exp_xv)=f+df(v)+\tfrac12\mathrm{Hess}f(v,v)+\dots$; the Gaussian second moments are $\delta_{ij}/(2t)$.

*Converse.* For $a_t>0$ with $\int p\,a_t\,e^{-t\delta(\cdot,y)}dV\le1$ for all $y$, $I(U;A)\ge\mathbb E\ln a_t(U)-tD$. Choose $a_t=(t/\pi)^{m/2}p^{-1}\exp[\mathrm{Scal}/(12t)-c/t^2]$. The constraint integral is $e^{-c/t^2}(t/\pi)^{m/2}T_t[e^{\mathrm{Scal}/(12t)}]=e^{-c/t^2}(1+O(t^{-2}))$: the term $+\mathrm{Scal}/(12t)$ of the function cancels the term $-\mathrm{Scal}/(12t)$ of the kernel. Hence

$$
R(D)\ge h_g(p)+\frac m2\ln\frac t\pi+\frac{\overline{\mathrm{Scal}}}{12t}-tD-\frac c{t^2},
$$

and $t=m/(2D)$ gives $h_g(p)-\tfrac m2\ln\frac{2\pi eD}m+\frac{\overline{\mathrm{Scal}}}{6m}D-O(D^2)$.

*Achievability.* Source-weighted kernel $p(y)e^{-t\delta}/T_tp(x)$. Since $\int\Delta_gp\,dV=0$, $D_t=\frac m{2t}-\frac{\overline{\mathrm{Scal}}}{12t^2}+O(t^{-3})$; the output density is $p(1-\Delta_gp/(4tp)+O(t^{-2}))$; and $I_t=h_g(p)+\tfrac m2\ln\frac t\pi+\frac{\overline{\mathrm{Scal}}}{12t}-tD_t+O(t^{-2})$. Eliminating $t$ gives the same expansion from above.

**(iv) Homogeneous spaces.** If isometries act transitively and $p$ is uniform, $Z_s=\tilde Z_s$ is constant and the output of the Gibbs kernel is uniform, so both bounds of (i) equal $\ln V_g-\ln Z_s-sD_s$.

**数值核对（精确求积）。** 归一化余量 $(R-R_{\rm lead})/D$：$S^2$ 在 $D\approx10^{-1},10^{-2},10^{-3}$ 为 0.16772、0.16679、0.16668（理论 1/6）；$S^3$ 为 0.33515、0.33351、0.33335（理论 1/3）；$S^1$、$T^2$ 为 0。斜率 0.500、1.000、1.500、1.000。

**度量更换。** $a\,g\preceq g'\preceq b\,g$ 时对一切 $D$ 有 $R_g(D/a)\le R_{g'}(D)\le R_g(D/b)$；$h_{g'}(p)=h_g(p)+\tfrac12\mathbb E\ln\det(g^{-1}g')$，故高分辨率下差值趋于该量。

**源侧"语义注水"的含义。** 最优测试信道在 $g$ 正交坐标中各向同性，图卡坐标下误差协方差为 $(D/m)G^{-1}$。度量判为不敏感的方向自动得到更大误差，这一分配完全由度量决定，不涉及信道增益。

## 3 Knowledge-conditioned RDF

$K$ finite, known to both terminals, selecting $(\mathcal M_k,g_k,\mu_k)$.

**知识库的两部分。** 知识库记作 $\mathcal K=(\theta,K)$。参数知识 $\theta$（预训练参考模型的权重、预测头、锚点字典）对每个样本相同，$I(U;\theta)=0$，它的作用是定义 $(\mathcal M,g,\mu)$ 本身；下面各条对固定的 $\theta$ 成立。情境知识 $K$ 是随样本变化、双方已知的随机变量，即本节的 $K$。两个 $\theta$ 对应两个信源，只能在两者之外定义的判据（如标签损失）下比较。

1. $R_K(D)=\min_{\sum p_kD_k\le D}\sum_kp_kR_k(D_k)$ (condition a feasible kernel; combine near-optimal kernels). Inserting the leading law and solving $m_k/(2D_k)=\nu$: $D_k=Dm_k/\bar m$ and $R_K(D)=\sum_kp_kh_{g_k}(p_k)-\tfrac{\bar m}2\ln\frac{2\pi eD}{\bar m}+o(1)$.
2. Common support, baseline $(\mathcal M,g_0,\bar\mu)$: the gain tends to $h_{g_0}(U)-\sum_kp_kh_{g_k}(p_k)=I(U;K)-\tfrac12\mathbb E\ln\det(g_0^{-1}g_K)$.
3. $g_k=g_0$: $0\le\Gamma_K(D)\le I(U;K)$ for every $D$ (a kernel ignoring $K$ has $I(U;A|K)=I(U;A)-I(A;K)$; conversely $I(U;A)\le I(U;K)+I(U;A|K)$). $g_k\preceq g_0$: $\Gamma_K\ge0$.
4. Nested knowledge with a criterion that depends on $K_1=\varphi(K_2)$ only: $R_{K_2}\le R_{K_1}\le R_{K_2}+I(U;K_2|K_1)$.
5. Mismatch $a\,g_{\rm t}\preceq g_{\rm r}\preceq b\,g_{\rm t}$: excess at most $R_{g_{\rm t}}(D/b)-R_{g_{\rm t}}(D/a)\to\tfrac m2\ln(b/a)$.
6. Strata of different dimension: a transmitter without $K$ pays at most $H(K)$ more, because the stratum can be described.

**数值核对（圆，720 态，$K=\pm1$ 等概率）。** 预测值与 $D=0.01$ 处的计算值（bit）：

| 知识改变的对象 | $I(U;K)$ | $\mathbb E\log_2w_K$ | 预测增益 | 计算增益 |
|---|---|---|---|---|
| 仅源分布 | 0.25515 | 0 | 0.25515 | 0.25515 |
| 仅度量（均匀源，$w=1+0.5\cos\theta$） | 0 | −0.10003 | 0.10003 | 0.10003 |
| 在高概率处变粗（$w_K=1-0.5K\cos\theta$） | 0.25515 | −0.40929 | 0.66444 | 0.66444 |
| 在高概率处变细（$w_K=1.25(1+0.5K\cos\theta)$） | 0.25515 | 0.53115 | −0.27600 | −0.27600 |

最后一行是"知识越多码率越低"的反例：知识携带 0.255 bit 的源信息，却因把度量变细而多付 0.276 bit。

**真实特征（固定判据，第 3、4 条）。** 两个骨干、8,192 个测试表示、256 个锚点，$K$ 取粗类或细类标签。无知识码率为 4 bit 时，粗类节省 [2.10, 2.16] bit（平方弦长）和 [2.73, 2.78] bit（后验 KL），均满足 $0\le\Gamma\le H(K)=4.32$ bit 和嵌套界。度量项未检验。见 [EXPERIMENTS.md](../EXPERIMENTS.md)。

## 4 Task-oriented RDF and the joint criterion

1. Sufficiency: the source may be replaced by $Q=q(U)$ on $\mathcal S=q(\mathcal M)$.
2. Task-rank law (proof in [structure_and_limits.md](structure_and_limits.md), Section 5): $R_T(D)=h_F-\tfrac r2\ln\frac{4\pi eD}r+o(1)$ for object and predictive receivers. It is the leading law for the semantic metric $\tfrac12G_q$ on $\mathcal S$.
3. Calibrated ordering: $\delta_T\le L_T^2d_{g_0}^2$ gives $R_T(D)\le R_0(D/L_T^2)$ for every $D$ (object receiver), and the difference is $\tfrac{m-r}2\ln(1/D)+O(1)$.
4. Joint cost: $\max\{R_T(D),R_{g_k}(D/\eta)\}\le R_\eta(D)\le R_{g_k}(D/(L_T^2+\eta))$ for every $D$; for fixed $\eta$, $R_\eta(D)=h_{g_\eta}(p)-\tfrac m2\ln\frac{2\pi eD}m+o(1)$ with $h_{g_\eta}=h_{g_k}+\tfrac{m-r}2\ln\eta+\tfrac12\mathbb E\sum_{i\le r}\ln(\eta+\lambda_i/2)$.

**数值核对（平环面，后验只依赖 $\theta$，$m=2$，$r=1$，$\eta=0.01$，$L_T^2=0.1286$）。**

| $D$（nat） | 任务 | 任务渐近线 | 联合 | 联合渐近线 | 校准广义 |
|---|---|---|---|---|---|
| $10^{-2}$ | 2.2224 | 2.2488 | 3.9148 | 3.9249 | 6.0013 |
| $10^{-3}$ | 3.9072 | 3.9097 | 7.2458 | 7.2468 | 9.3233 |
| $3\times10^{-4}$ | 4.7775 | 4.7782 | 8.9835 | 8.9838 | 11.060 |

$D$ 每缩小十倍，任务码率增加 1.68 bit，联合与广义增加 3.33、3.32 bit，对应 $\tfrac r2\log_210=1.66$ 与 $\tfrac m2\log_210=3.32$。
