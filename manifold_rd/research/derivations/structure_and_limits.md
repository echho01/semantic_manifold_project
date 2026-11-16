# 损失预算、对象代价的几何常数、任务秩律与信道极限

本文件给出主稿新增结论的完整推导：五项损失恒等式、凸性刻画、reach 上界、混合曲率局部常数、任务秩高分辨率律及其推论、链路逆定理与模拟映射的功率分配。数值核对见 [theory_checks.py](../theory_checks.py) 与 [测试](../tests/test_theory.py)。

[研究索引](../README.md) · [同码率分离下界](receiver_gap.md) · [连续流形的 Gibbs 界](geometric_bounds.md) · [计算结果](../results/theory_checks.json)

## 1 Loss budget

Let $Y$ be the task variable, $X$ the observation, $U=E_0(X,C)$, $J$ the encoder message, $V$ the channel output and $B\in\Delta$ the receiver's prediction. Assume the chain $Y-(X,C)-(U,C)-(J,C)-(V,C)-(B,C)$: each variable is conditionally independent of all earlier ones given its predecessor. Write $Q=P_{Y|U,C}$, $m_J=\mathbb E[Q\mid J,C]=P_{Y|J,C}$ and $m_V=\mathbb E[Q\mid V,C]=P_{Y|V,C}$.

**Theorem.** In bits,

$$
\mathbb E[-\log_2B(Y)]=H(Y|X,C)+I(Y;X|U,C)+I(Y;U|J,C)+I(Y;J|V,C)+\mathbb E D_{\mathrm{KL},2}(m_V\Vert B).
$$

**Proof.** Conditional on $(V,C,B)$ the law of $Y$ is $m_V$, since $B$ is independent of $Y$ given $(V,C)$. Therefore $\mathbb E[-\log_2B(Y)\mid V,C,B]=H(m_V)+D_{\mathrm{KL},2}(m_V\Vert B)$ and the left side equals $H(Y|V,C)+\mathbb ED_{\mathrm{KL},2}(m_V\Vert B)$. The chain gives $H(Y|J,V,C)=H(Y|J,C)$, $H(Y|U,J,C)=H(Y|U,C)$, $H(Y|X,U,C)=H(Y|X,C)$. Hence

$$
H(Y|V,C)=\underbrace{H(Y|V,C)-H(Y|J,C)}_{I(Y;J|V,C)}+\underbrace{H(Y|J,C)-H(Y|U,C)}_{I(Y;U|J,C)}+\underbrace{H(Y|U,C)-H(Y|X,C)}_{I(Y;X|U,C)}+H(Y|X,C).
$$

Each bracket is a conditional mutual information because the larger conditioning set leaves the entropy unchanged. Writing a conditional mutual information as the mean divergence between the two conditional laws of $Y$ gives $I(Y;U|J,C)=\mathbb ED_{\mathrm{KL},2}(Q\Vert m_J)$ and $I(Y;J|V,C)=\mathbb ED_{\mathrm{KL},2}(m_J\Vert m_V)$. $\square$

The reducible loss is the chain $Q\to m_J\to m_V\to B$. With a noiseless pipe ($V=J$) and a free decoder ($B=m_J$) only the compression term remains, and minimizing it at information $I(U;J|C)\le R$ is the information bottleneck problem $I(Y;U|C)-\max I(Y;J|C)$. A receiver that holds a table $\tilde b_V$ different from $m_V$ (other head, other knowledge base, stale anchors) pays exactly $\mathbb ED_{\mathrm{KL}}(m_V\Vert\tilde b_V)$; the object constraint is the special case in which the table must take values in $\mathcal S$.

**数值核对。** 64 个观测态、32 个表示、5 位索引、三种接收规则和四种链路共 12 个配置，五项之和与逐项枚举的对数损失相差不超过 $10^{-15}$ bit。另对 20 组随机有限链做了同一检验。

## 2 Convexity characterization

Let $\mathcal S$ be compact inside the open simplex. If $\mathcal S$ is convex, every $m_J\in\mathcal S$ and the two receivers coincide at every rate. If $\mathcal S$ is not convex, choose $a,b\in\mathcal S$ and $\lambda\in(0,1)$ with $m=\lambda a+(1-\lambda)b\notin\mathcal S$ and the two-point source with weights $(\lambda,1-\lambda)$. At zero rate the expected loss is linear in the law of the (independent) action, so a constant action is optimal. The free optimum is $m$. The object receiver pays in addition $\min_{s\in\mathcal S}D_{\mathrm{KL}}(m\Vert s)$, which is positive because $s\mapsto D_{\mathrm{KL}}(m\Vert s)$ is continuous on the compact set and vanishes only at $s=m$. For a finite source both distortion–rate functions are continuous at $R=0$, so the strict gap persists on an interval of rates.

## 3 Reach upper bound

**Lemma (centroid defect).** Let $\mathcal S\subset\mathbb R^N$ be closed with reach $\tau>0$ and let $Q\in\mathcal S$ have mean $m$ and variance $v<\tau^2$. Then $\operatorname{dist}(m,\mathcal S)\le\tau-\sqrt{\tau^2-v}\le v/\tau$.

**Proof.** Let $a$ be a nearest point, $\rho=\|m-a\|$, and $n=(m-a)/\rho$ when $\rho>0$. Because $a$ is nearest, $\rho\le\|m-Q\|$ pointwise and $\rho^2\le v<\tau^2$. For $0<t<\tau$ the point $a+tn$ lies on the segment from $a$ towards $m$ extended inside the reach, and its unique nearest point is $a$ [`federer1959curvature`, Theorem 4.8]. So no point of $\mathcal S$ lies in the open ball $B(a+tn,t)$: $\|Q-a-tn\|^2\ge t^2$, i.e. $\langle Q-a,n\rangle\le\|Q-a\|^2/(2t)$. Let $t\uparrow\tau$ and take expectations, using $\mathbb E\|Q-a\|^2=v+\rho^2$:

$$
\rho\le\frac{v+\rho^2}{2\tau}\quad\Longleftrightarrow\quad\rho^2-2\tau\rho+v\ge0 .
$$

The roots are $\tau\pm\sqrt{\tau^2-v}$. Since $\rho<\tau$ only $\rho\le\tau-\sqrt{\tau^2-v}=v/(\tau+\sqrt{\tau^2-v})$ is possible. $\square$

On a circle of radius $\tau$ the bound holds with equality, because $\|m-c\|^2=\tau^2-v$.

**Theorem (rounding).** Let all coordinates on $\mathcal S$ be at least $\eta$ and $d_J<\tau^2$ almost surely. Pinsker and $\|z\|_2^2\le\|z\|_1^2/2$ for zero-sum $z$ give $v_J\le d_J$. The lemma supplies $b_J\in\mathcal S$ with $\|m_J-b_J\|\le v_J/\tau$, and $D_{\mathrm{KL}}(m\Vert b)\le\chi^2(m\Vert b)\le\|m-b\|^2/\eta$. Hence $D_{\rm obj}\le\mathbb Ed_J+\mathbb E[d_J^2]/(\eta\tau^2)$ at the information of the given encoder.

## 4 Local constant and mixture curvature

Let $q(\xi)$ parametrize an $r$-dimensional $C^3$ submanifold $\mathcal S$ of the open simplex. Put $G_{ab}=\langle\partial_aq,\partial_bq\rangle_s$ with $\langle z,w\rangle_s=\sum_iz_iw_i/s_i$, and let $\mathrm{I\!I}_{ab}=(\partial_a\partial_bq)^\perp$ be the Fisher-normal part of the ordinary second derivative. A change of parameters adds to $\partial_a\partial_bq$ a tangent vector, so $\mathrm{I\!I}$ is a tensor: it is the embedding curvature for the mixture connection.

Consider a law supported in a Fisher ball of radius $\epsilon$, with coordinate mean $\bar\xi$, base point $\bar s=q(\bar\xi)$ and covariance $\Sigma^{ab}$.

1. $m=\mathbb Eq(\xi)=\bar s+\tfrac12\Sigma^{ab}\partial_a\partial_bq(\bar\xi)+O(\epsilon^3)$.
2. $D_{\mathrm{KL}}(q(\xi)\Vert m)=\tfrac12\|q(\xi)-m\|_{\bar s}^2+O(\epsilon^3)$ and $q(\xi)-m=(\xi-\bar\xi)^a\partial_aq+O(\epsilon^2)$, so $d=\tfrac12G_{ab}\Sigma^{ab}+O(\epsilon^3)$.
3. For $b=q(\bar\xi+\zeta)$ with $\zeta=O(\epsilon^2)$: $m-b=\tfrac12\Sigma^{ab}\partial_a\partial_bq-\zeta^a\partial_aq+O(\epsilon^3)$ and $D_{\mathrm{KL}}(m\Vert b)=\tfrac12\|m-b\|_{\bar s}^2+O(\epsilon^6)$. Minimizing the quadratic form over $\zeta$ removes the tangential component:

$$
\min_{b\in\mathcal S}D_{\mathrm{KL}}(m\Vert b)=\tfrac18\big\|\Sigma^{ab}\mathrm{I\!I}_{ab}\big\|_{\bar s}^2+O(\epsilon^5).
$$

Points $b$ at Fisher distance of order one from $\bar s$ cost a positive constant, so the minimum over all of $\mathcal S$ is the local one for small $\epsilon$.

For a curve with Fisher-unit tangent, $\Sigma=\sigma^2$, $d=\sigma^2/2$ and the price is $\tfrac18\sigma^4\kappa_m^2=\tfrac12\kappa_m^2d^2$ with $\kappa_m=\|\mathrm{I\!I}(e,e)\|$. For an isotropic cell $\Sigma^{ab}=(2d/r)G^{ab}$ the price is $d^2\|\mathbf H_m\|^2/(2r^2)$, $\mathbf H_m=G^{ab}\mathrm{I\!I}_{ab}$.

**后验圆。** $q_i(\theta)=1/3+a\cos(\theta-2\pi i/3)$，$q'=-a\sin(\cdot)$，$q''=-(q-c)$。平方 Fisher 速度 $g=\sum_i(q_i')^2/q_i$，混合曲率 $\kappa_m=\|q''-\langle q'',q'\rangle q'/g\|/g$。$a=0.2$ 时 $\kappa_m^2\in[2.4306,18.056]$，均匀角度下均值 7.6073。两点小胞元的比值在半宽 0.05 rad 时与 $\kappa_m^2/2$ 相差不到 0.2%。720 态稠密计算中，等信息间隙除以 $D_{\rm pred}^2$（nat）随码率依次为 5.36、4.42、3.90、3.81，趋近 $\tfrac12\mathbb E\kappa_m^2=3.804$。这一全局一致性是数值证据；严格证明需要说明两种最优编码器的胞元失真渐近相等。

**界的松紧。** 分离定理的系数 $\gamma=\eta^2/(4r_\circ^2)=0.074$，reach 上界的系数 $1/(\eta r_\circ^2)=125$，真实局部常数约 3.8。两条不等式来自 Pinsker 与 $\chi^2$ 界，对一般集合成立，因此常数保守。

## 5 Task-rank law

Assumptions: $q:\mathcal M\to\Delta$ smooth, positive, of constant rank $r$; $\mathcal S=q(\mathcal M)$ a compact embedded $C^3$ submanifold of dimension $r$; $P_Q$ has a continuous positive density $p_{\mathcal S}$ relative to the Fisher–Rao volume $V_F$.

**Local form of KL.** Let $\iota(q)=2\sqrt q$, an isometry from the open simplex with the Fisher–Rao metric onto an open subset of the sphere of radius 2. For $q\in\mathcal S$, $b\in\Delta$ put $x_i=\sqrt{b_i/q_i}-1$. Then $\sum_iq_i(1+x_i)^2=1$ gives $-2\sum_iq_ix_i=\sum_iq_ix_i^2$, and $\|\iota(q)-\iota(b)\|^2=4\sum_iq_ix_i^2$. With $-2\ln(1+x)=-2x+x^2+\varrho(x)$, $|\varrho(x)|\le\tfrac43|x|^3$ for $|x|\le\tfrac12$,

$$
D_{\mathrm{KL}}(q\Vert b)=\tfrac12\|\iota(q)-\iota(b)\|^2+\sum_iq_i\varrho(x_i),\qquad\Big|\sum_iq_i\varrho(x_i)\Big|\le\frac{\|\iota(q)-\iota(b)\|^3}{6\sqrt\eta},
$$

valid when $\|\iota(q)-\iota(b)\|\le\sqrt\eta$ because $|x_i|\le\|\iota(q)-\iota(b)\|/(2\sqrt\eta)$. Globally $D_{\mathrm{KL}}(q\Vert b)\ge-2\ln\sum_i\sqrt{q_ib_i}\ge\tfrac14\|\iota(q)-\iota(b)\|^2$.

**Achievability (object alphabet).** On $\mathcal S\times\mathcal S$ the cost is $\tfrac12\varrho_F^2(1+O(\varrho_F))$. The kernel $e^{-sD_{\mathrm{KL}}(q\Vert a)}dV_F(a)/\tilde Z_s(q)$ has $\tilde Z_s=(2\pi/s)^{r/2}(1+o(1))$, $sD_s\to r/2$, and its output density converges uniformly to $p_{\mathcal S}$. Its information is $h_F-\tfrac r2-\tfrac r2\ln(2\pi/s)+o(1)$; with $s=r/(2D_s)(1+o(1))$ this is $h_F-\tfrac r2\ln(4\pi eD_s/r)+o(1)$.

**Converse (predictive alphabet).** For any kernel into $\Delta$, $I\ge h_F-sD-\ln\sup_{b\in\Delta}Z_s(b)$ with $Z_s(b)=\int_{\mathcal S}e^{-sD_{\mathrm{KL}}(q\Vert b)}dV_F(q)$. Let $\mathcal S'=\iota(\mathcal S)$ with reach $\tau'$, $y=\iota(b)$, $a'$ a nearest point of $\mathcal S'$ and $\rho=\|y-a'\|$.

- If $\rho\ge\epsilon_1$ then $D_{\mathrm{KL}}\ge\epsilon_1^2/4$ on $\mathcal S$ and $Z_s(b)\le e^{-s\epsilon_1^2/4}V_F(\mathcal S)$.
- If $\rho<\epsilon_1<\tau'/2$, the ball argument of Section 3 gives $\langle x-a',y-a'\rangle\le\rho\|x-a'\|^2/(2\tau')$ for $x\in\mathcal S'$, hence $\|x-y\|^2\ge\rho^2+(1-\rho/\tau')\|x-a'\|^2$. On $\|x-y\|\le s^{-1/4}$ the local form gives $D_{\mathrm{KL}}\ge\tfrac12(1-O(s^{-1/4}))\|x-y\|^2$; elsewhere $D_{\mathrm{KL}}\ge\tfrac14s^{-1/2}$. The chordal Laplace integral on $\mathcal S'$ then yields

$$
Z_s(b)\le\Big(\frac{2\pi}{s}\Big)^{r/2}\frac{e^{-s\rho^2/3}}{(1-\rho/\tau')^{r/2}}(1+o(1))+e^{-\sqrt s/4}V_F(\mathcal S).
$$

The maximum over $\rho\in[0,\tau'/2]$ of $-s\rho^2/3-\tfrac r2\ln(1-\rho/\tau')$ is $O(1/s)$. So $\sup_bZ_s(b)\le(2\pi/s)^{r/2}(1+o(1))$ and $s=r/(2D)$ gives $R_{\rm pred}(D)\ge h_F-\tfrac r2\ln(4\pi eD/r)-o(1)$. Since $R_{\rm pred}\le R_{\rm obj}$, both equal $h_F-\tfrac r2\ln(4\pi eD/r)+o(1)$.

**数值核对。** 后验圆 $r=1$，$h_F=\ln2\pi+\tfrac12\mathbb E\ln g=1.02165$ nat。渐近线减去 720 态对象接收机的实际码率，在 0.95、2.08、3.41、4.92 bit 处分别为 0.128、0.032、0.005、0.0006 bit。

**推论。**

- 度量更换：$dV_{g'}=\sqrt{\det(g^{-1}g')}\,dV_g$，故 $h_{g'}(p)=h_g(p)+\tfrac12\mathbb E\ln\det(g^{-1}g')$，高分辨率下两种平方测地代价的码率差趋于该量。对 $g'=G_q+\eta g_0$，$\eta\downarrow0$ 时为 $\tfrac{m-r}2\ln\eta+O(1)$，适用范围是 $D\ll\eta$。
- 共享知识：上下文只改变同一 $(\mathcal M,g)$ 上的源分布时，$R(D)-R(D|C)\to h_g(U)-h_g(U|C)=I(U;C)$。24 态圆例中 $I(U;K)=0.25515$ bit，$\beta=8,16,32,64$ 的拉格朗日差为 0.24816、0.25481、0.25515、0.25514 bit。
- 维数随上下文变化：最小化 $\sum_cp_c[h_c-\tfrac{m_c}2\ln(2\pi eD_c/m_c)]$，约束 $\sum_cp_cD_c=D$，驻点条件 $m_c/(2D_c)=\nu$ 给出 $D_c=Dm_c/\bar m$ 与 $\sum_cp_ch_c-\tfrac{\bar m}2\ln(2\pi eD/\bar m)$。由于 $D_c\le D/p_c\to0$，各项余量一致为 $o(1)$。

## 6 Link converse, outage and exponent

With $\kappa$ complex uses per observation and a memoryless (or ergodic, receiver-CSI) channel of capacity $C_{\rm ch}$ bits per use, data processing and convexity of the RDF give $R_\delta(D|C)\le\kappa C_{\rm ch}$ for every block length, and separation achieves every $D$ with $R_\delta(D|C)<\kappa C_{\rm ch}$. Under quasi-static fading apply the converse conditionally on $h$: $R_\delta(D(h)|C)\le\kappa\log_2(1+|h|^2\mathsf{SNR})$. If the right side is below $R_\delta(D|C)$ then $D(h)>D$, so for Rayleigh fading

$$
\Pr\{D(h)>D\}\ge\Pr\{|h|^2<(2^{R_\delta(D|C)/\kappa}-1)/\mathsf{SNR}\}=1-\exp\!\big[-(2^{R_\delta(D|C)/\kappa}-1)/\mathsf{SNR}\big].
$$

Equating $h_F-\tfrac r2\ln(4\pi eD/r)+o(1)$ with $\kappa\ln(1+\mathsf{SNR})$ gives $D^\star=\tfrac r{4\pi e}e^{2h_F/r}(1+\mathsf{SNR})^{-2\kappa/r}(1+o(1))$. A fade changes the capacity term only; the source RDF is unchanged.

**Certified floor on a finite alphabet.** From $I_2+\beta D\ge L_\beta$ and $I_2\le\kappa C_{\rm ch}$: $D\ge\max_\beta(L_\beta-\kappa C_{\rm ch})/\beta$. For the two-use Gaussian link with gains 0.5, 1, 2 the floors are 0.10905, 0.03316, 0.00396 bit against measured 0.12823, 0.04561, 0.01034 bit.

## 7 Aligned analog allocation

Whitened task coordinates $w=G^{1/2}z$ with eigen-variances $\lambda_i$; real modes $y_i=h_ix_i+n_i$, $n_i\sim\mathcal N(0,\sigma^2)$; aligned map $x_{\pi(i)}=a_i\langle\nu_i,w\rangle$; unbiased estimate $y_{\pi(i)}/(h_{\pi(i)}a_i)$. Distortion $\sigma^2\sum_i(h_{\pi(i)}a_i)^{-2}$, power $\sum_ia_i^2\lambda_i\le P$. Cauchy–Schwarz:

$$
\Big(\sum_ia_i^2\lambda_i\Big)\Big(\sum_i\frac1{h_{\pi(i)}^2a_i^2}\Big)\ge\Big(\sum_i\frac{\sqrt{\lambda_i}}{h_{\pi(i)}}\Big)^2,
$$

with equality iff $a_i^2\propto1/(\sqrt{\lambda_i}h_{\pi(i)})$. The power on mode $\pi(i)$ is $a_i^2\lambda_i\propto\sqrt{\lambda_i}/h_{\pi(i)}$. The rearrangement inequality minimizes $\sum_i\sqrt{\lambda_i}\cdot h_{\pi(i)}^{-1}$ by pairing the largest $\lambda_i$ with the smallest $h^{-1}$. With equal $\lambda_i$ this is the allocation $p_j\propto\sqrt{c_j}$, $c_j=\sigma^2/h_j^2$, of the proposal slides.

The same expression is $\operatorname{tr}(G\,G_{\rm obs}^{-1})$ with $G_{\rm obs}=J_f^\top H^\top HJ_f/\sigma^2$, the Cramér–Rao bound for unbiased estimators of $z$ from $y=HJ_fz+n$. Capacity-achieving water-filling instead gives mode $i$ the power $(\mu-\sigma^2/h_i^2)_+$. The two rules move in opposite directions as $h_i$ decreases, and they answer different questions: the analog map must deliver every task direction in one shot, while the separated scheme maximizes total information and leaves the choice of directions to the source code.
