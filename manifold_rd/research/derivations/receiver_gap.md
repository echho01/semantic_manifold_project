# 预测后验与对象重构的率失真差距

给定消息先优化预测，再限制到可实现的对象集合。圆周上的定量界允许两种接收机分别优化编码器。

[研究索引](../README.md) · [符号和引用所用书目](../../ref.bib)

## Posterior reduction and the price of an object-valued receiver

Write $Q=P(Y_t=\cdot\mid U,C)$, suppress $C$ below, and measure KL in nats. Let $J$ be any encoder message. For every decoder $a(J)$, with $m_J=\mathbb{E}[Q\mid J]$, 

$$
\mathbb{E}D_{\mathrm{KL}}(Q\Vert a(J))=\mathbb{E}D_{\mathrm{KL}}(Q\Vert m_J)+\mathbb{E}D_{\mathrm{KL}}(m_J\Vert a(J)).
$$

Subtract the first term and condition on $J$: the remaining sum is $\sum_i m_{J,i}\ln(m_{J,i}/a_i(J))$. This is the Bregman centroid identity [`JMLR:v6:banerjee05b`]. For distortion depending on the source only via $Q$, average any $P_{J|U}$ over $U|Q$; it preserves $(Q,J)$ and costs $I(Q;J)\le I(U;J)$. Conversely one may compute $Q$ at the sender. If outputs must correspond to objects, the second term is at least the directed projection cost to the realizable set. At zero rate this compares the mean posterior with that set.

**Proposition (A finite-rate circle gap).** Let $Q_i(\theta)=1/3+a\cos(\theta-2\pi i/3)$ for $i=0,1,2$ and $0<a<1/3$. Its Euclidean mixture plane has center $c=(1/3,1/3,1/3)$, radius $r=\sqrt{3/2}\,a$, and all components are at least $\eta=1/3-a$. For any distribution on this circle, let $D_{\rm circ}(R)$ and $D_{\rm pred}(R)$ minimize expected KL at the same information budget, with circle and unrestricted predictive outputs respectively. Then 

$$
D_{\rm circ}(R)-D_{\rm pred}(R)\ \ge\ \frac{\eta^2}{4r^2}D_{\rm pred}(R)^2 .
$$

**Proof.** For each message let $m=\mathbb{E}[Q|J]$, $v=\mathbb{E}[\|Q-m\|_2^2|J]$, $d=\mathbb{E}[D_{\mathrm{KL}}(Q\Vert m)|J]$. By expanding around $c$, $v=r^2-\|m-c\|^2$. Thus every circle point has distance from $m$ at least $r-\|m-c\|\ge v/(2r)$. For differences of probability vectors, $\|z\|_1^2\ge2\|z\|_2^2$ by splitting equal positive and negative mass. Pinsker therefore implies $D_{\mathrm{KL}}(m\Vert a(J))\ge v^2/(4r^2)$. Next $\ln t\le t-1$ gives $D_{\mathrm{KL}}(Q\Vert m)\le\sum_i(Q_i-m_i)^2/m_i$, hence $v\ge\eta d$. Average the centroid identity and apply Jensen: the constrained distortion for this encoder is at least $\mathbb{E}d+\eta^2(\mathbb{E}d)^2/(4r^2)$. Its unrestricted receiver is feasible by data processing and achieves $\mathbb{E}d\ge D_{\rm pred}(R)$. Take the infimum over constrained encoders. For distortion in bits multiply the coefficient by $\ln 2$.

This is a one-sided inequality for a specified family. It does not establish a matching order for general posterior manifolds. A learned softmax head yields a measured proxy $Q_0$; the preceding Bayes claims concern the true $Q$. Report calibration, label loss and the head's clean accuracy alongside proxy KL.

## 强凸层集上的推广与反方向的构造

设可实现后验集合 $\mathcal S$ 的每个坐标至少为 $\eta>0$，在其凸包上函数 $\varphi$ 关于混合平面欧氏距离是 $\kappa$ 强凸、$L$ Lipschitz，且在 $\mathcal S$ 上恒为零。对任意消息 $J$，记条件均值为 $m$，参照上文令 $v$ 为条件方差，$d$ 为相对质心的条件 KL（nat）。在凸函数 $\varphi(q)-\kappa\|q\|^2/2$ 上用 Jensen 不等式，

$$
-\varphi(m)\ge\frac{\kappa}{2}\left(\mathbb E[\|Q\|^2\mid J]-\|m\|^2\right)=\frac{\kappa v}{2}.
$$

对任意可实现重构 $b$，$L\|b-m\|\ge|\varphi(b)-\varphi(m)|$。由概率差向量的零和性质与 Pinsker 有 $D_{\mathrm{KL}}(m\Vert b)\ge\|m-b\|_2^2$。又 $v\ge\eta d$，故投影代价至少为 $\gamma d^2$，其中 $\gamma=(\kappa\eta/(2L))^2$。条件平均后再对编码器取下确界，得到主稿中的同码率下界。圆的定义函数取 $\varphi(q)=(\|q-c\|^2-r^2)/2$，在圆盘上可取 $\kappa=1,L=r$。

公共上下文对应不同的 $\gamma_c$ 时，写 $d_c$ 为某编码器在上下文 $c$ 的质心平均失真。Cauchy–Schwarz 给出 $(\mathbb E d_C)^2\le\mathbb E[\gamma_C d_C^2]\mathbb E[1/\gamma_C]$。当右侧第二因子有限，条件率失真使用系数 $(\mathbb E[1/\gamma_C])^{-1}$。

在圆周上，将 $m$ 径向投影到 $b$（$m=c$ 时任取一点），则 $\|b-m\|=v/(r+\|m-c\|)\le v/r$。由对数上界推出 $D_{\mathrm{KL}}(m\Vert b)\le\sum_i(m_i-b_i)^2/b_i\le v^2/(\eta r^2)$。又 $v\le d$，同一编码器因此可达到 $\mathbb E d+\mathbb E[d^2]/(\eta r^2)$；该结论使用消息条件失真的二阶矩。要证明与下界同阶还须控制该二阶矩。

## 均匀圆的零率精确值

均匀角度下质心为 $c=(1,1,1)/3$。由于 $D_{\mathrm{KL},2}(c\Vert q)=-\tfrac13\log_2(27\prod_i q_i)$，取最近的正向 KL 投影等价于最大化三个坐标的乘积。三相余弦的和为零、两两乘积之和为 $-3/4$、三重积为 $\cos(3\theta)/4$。因此 $27\prod_i q_i=1-27a^2/4+27a^3\cos(3\theta)/4$。最大值取在 $\cos(3\theta)=1$，代入即为主稿的零率精确代价。[数值核对](../results/toy_zero_rate.csv)
