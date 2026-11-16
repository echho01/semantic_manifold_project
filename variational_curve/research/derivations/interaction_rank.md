# 交互残差、到 Blahut 最优的距离与交互秩

固定评价矩阵后，曲面变化只能通过编码 logit 中无法被行、列加法项吸收的部分影响率失真。本文件给出该结论的三个定量形式：相对 Blahut 最优的精确超额恒等式与二次上界，点积型编码头能否实现最优信道的秩判据，以及单一欧氏曲面的负型条件。

[研究索引](../README.md) · [质量与几何的可识别性](encoder_identifiability.md) · [数值实现](../experiments/interaction_rank.py) · [随机实例检验](../tests/test_interaction.py) · [圆周结果](../results/interaction_rank.json)

## 记号

源权重 $p_i>0$（$i\le n$），$M$ 个锚点，评价矩阵 $C^g\in\mathbb R_{\ge0}^{n\times M}$ 固定。编码器族为

$$
W_{ij}(C,m,s)=\frac{m_je^{-sC_{ij}}}{\sum_k m_ke^{-sC_{ik}}},\qquad s>0,\ m_j>0 .
$$

$L(W)=I_2(W)+\lambda D_g(W)$ 以 bit 计，$F(\lambda)=\min_W L(W)$ 为有限 Blahut 值，最优信道写作 $W^\star_{ij}=\pi^\star_j2^{-\lambda C^g_{ij}}/Z_i$，其中 $\pi^\star$ 是 $W^\star$ 自身的输出边缘，$Z_i=\sum_j\pi^\star_j2^{-\lambda C^g_{ij}}$。以 nat 计的斜率为 $\beta=\lambda\ln2$，故 $W^\star=W(C^g,\pi^\star,\beta)$。可加分离矩阵的线性空间记为 $\mathcal S_{\rm sep}=\{a_i+b_j\}$，$H_n=I_n-\mathbf 1\mathbf 1^\top/n$。双中心化 $T\mapsto H_nTH_M$ 是等权 Frobenius 内积下到 $\mathcal S_{\rm sep}^\perp$ 的正交投影，其核恰为 $\mathcal S_{\rm sep}$。

交互上确界范数定义为

$$
\epsilon(\Delta)=\min_{a,b}\max_{ij}|\Delta_{ij}-a_i-b_j| .
$$

行、列均值是可行的 $(a,b)$，所以 $\epsilon(\Delta)\le\max_{ij}|(H_n\Delta H_M)_{ij}|$。精确值是一个线性规划，实现见 `interaction_sup`。

## 一　相对 Blahut 最优的超额恒等式

**Proposition (Excess identity).** Let $V$ be any channel supported on the support of $\pi^\star$, with output marginal $\bar V$. Then

$$
L(V)-F(\lambda)=\sum_ip_i D_{\mathrm{KL},2}(V_i\Vert W^\star_i)-D_{\mathrm{KL},2}(\bar V\Vert\pi^\star).
$$

**Proof.** 对任意辅助边缘 $\pi$，$\sum_ip_iD_{\mathrm{KL},2}(V_i\Vert\pi)=I_2(V)+D_{\mathrm{KL},2}(\bar V\Vert\pi)$。取 $\pi=\pi^\star$。逐行有

$$
D_{\mathrm{KL},2}(V_i\Vert\pi^\star)+\lambda\sum_jV_{ij}C^g_{ij}
=\sum_jV_{ij}\log_2\frac{V_{ij}}{\pi^\star_j2^{-\lambda C^g_{ij}}}
=D_{\mathrm{KL},2}(V_i\Vert W^\star_i)-\log_2Z_i .
$$

把 $V=W^\star$ 代入同一式并利用 $\bar W^\star=\pi^\star$ 得 $F(\lambda)=-\sum_ip_i\log_2Z_i$。对 $i$ 加权求和并移项即得结论。支撑条件保证右侧两项有限。$\blacksquare$

右侧第一项按编码器自身的行概率加权。第二项非负，所以 $L(V)-F(\lambda)\le\sum_ip_iD_{\mathrm{KL},2}(V_i\Vert W^\star_i)$ 对任意 $V$ 成立（右侧可为无穷）。这是 Blahut 交替最小化的几何恒等式 [`csiszar1974`]，此处用于比较受限编码器与最优信道。

## 二　交互残差的二次上界

**Proposition (Quadratic bound).** Suppose the coding cost and scale satisfy

$$
sC^h-\beta C^g=a_i+b_j+E_{ij},\qquad |E_{ij}|\le\epsilon .
$$

Then the masses $m_j\propto\pi^\star_je^{b_j}$ give $L(W(C^h,m,s))-F(\lambda)\le\epsilon^2/(2\ln2)$ bits. Consequently $\inf_mL(W(C^h,m,s))-F(\lambda)\le\epsilon(sC^h-\beta C^g)^2/(2\ln2)$.

**Proof.** 记 $V=W(C^h,m,s)$。代入得

$$
V_{ij}=\frac{\pi^\star_je^{b_j}e^{-\beta C^g_{ij}-a_i-b_j-E_{ij}}}{\sum_k\pi^\star_ke^{b_k}e^{-\beta C^g_{ik}-a_i-b_k-E_{ik}}}
=\frac{W^\star_{ij}e^{-E_{ij}}}{\sum_kW^\star_{ik}e^{-E_{ik}}}.
$$

即 $V_i$ 是 $W^\star_i$ 以函数 $g=-E_{i\cdot}$ 作指数倾斜的结果。令 $P_\theta\propto W^\star_ie^{\theta g}$，$\Lambda(\theta)=\ln\mathbb E_{W^\star_i}e^{\theta g}$。则 $D_{\rm KL}(P_\theta\Vert P_0)=\theta\Lambda'(\theta)-\Lambda(\theta)=\int_0^\theta u\Lambda''(u)\,du$，且 $\Lambda''(u)=\operatorname{Var}_{P_u}(g)\le\epsilon^2$，因为取值范围不超过 $2\epsilon$ 的随机变量方差至多 $\epsilon^2$。取 $\theta=1$ 得 $D_{\rm KL}(V_i\Vert W^\star_i)\le\epsilon^2/2$ nat。代入超额恒等式并丢弃非负的边缘项。当 $\pi^\star$ 有零分量时，取正质量序列逼近 $m$，由连续性得到下确界形式。$\blacksquare$

三点说明。其一，$C^h=C^g$、$s=\beta$ 时 $\epsilon=0$，仅学质量即达到 Blahut 最优。其二，任何曲面都不低于 $F(\lambda)$，所以在同一曲面 $C^h$ 上“改曲面相对只学质量”的 Lagrangian 收益不超过 $\min_s\epsilon(sC^h-\beta C^g)^2/(2\ln2)$。其三，上确界范数是最坏情况。恒等式第一项才是精确量，当大残差落在 $W^\star$ 概率很小的配对上时二者可差几个数量级，见第六节数值。

### 两个任意代价之间的比值型比较

上面的结论以最优信道为参照。若只想比较两个都不最优的代价，可用下述较粗的结论。设 $C'=C+a_i+b_j+E_{ij}$，$|E_{ij}|\le\epsilon$，$\rho=2s\epsilon$。对任意正质量 $m$ 取 $m'_j\propto m_je^{sb_j}$，则 $W'=W(C',m',s)$ 与 $W=W(C,m,s)$ 满足 $W'_{ij}=W_{ij}e^{-sE_{ij}}/\mathbb E_{W_i}e^{-sE_{i\cdot}}$，故

$$
e^{-\rho}\le W'_{ij}/W_{ij}\le e^{\rho},\qquad D_g(W')\le e^{\rho}D_g(W),
$$

且以 nat 计

$$
I(W')\le I(W)+(e^{\rho}-1)\bigl(I(W)+2/e\bigr)+\min\{\rho,(e^{\rho}-1)^2\}.
$$

证明：$I(W')\le\sum_ip_iD_{\rm KL}(W'_i\Vert\bar W)$。拆成 $D_{\rm KL}(W'_i\Vert W_i)+\mathbb E_{W'_i}\iota_i$，$\iota_{ij}=\ln(W_{ij}/\bar W_j)$。令 $t=W'/W$，由 $t\ln t\le(t-1)+(t-1)^2$ 和 $\mathbb E_{W_i}t=1$ 得第一项不超过 $(e^\rho-1)^2$，同时不超过 $\max\ln t\le\rho$。第二项等于 $\mathbb E_{W_i}\iota+\mathbb E_{W_i}[(t-1)\iota]\le D_{\rm KL}(W_i\Vert\bar W)+(e^\rho-1)\mathbb E_{W_i}|\iota|$。又 $\mathbb E_{W_i}\iota^-=\sum_{j:W<\bar W}\bar W_j\,(r\ln(1/r))_{r=W_{ij}/\bar W_j}\le1/e$，故 $\mathbb E_{W_i}|\iota|\le D_{\rm KL}(W_i\Vert\bar W)+2/e$。合并即得。于是 $L(W')\le e^\rho L(W)+[(e^\rho-1)\,2/e+\min\{\rho,(e^\rho-1)^2\}]/\ln2$。该式对称地适用于两个方向，说明只学质量的最优值在交互残差下是连续的。

## 三　交互秩

点积型编码头的 logit 为 $\langle x_i,y_j\rangle+\beta_j$，$x_i,y_j\in\mathbb R^d$。由 $\|x_i-y_j\|^2=\|x_i\|^2+\|y_j\|^2-2\langle x_i,y_j\rangle$，它与双线性代价 $C^h_{ij}=\|x_i-y_j\|^2$ 的 Gibbs 族相同，源状态和锚点可以使用不同的映射。令 $T=H_nC^gH_M$，$r_{\rm int}=\operatorname{rank}T$。

**Theorem (Interaction rank).** Let $W^\star$ attain $F(\lambda)$ with $\pi^\star>0$. There exist $x_i,y_j\in\mathbb R^d$, $s>0$ and positive masses with $W(C^h,m,s)=W^\star$ if and only if $r_{\rm int}\le d$. When $\pi^\star$ has zero entries the statement holds with $C^g$ restricted to the supported columns and masses on that support.

**Proof.** $W^\star=W(C^g,\pi^\star,\beta)$。对缩放后的代价 $sC^h$ 与 $\beta C^g$（取公共尺度 1）应用吸收命题：存在 $m$ 使 $W(C^h,m,s)=W^\star$ 当且仅当 $sC^h-\beta C^g\in\mathcal S_{\rm sep}$，即 $H_n(sC^h)H_M=\beta T$。范数项属于 $\mathcal S_{\rm sep}$，故左侧为 $-2sH_nXY^\top H_M$，秩至多为 $d$，必要性得证。反之若 $\operatorname{rank}T\le d$，写 $T=UV^\top$，$U\in\mathbb R^{n\times d}$，$V\in\mathbb R^{M\times d}$，取 $X=-(\beta/(2s))U$，$Y=V$。因 $T$ 已双中心化，$H_nXY^\top H_M=-(\beta/(2s))T$，条件成立，质量取 $m_j\propto\pi^\star_je^{b_j}$，其中 $b_j$ 为可分离部分的列项。零质量情形把吸收命题限制到共同支撑列。$\blacksquare$

推论与说明：

- 嵌入只依赖 $\beta/s$。令 $s\propto\lambda$，同一组 $(X,Y)$ 对所有斜率有效。
- 族内 $\ln W_{ij}=\ln m_j-s\|x_i\|^2-s\|y_j\|^2+2s\langle x_i,y_j\rangle-\ln Z_i$ 的矩阵秩至多为 $d+2$。这与语言模型的 softmax 瓶颈 [`yang2018breaking`] 属于同一类秩障碍。这里的对象是率失真最优测试信道，判据只依赖评价矩阵。
- 若 $C^g_{ij}=\|u_i-v_j\|^2$，$u_i,v_j\in\mathbb R^{d_0}$，则 $T=-2H_nUV^\top H_M$，$r_{\rm int}\le d_0$。
- 最优边缘可能不唯一；定理对每一个最优对 $(W^\star,\pi^\star)$ 分别成立。

## 四　单一欧氏曲面与负型条件

**Corollary (Single surface).** Let states and anchors be the same $n$ points and share one map, $C^h_{ij}=\|F_i-F_j\|^2$. Let $C^g$ be symmetric with zero diagonal and $\pi^\star>0$. A surface in $\mathbb R^d$ with some $s>0$ and positive masses realizes $W^\star$ if and only if $-H_nC^gH_n\succeq0$ with rank at most $d$.

**Proof.** 必要性：$sC^h-\beta C^g=a_i+b_j$ 的左侧对称且对角为零。对角给出 $a_i+b_i=0$，对称给出 $a_i-b_i=a_j-b_j$，于是 $a_i\equiv c$、$b_i\equiv-c$，分离项恒为零。因此 $(\beta/s)C^g$ 必须是 $\mathbb R^d$ 中点集的平方欧氏距离矩阵。由 Schoenberg 定理 [`schoenberg1935`]，对称零对角矩阵 $D$ 是 $\mathbb R^d$ 中的平方距离矩阵当且仅当 $-\tfrac12H_nDH_n\succeq0$ 且秩至多为 $d$。充分性：取实现 $(\beta/s)C^g$ 的点集为 $F_i$，质量取 $\pi^\star$。$\blacksquare$

圆周例。连续圆上 $\delta^2=\pi^2/3+4\sum_{k\ge1}(-1)^k\cos(k\delta)/k^2$，$|\delta|\le\pi$。系数符号交替，$-H C^gH$ 既有正特征值也有负特征值，平方弧长不是负型。所以任何维数的单一欧氏曲面都无法精确实现圆周评价下的最优编码器。弧长本身满足 $|\delta|=\pi/2-(4/\pi)\sum_{k\,\mathrm{odd}}\cos(k\delta)/k^2$，为负型，弦长平方 $2-2\cos\delta$ 亦然。离散 96 点网格上的特征值核对见测试 `test_arc_is_not_of_negative_type_but_chord_is`。该推论只排除精确最优，匹配失真下的差距大小须数值计算。

## 五　谱曲面

取 $T$ 的前 $d$ 个奇异三元组得 $T_d$，令 $C_d=C^g-(T-T_d)$。$C_d$ 与秩 $d$ 的双线性代价只差一个可分离矩阵，因此 $W(C_d,m,\beta)$ 可由 $\mathbb R^d$ 中的点积编码头实现。残差 $T-T_d$ 已双中心化，其最大元不超过谱范数 $\sigma_{d+1}(T)$（Eckart–Young [`eckart1936`]）。由第二节，

$$
\inf_mL(W(C_d,m,\beta))-F(\lambda)\le\frac{(\beta\,\epsilon(T-T_d))^2}{2\ln2}\le\frac{(\beta\,\sigma_{d+1}(T))^2}{2\ln2}\ \text{bit}.
$$

谱截断最小化未加权的 Frobenius 残差。超额恒等式中起作用的是按 $p_iW^\star_{ij}$ 加权的散度，所以截断维数增大时匹配失真码率不必单调下降。

## 六　96 点圆周的数值核对（2026-11-16）

设置与主稿相同：$\theta_i=2\pi(i+1/2)/96$，$p_i=(1+0.85\cos\theta_i)/96$，锚点为同一组点，评价为平方最短弧。斜率网格为 `geomspace(0.25, 32, 25)`（nat），匹配失真值在 $D_g=0.30$ 处对各族自身的 $(D_g,I_2)$ 点线性插值，与早期实验约定一致。质量和尺度用 L-BFGS-B 最小化真实 $I+\beta D_g$，梯度即主稿的协方差一阶变分。

| 编码曲面 | $D_g=0.30$ 处码率 (bit) | 网格斜率 $\beta=1.5422$ 处超额 (bit) | 加权散度上界 (bit) | 上确界范数界 (bit) | $\max|T-T_d|$ | LP 精确 $\epsilon(T-T_d)$ |
|---|---|---|---|---|---|---|
| 有限 Blahut | 1.183191 | 0 | 0 | 0 | — | — |
| 弧长代价，均匀质量 | 1.264985 | — | — | — | — | — |
| 弧长代价，学习质量 | 1.183191 | — | — | — | — | — |
| 单位圆弦长，均匀质量 | 1.272414 | — | — | — | — | — |
| 单位圆弦长，学习质量 | 1.189914 | — | — | — | — | — |
| 谱曲面 $d=2$ | 1.189897 | 0.103462 | 0.106071 | 11.3989 | 2.5776 | 1.7759 |
| 谱曲面 $d=4$ | 1.210406 | 0.070200 | 0.071246 | 4.2622 | 1.5762 | 1.0723 |
| 谱曲面 $d=8$ | 1.187964 | 0.008596 | 0.008606 | 1.3252 | 0.8789 | 0.5944 |
| 谱曲面 $d=16$ | 1.183713 | 0.000743 | 0.000743 | 0.3596 | 0.4578 | 0.3106 |
| 谱曲面 $d=32$ | 1.183265 | 0.000062 | 0.000062 | 0.0818 | 0.2183 | 0.1544 |

“超额”一列使用匹配尺度 $s=\beta$ 与中心化给出的质量 $m_j\propto\pi^\star_je^{b_j}$，对应第二节的构造；“码率”一列另行优化质量与尺度并在匹配失真处比较。数值恒等式误差低于 $2\times10^{-9}$ bit。其余观察：

- $T$ 的数值秩为 95，奇异值成对出现：192.07、48.07、21.40、12.07、7.75、5.40……。$-T$ 的最小特征值为 $-48.07$，弦长平方对应矩阵的最小特征值在舍入误差内为零。
- 该斜率处 $\pi^\star$ 的 96 个分量全部为正，$D_g=0.3242$；Blahut 对偶间隙不超过 $4.7\times10^{-5}$ nat。
- 学习质量的单位圆弦长曲面比 Blahut 高 0.0067 bit，与 $d=2$ 谱曲面（同为一次谐波）一致到 $2\times10^{-5}$ bit。
- $d=4$ 劣于 $d=2$。偶数个谐波的截断在零弧长处二阶导为零（数值上 $T_4$ 沿对角的二阶差分为 0，$d=2,6$ 时为真值的 2 倍），而编码器的概率集中在小弧长配对。另算得 $d=6,10,12,24$ 的码率为 1.193732、1.185458、1.183575、1.183323。
- 上确界范数界在 logit 残差小于 1 之后才有信息量；加权散度与真实超额在所有维数上几乎重合。

~~~sh
python tests/test_interaction.py
python experiments/interaction_rank.py --output results/interaction_rank.json
~~~

## 适用范围

结论针对固定的有限源、锚点和评价矩阵。连续流形上的类比需要核算子的谱和一致余项控制，此处未证。交互秩刻画“能否精确实现最优测试信道”，受限族在匹配失真下损失多少码率仍由数值比较回答。真实特征云上的交互秩与残差尚未测量。
