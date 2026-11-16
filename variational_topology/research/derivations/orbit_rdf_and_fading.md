# 对称轨道上的关系率失真与衰落信道下的询问下界

持久图在读出图的自同构下不变，带位置身份的连通回答则会改变。本文给出轨道信源的关系率失真闭式解，并把多询问的信息下界推广到信道使用次数、遍历衰落和中断三种形式。

[研究索引](../README.md) · [连通风险与双模板基准](connectivity.md) · [多询问的信息量下界](query_converse.md) · [轨道实验实现](../experiments/orbit_rdf.py) · [数值核验](../tests/test_orbit.py) · [320 张字形的实测](../results/orbit_rdf.json)

## 记号

$G=(V,E)$ 为固定的有限连通读出图，$f:V\to[0,1]$ 为顶点分数，$m_f(a,b)=\min_{\pi:a\leadsto b}\max_{v\in\pi}f(v)$ 为合并层次，$B_f(a,b,t)=\mathbf 1\{m_f(a,b)\le t\}$ 为询问 $L=(a,b,t)$ 的回答。询问分布 $\nu$ 与信源独立，$d_\nu(f,f')=\mathbb E_\nu|B_f(L)-B_{f'}(L)|$。对 $\sigma\in\operatorname{Aut}(G)$ 记 $(f\circ\sigma)(v)=f(\sigma v)$，于是 $(f\circ\sigma)\circ\tau=f\circ(\sigma\tau)$，这是右作用。

## 一　合并函数的结构与持久图丢失的信息

**命题 1（合并函数）。** $m_f$ 对称，$m_f(a,a)=f(a)\le m_f(a,b)$，且对任意 $a,b,c$ 有强三角不等式 $m_f(a,c)\le\max\{m_f(a,b),m_f(b,c)\}$。$f$ 的下水平集 $H_0$ 持久图是数组 $m_f$ 的函数，且该函数在用任意置换 $\varsigma\in S_V$ 同时重标两个自变量时不变。

**证明。** 路径含两个端点，所以 $m_f(a,b)\ge\max\{f(a),f(b)\}$，取平凡路径得 $m_f(a,a)=f(a)$。把实现 $m_f(a,b)$ 与 $m_f(b,c)$ 的两条路径在 $b$ 处连接，得到一条从 $a$ 到 $c$ 的路径，其最大分数为二者较大值，取最小即得强三角不等式。阈值 $t$ 下的活动顶点集为 $\{a:m_f(a,a)\le t\}$，其上的连通分量是关系 $m_f(a,b)\le t$ 的等价类（强三角不等式给出传递性）。因此 $m_f$ 决定了各阈值下的分量划分及其包含关系，即合并树。按长者规则，$H_0$ 持久图由合并树各分支的出生值和合并值组成，是合并树的函数。置换 $\varsigma$ 只改变叶子的名称，出生值与合并值的多重集不变。∎

**命题 2（自同构不变性）。** 对 $\sigma\in\operatorname{Aut}(G)$，

$$
m_{f\circ\sigma}(a,b)=m_f(\sigma a,\sigma b),\qquad B_{f\circ\sigma}(a,b,t)=B_f(\sigma a,\sigma b,t).
$$

若 $\sigma$ 还作用在由 $G$ 构造的复形上（$G$ 的团复形，或 $\sigma$ 为网格对称时的像素立方复形），则 $f\circ\sigma$ 与 $f$ 的所有未标记持久图在每个同调维数上相同。

**证明。** 自同构把 $a$ 到 $b$ 的路径一一映为 $\sigma a$ 到 $\sigma b$ 的路径，$\max_{v\in\pi}(f\circ\sigma)(v)=\max_{u\in\sigma\pi}f(u)$。两边对路径取最小得第一式，第二式随之成立。对每个阈值 $t$，$\sigma$ 把 $f\circ\sigma$ 的下水平子复形同构地映到 $f$ 的下水平子复形，并与阈值间的包含映射交换。两个滤过因此同构，持久模同构，持久图相同。∎

**推论。** 任何只依赖未标记持久图的损失在 $f$ 与 $f\circ\sigma$ 之间为零。关系失真为

$$
d_\nu(f,f\circ\sigma)=\mathbb E_\nu\,|\mathbf 1\{m_f(a,b)\le t\}-\mathbf 1\{m_f(\sigma a,\sigma b)\le t\}|,
$$

阈值连续且与端点独立时等于 $\mathbb E|F(m_f(a,b))-F(m_f(\sigma a,\sigma b))|$。它为零当且仅当合并函数在 $\nu$ 的支撑上几乎处处 $\sigma$ 不变。$H=W$ 的四邻接网格的自同构群是 8 阶二面体群，主稿的水平镜像是其中一个元素。

阈值积分风险还可按错误类型拆分：$|F(m_f)-F(m_{\hat f})|=[F(m_{\hat f})-F(m_f)]_++[F(m_f)-F(m_{\hat f})]_+$。$m_{\hat f}>m_f$ 时在两层之间的阈值上原图连通而接收图不连通，故第一项是误断开的概率，第二项是误连接的概率。训练时可给两项不同权重。

## 二　轨道信源的关系率失真

**设定。** $\Gamma\le\operatorname{Aut}(G)$，模板 $x$ 的轨道 $\mathcal O=\{x\circ\gamma:\gamma\in\Gamma\}$，$K=|\mathcal O|=|\Gamma|/|\mathrm{Stab}(x)|$。$\nu$ 与信源独立，并在 $(a,b,t)\mapsto(\gamma a,\gamma b,t)$ 下不变（例如有放回的均匀有序点对乘以任意阈值分布）。信源 $X$ 在 $\mathcal O$ 上均匀，合法重构限定在 $\mathcal O$ 内。记 $\varrho_1,\dots,\varrho_K$ 为 $\{d_\nu(x,x')\}_{x'\in\mathcal O}$。

**引理（平衡失真矩阵）。** 对任意分数 $g,g'$ 和 $\tau\in\Gamma$，$d_\nu(g\circ\tau,g'\circ\tau)=d_\nu(g,g')$。因此 $\mathcal O$ 上的失真矩阵每一行、每一列都是 $(\varrho_i)$ 的一个排列。

**证明。** 由命题 2，$d_\nu(g\circ\tau,g'\circ\tau)=\mathbb E_\nu|B_g(\tau a,\tau b,t)-B_{g'}(\tau a,\tau b,t)|$。$\nu$ 在 $\tau$ 下不变，换元 $(a,b)\mapsto(\tau a,\tau b)$ 后即为 $d_\nu(g,g')$。$\Gamma$ 以右乘可迁地置换 $\mathcal O$：取 $y=x\circ\tau$，则 $\{d_\nu(y,y')\}_{y'\in\mathcal O}=\{d_\nu(x\circ\tau,x'\circ\tau)\}_{x'\in\mathcal O}=\{d_\nu(x,x')\}_{x'\in\mathcal O}$。$d_\nu$ 对称，故列与行相同。∎

**定理 1（轨道率失真）。** 令 $\delta_{\max}=K^{-1}\sum_i\varrho_i$。关系失真不超过 $\delta$ 所需的最小信息（bit）为

$$
R_Q(\delta)=\log_2K-H_2(P_\beta),\quad 0\le\delta<\delta_{\max};\qquad R_Q(\delta)=0,\quad\delta\ge\delta_{\max},
$$

其中 $P_\beta(i)\propto2^{-\beta\varrho_i}$，$\beta\in(0,\infty]$ 由 $\sum_iP_\beta(i)\varrho_i=\delta$ 确定（$\beta=\infty$ 表示集中在 $\varrho_i=0$ 的指标上的均匀分布）。同时，任何只依赖未标记持久图的失真在 $\mathcal O\times\mathcal O$ 上恒为零，其率失真函数恒为零。

**证明。** 定义 $\varphi(\delta)=\max\{H_2(P):\sum_iP(i)\varrho_i\le\delta\}$。熵是凹函数、约束是线性的，所以 $\varphi$ 凹且单调不减。$\delta\ge\delta_{\max}$ 时均匀分布可行，$\varphi=\log_2K$。$\delta<\delta_{\max}$ 时由拉格朗日条件最大值在 Gibbs 分布 $P_\beta$ 处取得。

*逆定理。* 取任意满足 $\mathbb E d_\nu(X,\hat X)\le\delta$ 的测试信道。给定 $\hat X=\hat x$，记条件失真为 $\delta_{\hat x}=\sum_xP(x\mid\hat x)d_\nu(x,\hat x)$。由引理，列 $\{d_\nu(x,\hat x)\}_x$ 是 $(\varrho_i)$ 的排列，因此 $H_2(X\mid\hat X=\hat x)\le\varphi(\delta_{\hat x})$。对 $\hat X$ 取期望并用 Jensen 不等式和单调性，$H_2(X\mid\hat X)\le\varphi(\mathbb E\delta_{\hat X})\le\varphi(\delta)$。于是 $I_2(X;\hat X)=\log_2K-H_2(X\mid\hat X)\ge\log_2K-\varphi(\delta)$。

*可达性。* 令 $\hat X$ 在 $\mathcal O$ 上均匀，取反向信道 $P(x\mid\hat x)=2^{-\beta d_\nu(x,\hat x)}/Z$。列是同一组数的排列，故归一化常数 $Z=\sum_i2^{-\beta\varrho_i}$ 与 $\hat x$ 无关。$X$ 的边缘为 $K^{-1}\sum_{\hat x}2^{-\beta d_\nu(x,\hat x)}/Z=K^{-1}$，因为行也是同一组数的排列。所以 $X$ 确为均匀信源，平均失真为 $\sum_iP_\beta(i)\varrho_i=\delta$，$I_2(X;\hat X)=\log_2K-H_2(P_\beta)$。$\delta\ge\delta_{\max}$ 时任一常量重构的失真为行均值 $\delta_{\max}$，零信息可达。

*持久图部分。* 命题 2 表明轨道上任意两元素的未标记持久图相同，常量重构即达到零失真。∎

证明只用到行列的排列性质，对非交换群和非平凡稳定子同样成立，不需要构造群差变量。

**推论（双模板）。** $K=2$、$\varrho=(0,\alpha)$、$\alpha>0$ 时 $P_\beta=(1-p,p)$，$\delta=\alpha p$，$\delta_{\max}=\alpha/2$，得 $R_Q(\delta)=1-h_2(\delta/\alpha)$（$0\le\delta<\alpha/2$），与[双模板基准](connectivity.md)一致。

**零失真端点。** 若 $\nu$ 对所有顶点对（含 $a=b$）和阈值有满支撑，则 $d_\nu(x,x')=0$ 当且仅当 $x=x'$，因为 $B_x(a,a,t)=\mathbf 1\{f_x(a)\le t\}$ 已确定 $f_x$。此时 $R_Q(0)=\log_2K$，稳定子平凡时等于 $\log_2|\Gamma|$，方形网格上最多为 3 bit。

**适用范围。** 定理固定了重构字母表和与信源独立的询问分布。允许轨道外的重构时零率失真可以更小。按原图前景采样的询问会透露信源信息，应改用[条件信息界](query_converse.md)。

## 三　信道使用次数、遍历衰落与中断

沿用[多询问下界](query_converse.md)的设定：一次传输，询问 $L_{1:d}$，回答 $B_j=B_X(L_j)$，误判率 $e_j=\Pr(B_j\ne B_{\hat X}(L_j))$，共享上下文 $C$，接收量 $V$，$H_Q=H_2(B_{1:d}\mid L_{1:d},C)$。给定 $(X,C)$ 时询问采样与信道噪声独立。主定理为 $[H_Q-\sum_jh_2(e_j)]_+\le I_2(X;V\mid C)$。记 $\Delta=H_Q-\sum_jh_2(e_j)$，信道为 $Y=hZ+N$，$N\sim\mathcal{CN}(0,\sigma^2I_k)$，$\mathbb E\|Z\|^2\le k$。

**推论 1（信道使用次数）。** AWGN（$h=1$）下 $I_2(X;V\mid C)\le k\log_2(1+1/\sigma^2)$，故

$$
k\ \ge\ \frac{[\Delta]_+}{\log_2(1+1/\sigma^2)} .
$$

**推论 2（相干块衰落，遍历形式）。** 设每次传输一个增益 $h$，接收端已知、发送端未知，且 $h$ 与 $(X,L_{1:d},C)$ 独立。则对平均误判率 $e_j$，

$$
[\Delta]_+\ \le\ k\,\mathbb E_h\log_2\!\left(1+|h|^2/\sigma^2\right).
$$

Rayleigh 衰落（$|h|^2\sim\mathrm{Exp}(1)$）时右端为 $k\,e^{\sigma^2}E_1(\sigma^2)/\ln2$，$E_1$ 为指数积分。

**证明。** 此时 $V=(Y,h)$，主定理直接适用。$I(X;V\mid C)=I(X;h\mid C)+I(X;Y\mid h,C)=I(X;Y\mid h,C)$。发送端不知道 $h$，故 $Z$ 与 $h$ 独立，给定 $h$ 时能量约束不变，条件信道是增益 $|h|^2$ 的复高斯信道，$I_2(X;Y\mid h=h_0,C)\le k\log_2(1+|h_0|^2/\sigma^2)$。对 $h$ 取期望即得。Rayleigh 闭式来自 $\int_0^\infty\ln(1+g/\sigma^2)e^{-g}\,dg=e^{\sigma^2}E_1(\sigma^2)$。∎

该式对单次块衰落传输的平均误判率成立，不需要跨衰落块编码。发送端有信道信息并做功率自适应时，右端应换成相应的容量表达式。

**推论 3（中断形式）。** 设目标 $\epsilon_j\le1/2$ 且 $R_\epsilon=H_Q-\sum_jh_2(\epsilon_j)>0$。记 $e_j(h)=\Pr(B_j\ne\hat B_j\mid h)$。则

$$
\Pr_h\{\exists j:\ e_j(h)>\epsilon_j\}\ \ge\ \Pr_h\{k\log_2(1+|h|^2/\sigma^2)<R_\epsilon\},
$$

Rayleigh 衰落时右端为 $1-\exp[-\sigma^2(2^{R_\epsilon/k}-1)]$。

**证明。** 在给定 $h=h_0$ 的条件概率空间上应用主定理。$h$ 与 $(X,L_{1:d},C)$ 独立，故 $H_2(B_{1:d}\mid L_{1:d},C,h=h_0)=H_Q$，于是 $H_Q-\sum_jh_2(e_j(h_0))\le k\log_2(1+|h_0|^2/\sigma^2)$。若所有 $e_j(h_0)\le\epsilon_j\le1/2$，由 $h_2$ 在 $[0,1/2]$ 上单调增得 $R_\epsilon\le k\log_2(1+|h_0|^2/\sigma^2)$。取逆否命题并对 $h$ 求概率。∎

**何时取等。** 若回答是相互独立的等概率比特，且 $V$ 是各比特经过交叉概率为 $e_j$ 的独立二元对称信道的输出，则 $H_Q=d$，$I_2(X;V\mid C)=\sum_j(1-h_2(e_j))$，主定理取等。回答相关时联合条件熵只把共同信息计一次。

## 四　数值记录

- [单元测试](../tests/test_orbit.py)：在 $4\times4$ 随机灰度分数上用 Floyd–Warshall 型递推算出全部 $m_f$，核对强三角不等式、8 个二面体元素下的 $m_{f\circ\sigma}(a,b)=m_f(\sigma a,\sigma b)$ 以及 $H_0$ 持久图不变；核对灰度与二值轨道失真矩阵的行排列性质；闭式解与 Blahut–Arimoto 解在多个斜率上一致；双模板特例与 Rayleigh 遍历闭式分别与二元熵公式和数值积分一致。
- 合成模板：200 个 $5\times5$ 随机二值模板（密度 0.45，种子 20261115），全部满足行排列性质，闭式与 BA 的最大差为 $8.9\times10^{-16}$ bit。
- Omniglot：与镜像实验相同的 320 张官方 evaluation 字形（ID 哈希一致，镜像 $\alpha$ 与[既有记录](../results/reflection_queries.json)逐图差为 0）。全部 320 张在 8 个网格对称下稳定子平凡，轨道大小为 8。全像素询问下 $\delta_{\max}$ 均值 0.011541、中位数 0.010490；每张字形与其他轨道元素的最小、最大误判率的均值为 0.011312 和 0.014007。$\delta/\delta_{\max}=0,0.1,0.25,0.5$ 处的平均码率为 3.0000、2.3124、1.6114、0.7698 bit/字形。镜像二元轨道在相同比例处为 1、0.7136、0.4564、0.1887 bit，只取决于比例。
- 信道公式求值：0、5、10、15、20 dB 下，AWGN 每次复使用至多 1.0000、2.0574、3.4594、5.0278、6.6582 bit，Rayleigh 遍历值为 0.8603、1.7160、2.9065、4.3302、5.8840 bit。这些是公式求值，不是链路测量。

~~~sh
python tests/test_orbit.py
python experiments/orbit_rdf.py --data data/omniglot --reference results/reflection_queries.json --output runs/orbit_rdf.json --figure ../figures/orbit_rdf.png
~~~
