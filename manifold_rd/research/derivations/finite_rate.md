# 有限字母表的率失真证书与物理链路

从变分目标推出 BA 更新、对偶可行解和时分上界。互信息、固定长度索引和信道使用次数分别记账。

[研究索引](../README.md) · [符号和引用所用书目](../../ref.bib)

## Finite certificates and the physical link

For finite source mass $p_i$, fixed costs $c_{ij}$ and slope $\beta\ge0$, minimize $I_2+\beta D$. Introduce a trial prior $r$; $\sum_i p_iD_{\mathrm{KL},2}(W_i\Vert r)= I_2+D_{\mathrm{KL},2}(\bar W\Vert r)$, $\bar W_j=\sum_i p_iW_{ij}$. The row minimizer and marginal update are 

$$
Z_i=\sum_j r_j2^{-\beta c_{ij}},\quad W_{ij}=r_j2^{-\beta c_{ij}}/Z_i,\quad r'_j=\sum_i p_i W_{ij}.
$$

These give a feasible upper value $U=I_2(W)+\beta D(W)$ [`blahut1972computation`]. Put $G_j=\sum_i p_i2^{-\beta c_{ij}}/Z_i$ and $\alpha_i=-\log_2 Z_i-\log_2\max_jG_j$. Then $\sum_i p_i2^{\alpha_i-\beta c_{ij}}\le1$ for every $j$. For any channel, apply nonnegativity of relative entropy conditionally given $j$ to the conditional source distribution against the normalized weights $p_i2^{\alpha_i-\beta c_{ij}}$. Summation gives $I_2+\beta D\ge L=\sum_i p_i\alpha_i$. Thus $L-\beta D_0$ is a rate lower bound at target $D_0$; interpolating achievable $(D,I)$ by independent time sharing gives an upper bound. Weight these bounds by $p_c$ at a common slope for shared context. They certify only the declared finite alphabets.

When transmission uses $k$ channel uses per source sample, data processing gives $R_\delta(D)\le k C_{\rm ch}$ for an admissible scheme under the corresponding memoryless channel assumptions. A scalar real AWGN use with gain $h$ has capacity $\tfrac12\log_2(1+h^2P/\sigma^2)$. A small $h$ can cause outage; it does not change a fixed source RDF. Conversely, the index experiment samples $J\sim W(\cdot|U)$, sends $\lceil\log_2 M\rceil$ uncoded BPSK symbols, and decodes an index before evaluating loss. With unit-energy BPSK and $E_b/N_0=\gamma$, hard detection in AWGN has error $\Phi(-\sqrt{2\gamma})$ per bit. In coherent block Rayleigh the same gain affects all bits of an index. Record mutual information, actual uses, bit error and index block error separately. Codebook and model are pre-shared.

## 同一有限成本矩阵对源权重的敏感性

令 $M\ge2$ 为重构个数，$0\le c_{ij}\le c_{\max}$，$d_{\mathrm{TV}}(p,\tilde p)\le\epsilon\le1-1/M$。固定 $W$，每行平均失真都在 $[0,c_{\max}]$，因此失真变化至多 $u=\epsilon c_{\max}$。输出边缘的总变差也至多 $\epsilon$。对两个边缘作最大耦合，错误概率为实际总变差 $t$，用错误指示并条件化得单向条件熵至多 $h_2(t)+t\log_2(M-1)$。该式在所设区间单调。行熵的范围为 $[0,\log_2 M]$，更换输入权重使其均值变化至多 $\epsilon\log_2 M$。所以同一 $W$ 的信息差不超过 $v=h_2(\epsilon)+2\epsilon\log_2 M$。

把 $\tilde p$ 下失真 $\le D-u$ 的任意信道移到 $p$ 下，给出 $R_p(D)\le R_{\tilde p}(D-u)+v$；反向移动 $p$ 下失真 $\le D$ 的信道给出 $R_p(D)\ge R_{\tilde p}(D+u)-v$。适用于偏移后的可行目标。该敏感性估计假定两次计算使用同一合法字母表和成本矩阵。替换数据集、锚点或任务头时需要另证。
