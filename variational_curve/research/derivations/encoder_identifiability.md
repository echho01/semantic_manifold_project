# 哪些曲面变化能与码字质量的变化区分开

[研究索引](../README.md) · [连续曲面的一阶变分](surface_variation.md) · [检验实现](../experiments/measure_deformation.py) · [2,048 个测试样本的实测](../results/encoder_contrasts.json)

## 有限模型与可识别条件

评价矩阵 $C^g$、输入概率 $p_i$、锚点和锚点身份固定。编码器另有代价矩阵 $C$、温度参数 $s>0$ 和严格正的列质量 $m$。自然指数定义随机行

$$
W_{ij}(C,m)=\frac{m_j e^{-sC_{ij}}}{\sum_k m_k e^{-sC_{ik}}}.
$$

令 $\Delta=C'-C$。问题是能否找到某组 $m'$ 使所有行和所有列的实际分配概率保持相同。相同行概率等价于任意两列的对数优势相同，因此

$$
s(\Delta_{ij}-\Delta_{ik})=
\log(m'_j/m_j)-\log(m'_k/m_k).
$$

右侧对输入行 $i$ 恒定。取行 $i'$ 与之相减，可得所有行列交叉差为零。反过来，任选 $(i_0,k_0)$，定义 $a_i=\Delta_{ik_0}$ 和 $b_j=\Delta_{i_0j}-\Delta_{i_0k_0}$。交叉差为零给出 $\Delta_{ij}=a_i+b_j$。取 $m'_j=m_j e^{s b_j}/\sum_\ell m_\ell e^{s b_\ell}$，列项抵消，行项被行归一化消去，原信道即被恢复。这两个方向都要求比较同一列身份。质量有零分量时先限制到共同的有效支撑，两个模型支撑不同则直接可区分。温度可变时对 $s'C'-sC$ 应用同一判据。

## 实测量

令 $H_n=I_n-\mathbf1\mathbf1^\top/n$，$H_M$ 同理，并取 $T=s(C'-C)$。行、列加法项构成一个线性空间。其在等权 Frobenius 内积下的正交补投影为 $H_n T H_M$。左、右两次中心化使每一列和每一行的均值都为零。判据可以用一个标量记录

$$
e_{\rm int}=\sqrt{\frac1{nM}\sum_{ij}(T_{ij}-\bar T_{i\cdot}-\bar T_{\cdot j}+\bar T_{\cdot\cdot})^2}.
$$

理论上 $e_{\rm int}=0$ 当且仅当所有交叉差都为零。写成均值避免显式构造中心化矩阵。指标度量 logit 的变化；它没有信息率单位。对同一固定矩阵，实测另外计算 $I_2(W)$、$D_g(W)$ 和 $I_2+\lambda D_g$。比较 $R(D_g)$ 另需变化斜率并匹配失真。

## 对参数的导数

设 $S_{ij}=\delta\log m_j-s\,\delta C_{ij}-C_{ij}\delta s$。微分行归一化得 $\delta W_{ij}=W_{ij}(S_{ij}-\sum_k W_{ik}S_{ik})$。写出互信息的两个熵项，$\sum_{ij}p_i\delta W_{ij}=\sum_j\delta\bar W_j=0$，得到

$$
\delta (I_2+\lambda D_g+\gamma B)=
\sum_i p_i\operatorname{Cov}_{W_i}\!\left(\log_2\frac{W_{ij}}{\bar W_j}+\lambda C^g_{ij},\, S_{ij}\right)+\gamma\delta B .
$$

以 $C^h_{ij}=\|F(z_i)-F(z_j)\|^2$ 为例，其变分为 $2\langle F(z_i)-F(z_j),\delta F(z_i)-\delta F(z_j)\rangle$。更新后评价仍使用同一 $C^g$。对任意辅助边缘 $\pi$，$\sum_i p_iD_{\mathrm{KL},2}(W_i\Vert\pi)=I_2(W)+D_{\mathrm{KL},2}(\bar W\Vert\pi)$。因此小批量中可联合训练 $\pi$，测试时用真实全测试集边缘。在 $C=C^g,\ s=\lambda\ln2$ 且允许任意质量和边界极限时，此恒等式和逐行最小化给出 BA 最优点。

## 有伸缩界的实现

令 $z=(2\sqrt q,\sqrt\eta\,u)$，其中 $u$ 已单位化。第一块对 $u$ 的微分满足 $\sum_\ell(dq_\ell)^2/q_\ell$，因此弦长在局部产生 Fisher 拉回和 $\eta$ 加权的特征度量。映射 $F(z)=z+\alpha B\tanh(Az)$ 满足 $\|A\|_{\rm op},\|B\|_{\rm op}\le1$，$0\le\alpha<1$。由正、反三角不等式得

$$
(1-\alpha)\|z-v\|\le\|F(z)-F(v)\|\le(1+\alpha)\|z-v\|.
$$

平方后给出有限点对代价界，取微分还得到相应的局部度量界。实现将矩阵的 Frobenius 范数截到 1，足以满足上述条件。邻锚 $v$ 的二阶差分使用环境弦中点 $(z+v)/2$，故该正则项惩罚映射延拓的非线性。把它解释成曲面的内在曲率还需额外假设。
