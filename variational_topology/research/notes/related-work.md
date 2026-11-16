# Paper 2：相关工作与贡献定位

核查日期：2026-10-20（UTC）。本记录服务于 `topology` 分支的 **Topological Semantic Communication: Rate, Noise, and Knowledge Synchronization**。正式书目信息见 [现用书目](../../ref.bib)，原始来源、DOI、出版版本及当前 65 条书目/30 个 Paper 2 引用的核验见 [paper2-verification.md](paper2-verification.md)。本稿定理及其成立条件见 [证明索引](../derivations/coverage_and_dictionary.md)。

## 1. 本稿的问题与比较单位

本稿在选定表示、滤流和任务后，研究三类具体通信问题：共享字典需要何种分辨率才能认证源结构；接收端怎样校准加噪特征的滤流尺度并控制整块结构错误；接收端已有结构参考时，需要多少比特完成同步，以及保留任务分支身份增加多少信息。

相关工作据此按**研究对象、已有侧信息、优化资源和恢复事件**比较。已有文献已经研究 persistence diagram 的稳定性、量化、sketching、信息论描述和高维噪声下的失效。本稿的贡献定位建立在下面可核查的模型差别和具体定理上。

| 本稿操作模型 | 优化或评价的量 | 对应结果 |
|---|---|---|
| 带覆盖和 reach 承诺的共享字典 | 指定均匀重构误差下的每符号索引长度；持久像认证 | `thm:threshold`、`prop:identifiability` |
| 已知噪声能量的模拟特征块 | 参考滤流误差、去偏误差、接触与真分区的整块可靠性 | `thm:reparam`、`cor:debias`、`prop:tight`、`cor:price` |
| 固定解释的两个 KB，或局部 diagram 更新 | 稠密字典失配极限；满足最坏情况误差的固定消息长度 | `thm:floor`、`thm:sync` |
| 标量分数或携带 Reeb 边身份的任务 | 小失真 Shannon 率及分支条件熵开销 | `thm:reeb`、`prop:screen` |

## 2. Diagram 估计、压缩与信息论

| 原始工作 | 已有结果及其资源模型 | 本稿的对应结果与差别 |
|---|---|---|
| Divol–Lacombe, [Estimation and Quantization of Expected Persistence Diagrams](https://proceedings.mlr.press/v139/divol21a.html), ICML 2021 | 研究随机 diagram 的 expected persistence diagram，这是一种统计汇总；估计其分布相关对象，并以小支撑测度进行量化，给出收敛率。 | `thm:sync` 的对象是已知分离参考附近的有限 diagram。独立端点更新形成明确局部源类，目标是在最坏情况下用固定二进制报文达到指定 bottleneck 误差。 |
| Sheehy–Sheth, [Sketching Persistence Diagrams](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SoCG.2021.57), SoCG 2021 | 给出不同加权点数的近似序列，每个规模具有 2-approximation 保证；整个序列可紧凑存储，并加速 diagram 距离计算。 | 本稿将端点格子索引联合编码为整数，证明精确位数 $\lceil\log_2(m^{2K})\rceil$。比较依据是消息字母表和误差，而不是仅比较输出 diagram 的点数。 |
| Perumalsamy–Gunalan–Thinakaran, [An Information Theory of Persistent Homology: Entropy, the Data Processing Inequality, and Rate–Distortion Bounds for Topological Features](https://www.mdpi.com/2227-7390/14/8/1385), Mathematics 2026 | 已直接提出 persistent homology 的熵、互信息与率失真描述。 | 本稿规定有限定长的更新协议、指定特征信道和带有限跨边错误代价的任务源，以实际构造与 converse 确定各自通信成本。 |

因此，Introduction 同时引用 `divol2021quantization`、`sheehy2021sketching`、`perumalsamy2026information`，再给出本稿的局部源类、侧信息和精确消息长度。描述长度的简单上界 `prop:desc` 与局部最优长度 `thm:sync` 分别陈述；essential bars、特征对应和端点有效性都有明确约定。

局部同步源类允许 $K$ 个分离条的 $2K$ 个端点独立更新。这一 diagram 源类本身足以定义和证明码长问题；若改成固定点数欧氏 Rips 点云产生的 diagram，需要重新研究可实现的端点变化范围。

## 3. 高维噪声与几何推断

| 原始工作 | 已有结果与假设方向 | 本稿的对应推进 |
|---|---|---|
| Damrich–Berens–Kobak, [Persistent Homology for High-dimensional Data Based on Spectral Methods](https://proceedings.neurips.cc/paper_files/paper/2024/hash/4a32a646254d2e37fc74a38d65796552-Abstract-Conference.html), NeurIPS 2024 | 分析高维噪声对常规 persistence 的影响，并通过 kNN 图的 diffusion/effective-resistance 距离改善拓扑恢复。 | 本稿固定各向同性特征信道，构造保留干净距离的噪声参考滤流，证明有限样本误差与去偏预算，并分析明确阈值图接收机的分区事件。 |
| Hiraoka–Imoto–Kanazawa–Liu, [Curse of Dimensionality on Persistence Diagrams](https://www.aimsciences.org/article/doi/10.3934/fods.2025008), Foundations of Data Science 2026 | 给出 HDLSS 条件下 diagram 不一致的渐近分析及 normalized PCA 的理论改善。 | 本稿保留 $n,D_a,\delta$ 参数，区分噪声偏移和随机距离涨落；精确的单对非中心卡方分布及保留依赖关系的两侧界连接有限样本与渐近尺度。 |
| Chazal–Cohen-Steiner–Mérigot, [Geometric Inference for Probability Measures](https://link.springer.com/article/10.1007/s10208-011-9098-0), FoCM 2011 | 以概率测度和 distance-to-measure 等对象进行更稳健的几何推断。 | 本稿明确自己分析的原始/去偏阈值图及已知中心接收机。利用额外分布或全局结构的接收机需要以其自身定义比较。 |

比较高维结论时先统一噪声单位。本稿使用 $\mathbb E\|N_i\|^2=\sigma^2$；固定每坐标方差会使总能量随 $D_a$ 变化，固定总能量则会使每坐标方差随 $D_a$ 缩小。

本稿的接触结果为

$$
1-(1-p_\theta)^m\le\Pr(r_{\rm out}\le\theta)\le\min\{1,m^2p_\theta\}.
$$

单对概率 $p_\theta$ 由非中心卡方 CDF 精确计算。选取 $m$ 对互不相交的噪声向量证明下界，对全部 $m^2$ 对取 union bound 证明上界。在联合极限下得到常数 $2^{-5/4}$ 与 $1/2$ 的 bracket。

这里跨簇接触、整个图在零半径连通、阈值分量数正确和真分区恢复分别定义。`cor:price` 对已知中心与阈值图统一采用整块错误概率，噪声尺度之比为 $(D_a/\log n)^{1/4}$。不同事件或逐点/整块可靠性的混合比较会改变对数因子。

## 4. 几何恢复与稳定性的基础地位

| 原始工作 | 本稿使用的基础 | 本稿自行证明的部分 |
|---|---|---|
| Federer, [Curvature Measures](https://doi.org/10.1090/S0002-9947-1959-0110078-1), 1959；Niyogi–Smale–Weinberger, [Finding the Homology of Submanifolds with High Confidence from Random Samples](https://link.springer.com/article/10.1007/s00454-008-9053-2), 2008 | reach 以下的唯一最近点几何、同调恢复和局部体积控制。 | 选用 $\eta=r/6$，以严格覆盖/扰动预算认证 Rips 持久像；将最小覆盖数对应到均匀重构约束下的精确索引长度。 |
| Chazal–Oudot, [Towards Persistence-Based Reconstruction in Euclidean Spaces](https://doi.org/10.1145/1377676.1377719), 2008 | 用持久像完成几何重构的路线。 | 通过 Čech–Rips 两条包含链给出本稿尺度常数，附录完整证明扰动后的秩夹逼。 |
| Cohen-Steiner–Edelsbrunner–Harer, [Stability of Persistence Diagrams](https://link.springer.com/article/10.1007/s00454-006-1276-5), 2007；Chazal–de Silva–Oudot, [Persistence Stability for Geometric Complexes](https://link.springer.com/article/10.1007/s10711-013-9937-z), 2014 | 函数和几何复形的 persistence 稳定性。 | 具体特征信道的平方距离浓度、顶点保持零的参考滤流、去偏与方差误差预算，以及固定 KB 解释的有限分辨率区间。 |
| Lesnick, [The Theory of the Interleaving Distance on Multidimensional Persistence Modules](https://link.springer.com/article/10.1007/s10208-015-9255-y), 2015 | 普适函数稳定性与 interleaving 的最大性。 | 将该既有结果用于选择结构失真。主文明确素域和逐点有限维模块，未把普适性包装为本稿的新定理。 |

`thm:threshold` 的 converse 来自“均匀最近原子误差”等价于覆盖条件。`prop:identifiability` 则构造相同索引法则下不同支持拓扑，说明为何必须有覆盖承诺。这两个结论各自回答已定义的操作问题，不借拓扑类型计数混合每符号码率和整条 descriptor 的长度。

## 5. 学习式语义通信与知识同步

| 原始工作 | 与本稿的联系 | 引用的准确作用 |
|---|---|---|
| Xie et al., [Deep Learning Enabled Semantic Communication Systems](https://doi.org/10.1109/TSP.2021.3071210), 2021；Gündüz et al., [Beyond Transmitting Bits](https://doi.org/10.1109/JSAC.2022.3223408), 2023；Shao–Mao–Zhang, [Learning Task-Oriented Communication for Edge Inference](https://doi.org/10.1109/JSAC.2021.3126087), 2022 | 从任务和表示规定通信目标。 | 建立工程问题的背景；本稿通过固定失真、源类和接收机分析通信成本。 |
| Moor et al., [Topological Autoencoders](https://proceedings.mlr.press/v119/moor20a.html), ICML 2020 | 将拓扑结构用于表示学习。 | 说明结构约束如何进入学习目标。 |
| Erak et al., [Topology-Preserving Deep Joint Source-Channel Coding for Semantic Communication](https://arxiv.org/abs/2603.17126), 2026 预印本 | 在端到端 DeepJSCC 训练中加入图像和潜特征的 persistence 正则项。 | 动机与方法家族引用。作者公开记录写为 submitted，书目保留预印本类型。本文指定的特征信道并非其神经编码器的性能极限。 |
| Wang–Guo, [Transceiver Cooperative Learning-aided Semantic Communications Against Mismatched Background Knowledge Bases](https://arxiv.org/abs/2301.03133), 2023 预印本 | 通过协作学习和参数交换缓解 KB 不一致。 | KB 不匹配的应用背景。本稿定量分析固定解释规则下的极限，并研究显式结构更新消息。 |

`thm:floor` 给出 $|d_{\rm topo}(K)-\delta_{\rm KB}|\le\varepsilon_R(K)$。因此稠密字典极限与有限 $K$ 的误差分开报告，有限曲线允许从下方接近极限。额外配准和同步会改变解释规则，应按新的接收机模型评价。

本分支的实证部分使用已运行的规则几何构造、Monte Carlo 信道块和同条目的有限数据表示。相应数据和脚本给出来源、随机种子、条件、事件定义和原始结果。将来若加入神经编解码器的性能对比，应共同约定信道使用数、功率、源数据和结构恢复事件；当前理论可以给出校准与精度参照。

## 6. Reeb 任务与 Shannon 率失真

Dey–Wang 的 [Reeb Graphs: Approximation and Persistence](https://link.springer.com/article/10.1007/s00454-012-9463-z) 和 Kaluba–Marzantowicz–Silva 的 [On Representation of the Reeb Graph as a Sub-Complex of Manifold](https://www.tmna.ncu.pl/static/published/2015/v45n1-15.pdf) 提供 Reeb 图与其商映射的结构依据。Linder–Zamir 的 [1994 年 Shannon 下界研究](https://doi.org/10.1109/18.340474) 及 Koch 的 [2016 年定理](https://www.tsc.uc3m.es/~koch/files/TIT.2016.2604254.pdf) 支撑本稿每条边上有界、有限微分熵分数源的标量高分辨率公式。

本稿把输出规范分成标量 $T=f(S)$ 与分支输出 $(E,T)$。对同边平方误差、跨边固定正代价，完整 converse 允许 $O(D)$ 的分支错误而保持常数项，得到

$$
R_{\rm task}(D)-R_f(D)\longrightarrow H(E\mid T).
$$

这个开销来自输出要求保留分支身份。若任务仅需分数，分支不确定性无需额外传送。对熵编码量化器，$H(E,Q_\Delta(T))-H(Q_\Delta(T))=H(E\mid Q_\Delta(T))$ 给出有限分辨率恒等式，可用于验证上述开销；标量量化器的绝对率与 Shannon 最优率另有区别。

`prop:screen` 的率失真充分性来自商变量与数据处理；同调满射在 Reeb/Morse 条件下引用原始定理。persistence diagram 提供稳定摘要，完整 Reeb 分支及其概率仍由任务模型规定。

## 7. 当前稿件的贡献表述

可以由现稿完整证明和已核验近邻支持的表述是：

1. **认证字典分辨率。** 给出均匀重构目标下的精确索引长度及可计算 Rips 持久像证书，并构造缺少覆盖承诺时同索引法则、不同源拓扑的例子。
2. **校准后的结构可靠性。** 对明确的 Gaussian 特征信道给出噪声参考与去偏误差；处理重复点；用有限非中心卡方概率和两侧界分析接触，并以相同整块可靠性比较接收机。
3. **知识同步与任务缩减。** 在固定解释规则下估计结构失配极限；对分离局部 diagram 源类给出最优定长更新码；在有限跨边错误代价下证明分支条件熵开销。

这些表述的核心是可验证的条件、实际消息和匹配证明。文献检索结果的数量不承担原创性论证的作用；研究定位以第 2–6 节的近邻映射为依据。

## 8. Paper 1 的独立记录

Paper 1 的文献与已完成审阅对应 main 提交 [`8aad9ea2f51cabb1cd2476dd633582aa9c736e5d`](https://github.com/echho01/semantic-manifold/tree/8aad9ea2f51cabb1cd2476dd633582aa9c736e5d)：

- [相关工作与定位](https://github.com/echho01/semantic-manifold/blob/8aad9ea2f51cabb1cd2476dd633582aa9c736e5d/refs/related-work.md)
- [原始来源核验](https://github.com/echho01/semantic-manifold/blob/8aad9ea2f51cabb1cd2476dd633582aa9c736e5d/refs/verification.md)
- [该提交的完整书目](https://github.com/echho01/semantic-manifold/blob/8aad9ea2f51cabb1cd2476dd633582aa9c736e5d/refs/references.bib)

本分支记录聚焦 Paper 2；两个分支的论文版本和已完成文献记录通过上述固定提交链接对应。
