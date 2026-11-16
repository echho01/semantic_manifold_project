# Paper 2：文献来源与引用核验

核验日期：2026-10-20（UTC）。本记录对应 `topology` 分支的 **Topological Semantic Communication: Rate, Noise, and Knowledge Synchronization**。书目数据位于 [现用书目](../../ref.bib)，近邻比较见 [related-work.md](related-work.md)，本稿数学结果的适用条件与证明位置见 [T2-proof-notes.md](../derivations/coverage_and_dictionary.md)。

## 1. 核验范围与书目快照

核验分为三层：核对作者、题名、出版类型与 DOI 等身份信息；核对被引用内容能否支持主文的相邻论断；检查 BibTeX 结构、重复记录和 LaTeX 引用解析。书目存在且元数据正确，并不自动意味着本稿的定理已由该文献证明。下表与第 2–5 节分别记录这些层次。

| 检查时点 | 共享书目条数 | Paper 1 引用 key 数 | Paper 2 引用 key 数 | 缺失 key | 重复 key / DOI / 规范化题名 |
|---|---:|---:|---:|---:|---|
| Paper 2 首轮书目清理完成 | 64 | 22 | 32 | 0 | 0 / 0 / 0 |
| 替换旧实验协议并加入 scikit-learn 来源后 | **65** | **22** | **30** | **0** | **0 / 0 / 0** |

首次清理把 main 已核验的 **34 条同 key 元数据**合并回当前分支，并补入 **10 条 main 中已核验但本分支缺少的记录**。随后加入 Koch、Kaluba、Damrich、Hiraoka 四条直接支持 Paper 2 的记录，形成 64 条快照。新的有限数据实验增加 `pedregosa2011scikit`，总数成为 65。Paper 2 引用数减少是因为旧 CLIP/DINOv2/Sentence-BERT 实验协议被已经运行的有限数据实验替换。

补齐的核心 keys 包括 `chazal2014persistence`、`divol2021quantization`、`sheehy2021sketching`、`perumalsamy2026information`、`koch2016shannon`、`kaluba2015reeb`、`damrich2024persistent` 和 `hiraoka2026curse`。Chazal 2014 仅保留 `chazal2014persistence` 这一记录；Divol 和 Shao 的现用 keys 分别为 `divol2021quantization` 和 `shao2022taskoriented`。

本轮为 Lesnick、Chazal 2008/2011、Boissonnat、Petersen、Dey、Laurent–Massart、Linder–Zamir、Yu 补齐 DOI；之后补入 Niyogi–Smale–Weinberger 的正式刊版 DOI。共享书目中的未引用条目继续保留，以兼容系列论文。**65 条**是当前记录总数；本记录的逐项引用内容核验围绕下列 **30 条 Paper 2 实际引用**展开。

## 2. 几何与 persistence 的原始依据

| BibTeX key | 经核验的原始文献与主来源 | 本稿引用内容及条件 |
|---|---|---|
| `federer1959curvature` | H. Federer, *Curvature Measures*, Transactions of the AMS 93(3):418–491, 1959。[DOI / AMS](https://doi.org/10.1090/S0002-9947-1959-0110078-1) | reach 与最近点几何。`lem:intrinsic` 在 reach 以下使用唯一最近点及形变收缩。 |
| `niyogi2008finding` | P. Niyogi, S. Smale, S. Weinberger, *Finding the Homology of Submanifolds with High Confidence from Random Samples*, DCG 39:419–441, 2008。[出版社](https://link.springer.com/article/10.1007/s00454-008-9053-2)；DOI `10.1007/s00454-008-9053-2` | 流形同调恢复及局部体积下界；本稿明确给定覆盖与 reach 承诺，再推导索引长度。使用正式刊版 DOI，区分作者旧稿上出现的早期编号。 |
| `cohensteiner2007stability` | D. Cohen-Steiner, H. Edelsbrunner, J. Harer, *Stability of Persistence Diagrams*, DCG 37(1):103–120, 2007。[出版社](https://link.springer.com/article/10.1007/s00454-006-1276-5)；DOI `10.1007/s00454-006-1276-5` | 同一空间上的函数 sup-norm 稳定性。用于单形滤值变化和 sublevel 任务摘要。 |
| `chazal2014persistence` | F. Chazal, V. de Silva, S. Oudot, *Persistence Stability for Geometric Complexes*, Geometriae Dedicata 173:193–214, 2014。[出版社](https://link.springer.com/article/10.1007/s10711-013-9937-z)；DOI `10.1007/s10711-013-9937-z` | 几何复形的稳定性及紧支持的 q-tame 框架。本稿采用“边距的一半”参数，使 Rips 的 Hausdorff 稳定性常数为 1。 |
| `chazal2016structure` | F. Chazal, V. de Silva, M. Glisse, S. Oudot, *The Structure and Stability of Persistence Modules*, SpringerBriefs in Mathematics, 2016。[出版社](https://link.springer.com/book/10.1007/978-3-319-42545-0)；DOI `10.1007/978-3-319-42545-0` | 一参数 persistence module 的结构、interleaving 与 diagram 的等距及稳定性。 |
| `edelsbrunner2010computational` | H. Edelsbrunner, J. Harer, *Computational Topology: An Introduction*, AMS, 2010。[DOI / AMS](https://doi.org/10.1090/mbk/069) | Čech、Rips、nerve 与计算拓扑基础；本稿始终显式说明半径约定。 |
| `lesnick2015theory` | M. Lesnick, *The Theory of the Interleaving Distance on Multidimensional Persistence Modules*, FoCM 15(3):613–650, 2015。[出版社](https://link.springer.com/article/10.1007/s10208-015-9255-y)；[作者稿](https://arxiv.org/abs/1106.5305)；DOI `10.1007/s10208-015-9255-y` | `thm:canon` 调用原有普适性：素域、逐点有限维、一参数模块，以及对所有共同空间函数的普适稳定性。以正同调次数的几何提升得到最大性。 |
| `chazal2008towards` | F. Chazal, S. Y. Oudot, *Towards Persistence-Based Reconstruction in Euclidean Spaces*, SoCG 2008:232–241。[DOI / ACM](https://doi.org/10.1145/1377676.1377719)；[作者机构记录](https://researchportal.ip-paris.fr/en/publications/towards-persistence-based-reconstruction-in-euclidean-spaces/) | 持久像型几何重构。当前论文自行给出两个 Rips 尺度与严格覆盖、扰动常数的完整证明。 |
| `boissonnat2019reach` | J.-D. Boissonnat, A. Lieutier, M. Wintraecken, *The Reach, Metric Distortion, Geodesic Convexity and the Variation of Tangent Spaces*, JACT 3(1–2):29–58, 2019。[出版社](https://link.springer.com/article/10.1007/s41468-019-00029-8)；DOI `10.1007/s41468-019-00029-8` | reach 控制内外距离，用于覆盖球的体积上界。 |
| `petersen2016riemannian` | P. Petersen, *Riemannian Geometry*, 第三版，GTM 171, Springer, 2016。[出版社](https://link.springer.com/book/10.1007/978-3-319-26654-1)；DOI `10.1007/978-3-319-26654-1` | Bishop 体积比较；与 Gauss 方程和 reach 距离界一起支撑 `lem:cover` 中安全常数 $1.1$。 |

这些结果构成本稿的数学基础。覆盖数的码长解释、具体噪声参考和局部同步源类在本稿中另外定义；经典稳定性或同调恢复定理本身不直接给出所有通信任务的最小位数。

## 3. 概率、率失真与 Reeb 商

| BibTeX key | 经核验的原始文献与主来源 | 本稿引用内容及条件 |
|---|---|---|
| `cover2006elements` | T. M. Cover, J. A. Thomas, *Elements of Information Theory*, 第二版，Wiley, 2006。[DOI / Wiley](https://doi.org/10.1002/047174882X) | 数据处理、Fano、不确定性分解及源编码。概念率依赖概念概率，整块与逐符号可靠性各有定义。 |
| `laurent2000adaptive` | B. Laurent, P. Massart, *Adaptive Estimation of a Quadratic Functional by Model Selection*, Annals of Statistics 28(5):1302–1338, 2000。[Project Euclid](https://projecteuclid.org/journals/annals-of-statistics/volume-28/issue-5/Adaptive-estimation-of-a-quadratic-functional-by-model-selection/10.1214/aos/1015957395.full)；DOI `10.1214/aos/1015957395` | 原文 Lemma 1 的卡方两侧尾界。跨项另用 Gaussian 尾界，两者经 union bound 给出有限样本的平方距离误差。 |
| `linder1994asymptotic` | T. Linder, R. Zamir, *On the Asymptotic Tightness of the Shannon Lower Bound*, IEEE TIT 40(6):2026–2031, 1994。[DOI / IEEE](https://doi.org/10.1109/18.340474) | 高分辨率标量率失真的经典依据；本稿明确分支条件分数源及平方失真。 |
| `koch2016shannon` | T. Koch, *The Shannon Lower Bound Is Asymptotically Tight*, IEEE TIT 62(11):6155–6161, 2016。[作者全文](https://www.tsc.uc3m.es/~koch/files/TIT.2016.2604254.pdf)；[作者机构记录](https://researchportal.uc3m.es/display/act453052)；DOI `10.1109/TIT.2016.2604254` | 有限微分熵且整数部分熵有限时，norm-based difference distortion 的 Shannon 下界渐近紧；本稿的有界分数支撑满足整数部分条件，允许可积的临界值密度奇异性。 |
| `dey2013reeb` | T. K. Dey, Y. Wang, *Reeb Graphs: Approximation and Persistence*, DCG 49(1):46–73, 2013。[出版社](https://link.springer.com/article/10.1007/s00454-012-9463-z)；DOI `10.1007/s00454-012-9463-z` | Reeb 图、近似和 persistence 的背景。主文把标量分数和带边身份的商变量区分开。 |
| `kaluba2015reeb` | M. Kaluba, W. Marzantowicz, N. Silva, *On Representation of the Reeb Graph as a Sub-Complex of Manifold*, TMNA 45(1):287–307, 2015。[原始期刊 PDF](https://www.tmna.ncu.pl/static/published/2015/v45n1-15.pdf)；[Project Euclid / DOI](https://projecteuclid.org/journals/topological-methods-in-nonlinear-analysis/volume-45/issue-1/On-representation-of-the-Reeb-graph-as-a-sub-complex/10.12775/TMNA.2015.015.full)；DOI `10.12775/TMNA.2015.015` | 原文构造的嵌入复形与 Reeb 投影复合诱导基本群同构，故投影诱导满射。原文假设闭、连通、维度至少二及孤立临界点；本稿逐连通分量应用，维度一另由商为同胚处理。 |

Kaluba 的页码以原始 21 页 PDF 和 Project Euclid 的 287–307 为准；个别自动引文导出中的不同结束页不覆盖原始刊版。Koch 的年份使用正式 IEEE 刊版 2016，未把作者稿日期当作出版年。

## 4. 直接近邻：diagram 压缩及高维噪声

| BibTeX key | 核验后的正式记录与主来源 | 已有研究对象 | 本稿中可核实的差异 |
|---|---|---|---|
| `divol2021quantization` | V. Divol, T. Lacombe, *Estimation and Quantization of Expected Persistence Diagrams*, ICML 2021, PMLR 139:2760–2770。[官方论文页及 BibTeX](https://proceedings.mlr.press/v139/divol21a.html) | 随机 diagram 的 expected persistence diagram 及其有限支撑量化，含估计和近似的收敛结果。 | 本稿的同步源类是分离参考附近的有限 diagram，资源是保证全部源元素达到误差要求的固定消息长度。 |
| `sheehy2021sketching` | D. R. Sheehy, S. Sheth, *Sketching Persistence Diagrams*, SoCG 2021, LIPIcs 189:57:1–57:15。[官方页及 BibTeX](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SoCG.2021.57)；DOI `10.4230/LIPIcs.SoCG.2021.57` | 按不同加权点数构造一系列 diagram sketches，给出 2-approximation、紧凑存储及距离计算收益。 | 本稿对已知局部侧信息的端点格子进行联合整数编码，证明整条报文的精确二进制长度。 |
| `perumalsamy2026information` | D. Perumalsamy, C. Gunalan, R. Thinakaran, *An Information Theory of Persistent Homology: Entropy, the Data Processing Inequality, and Rate–Distortion Bounds for Topological Features*, Mathematics 14(8):1385, 2026。[出版社](https://www.mdpi.com/2227-7390/14/8/1385)；DOI `10.3390/math14081385` | 直接研究 persistent homology 的信息论描述，包括 diagram 层面的熵、互信息和率失真。 | 本稿给出特定接收机与源类：有限定长的同步码、加噪特征滤流和具有分支错误代价的任务源。元数据和摘要支持此定位，不由本记录宣称重证该文全部定理。 |
| `damrich2024persistent` | S. Damrich, P. Berens, D. Kobak, *Persistent Homology for High-dimensional Data Based on Spectral Methods*, NeurIPS 37:41954–42014, 2024。[官方论文页](https://proceedings.neurips.cc/paper_files/paper/2024/hash/4a32a646254d2e37fc74a38d65796552-Abstract-Conference.html)；[全文](https://proceedings.neurips.cc/paper_files/paper/2024/file/4a32a646254d2e37fc74a38d65796552-Paper-Conference.pdf)；DOI `10.52202/079017-1328` | 高维噪声下的拓扑恢复；利用 kNN 图的 diffusion 和 effective-resistance 距离改善鲁棒性。 | 本稿分析明确的各向同性特征信道，保留干净距离构造参考滤流，给出尺度校准、有限样本误差和指定阈值图的接触界。 |
| `hiraoka2026curse` | Y. Hiraoka, Y. Imoto, S. Kanazawa, E. Liu, *Curse of Dimensionality on Persistence Diagrams*, Foundations of Data Science 8:236–286, 2026。[出版社](https://www.aimsciences.org/article/doi/10.3934/fods.2025008)；[作者稿](https://arxiv.org/abs/2404.18194)；DOI `10.3934/fods.2025008` | HDLSS 条件下 persistence diagrams 的不一致与 normalized PCA 的理论改善。 | 本稿保留有限 $n,D_a,\delta$ 参数，并明确总噪声能量、重复点、接触与分区事件、整块可靠性比较。 |

### 重要版本与模型区别

- **NeurIPS 卷号。** Damrich 的官方 proceedings 和 BibTeX 为卷 **37**。其他论文对它的引用不是其卷号的最终依据。
- **Hiraoka 年份。** 出版社记载 2025-06-10 early access；正式卷为 **2026 年卷 8**，页码 236–286，没有列 issue，因此书目不补造 issue。
- **噪声归一化。** 本稿 $\sigma^2$ 为特征向量总噪声能量。比较高维结果时，先把已有工作的每坐标噪声方差转换到相同约定；固定每坐标方差与固定总能量是不同的联合极限。
- **普适性与近邻定位。** 稳定性、bottleneck 最大性、diagram 量化和高维失效均已有文献。当前贡献通过源类、消息、接收机和可靠性事件的具体区别表达，见 [related-work.md](related-work.md)。

## 5. 通信动机、替代方法与软件

| BibTeX key | 核验的记录与主来源 | 在 Paper 2 中的作用 |
|---|---|---|
| `xie2021deepsc` | H. Xie, Z. Qin, G. Y. Li, B.-H. Juang, *Deep Learning Enabled Semantic Communication Systems*, IEEE TSP 69:2663–2675, 2021。[DOI / IEEE](https://doi.org/10.1109/TSP.2021.3071210) | 学习式语义通信系统的背景。 |
| `gunduz2023beyond` | D. Gündüz et al., *Beyond Transmitting Bits: Context, Semantics, and Task-Oriented Communications*, IEEE JSAC 41(1):5–41, 2023。[DOI / IEEE](https://doi.org/10.1109/JSAC.2022.3223408) | 用上下文、语义及任务规定通信目标的背景综述。 |
| `shao2022taskoriented` | J. Shao, Y. Mao, J. Zhang, *Learning Task-Oriented Communication for Edge Inference: An Information Bottleneck Approach*, IEEE JSAC 40(1):197–211, 2022。[DOI / IEEE](https://doi.org/10.1109/JSAC.2021.3126087) | 任务导向学习与 information bottleneck 动机；本稿另行定义标量与 Reeb 分支任务的失真。 |
| `wang2023tclsc` | Y. Wang, S. Guo, *Transceiver Cooperative Learning-aided Semantic Communications Against Mismatched Background Knowledge Bases*, arXiv:2301.03133, 2023。[作者发布的 arXiv 记录](https://arxiv.org/abs/2301.03133) | KB 不一致及协作参数交换的工程动机。引用的出版类型为预印本。本稿的失配量受固定解释规则约束，额外配准会改变模型。 |
| `moor2020topological` | M. Moor, M. Horn, B. Rieck, K. Borgwardt, *Topological Autoencoders*, ICML 2020, PMLR 119:7045–7054。[官方论文页](https://proceedings.mlr.press/v119/moor20a.html) | 表示学习中的拓扑损失与结构保持动机。 |
| `erak2026topojscc` | O. Erak, O. Alhussein, F. Fang, S. Muhaidat, *Topology-Preserving Deep Joint Source-Channel Coding for Semantic Communication*, arXiv:2603.17126, 2026。[作者发布的 arXiv 记录](https://arxiv.org/abs/2603.17126) | 学习式 TopoJSCC 使用图像与潜特征 persistence 正则项。当前公开记录为预印本，未据“submitted”声称期刊录用。本稿的指定特征信道结果不作为该神经编码器的最优性极限。 |
| `chazal2011geometric` | F. Chazal, D. Cohen-Steiner, Q. Mérigot, *Geometric Inference for Probability Measures*, FoCM 11(6):733–751, 2011。[出版社](https://link.springer.com/article/10.1007/s10208-011-9098-0)；DOI `10.1007/s10208-011-9098-0` | distance-to-measure 等使用分布信息的几何推断方法；说明不同接收机可以利用不同信息。 |
| `bauer2021ripser` | U. Bauer, *Ripser: Efficient Computation of Vietoris–Rips Persistence Barcodes*, JACT 5(3):391–423, 2021。[出版社](https://link.springer.com/article/10.1007/s41468-021-00071-5)；DOI `10.1007/s41468-021-00071-5` | Rips persistence 的计算来源。实现输出须转换为本文半径约定。 |
| `pedregosa2011scikit` | F. Pedregosa et al., *Scikit-learn: Machine Learning in Python*, JMLR 12(85):2825–2830, 2011。[JMLR 官方页](https://jmlr.org/papers/v12/pedregosa11a.html) | 已运行的有限数据实验的软件来源；软件引用与正 reach 等理论假设分别核验。 |

## 6. 可复核的内容审计结论

| 主稿论断 | 来源能支持的内容 | 本稿已明确的边界或新增证明 |
|---|---|---|
| bottleneck 是普适稳定失真的最大者 | Lesnick 的几何提升、interleaving 普适性；一参数等距 | 主稿注明素域、逐点有限维及普适函数稳定性；没有把这个既有定理列为原创贡献。 |
| 字典分辨率认证同调 | reach 下恢复、持久像推断、覆盖数体积界 | 主文与 Appendix A 证明指定 Rips 尺度、严格噪声预算、均匀重构意义下的精确索引率，以及没有覆盖承诺时的索引不识别性。 |
| 高维噪声影响 persistence | Damrich、Hiraoka 已分析失效与缓解 | 主文与 Appendix B 给出具体噪声参考、有限样本边误差和保持跨对依赖的接触 bracket。 |
| diagram 压缩与率失真 | Divol、Sheehy、Perumalsamy 已覆盖相邻对象 | 本稿局部同步定理规定侧信息、可独立更新的端点源类，并构造达到精确最坏情况长度的固定消息。 |
| KB 差异形成极限 | persistence 稳定性可估计几何差异；Wang–Guo 说明 KB 不一致的应用背景 | 在固定解释规则下证明双侧有限分辨率区间及稠密字典极限，允许有限字典的误差低于极限。 |
| 分支任务有额外信息成本 | Reeb 商结构；标量 Shannon 下界渐近紧 | Appendix C 处理有限跨边错误代价，保留常数项并得到 $H(E\mid T)$；diagram 本身不编码所有 Reeb 分支身份与概率。 |

## 7. 引用检查的复现方法

在仓库根目录运行以下只读检查可重数引用 key，并验证每个引用存在于书目中：

```python
from pathlib import Path
import re
from collections import Counter

bib = Path('refs/references.bib').read_text()
keys = re.findall(r'@\w+\{([^,]+),', bib)
assert all(n == 1 for n in Counter(keys).values())
known = set(keys)
print('bibliography entries:', len(keys))

for directory in ('paper1-curvature', 'paper2-topology'):
    cited = set()
    for path in Path(directory).rglob('*.tex'):
        text = re.sub(r'(?<!\\)%.*', '', path.read_text())
        pattern = r'\\cite\w*\*?(?:\[[^\]]*\])*\{([^}]+)\}'
        for match in re.finditer(pattern, text):
            cited.update(key.strip() for key in match.group(1).split(','))
    missing = sorted(cited - known)
    print(directory, len(cited), 'missing:', missing)
    assert not missing
```

文献用途分别列于上表。当前主稿与各推导共享工作目录中的 ref.bib。

Paper 1 的已完成来源核验独立保留在 main 提交 [`8aad9ea2f51cabb1cd2476dd633582aa9c736e5d`](https://github.com/echho01/semantic-manifold/tree/8aad9ea2f51cabb1cd2476dd633582aa9c736e5d) 的 [refs/verification.md](https://github.com/echho01/semantic-manifold/blob/8aad9ea2f51cabb1cd2476dd633582aa9c736e5d/refs/verification.md)。
