# 拓扑：论文借鉴、数据与判据

| 已核对的材料 | 采用之处 | 区别及应报告的对照 |
|---|---|---|
| [TopoJSCC，IEEE Communications Letters](https://doi.org/10.1109/LCOMM.2026.3717806)；[预印本](https://arxiv.org/abs/2603.17126) | 单图→有限个归一化复符号→AWGN/块 Rayleigh→重构；Omniglot 和 DeepGlobe 道路掩膜；图像 cubical PH、批次收发潜变量 Rips PH、H0/H1 的 W2 | 实现了独立的 topojscc 对照。采用连续临界值、可调 PH 子批量、显式瓶颈、衰落接收机 MMSE 均衡，默认混合 SNR 训练。需要复核原表时分信道/SNR 独立训练，报告差异；不复制论文图表数字 |
| [TopoAE，ICML 2020](https://proceedings.mlr.press/v119/moor20a.html) | 对参考与接收空间分别选 MST 边，再比较两边相同端点的边长 | 对照 topoae 使用批次 H0、归一化尺度和源图像/接收码对应；不是原作者网络；不把随机批次中的环说成单幅图像的环 |
| [MALIS，NeurIPS 2009](https://arxiv.org/abs/0911.5372) | 用路径瓶颈表达有身份的点对连通，最小生成树定位临界边/点 | 该图论原语已有先例。本项目优化有噪通信后的阈值积分询问误判，并测等带宽下的代价；固定询问字典给出资源下界 |
| [Topological Flow Matching，ICLR 2026](https://proceedings.iclr.cc/paper_files/paper/2026/hash/238e6167319ad5da33c7b38594a1edb1-Abstract-Conference.html) | 固定图或复形上的 Laplacian 可以编码已知领域结构 | 留作已知网格先验的扩展；没有把流匹配当作本项目已实现的算法 |
| [TopoCode](https://arxiv.org/abs/2411.12825)、[clDice，CVPR 2021](https://openaccess.thecvf.com/content/CVPR2021/html/Shit_clDice_-_A_Novel_Topology-Preserving_Loss_Function_for_Tubular_Structure_CVPR_2021_paper.html) | 前者用于评估额外结构信息的代价；后者提供道路中心线保真的感知对照 | 当前主矩阵不实现 BPG+LDPC、带边信息的 TopoCode 或 clDice；列作扩展基线，不能称为已复现结果 |

指定的 [TopoJSCC 本地全文](../../refs/topology/2026_CL_Topology_Preserving_Deep_Joint_Source_Channel_Coding_for_Semantic_Communication.pdf) 与[拓扑文献目录](../../refs/topology)保留在公共资料库。

**数据。** [Omniglot 官方源](https://github.com/brendenlake/omniglot)：按提交和 SHA-256 下载，黑笔画转亮前景，整字母表划分。排序后每第 5 个 background 字母表用于验证，其余训练，evaluation 字母表只测试。[DeepGlobe 原始挑战](https://openaccess.thecvf.com/content_cvpr_2018_workshops/w4/html/Demir_DeepGlobe_2018_A_CVPR_2018_paper.html)：用带标签的完整道路二值掩膜。代码按原始 tile ID 的稳定哈希划分 70/15/15，禁止先切 patch，拒绝重复 ID；这是自行留出的测试集，不是官方盲测榜单。默认 32×32 为成本可控的入口；道路正式结论需在更高分辨率及同一原图分组裁块上复查细路是否被缩没，含地理位置时优先按地理区块留出。裁块和地理分组方案属于后续扩展。此次未取得或运行 DeepGlobe 全量训练。

**问题、单位与系统。** 同一图像的四邻接顶点滤过，回答两点是否在同一阈值下连通；顶点身份是固定坐标。选择无穷范数地面距离的 W2、有限 H0/H1 图，并可与该查询发生分歧。训练时 query 只需最小生成树瓶颈，proposed 再比较源/接收表示间的临界边；通道回传梯度。对既定复使用次数 k 和平均功率 1 学习，`--groups 4 --power learned` 要和相同分组的 `--power uniform` 比，必须同时保留 MSE+可学功率对照。固定高度/宽度下实际 CBR=k/(H×W)，以四舍五入后的 k 计；这里的 SNR 是复符号 Es/N0，不是工作一、二的 Eb/N0。

**矩阵与判断。** 固定 rho=.4 扫 0/5/10/15/20 dB，两信道均测；15 dB 另扫 rho=.05/.10/.20/.30/.40/.50。基线 mse / topojscc（可分别置零 `--lambda-img` 或 `--lambda-lat`）/ topoae / query / proposed。默认训练阶段共享信道/SNR 分布；如需匹配原论文带宽曲线，可显式 `--channels awgn --snr 15` 训练再评估。正式矩阵至少三个训练种子、同一测试图像和噪声种子；在验证集挑选权重，在相同 CBR、相近 PSNR 下比较询问风险。指标：全图像素 MSE 的 PSNR、图像 H0/H1 W2、混合端点风险和均匀端点风险、实际 k、分组功率。记录小批量 PH 的规模、耗时、显存；小样本请不要把不同噪声重复当成独立训练。`*_paired.jsonl.gz` 保留图像、条件、重复和逐图指标，可按字母表或道路原图组做配对自助区间。如果只在训练端点分布获益，或在相同 PSNR 下没有改进，应如实限定结论。

**把“误判率”换算为资源下界。** [推导](derivations/query_converse.md)与 [query_bound.py](experiments/query_bound.py)对应。为所有测试图使用同一批事先选好的询问，联合计算回答向量熵，扣除各询问二元误差熵，给出 bit/图下界和复 AWGN 使用次数下界。可传 `--dictionary queries.json`，格式为二维数组，例如 16×16 图像上的 `[[34,47,0.35],[114,117,0.6]]`。需在校验集确定询问后冻结。有限样本置信下界假定独立图像，字母表内相关时仅作描述或另做聚类留出；小样本为 0 是有效但不紧的下界。一般字形集合上的短验证只得到零下界，记录已删除；非平凡的验证在对称轨道信源上完成（见后文）。

**神经编解码器的正式流程。** 见后文“神经编解码器的收敛实验”和 [final_protocol.sh](experiments/final_protocol.sh)。早期 12 张图、1 轮的短验证记录带有 Sigmoid 输出层饱和的缺陷，已于 2026-11-19 连同 `matrix.json`、`run_matrix.py` 从仓库删除，可在 Git 历史的初始提交中找到。

~~~sh
python experiments/fetch_omniglot.py --data data/omniglot
bash experiments/final_protocol.sh
python experiments/summarize_final.py runs/final
~~~

## 镜像对照的已测结果

从 Omniglot 的 20 个官方 evaluation 字母表各抽 16 张，种子 20261115，最近邻缩放为 32×32，黑色笔画阈值 128。逐图统计两套四邻接分量的交叉计数，枚举意义下包含全部有放回的有序点对。镜像保持 cubical 滤过逐级同构。均匀端点的逐图平均误判为 0.012096，前景端点为 0.778214。误判在全部 320 张图上均为正。这是坐标身份的受控检验，与无线重构性能评测分开记录。[精确推导](derivations/connectivity.md) · [脚本](experiments/reflection_queries.py) · [逐图原始记录](results/reflection_queries.json)

~~~sh
python experiments/reflection_queries.py --data data/omniglot --output runs/reflection.json --figure ../figures/reflection.png
~~~

## 对称轨道与信道下界的补充结果（2026-11-16）

主稿目标为 IEEE GLOBECOM 2027（会议双栏，含参考文献 6 页）。新增内容不涉及神经网络训练。

- 沿用镜像实验的 320 张字形和种子，ID 哈希一致，镜像误判率与[既有记录](results/reflection_queries.json)逐图相同。全部字形在方形网格的 8 个对称下稳定子平凡。按[轨道率失真定理](derivations/orbit_rdf_and_fading.md)计算：全像素询问下 $\delta_{\max}$ 均值 0.011541，在 $\delta/\delta_{\max}=0,0.1,0.25,0.5$ 处平均需要 3.0000、2.3124、1.6114、0.7698 bit/字形，而任何持久图失真的码率为零。[逐图记录与曲线](results/orbit_rdf.json)
- 200 个 $5\times5$ 随机模板上闭式解与 Blahut–Arimoto 解最大相差 $8.9\times10^{-16}$ bit。[单元测试](tests/test_orbit.py)只依赖 NumPy 和 SciPy。
- 信道使用次数下界、Rayleigh 遍历形式与中断形式是公式求值。它们给出任何编解码器都需满足的必要条件，不代表已训练模型的性能。
- 轨道定理把重构限定在轨道内并要求询问分布与信源独立。按前景采样的询问与轨道外的重构不在该定理范围内。

~~~sh
python tests/test_orbit.py
python experiments/orbit_rdf.py --data data/omniglot --reference results/reflection_queries.json --output runs/orbit_rdf.json --figure ../figures/orbit_rdf.png
~~~

## 关系泛函的一阶变分与变分分配（2026-11-16）

主稿改以关系泛函及其一阶变分为主线，[推导](derivations/relation_variation.md)。新增计算不涉及神经网络训练，依赖 NumPy、SciPy、Pillow（绘图另需 Matplotlib），全量约一分钟。

- **合成核对。** 100 个 $6\times6$ 随机单射场上，并查集载荷公式与逐对枚举一致，$\|\delta\|_\infty<\gamma_f/2$ 时 $\mathcal R[f+\delta;f]=\sum_vc_f(v)|\delta(v)|$ 的误差在 $10^{-16}$ 量级；超过半间隔后公式可以失效。[测试](tests/test_variation.py)另验证带并列值的方向导数、二值翻转余项界和三种分配的最优性。
- **影响图（真实 Omniglot，320 张，与镜像实验同一批，逐图镜像误判率差为 0）。** 逐像素翻转并精确计数。总影响均值 0.08265、中位数 0.05887；前景占像素 8.43%，承担影响的 70.83%；影响最大的 1% 像素占 30.94%，前 5% 占 65.59%。
- **一阶律。** 全部像素以概率 $\varepsilon$ 独立翻转，实测风险与 $\varepsilon\sum_v\mathrm{Inf}_x(v)$ 之比在 $\varepsilon=10^{-4},10^{-3},10^{-2}$ 为 0.955、0.992、0.951。前两者的配对差区间含零，第三者不含零。最坏情形余项 $S^2/2$ 在 $32\times32$ 图上没有信息量。
- **未编码 BPSK 的变分功率。** 平均每像素能量相同，每图 100 次公共噪声。4、6、8 dB 处关系风险降到均匀功率的 4.3%、3.6%、4.5%，像素误码率升到 2.33、2.29、2.25 倍。发端用自己的图像计算分配，接收端符号判决不需要知道分配。这是对读出的受控实验，不构成对神经编解码器的排名。

主稿为腾出篇幅删去的内容仍保留在研究目录：digits 字典表（保留首尾两个数）和分支熵诊断、轨道码率曲线图 [orbit_rdf.png](../figures/orbit_rdf.png)、双失真源编码核的公式、两模式实验的充分半径、询问逆定理短程核对的零下界说明。

~~~sh
python tests/test_variation.py
python experiments/relation_variation.py --data data/omniglot --reference results/reflection_queries.json --output runs/relation_variation.json --maps runs/influence_maps.npz --figure ../figures/influence.png
~~~

## 针对外部审阅意见的修补与神经编解码器的收敛实验（2026-11-18）

| 审阅意见 | 判断 | 处理 |
|---|---|---|
| 定理 3 证明有“重要推导错误” | 不成立，审阅意见自己逐步检查后也确认证明严密 | 不改 |
| 命题 4(a) 的 $c^{2/3}$ 模拟分配在实际信噪比下超出适用域 | 成立 | 实测后删除该规则（下）；同时删去没有图像实验支撑的方差型规则，命题只保留 BPSK 分配 |
| BPSK 分配依赖逐图计算，功率图泄露结构 | 事实成立，不是错误 | 新增与图像无关的共享功率图对照（下） |
| 定理 3 的真实数据下界为 0 | 成立 | 改在对称轨道信源上验证，得到非平凡下界（下） |
| 神经编解码器没有收敛结果，系统模型与实验脱节 | 成立 | 修复训练缺陷后跑到收敛并做三种子配对评估（下）；系统模型改为一条一般链路的两个实例（逐位置链路、学习编解码器） |
| 关系损失每图只采 32 对，梯度稀疏 | 是方差问题 | 微调时用 32 张图、128 对端点；更宽的设置（128 张图）结果相近 |

**训练缺陷。** 原解码器输出层是 Sigmoid。稀疏字形加像素 MSE 下，输出在约 20 步内被推到 0，梯度降到 $10^{-21}$，此后不再学习。此前五种方法 PSNR 都是 7.14 dB、四组学习率和宽度得到完全相同的 MSE，都是这个原因。现改为线性输出，评估时截断到 $[0,1]$（`--activation linear`，旧检查点按 sigmoid 读取）。修复后 MSE 基线在验证集 10 dB 下 PSNR 约 28 dB。

**三项受控对照**（`experiments/allocation_controls.py`，[数值](results/allocation_controls.json)）：

- 模拟规则：40 个 $16\times16$ 灰度字形场，精确全点对风险，公共噪声。分配/均匀风险比在噪声标准差 $10^{-6}$ 时为 0.76（线性化预测 0.73），在 $10^{-3}$、$10^{-2}$、0.03、0.1 时为 1.84、1.08、1.06、1.12，0.3 时为 0.96。噪声一旦打乱排序，该规则比均匀功率更差。
- 共享功率图：用 200 张训练字母表字形的平均影响图给所有测试字形分配功率，4/6/8 dB 下风险为均匀功率的 46%、40%、35%；逐图计算为 4.4%、3.7%、3.1%。
- 轨道信源上的逆定理：320 个字形各取 8 个对称像作等概信源，贪心选 3 个固定询问，$H_Q$ 平均 2.48 bit。3 个 BPSK 符号在 0/3/6/10 dB 下，下界为 0.98、1.83、2.37、2.48 bit，链路实际承载 2.16、2.74、2.97、3.00 bit，无一违反。

**神经编解码器。** GPU 服务器，[流程脚本](experiments/final_protocol.sh)，[汇总](results/neural_final/summary.json)。Omniglot 训练字母表 15,700 张，$32\times32$，$k=410$，hidden 64，四组等功率。每个种子先用 MSE 训练 100 轮（lr $3\times10^{-4}$），再从同一权重微调 20 轮（lr $10^{-4}$）得到四臂：MSE、关系项（$\lambda_Q=0.01$）、关系项加临界边项（$\lambda_M=0.02$）、TopoJSCC 式损失。三个种子，2,048 张官方测试字形，公共噪声，每图 64 对固定端点。关系风险相对 MSE 臂的逐图配对变化（%，95% 区间）：

| 信道 | MSE：PSNR / 风险 | 关系项：PSNR / 风险 | 变化 |
|---|---|---|---|
| AWGN 0 dB | 24.47 / 0.0550 | 24.43 / 0.0542 | −1.4 [−1.6, −1.3] |
| AWGN 10 dB | 27.76 / 0.0369 | 27.68 / 0.0365 | −1.1 [−1.3, −0.8] |
| AWGN 20 dB | 28.25 / 0.0348 | 28.17 / 0.0345 | −0.9 [−1.1, −0.6] |
| Rayleigh 0 dB | 20.89 / 0.0825 | 20.87 / 0.0807 | −2.1 [−2.3, −2.0] |
| Rayleigh 10 dB | 26.02 / 0.0435 | 25.97 / 0.0428 | −1.6 [−1.8, −1.4] |
| Rayleigh 20 dB | 27.92 / 0.0358 | 27.84 / 0.0355 | −1.0 [−1.3, −0.8] |

- 关系项的增益为 1–2%，三个种子中两个下降、一个基本不变，种子间差异与效应同量级。
- 加临界边项后变化为 −0.4% 到 +2.1%，PSNR 再低约 0.1 dB，没有一致增益。
- TopoJSCC 式损失在默认权重下与 MSE 相差不超过 0.1%。
- 从头训练时增大 $\lambda_Q$ 会同时降低 PSNR 和升高风险：验证集 10 dB 下，$\lambda_Q=0.01/0.03/0.1/0.3$ 的 PSNR 为 27.0/25.2/22.8/20.0 dB，风险为 0.038/0.045/0.056/0.064（MSE 为 28.0 dB、0.036）。

结论：在这个学习编解码器里，关系损失没有带来与逐位置功率分配可比的保护。主稿据此调整：系统模型不再以卷积编解码器为主体，摘要和结论如实报告这一结果，并把“如何把逐位置保护带进压缩表示”列为未解决问题。未做的部分：带宽比扫描、可学习功率分组、道路数据。

## 第二轮审阅后的修补（2026-11-19）

- **打乱位置对照。** 把逐图功率值随机置换到其他像素，4/6/8 dB 下关系风险为均匀功率的 2.45、2.72、2.84 倍（逐图分配为 0.043、0.035、0.031，共享图为 0.46、0.38、0.37；本次重跑的噪声抽样与上次不同，共享图数值略有变化）。增益来自功率放在哪里，不是来自功率不均匀本身。
- 功率分配命题限定为内点极小元。
- TopoJSCC 式对照的 cubical 复形把像素当顶维方格，对角像素相连；查询用四邻接。因篇幅限制，这一点只记录在此，主稿未写。该对照与像素损失无差别，未再加同邻接的 $H_0$ 对照。
- 已删除 `experiments/matrix.json`、`run_matrix.py`、`results/omniglot_smoke_*`、`query_bound_smoke.json` 和 `smoke_provenance.json`。
