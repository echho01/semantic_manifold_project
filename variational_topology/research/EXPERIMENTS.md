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

**把“误判率”换算为资源下界。** [推导](derivations/query_converse.md)与 [query_bound.py](experiments/query_bound.py)对应。为所有测试图使用同一批事先选好的询问，联合计算回答向量熵，扣除各询问二元误差熵，给出 bit/图下界和复 AWGN 使用次数下界。可传 `--dictionary queries.json`，格式为二维数组，例如 16×16 图像上的 `[[34,47,0.35],[114,117,0.6]]`。需在校验集确定询问后冻结。有限样本置信下界假定独立图像，字母表内相关时仅作描述或另做聚类留出；小样本为 0 是有效但不紧的下界。[短验证结果](results/query_bound_smoke.json)不支持非平凡的正式资源主张。

**已执行。** 已在官方 Omniglot 上对全部五种模式训练、保存、独立测试（每种 12 张训练、8 张留出、1 epoch、两信道两 SNR，16×16）。[新方法的原始记录](results/omniglot_smoke_proposed.json)、[TopoJSCC 对照](results/omniglot_smoke_topojscc.json)。这种规模用于检查数据、梯度、功率和评价流程，不排名。已有 topology 分支的 [digits 诊断](code/results/topology/digits_summary.json) 与新图像任务是两种不同证据。旧环/双簇 24 dB 两策略的 102/350 与 53/350 仅适用于原型、增益及解码器固定的那组模拟。

从本项目根目录运行；先打开生成的清单确认工作量。`plan` 只列配置；`train` 保存检查点；`validation` 与 `test` 各自读取已训练的检查点，写独立 JSON。复制 `experiments/matrix.json` 后可用 `--spec my-matrix.json` 只跑选定设置。全面板作业较多，请先跑 pilot；不要用正式测试集挑超参数。
~~~sh
python experiments/fetch_omniglot.py --data data/omniglot
python experiments/run_matrix.py plan --profile full --data data/omniglot --out runs/matrix
python experiments/run_matrix.py train --profile pilot --data data/omniglot --out runs/matrix
python experiments/run_matrix.py validation --profile pilot --data data/omniglot --out runs/matrix
python experiments/run_matrix.py test --profile pilot --data data/omniglot --out runs/matrix
~~~
道路数据正式预案：复制 matrix.json，将 full.train.dataset 改为 roads，`--data` 指向掩膜根目录、`--out` 指向单独的 roads 结果目录。

已有结果做配对检验的命令（差值方向是第一个减第二个，区间以这两份固定检查点为条件）：

~~~sh
python experiments/compare.py results/omniglot_smoke_proposed.json results/omniglot_smoke_topojscc.json --channel awgn --snr 0 --metric hybrid_query
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
