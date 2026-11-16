# 工作三研究索引 · 保持位置连通关系的拓扑学习

[Conference 主稿](../main.tex) · [论文 PDF](../build/main.pdf) · [书目](../ref.bib) · [信道、对照方法与实验矩阵](EXPERIMENTS.md)

接收机重构普通图像或掩膜，评估器在相同位置询问能否连通。一次传输可以支持多项询问，估计所需资源时要共同计算这些答案的信息量。

| 数学推导 | 与当前工作的关系 |
|---|---|
| [连通风险、MST、双模板精确 RDF](derivations/connectivity.md) | 从固定端点和阈值分布推导接收误判率、树上的梯度、镜像实验的全对计数 |
| [多询问的信息量下界](derivations/query_converse.md) | 条件化询问身份，按一张图像的一次传输计费 |
| [对称轨道率失真与衰落下界](derivations/orbit_rdf_and_fading.md) | 持久图在读出图自同构下不变、轨道信源的关系率失真闭式解、信道使用次数与遍历衰落和中断形式的下界 |
| [覆盖、可识别性与有限字典](derivations/coverage_and_dictionary.md) | **origin/topology** 分支（fe42dd32）的 reach、Rips 滤过约定与字典恢复前提 |
| [点云噪声与分区恢复](derivations/noise_and_partition.md) | 同一分支的重参数化、去偏、分量碎裂和跨分量接触 |
| [字典同步、知识失配与 Reeb 分支](derivations/synchronization_and_reeb.md) | 同一分支的有限覆盖误差、最短更新和任务分支条件熵 |
| [裕量、描述符精度与资源](derivations/coding_certificates.md) | main 分支（8aad9ea2）的结构证书、描述符覆盖与打包、精度分配 |
| [双失真与功率变分](derivations/structural_variation.md) | 早期模型的乘积源特例、无线可靠性与充分功率条件 |

| 实验 | 入口与结果 |
|---|---|
| 320 张真实 Omniglot 字形和镜像，按两种端点分布精确计数 | [实现](experiments/reflection_queries.py) · [逐图数据](results/reflection_queries.json) |
| 同一批 320 张字形的 8 元二面体轨道，闭式码率与 BA 核对 | [实现](experiments/orbit_rdf.py) · [核验](tests/test_orbit.py) · [逐图数据与信道公式求值](results/orbit_rdf.json) |
| digits 成对表示、有限字典与知识失配 | [实现](code/embeddings/run_topology.py) · [digits 记录](code/results/topology/digits_summary.json) · [字典、信道和同步实验说明](notes/dictionary_experiments.md) |
| 环与双簇特例、350 对试验的无线分配 | [复现](experiments/reproduce_early.py) · [功率、误判与置信区间](code/results/early_topology/reproduced.json) |
| 五种图像编码目标、AWGN、相干块衰落与信息下界 | [训练和测试](experiments/topology_jscc.py) · [拓扑损失](experiments/structure.py) · [信息界估计](experiments/query_bound.py) · [短程运行范围](results/smoke_provenance.json) |

以下从本目录执行。研究依赖包括 PyTorch、SciPy 和 GUDHI。镜像实验只需 NumPy、Pillow；轨道实验和 `tests/test_orbit.py` 另需 SciPy；重绘论文图另需 Matplotlib。
~~~sh
python -m pip install -r requirements.txt
python tests/test_connectivity.py
python tests/test_orbit.py
python experiments/fetch_omniglot.py --data data/omniglot
python experiments/reflection_queries.py --data data/omniglot --output runs/reflection.json --figure ../figures/reflection.png
python experiments/orbit_rdf.py --data data/omniglot --reference results/reflection_queries.json --output runs/orbit_rdf.json --figure ../figures/orbit_rdf.png
python experiments/topology_jscc.py train --data data/omniglot --dataset omniglot --device cuda --mode proposed --rho .4 --groups 4 --power uniform --checkpoint runs/topo.pt
python experiments/topology_jscc.py test --data data/omniglot --device cuda --split test --snr 0 5 10 15 20 --checkpoint runs/topo.pt --output runs/test.json
~~~

比较真实无线图像性能的多种子训练仍待完成。镜像、有限字典与已知原型信道各有自己的评价对象和资源单位。
