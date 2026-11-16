# 拓扑字典、信道与分支实验记录

[研究索引](../README.md) · [论文实验协议](../EXPERIMENTS.md)


本页对应 `topology` 分支的当前论文及实际执行记录。完整结果保存在
[`code/results/topology/results.json`](../code/results/topology/results.json)，
四张论文图由同一次运行的测量值生成。以下数值来自这些记录，预训练模型实验单独列为扩展协议。

## 1. 完整 CPU 复现

在仓库根目录安装依赖：

```bash
python -m pip install -r requirements.txt
```

随后用一个 shell 命令执行跨库核验、14 项回归、全部六阶段实验和绘图：

```bash
python code/tests/run_topology_checks.py && \
python code/synthetic/run_topology_suite.py \
  --output code/results/topology --reparam-trials 20 --channel-trials 120
```

此流程不需要 GPU、预训练模型或外部数据下载。Digits 数据由 scikit-learn 随包提供。
`bash code/synthetic/run_paper2.sh` 也可运行全部实验，脚本先执行 TDA 跨库自检；
完整的 14 项回归由上面的 `run_topology_checks.py` 执行。
单阶段可用 `--stages dictionary,reparam` 等逗号列表；单阶段运行只把实际执行的阶段写入
`results.json`，完整投稿复现应使用全部六阶段的默认列表。

当前记录使用 Python 3.12.14、NumPy 2.3.5、SciPy 1.17.0、Ripser 0.6.15、
Persim 0.4.1、GUDHI 3.13.0、scikit-learn 1.8.0、Matplotlib 3.10.8。
版本、配置、运行时间及关键源文件 SHA-256 均由驱动自动写入 `results.json`。

### 数值约定

- 所有 Rips 半径均为最大边长的一半；`tda.rips` 已转换 Ripser 的输出，调用方不应再除以二。
- `sigma` 是总噪声幅度：每坐标方差为 `sigma**2 / D_a`，每个向量的期望噪声能量为 `sigma**2`。
- 高维合成实验在每个重复中固定一组标准噪声，贯穿同一维度的整个 `sigma` 扫描。
- 重参数化图的误差条为均值加减 `1.96 × Monte Carlo 标准误`；事件概率报告 95% Wilson 区间。
  置信度针对单个固定条件，没有宣称整张参数扫描同时覆盖。
- 概率定理允许指定的小概率失败；“本次未观测到超界”是运行结果，不等于概率为一的保证。

## 2. 六阶段设置与实际结果

| 阶段 | 本次规模 | 实际结果 | 原始记录 |
|---|---|---|---|
| `dictionary` | 两圆与 Clifford 环面，共 22 个规则字典 | 7 个满足覆盖证书的配置全部恢复目标 Betti 数 | [dictionary.csv](../code/results/topology/dictionary.csv) |
| `reparam` | 4 个环境维数 × 3 个噪声幅度 × 20 次重复；另 100 次校准测量 | 240 次重参数化、240 次去偏及 100 次校准误差均未超过相应界 | [reparameterization.csv](../code/results/topology/reparameterization.csv)、[calibration.csv](../code/results/topology/calibration.csv) |
| `channel` | 每概念 64 点，4 个维数，每维 120 条噪声路径 | 保存 36,000 行原始测量、480 条首次事件记录和 1,500 行概率汇总 | [channel_raw.csv](../code/results/topology/channel_raw.csv)、[channel_probabilities.csv](../code/results/topology/channel_probabilities.csv) |
| `mismatch` | 两种变形 × 7 个字典大小；500 次实际同步编码与解码 | 14 条曲线测量全部满足双侧带；500 次码长及重构检查全部通过 | [mismatch.csv](../code/results/topology/mismatch.csv)、[sync_codes.csv](../code/results/topology/sync_codes.csv) |
| `task` | 500,000 个表面积均匀环面样本，6 个量化精度；20,000 点诊断云 | 解析表达式的数值积分给出 0.3537667163 bit 的分支开销，最细经验量化给出 0.3541681590 bit | [reeb_quantization.csv](../code/results/topology/reeb_quantization.csv)、[task_summary.json](../code/results/topology/task_summary.json) |
| `digits` | 1,797 个图像；1,078/719 内部分层划分；64 维像素与 16 维 PCA | 对齐的 200 个留出条目上，有限参考图差异为 0.0555076748；10 个字典测量全部满足双侧带 | [digits_report.json](../code/results/topology/digits_report.json) |

### 2.1 字典覆盖与持久像恢复

两圆为半径 1、圆心 `(-1.5, 0)` 和 `(1.5, 0)` 的不交圆，整体 reach 为 0.5。
每圆取 `m` 个等角点，已知覆盖半径为 `2 sin(pi/(2m))`。
Clifford 环面为 `S1 × S1` 在四维中的标准嵌入，reach 为 1；`m × m` 网格的覆盖半径为
`2 sqrt(2) sin(pi/(2m))`。接收机使用 `eta = reach/6` 与 `(2 + sqrt(2)) eta` 两个尺度的持久像。

本次扫描中，两圆首次在 76 个原子、7 bit 索引处满足证书并恢复 `(beta0, beta1) = (2, 2)`；
环面首次在 729 个原子、10 bit 索引处满足证书并恢复 `(1, 2, 1)`。
这些是被扫描规则构造的结果，不是最优覆盖数或所有拓扑推断算法的最低码率。
依据见 [dictionary_summary.json](../code/results/topology/dictionary_summary.json)。

### 2.2 高维重参数化与噪声校准

干净云固定为三谐波闭曲线的 80 个等间隔参数点；同一个六维云以正交标架嵌入
`D_a = 64, 256, 1024, 4096`。总噪声幅度取 `0.5, 1, 2`，每种配置 20 次重复，
理论参数 `delta = 0.05`。测量取 `H0`、`H1` 两个 bottleneck 距离的最大值。

参考滤过直接由非对角距离 `sqrt(d_ij**2 + 2*sigma**2)` 构建，保留重复干净样本产生的全部
零合并半径。它解决了直接变换普通 `H0` 持久图会漏掉零长度条的问题。
当 `sigma = 1`、维数从 64 增至 4096 时，重参数化参考误差均值由 0.16409 降至 0.02377，
去偏误差由 0.21544 降至 0.09293；直接比较干净与含噪图的误差均值则为 0.37665 和 0.36358。
这些数值描述固定曲线和本次样本量上的效应。

校准实验固定 `D_a = 1024`、`sigma = 1`，对方差引入 `-20%, -10%, 0, 10%, 20%` 的相对误差。
其界使用论文中的 `beta + 2*v`，其中 `v = abs(estimated_variance - true_variance)`。
具体配置与记录见 [reparam_summary.json](../code/results/topology/reparam_summary.json)。

### 2.3 组件计数、跨接触、真实分区与全图坍缩

两概念各有 64 个重复中心样本，中心间隔 `Delta = 1`，接收阈值 `theta = 0.25`。
扫描上述四个维数。每个维数的噪声网格合并了 17 个归一化水平 `0.15..0.75`、
41 个步长为 `0.005` 的归一化水平 `0.25..0.45`、16 个绝对水平 `0.30..0.45`，
以及 `0, 0.25, sqrt(2)/4, 0.5`；合并后每个维数有 75 个噪声水平。
同一重复共享标准噪声，因而首次事件不会因逐格重新抽样而人为改变。

| 记录字段 | 实际含义 |
|---|---|
| `raw_count` / `debiased_count` | 阈值图的组件数；正确数量本身不保证正确成员分组 |
| `raw_spurious` | 原始图在 `theta` 处的组件数超过 2 |
| `zero_cross_contact` | 至少一条去偏跨概念距离等于零 |
| `theta_cross_contact` | 最小跨概念半径 `r_out <= theta`，存在阈值内跨边 |
| `true_partition` / `partition_failure` | 分别是 `r_in <= theta < r_out` 及其否定，其中 `r_in` 为两个概念内部的最大 MST 连接半径 |
| `debiased_collapse` | 全部 128 个顶点在零半径处已经连通 |

`channel_probabilities.csv` 给零接触与阈值接触同时保存精确单对非中心卡方概率，以及
`1-(1-p)**m`、`min(1,m*m*p)` 两个严格界。配对之间共享噪声，因此没有把所有 `m**2` 对当成独立样本。

本次首次事件的经验中位数如下。归一化幅度为 `sigma / [Delta*(D_a/log(m))**0.25]`。

| `D_a` | 原始碎裂中位数：绝对 `sigma` | 零跨接触中位数：归一化 | 分区失败中位数：归一化 | 全图零坍缩右删失条数 / 120 |
|---|---:|---:|---:|---:|
| 64 | 0.37 | 0.3900 | 0.3150 | 9 |
| 256 | 0.37 | 0.4125 | 0.3450 | 3 |
| 1024 | 0.36 | 0.4350 | 0.3700 | 2 |
| 4096 | 0.36 | 0.4400 | 0.3800 | 9 |

首次事件按网格首次出现记录，不是连续参数的精确临界值。未在扫描内出现的路径按右删失计入经验分位数，
`null` 表示该分位数超过扫描上限；不应删除这些路径再计算“坍缩阈值”。完整记录见
[channel_onsets.csv](../code/results/topology/channel_onsets.csv) 和
[channel_summary.json](../code/results/topology/channel_summary.json)。

`channel_events.pdf` 的 (d) 面板固定 **128 个条目的 10% 块错误率目标**。
图接收机使用每个固定 `sigma` 的 `partition_failure` 频率，记录最后一个低于 10% 与第一个达到 10% 的相邻采样水平；
它与上表的“每条路径首次事件中位数”是不同统计量。
已知中心接收机使用精确块错误率 `1-(1-Q(Delta*sqrt(D_a)/(2*sigma)))**n`，
因此解析阈值为 `Delta*sqrt(D_a)/(2*Q_inverse(1-0.9**(1/n)))`，其中 `n=128`。
图接收机的交叉区间描述有限网格上的经验频率，不能当作连续噪声阈值的精确值或统计置信区间；
两端的噪声水平、失败频率、Wilson 区间及已知中心阈值保存在
[channel_block_reliability.csv](../code/results/topology/channel_block_reliability.csv)。
其中 `graph_last_below_sigma` 与 `graph_first_at_or_above_sigma` 是经验交叉的两端，
`known_center_block_sigma` 是相同块错误率目标下的解析阈值。

| `D_a` | 图接收机的经验交叉区间：`sigma` | 两端失败次数 / 120 | 已知中心的解析 `sigma` |
|---|---:|---:|---:|
| 64 | [0.54467, 0.55457] | 9/120 → 14/120 | 1.27077 |
| 256 | [0.85431, 0.86832] | 9/120 → 13/120 | 2.54153 |
| 1024 | [1.33692, 1.34682] | 10/120 → 12/120 | 5.08307 |
| 4096 | [1.93270, 1.96071] | 10/120 → 12/120 | 10.16614 |

有限实验中的接触数值不确立通用精确常数；理论保留已证明的渐近常数区间。

### 2.4 有限分辨率失配与实际同步消息

720 个配对条目给出一个圆形发送端参考，以及分裂成两圆、椭圆收缩两种接收端几何。
字典由发送端最远点采样确定，原子数为 `8, 16, 32, 64, 128, 256, 512`。
平面 alpha 滤过计算与偏移滤过对应的持久性；每次使用接收端的真实有限参考覆盖半径检查
`abs(distortion - reference_delta) <= cover_radius`。

14 个测量均满足此双侧带，其中 10 个低于对应有限参考图的失配距离。
因此失配定理给出的是指定解释规则的稠密字典极限及有限分辨率误差带，不保证任意有限码率曲线都在极限之上。
`reference_delta` 是 720 条目参考云的距离；它与连续支持的距离仍有采样差异。
六个椭圆变形幅度及圆退化为区间的独立计算也保存在
[mismatch_summary.json](../code/results/topology/mismatch_summary.json) 与
[deformation.csv](../code/results/topology/deformation.csv)。

同步阶段采用已共享身份、互相分离的局部特征图：`K = 1, 2, 6, 12`，`eps0 = 0.08`，
精度比 `eps0/eps1 = 2, 2.5, 4.5, 6.5, 9.5`，每种组合 25 次，包含偏移立方体的两个角点。
编码器生成真实二进制消息，再由解码器恢复所有端点。500 次检查同时验证 bottleneck 精度与
`ceil(log2(m**(2*K)))` 位的最优固定消息长度，其中 `m = ceil(eps0/eps1)`。
例如 `K=12, m=3` 时，联合整数编码使用 39 bit，逐坐标编码使用 48 bit。
这项最优性针对论文所列的局部对应关系和分离条件，不包含特征创建、删除或字典初始化成本。

### 2.5 分支任务的真实熵与图诊断

任务源为主半径 2、次半径 1 的旋转环面，按表面积均匀采样，得分为第一坐标。
真实 Reeb 分支由解析几何标注；先对方位角积分，再对另一角积分，得到
`H(E|T) = 0.3537667163287371` bit。
500,000 个样本上的中间分支区域概率为 0.354436，其 Wilson 区间为
`[0.3531112508, 0.3557629859]`。

对步长 `0.4, 0.2, 0.1, 0.05, 0.02, 0.01` 的得分量化，直接计算经验熵差
`H(E,Q(T)) - H(Q(T))`；最细步长得到 0.3541681590 bit。
这验证分支信息的额外代价，不把标量量化器的有限精度性能等同于 Shannon 率失真函数。

另用 20,000 点、3 个连接半径和 3 个分箱数测量 slab 分量熵与 Mapper 图。
其熵估计在 0.29965–0.33282 bit 间变化，反映有限样本与分箱选择。
**这些是诊断量，不是经过认证的 Reeb 图或真实条件熵**；也没有要求 Mapper 的 Betti 数小于某一阈值下的 Rips 条数。
参数敏感性见 [reeb_diagnostics.csv](../code/results/topology/reeb_diagnostics.csv)。

### 2.6 离线 Digits 表示实验

使用 scikit-learn `load_digits` 的全部 1,797 个 8×8 图像，将像素除以已知量程 16。
按标签分层、固定种子做 **1,078 个 reference / 719 个 evaluation 的内部划分**；它不是数据集另有的官方测试集。
`pixels64` 保留 64 个像素坐标；`pca16` 只在 reference ID 上拟合 16 维 PCA，保留训练方差约 85.03%。

归一化中心、距离尺度、类中心、阈值与线性得分探针全部在 reference split 上拟合。
200 个共同 evaluation ID 用于跨表示持久图与字典实验；分类曲线另从 evaluation split 每类抽 20 个样本。
输入缓存按唯一 `item_ids` 的交集排序对齐，随后逐 ID 检查标签一致性。
缺失 ID、重复 ID 或对齐后标签冲突会直接报错，不会按最短数组长度截断对齐。

在各自训练尺度归一化后，两个表示的有限参考图差异为 0.0555076748，双向共 10 个字典测量均满足误差带。
无噪声时，留出分类子集的类中心准确率为 89.0% 和 88.5%；`sigma=1` 时为 78.0% 和 55.5%。
同一无噪声分类子集的 H0 计数分别为 10 和 7，说明图组件数与标签数并不自动一致。
这是一组可重复的离线表示案例，不构成 CLIP、DINOv2 或 SBERT 的实测结果。

参考、评估与持久图子集的完整 ID、两个缓存的 SHA-256、归一化尺度及探针协议见
[digits_report.json](../code/results/topology/digits_report.json)。其中 `reeb` 段输出的 slab/Mapper 量同样仅作诊断。

## 3. 回归、来源与图表

[validation/results.json](../code/results/topology/validation/results.json) 和
[test_report.txt](../code/results/topology/validation/test_report.txt) 记录 14/14 项回归通过，并包含独立库对照：

| 独立核验 | 本次结果 |
|---|---:|
| 原生 bottleneck 与 Persim，300 对随机图 | 最大绝对差 0 |
| 原生 H0/MST 与 Ripser | 最大绝对差 `5.05e-8` |
| 平面 alpha 的 H0 与 MST | 最大绝对差 `6.94e-18` |
| 原生平面 alpha 与 GUDHI，20 个点云 | 最大 bottleneck 差 `1.33e-15` |

回归还覆盖重复点的噪声参考、零 MST 边、essential 类、退化区间、非二次幂混合进制消息、
无效码字、缓存 ID 重排与拒绝路径、训练集探针要求和 Mapper 覆盖条件。

固定基础种子为 `20261020`。重参数化噪声种子为 `SEED + 100000 + D_a*100 + trial`；
信道种子为 `SEED + 200000 + D_a*1000 + trial`；同步偏移用 `SEED + 300000`；
任务熵样本用 `SEED + 400000`，诊断云用其后一种子；Digits 划分使用基础种子。
CSV 保留每个噪声重复的种子，NPZ 保存固定参考云；主报告与验证报告分别保存其覆盖的关键代码 SHA-256。

| 论文矢量图 | 对应记录 |
|---|---|
| [dictionary_mismatch.pdf](../code/results/topology/dictionary_mismatch.pdf) | `dictionary.csv`、`mismatch.csv`、`sync_codes.csv` |
| [noise_reparameterization.pdf](../code/results/topology/noise_reparameterization.pdf) | `reparameterization.csv`、`calibration.csv` |
| [channel_events.pdf](../code/results/topology/channel_events.pdf) | `channel_probabilities.csv`、`channel_onsets.csv`、`channel_block_reliability.csv`、`channel_summary.json` |
| [task_digits.pdf](../code/results/topology/task_digits.pdf) | `reeb_quantization.csv`、`digits_report.json` |

## 4. 预训练表示扩展协议

CLIP、DINOv2 和 MiniLM/SBERT 的模型提取与缓存分析入口已经接入，**本次未执行模型下载或推理**。
按 [code/embeddings/README.md](../code/embeddings/README.md) 的命令，先对同一数据集与 split 提取带 ID 的缓存，
再运行五个缓存分析阶段 `floor,dictionary,reparam,cliff,reeb`。
所有训练预处理与接收机拟合继续使用单独 reference split；任务类别是否对应几何组件由测量判断。
扩展结果应独立报告其数据、模型、尺度与样本量，不能用本页的合成或 Digits 数值代替。
