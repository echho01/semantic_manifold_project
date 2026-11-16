# 表示缓存与拓扑实验

`run_topology.py` 在配对特征缓存上执行 `floor,dictionary,reparam,cliff,reeb` 五阶段分析。
它只读取数组，不加载预训练模型。当前已经执行的案例是离线 Digits 的 64 维像素与 16 维 PCA；
CLIP、DINOv2、MiniLM/SBERT 是可运行的扩展协议，本次没有执行其模型下载或推理。
完整合成实验与结果解释见 [Paper 2 实验说明](../../notes/dictionary_experiments.md)。

## 1. 复现已执行的 CPU 结果

在仓库根目录运行：

```bash
python -m pip install -r code/requirements.txt
python code/tests/run_topology_checks.py && \
python code/synthetic/run_topology_suite.py \
  --output code/results/topology --reparam-trials 20 --channel-trials 120
```

此命令执行 14 项回归、独立库对照、全部六阶段实验及绘图。只重做 Digits 可使用
`python code/synthetic/run_topology_suite.py --stages digits`；此时主报告仅记录该阶段。
结果在 [`../results/topology/`](../results/topology)，其中
[`digits_report.json`](../results/topology/digits_report.json) 给出模型数组、完整 ID、尺度、缓存 SHA-256 和测量值。

### 本次 Digits 的训练与评估协议

使用全部 1,797 个 8×8 图像。按标签分层、`seed=20261020` 做 1,078 个 reference 与 719 个 evaluation
的内部划分，没有把它称作官方测试集。像素除以已知量程 16；16 维 PCA 仅在 reference ID 上拟合。
PCA 保留约 85.03% 的 reference 方差。

所有数据依赖的预处理与接收机拟合使用 reference split：居中、两两距离中位数尺度、类中心、
H0 阈值及 ridge 得分探针。200 个共同 evaluation ID 用于 `H0/H1` 图与字典实验；
分类曲线另从 evaluation split 每类选 20 点，阈值由独立 reference 子集确定。

实际结果包括：有限参考图差异 0.0555076748，双向 10 个字典测量全部在覆盖误差带内。
无噪声类中心准确率为 89.0%（像素）和 88.5%（PCA），`sigma=1` 时为 78.0% 和 55.5%。
无噪声 H0 计数分别为 10 与 7；类别数与几何组件数是不同指标。
这些是单个固定划分下的离线结果，没有作为预训练模型表现或普遍优劣结论。

## 2. 缓存格式与 ID 对齐

输入位置为 `<paths.root>/<model>__<dataset>/arrays.npz`：

| 数组 | 要求 |
|---|---|
| `image` 或 `text` | 有限数值的二维数组，形状为 `[n_items, dimension]` |
| `item_ids` | 一维、逐条唯一、与特征行等长；不同模型必须指向同一数据集和 split 中的相同条目 |
| `labels` | 运行 `cliff` 或 `reeb` 时需要；长度与特征一致，同一 ID 在各缓存中的标签必须一致 |

驱动按唯一 ID 的交集排序，再重排每个缓存，最后逐 ID 比较标签。
缺失/重复 ID、非有限特征或标签冲突会报错。已有旧缓存若没有 `item_ids`，应重新运行提取阶段，
不能通过行号或截断到最短数组来猜测对应关系。`data_sources.py` 与 `run_all.py --stages extract`
会把原始条目身份写入缓存；数据集、split 与配置也应保持一致。

完成对齐后，驱动按固定种子生成不交的 reference/evaluation split；存在标签时使用分层划分。
每个表示分别在 reference split 上拟合均值与距离尺度，随后应用于 evaluation 数据。
报告中的 `normalization_reference_median` 保存使用的尺度，有限图距离以这些固定训练尺度为单位。
它不表示不同归一化规则会产生相同的跨表示距离。

## 3. 可选预训练表示协议

以下命令从 `code/embeddings/` 执行。首次提取需要相应模型权重与数据下载；模型提取可以使用 GPU，
缓存后的拓扑阶段只需 CPU。本次交付的结果中不包含以下运行。

```bash
cd code/embeddings
python -m pip install -r requirements.txt

python run_all.py --config config.yaml --model clip_vit_b32 \
  --dataset cifar10 --stages extract
python run_all.py --config config.yaml --model dinov2_s \
  --dataset cifar10 --stages extract
python run_topology.py --config config.yaml --dataset cifar10 \
  --modality image --models clip_vit_b32,dinov2_s

python run_all.py --config config.yaml --model minilm \
  --dataset agnews --stages extract
python run_topology.py --config config.yaml --dataset agnews \
  --models minilm --modality text --stages reparam,cliff,reeb
```

`config.yaml` 中 `minilm` 使用 `sentence-transformers/all-MiniLM-L6-v2`；
`clip_vit_b32` 使用 `openai/clip-vit-base-patch32`，`dinov2_s` 使用 `facebook/dinov2-small`。
图像比较显式选择 `--modality image`。CLIP 的类别文字模板包含标签信息，不应把它作为同一图像实验的独立无标签表示。

提取后的 `bash run_topology.sh config.yaml` 可按上述图像对与文本单模型顺序分析缓存。
输出为 `runs/topology__<dataset>/report.json`，也可用 `--tag` 指定不同输出名。
扩展运行应保存配置及实际模型/数据版本；现有报告自动保存输入缓存和关键分析代码的 SHA-256。

### 缓存分析阶段

| 阶段 | 测量及判读 |
|---|---|
| `floor` | 同一留出条目上的有限参考图距离矩阵；没有自动认定两个模型的连续支持已知 |
| `dictionary` | 在发送端选择同一组条目，在接收端重构，检查 `abs(d_topo-floor) <= eps_R` |
| `reparam` | 将含噪图与完整噪声参考滤过比较，再检查去偏误差；`sigma` 为总幅度 |
| `cliff` | 用 reference 数据选择阈值和类中心，在留出样本上分别报告组件计数及标签准确率 |
| `reeb` | reference 上拟合得分，再测量留出云的 slab 分量熵和 Mapper 图，作为有限样本诊断 |

`--stages` 接受这些名称的逗号列表。`floor` 和 `dictionary` 需要至少两个表示；
`cliff` 和 `reeb` 需要一致标签。`config.yaml` 的扩展协议基础种子为 0，
已执行的 Digits 和主合成套件基础种子为 20261020；报告记录各自真实配置。
配置中的 `noise_levels` 和 `cliff_noise_levels` 都是 `sigma` 幅度，实际总噪声能量为其平方。

## 4. 结果解释

正确组件数不保证正确条目分组。合成主套件另外保存 `r_in`、`r_out`、跨接触、真实分区与全图零坍缩，
缓存分析的组件数与分类准确率则是分别测量的统计量。合成图 `channel_events.pdf` 的 (d) 面板固定
128 个条目的 10% 块错误率目标，对照图接收机的经验频率交叉区间与已知中心接收机的解析阈值。
前者按每个固定 `sigma` 的 `partition_failure` 频率计算，不是逐路径首次事件的中位数；
有限采样网格给出的区间也不是连续阈值的精确值或统计置信区间。

有限字典的失配可以在参考图差异的两侧波动，使用双侧覆盖误差带判断，而不是要求永不低于它。
重参数化界带有失败概率，校验时应同时报告指定的概率参数、重复数、误差口径与观测到的超界次数。
噪声参考必须保留重复样本的零合并半径；原生实现对此有专门回归。

`reeb` 输出中的 `binned_component_entropy_bits` 是 slab 分量的经验条件熵。
为兼容旧调用保留的 `H_E_given_f_bits` 字段也具有相同诊断含义，不能解释为已认证的真实 Reeb 条件熵。
Mapper 的 Betti 数受覆盖、连接半径与样本量影响；没有以某个 Rips 持久条计数作为通用上界。
主套件的环面任务另有解析分支标签与积分真值，承担论文分支开销的定量验证。

## 5. 模块与检查

| 文件 | 作用 |
|---|---|
| `data_sources.py` | 数据加载与原始条目 ID，避免与 Hugging Face `datasets` 包重名 |
| `models.py`、`run_all.py --stages extract` | 可选预训练特征提取与带 ID 缓存 |
| `run_topology.py` | 严格 ID 对齐、reference/evaluation 划分、报告与来源信息 |
| `topology_exp.py` | 训练尺度、有限失配、字典、噪声、留出类中心与探针测量 |
| `selftest_topology.py` | 小型合成缓存流程自检 |
| `../tests/run_topology_checks.py` | 14 项回归和独立 TDA 库核验，保存正式验证记录 |

本次正式验证为 14/14 项回归通过；跨库对照包括 300 对随机持久图、MST/Ripser 与 20 个原生 alpha/GUDHI 点云。
完整误差与关键源文件哈希见
[`validation/results.json`](../results/topology/validation/results.json)。
同目录的 `geometry.py`、`rd.py`、`manifold_gain.py`、`task_mismatch.py` 及其旧自检属于保留的 Paper 1 分析路径，
不计入本次 Paper 2 的实验通过数。
