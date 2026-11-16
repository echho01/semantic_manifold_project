# 率失真：实验计划与现有证据

**最先检验的主张。** 固定任务和观察时，强制“解码为一个能实现的对象”是否付出可测代价？要分别比较操作定义、有限源上下界和实际数字链路。网格 BA 上界不能充当全体连续重构的下界；拟合后验与真实 Bayes 后验也不能混称。

| 实验 | 自变量与对照 | 判据、已有材料 |
|---|---|---|
| 后验圆、条件源 | 圆上重构/圆盘网格、均匀/共享 K、不同网格 | 固定 D 给上下界而非单线。已有 [有限上下界](results/finite_rd_bounds.csv)、[零率](results/toy_zero_rate.csv)、[高斯链路](results/gaussian_link.csv) |
| 真实特征，经验源编码 | 固定任务头；嵌套比较 M 个训练对象和这 M 个对象加 M 个相邻锚点后验中点；改变 beta | 训练集输出目标上下界及 gap；未参与训练的数据上用冻结先验重新计算条件分布、实际边缘、I 和 KL。测试 I 不是训练问题的下界。“有限中点”不是“完整单纯形” |
| 实际比特及无线链路 | 固定宽索引；无噪声、单位能量 BPSK+AWGN、每索引一个相干 Rayleigh 增益；相同常规比特映射，无纠错/重传 | 实际使用次数、BER、索引 BLER、KL、准确率、真实标签 NLL。固定 2M 输出比 M 多一位；在相同带宽比较 2M 个对象与 M+M，须两次使用相同底层训练排列；在相同信息率时单独比较曲线 |
| 任务与模型失配 | 同模型余弦近邻、同模型 KL 近邻、直接发标签决策；另用同一 ID 所代表的异模型锚点解码 | [既有双模型预实验](results/semasia_cifar100.json)提供跨模型真实负对照。直接发送 20 类判决要 5 位；无效字样计为错误/擦除，无免费回退，硬判决不报伪造 KL |

**真实点云上的几何对照。** 除后验 KL，训练入口另有 `--criterion task_chord` 与 `--criterion graph_geodesic --neighbors 12`，均固定对象码本。前者在冻结的 `(2√q, √η u)` 坐标用平方弦长；后者在训练点构建带权并集 kNN 图，计算锚点的平方最短路。测试图只能挂到训练顶点，不参与训练图拓扑或相互连边。报告构图 k、连通性、成本单位和对 k=8/12/24 的敏感性；断图直接报错。比较弦长和图路径、以及下游 KL/标签风险是否一致；平方图长没有 bit 单位。[已运行的图测地近似](results/semasia_smoke_geodesic.json)和[相同锚点弦长对照](results/semasia_smoke_task_chord.json)使用相同 ResNet-18 短验证数据。正式矩阵同时覆盖这些判据。参见 [连续几何上下界与小失真推导](derivations/geometric_bounds.md) 和 [图实现](experiments/reference_graph.py)。

**数据与划分。** [SEMASIA 原论文](https://arxiv.org/abs/2605.09485)、[官方 CIFAR-100](https://huggingface.co/datasets/spaicom-lab/semasia-cifar100)。固定 revision `ace0aa91fc0b79953ec1ac619a25072d3a5b3037`，SHA-256 见 [读取器](experiments/semasia_io.py)。ResNet-18 512 维、MobileNetV3-Small 1024 维；每个模型单独按训练 ID 的稳定随机排列取 8192 个拟合、2048 个温度校准、首 M 个拟合 ID 作锚；验证用官方训练集余下 ID，测试只从官方 test 取样。拟合线性头使用原标签，但编码决策看测试特征预测，不看真标签。默认任务是 20 个粗类，100 个细类作任务敏感性比较。ID 键带 split，禁止把不同模型向量按行号或不同维坐标逐维相减。码本、头、知识预先共享；若要更新，另计描述开销。

**正式表格。** 两个模型，训练种子 1115/1116/1117，探查码字数 M=16/64/256（做严格等宽对照时另加 32/128/512 对象码本），beta=0.25/0.5/1/2/4/8。状态：控制住合法输出和双边上下界可复核有限圆例，真实数据的全矩阵尚未执行。对每个模型在共同 D 区间比较上下界/可达率，意外饱和时报告有限字母表、估计维数对尺度/子样本的敏感性。对同一测试 ID 的差按图像配对重采样；跨训练种子单列区间。不要把多次噪声或相同图像当成独立源样本。报告训练 gap 是否达到容差与运行时间；只对训练的有限离散问题使用“证书”。

**这次执行的例子。** [GPU 测试 JSON](results/semasia_smoke.json)：真实 ResNet-18，768 拟合/192 校准/192 官方测试，16 个锚和 beta=3。训练 Lagrangian 对偶间隙约 0.00000683 bit；独立测试的经验信息率 2.739 bit/图，实际发送 4 个 BPSK 符号/图。二者用途不同。[第二模型与扩充字母表](results/semasia_smoke_augmented_mobile.json)用于检查另一条代码路径，两项改变同时发生，不能比较算法优劣。旧版较大预实验用 8192/2048/2048，256 锚点下的真实准确率应与更便宜的直接标签判决比较。JSON 与逐图 .npz 只作复现实例，不是正式多种子结论。

从本项目根目录运行；先打开生成的清单确认工作量。`plan` 只列配置；`train` 保存检查点；`validation` 与 `test` 各自读取已训练的检查点，写独立 JSON。复制 `experiments/matrix.json` 后可用 `--spec my-matrix.json` 只跑选定设置。全面板作业较多，请先跑 pilot；不要用正式测试集挑超参数。
~~~sh
python experiments/run_matrix.py plan --profile full --data data/semasia --out runs/matrix
python experiments/run_matrix.py train --profile pilot --device cuda --download --data data/semasia --out runs/matrix
python experiments/run_matrix.py validation --profile pilot --device cuda --data data/semasia --out runs/matrix
python experiments/run_matrix.py test --profile pilot --device cuda --download --data data/semasia --out runs/matrix
~~~

已有结果做配对检验的命令（差值方向是第一个减第二个，区间以这两份固定检查点为条件）：

~~~sh
python experiments/compare.py results/semasia_smoke.json results/semasia_smoke_task_chord.json --channel awgn --snr 0 --metric kl_bits
~~~

## 与稿件同步的重复计算（2026-11-16）

以 24 点网格重算 64 个操作点、8 个目标失真区间、5 个零率例和 5 个高斯增益。四份 CSV 与原记录逐项一致至浮点误差。重绘主稿零率图的命令为 `python toy_geometry_check.py --grid 24 --output runs/circle --figure ../figures/toy_posterior_manifold.png`。

## 结构结论的数值核验（2026-11-16）

`theory_checks.py` 只依赖 NumPy（绘图另需 Matplotlib），在 720 态网格上约需 12 分钟，主要耗时在按预测接收机的信息量对对象接收机做斜率二分。输出 [theory_checks.json](results/theory_checks.json) 与主稿图 3。

| 检验 | 设置 | 结果 |
|---|---|---|
| 五项损失恒等式 | 64 观测态并成 32 个表示，$\beta=8$ 的对象 BA 编码器（0.4707 bit），5 位 BPSK 索引，查表 / m-投影 / 后验均值三种接收 | 12 个配置下五项之和与直接枚举的对数损失相差 $<10^{-15}$ bit；无噪时对象约束 0.02398 bit；2 dB Gray 映射下信道项 0.00292、查表解码项 0.03680、投影解码项 0.02670 |
| 等信息间隙 | 720 态均匀圆，预测接收机用自由质心 Blahut 迭代，对象接收机在相同信息量处取值 | 间隙除以 $D_{\rm pred}^2$（nat）从 5.36 降到 3.81，$\tfrac12\mathbb E\kappa_m^2=3.804$；分离定理下界系数 0.074，reach 上界系数 125 |
| 任务秩渐近线 | $h_F=1.02165$ nat，$r=1$ | 渐近线与对象码率之差在 4.92 bit 处为 0.0006 bit |
| 容量处失真下限 | 24 态对象问题的对偶下界，容量取两次实高斯使用 | 增益 0.5、1、2 时下限 0.10905、0.03316、0.00396 bit，实测 0.12823、0.04561、0.01034 bit |
| 共享知识的节省 | 24 态圆、$K\in\{-1,1\}$ | $I(U;K)=0.25515$ bit，$\beta=16,32$ 的拉格朗日差 0.25481、0.25515 bit |

稠密网格的数值是可达点的计算，没有连续字母表的对偶证书；主稿表 I 的 24 态区间才是证书。等信息间隙趋于混合曲率常数是数值观察，单胞元的局部命题有证明。
