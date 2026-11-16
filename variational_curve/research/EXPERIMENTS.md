# 曲面：四组对照与 SEMASIA

**主要比较是 joint 对 density。** 相同的任务头、同一固定评价代价、锚点 ID、码本大小、核温度、发射能量和训练预算。在冻结后验的 Hellinger 坐标与原始单位特征的拼接空间学习低秩残差。它是一个有局部度量解释的环境弦长模型；不能把样本间的弦长叫作真实曲面测地线。待发表结果应报告哪些受限接口需要学曲面，以及哪些仅学概率就够了。

| 编码器 | 变动参数 | 检验什么 |
|---|---|---|
| fixed | 无编码器可学参数；输出边缘辅助变量单独估计 | 参考表示和均匀锚质量 |
| density | 每个锚的编码质量 | 仅改变各码字使用的相对意愿 |
| geometry | 有全局伸缩界的残差映射 | 学曲面能否在均匀质量时改变有意义的邻域 |
| joint | 上述两者 | 是否超出仅概率、仅曲面的效果 |
| 物理对照 | 最近余弦锚、最近后验 KL 锚、直接任务判决 | 实际的同样本固定宽度传输；直接决策按自身实际位数计 |

**来自两代方法的基准。** 旧 96 点非均匀圆在 D=0.30 处的真实插值表见 [原始摘要](code/early_studies/results/communication_final/summary.json)：仅学密度 1.18319194，联合 1.18319195 bit/状态，不能归因为曲率收益。新版固定流形的球面、曲面族、投影噪声与中断有 [逐条件数据和图](code/results/paper1)；保持原指标和自然对数/bit 换算。残差的二阶平滑罚只作离散正则项，不拿它估计 Willmore 能量。

**真实点云与通信系统。** 使用[官方数据](https://huggingface.co/datasets/spaicom-lab/semasia-cifar100)的 ResNet-18 和 MobileNetV3-Small；校验每个模型各 split 的 SHA-256 和 model_name。按排序后的图像 ID 同步抽取。拟合、温度校准、验证三者从官方 train 互斥抽取；最后才用官方 test。模型、单位归一化、参考损失权重 eta 固定，20 粗类为主、100 细类为敏感性任务。不同模型的特征维数不需相同，不把两者拼成一个物理流形。`experiments/semasia_curve.py` 优化互信息加固定失真；反向传播使用可优化的全局边缘上界，每个 epoch 与测试均用实际输出边缘重算 I。训练结束后冻结映射和锚质量。收发双方预先共享同一头和锚点；采样码字索引后真正生成比特、调制并硬判决。无信道、AWGN 和每索引一次相干 Rayleigh 的 SNR 标记均为 Eb/N0；一个比特占一个实 BPSK 使用。

**实验矩阵与判定。** 扫 M=16/64/256、lambda=1/3/8、两骨干、3 个训练种子。在一组源/任务内预先固定 eta，先在验证集为各方法公平扫描尺度 s（默认 5）和学习率，比较在相同固定失真的 I，另外比较同等实际位数/功率的性能。输出真实标签 NLL/准确率、后验 KL、固定几何失真、BER/BLER、各组参数量、总训练时间；消融禁用平滑或改变秩。预期可证伪：若 joint 对 density 的配对区间跨零或检索/直接判决更优，只报告该工作点和失败范围。`results/*_paired/` 下的图像 ID、每次噪声重复 KL/正确性可用于配对，先在每图内平均重复，按 ID 重采样，在独立训练种子上报告波动。更换评价准则单列，不能拿一个方法自己学出的距离评价所有方法。

**已执行的短验证。** 真实 ResNet-18 的四组都训练两轮、重新加载检查点，并在 192 张官方测试图上做干净、AWGN、Rayleigh 评估：见 [固定](results/semasia_smoke_fixed.json)、[仅概率](results/semasia_smoke_density.json)、[仅曲面](results/semasia_smoke_geometry.json)、[联合](results/semasia_smoke.json)。16 个锚，4 个真实使用；联合模式测试 I 约 3.448 bit/图，预期干净后验 KL 约 1.748 bit。不同损失、不同约束下不能拿这个数和工作一的 BA 比较优劣。全数据量、训练收敛和独立几何收益还未得到验证。

从本项目根目录运行；先打开生成的清单确认工作量。`plan` 只列配置；`train` 保存检查点；`validation` 与 `test` 各自读取已训练的检查点，写独立 JSON。复制 `experiments/matrix.json` 后可用 `--spec my-matrix.json` 只跑选定设置。全面板作业较多，请先跑 pilot；不要用正式测试集挑超参数。
~~~sh
python experiments/run_matrix.py plan --profile full --data data/semasia --out runs/matrix
python experiments/run_matrix.py train --profile pilot --device cuda --download --data data/semasia --out runs/matrix
python experiments/run_matrix.py validation --profile pilot --device cuda --data data/semasia --out runs/matrix
python experiments/run_matrix.py test --profile pilot --device cuda --download --data data/semasia --out runs/matrix
~~~

已有结果做配对检验的命令（差值方向是第一个减第二个，区间以这两份固定检查点为条件）：

~~~sh
python experiments/compare.py results/semasia_smoke.json results/semasia_smoke_density.json --channel clean --metric accuracy
~~~

## 真实特征上的几何可识别性实测

对四个保存的两轮训练检查点，在同一批 2,048 个官方 ResNet-18 测试 ID 上计算行列双中心化的 logit 变化，所有模式共享评价成本、16 个锚点、$s=5$ 和 $\lambda=3$。固定、密度、曲面、联合模式的 RMS 依次为 0、0、0.262074、0.261949。曲面模式相对固定模式减少 0.02302 bit/图的互信息，同时增加 0.00539 的固定失真。论文把这两个变化同时列出。等失真结论依赖正式斜率扫描与独立训练种子。[说明与命令](README.md) · [带 ID 和检查点哈希的测量](results/encoder_contrasts.json)

## 交互秩与受限曲面的圆周核对（2026-11-16）

在同一 96 点非均匀圆周、同一 25 个斜率和同一平方弧长评价下，新增三类受限编码曲面：单位圆弦长（单一欧氏曲面）、由双中心化评价矩阵截断奇异值分解得到的 $d$ 维点积曲面、以及原有弧长代价。质量与尺度以真实 $I+\beta D_g$ 为目标优化，在 $D_g=0.30$ 处对各自的可达点插值。Blahut 参考复现 1.183191 bit，均匀质量弧长核复现 1.264985 bit。学习质量的弦长曲面为 1.189914 bit；谱曲面在 $d=2,4,8,16,32$ 时为 1.189897、1.210406、1.187964、1.183713、1.183265 bit。$d=4$ 劣于 $d=2$，原因和超额恒等式的逐项数值见[推导第六节](derivations/interaction_rank.md)。这些是随机测试信道的信息率，未涉及固定长度索引。真实特征云的交互谱尚未计算。[结果](results/interaction_rank.json)

~~~sh
python tests/test_interaction.py
python experiments/interaction_rank.py --output results/interaction_rank.json --dimensions 2 4 6 8 10 12 16 24 32
~~~
