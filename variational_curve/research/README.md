# 工作二研究索引 · 编码曲面的变分泛函

[Letter 主稿（拟投 IEEE Communications Letters）](../main.tex) · [论文 PDF](../build/main.pdf) · [书目](../ref.bib) · [实验矩阵及运行记录](EXPERIMENTS.md)

本工作固定语义度量 $g$，把编码曲面 $F$（诱导度量 $h$）当作率失真问题的变量。主稿推出两个泛函：各向同性数字编码器的等距亏量 $\Lambda[h]$（含一阶变分、极小元和任务无关方向的各向异性推论），以及模拟映射的面积–reach–功率界。有限码本上的可吸收判据与超额恒等式是它的离散对应。交互秩定理、单一欧氏曲面推论和谱曲面表已从主稿移出，推导、结果和测试仍保留在本目录。

| 数学问题 | 内容与来源 |
|---|---|
| [变分曲面泛函：等距亏量、维数、面积与 reach](derivations/variational_surface.md) | 主稿定理 1–2、命题 1、推论 1 的完整证明，曲面能控制哪种曲率，两种分辨率区的数值核对 |
| [编码几何与码字质量的可识别性](derivations/encoder_identifiability.md) | 主稿命题、双中心化实测量、离散一阶变分与残差映射伸缩界 |
| [交互残差、超额恒等式与交互秩](derivations/interaction_rank.md) | （超额恒等式在主稿命题 2；交互秩部分不在主稿）相对 Blahut 最优的精确超额与 $\epsilon^2/(2\ln2)$ 上界、点积编码头的秩判据、单一欧氏曲面的负型条件、谱曲面 |
| [连续编码曲面与有限实现](derivations/surface_variation.md) | 早期变分模型中的材料导数、测度、似然接收机及任务失配 |
| [参考几何与率失真核展开](derivations/geometry_and_rate.md) | main 分支的内在与外在损失、曲率系数、匹配上下界 |
| [高码率码点分配](derivations/allocation.md) | 固定码本规模、熵约束及局部胞元近似 |
| [投影信道、任务约化与失配](derivations/projection_and_task.md) | 法管投影的四阶展开、分量误判、任务充分性和度量界 |

| 实验 | 可复现入口与证据 |
|---|---|
| 模拟映射：等功率等 reach 下的等距与压扩曲面 | [计算](experiments/analog_mapping.py) · [测试](tests/test_variational.py) · [结果](results/analog_mapping.json) · 主稿表 `tab:analog` |
| ResNet-18 特征：原始曲面、学习曲面与等距参照在相同失真处的超额码率（3 个种子） | [计算](experiments/raw_surface.py) · [结果](results/raw_surface.json) · 主稿表 `tab:real` |
| 720 态圆与环面：等距亏量、密度区与几何区、各向异性推论 | [计算](experiments/surface_functional.py) · [测试](tests/test_variational.py) · [结果](results/surface_functional.json) · [主稿图 1](../figures/surface_functional.png) |
| 96 点圆周的固定、曲面、密度、联合与 BA 对照 | [复现](experiments/reproduce_early.py) · [绘图](experiments/plot_circle_ablation.py) · [各斜率数据](code/results/early_surface/G5_geometry_density.csv) · [等失真结果](code/results/early_surface/reproduced.json) |
| 96 点圆周的交互谱、弦长曲面与 $d$ 维谱曲面的等失真码率 | [计算](experiments/interaction_rank.py) · [随机实例检验](tests/test_interaction.py) · [结果](results/interaction_rank.json) |
| 连续球面、旋转曲面和投影噪声 | [计算程序](code/synthetic/run_paper1.py) · [设置和数值误差](notes/fixed_geometry_experiments.md) · [原始数值](code/results/paper1/results.json) |
| SEMASIA 上的有限锚点编码与物理链路 | [曲面训练与测试](experiments/semasia_curve.py) · [独立的几何变化检验](experiments/measure_deformation.py) |

以下命令在本目录执行。共享缓存只存官方数据分片，工作一和工作二可独立运行。
~~~sh
python -m pip install -r requirements.txt
python tests/test_surface.py
python tests/test_interaction.py
python tests/test_variational.py
python experiments/surface_functional.py --output results/surface_functional.json --figure ../figures/surface_functional.png
python experiments/interaction_rank.py --output results/interaction_rank.json
python experiments/semasia_curve.py train --data data/semasia --download --device cuda --book 64 --mode joint --lam 3 --epochs 30 --checkpoint runs/joint.pt
python experiments/semasia_curve.py test --data data/semasia --download --device cuda --checkpoint runs/joint.pt --output runs/joint-test.json
python experiments/measure_deformation.py --data data/semasia --download --checkpoints results/checkpoints/curve-fixed.pt results/checkpoints/curve-density.pt results/checkpoints/curve-geometry.pt results/checkpoints/curve-joint.pt --output runs/contrasts.json
~~~

[实验协议](EXPERIMENTS.md)给出四种模式的批量执行方法、模型选择和多种子比较。圆周上密度模型已接近 BA，这对应亏量为零的低分辨率区；曲面度量与语义度量不成比例时，高分辨率码率损失趋于 $\Lambda[h]$，优化质量不能消除。把曲面限制为低维点积或单一欧氏弦长后才出现可测的码率代价。SEMASIA 上测得交叉差，当前单斜率数据还不能判定等失真收益。交互秩的两项检验只需 NumPy 与 SciPy。
