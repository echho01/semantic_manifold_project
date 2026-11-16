# 固定几何的数值实验与复现

[研究索引](../README.md) · [论文实验协议](../EXPERIMENTS.md)


本记录对应已运行的 CPU 实验：连续球面率失真积分、旋转曲面 BA、球面投影误差求积，以及两个同心圆的分量误判。统一入口为 [`run_paper1.py`](../code/synthetic/run_paper1.py)，完整记录为 [`results.json`](../code/results/paper1/results.json)。下文参数、数值和图表均来自该记录。

## 1. 复现命令与单位

在本工作的 research 目录执行，Python 版本不低于 3.10：

```bash
python -m pip install -r requirements.txt
python code/synthetic/run_paper1.py --output code/results/paper1
python code/tests/test_numerics.py -v
```

已运行环境为 Python 3.12.14、NumPy 2.3.5、SciPy 1.17.0，使用 CPU。JSON 保存软件版本、参数、源码 SHA256 和随机种子；CSV 保存供论文表格使用的数值；PDF、SVG 和 PNG 由同一绘图代码生成。

连续率失真率量使用 **nat/source symbol**，重构点约束在源流形上。密度相对于黎曼体积定义。`v` 表示每个环境坐标的噪声方差，`sigma_total² = D_a v`；中断曲线的横轴为每坐标标准差 `sqrt(v)`。球面半径和环面小半径取 1。

## 2. 连续球面率失真

对均匀源 $S^d(1)$，分别采用测地失真 $\theta^2$ 和弦长失真 $4\sin^2(\theta/2)$。旋转对称性给出

$$
Z(t)=|S^{d-1}|\int_0^\pi e^{-t\delta(\theta)}\sin^{d-1}\theta\,d\theta,
\quad D(t)=-\partial_t\ln Z(t),
$$
$$
R(t)=\ln\operatorname{Vol}(S^d)-\ln Z(t)-tD(t).
$$

求积先使用 $u=\sqrt t\,\theta$ 保持小失真积分窗口的分辨率。参数为 $d\in\{1,2,3,5,8\}$、$t\in\{500,1000,2000,4000\}$。令

$$
R_0(D)=h_g(p)-\frac d2\ln\frac{2\pi eD}{d},
$$

对 $[R(D)-R_0(D)]/D$ 关于 $D$ 作二次拟合，取零失真截距。测地与弦长系数的理论值分别为 $(d-1)/6$ 和 $(d-2)/8$；十种设置的最大绝对系数误差为 **$3.203\times10^{-10}$**。

原始数据：[`sphere_curve.csv`](../code/results/paper1/sphere_curve.csv)、[`sphere_coefficients.csv`](../code/results/paper1/sphere_coefficients.csv)。

## 3. 非恒定曲率与非均匀密度：旋转曲面 BA

环面小半径为 1，面积密度为 $p(\theta)\propto e^{k\cos\theta}$，四组参数为 $(R,k)=(2,0),(2,1),(2,-1),(3,0)$。另使用赤道半径 1.3、极半径 0.7 的均匀旋转椭球。

子午线使用 768 个求积点：环面在 $[0,2\pi)$ 上作周期求积；椭球在 $(0,\pi)$ 上作 Gauss–Legendre 求积，面积权重均为正。方位角积分解析完成：

$$
K_{ij}(t)=e^{-te_{ij}}I_{0e}(2tf_if_j),
\quad e_{ij}=(f_i-f_j)^2+(g_i-g_j)^2.
$$

其中 $I_{0e}(z)=e^{-z}I_0(z)$。BA 优化环质量，信息量按诱导输出边缘计算。环面的停止容差为 $10^{-12}$ nat，椭球为 $10^{-9}$ nat；最大最终拉格朗日对偶间隙为 **$9.963\times10^{-10}$ nat**。

环面使用 $t=100,200,400,800$，椭球使用 $t=200,400,800,1600$，采用与球面相同的残差外推。

| 源 | 理论系数 | 拟合系数 |
|---|---:|---:|
| 环面 $(R,k)=(2,0)$ | -0.07216878 | -0.07216878 |
| 环面 $(2,1)$ | -0.04575612 | -0.04575612 |
| 环面 $(2,-1)$ | -0.11662433 | -0.11662433 |
| 环面 $(3,0)$ | -0.06629126 | -0.06629126 |
| 椭球 $(1.3,0.7)$ | -0.05487468 | -0.05487467 |

理论系数为 $-\mathbb E_p[(\kappa_1-\kappa_2)^2]/16$。最大绝对拟合误差为 **$6.129\times10^{-9}$**。均匀 $R=2$ 环面在 $t=800$ 时把网格加密到 1536 点，率量变化为 $8.882\times10^{-15}$ nat，残差除以失真后的变化为 $1.421\times10^{-12}$。这些是求解间隙、网格变化和外推误差三种分别记录的量。

原始数据：[`surface_coefficients.csv`](../code/results/paper1/surface_coefficients.csv)；各曲线、迭代次数、输出边缘误差、网格加密结果见 JSON 的 `surfaces` 字段。

## 4. 球面投影与分量误判

球面投影用独立随机变量 $X\sim\mathrm{Gamma}(d/2,1)$、$Z\sim N(0,1)$ 表示角度：

$$
\Theta=\operatorname{atan2}(\sqrt{2vX},1+\sqrt v Z).
$$

每变量使用 80 个 Hermite／广义 Laguerre 节点，并检查加密到 120 个节点。对 $d=1,2,3,4,5$ 和 $v=0.001,0.002,0.004,0.008$，外推 $(\mathbb E\Theta^2-dv)/v^2$，其理论值为 $d(5-2d)/3$。最大绝对系数误差为 **$3.104\times10^{-6}$**；最大求积加密误差约为 $2.394\times10^{-16}$，单位为平方长度。原始数据见 [`projection_coefficients.csv`](../code/results/paper1/projection_coefficients.csv) 和 [`projection_curve.csv`](../code/results/paper1/projection_curve.csv)。

分量误判实验使用等概率、半径为 1 和 1.6 的同心圆，嵌入 $\mathbb R^{32}$。决策边界半径为 1.3，分量选择只依赖前两个坐标的接收半径。使用非中心卡方分布计算精确概率，再分别运行 400,000 个独立样本。每坐标噪声标准差为 `0.08, 0.10, 0.12, 0.15, 0.177, 0.20, 0.25`，对应种子 `20261019` 至 `20261025`。

| 每坐标标准差 | Monte Carlo 误判率 | 精确径向概率 |
|---:|---:|---:|
| 0.080 | 0.0000975 | 0.000090374 |
| 0.120 | 0.0063425 | 0.006355785 |
| 0.177 | 0.0460600 | 0.046216148 |
| 0.250 | 0.1188075 | 0.118475201 |

[`outage.csv`](../code/results/paper1/outage.csv) 同时保存失败次数、Wilson 95% 区间及 $Q(\Delta_0/(2\sqrt v))$ 径向近似。该近似与一般高斯范数安全上界分别解释；实验不把某一噪声尺度指定为所有流形的必然“悬崖”。

## 5. 论文图表与回归记录

| 图表／记录 | 文件 |
|---|---|
| RD 系数与非均匀曲面残差 | [`rd_validation.pdf`](../code/results/paper1/rd_validation.pdf)，同名 SVG/PNG |
| 投影二阶项与同心圆中断 | [`channel_validation.pdf`](../code/results/paper1/channel_validation.pdf)，同名 SVG/PNG |
| 全部参数、数值曲线及精度记录 | [`results.json`](../code/results/paper1/results.json) |
| 10 项 CPU 回归测试 | [`test_report.txt`](../code/results/validation/test_report.txt)、[`summary.json`](../code/results/validation/summary.json) |

回归测试覆盖连续系数、BA 间隙、高维局部几何、投影平移、监督探针秩、成对度量、量化器迁移及驱动分阶段恢复。表图所用结果不依赖预训练模型下载。

## 6. 可选扩展：真实嵌入诊断

[`code/embeddings`](../code/embeddings) 提供 CLIP、MiniLM/BERT、DINOv2 的提取与诊断入口，运行说明见其 README。该扩展尚未生成当前论文中的预训练编码器实验结果。

所有码本、探针和类别中心使用固定的 80/20 拟合／评估划分；不同编码器通过原始 item IDs 对齐。`rd` 测量固定码率 $\log_2 N$ 下的独立测试集环境 MSE，重构为环境空间质心；其拟合残差作为经验量化诊断，连续、流形内重构的 Shannon 系数另行报告。局部投影使用独立参考样本；局部曲率半径不作为全局 reach；类别文本提示保存为类别原型。实际运行后，驱动才会写入对应的 `arrays.npz` 和 `report.json`。
