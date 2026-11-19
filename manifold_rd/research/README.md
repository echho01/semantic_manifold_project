# 工作一研究索引 · 语义流形率失真

[Transactions 主稿（拟投 IEEE Transactions on Communications）](../main.tex) · [论文 PDF](../build/main.pdf) · [参考文献](../ref.bib) · [实验协议、运行矩阵与记录](EXPERIMENTS.md)

论文与仿真的分工（知识库 $\mathcal K=(\theta,K)$ 在主稿中如何定义、在实验中如何实例化）见[实验记录开头的对照表](EXPERIMENTS.md)。研究按“流形建模 → 语义失真 → 语义率失真”组织，给出广义、知识库条件和任务导向三层率失真及其相互关系。接收端必须输出可实现对象还是可以输出预测分布，是其中“合法输出”一层。主线的取舍和对 Word、开题 PPT 的逐条判断见[批判性对照](notes/06_framework_critique.md)。

| 推导 | 在论文和实验中的用途 |
|---|---|
| [流形建模、语义失真与三层语义率失真](derivations/framework.md) | 主线：语义度量与校准、信道诱导失真、广义率失真及曲率项、知识库条件率失真与增益分解、任务秩律与联合判据 |
| [任务、共享知识与局部度量](derivations/task_source.md) | 观测、后验、合法输出、Fisher 拉回和共享上下文的失真分配 |
| [后验质心、投影与同码率分离](derivations/receiver_gap.md) | 圆的零率精确值、正码率下界、强凸层集推广与径向投影上界 |
| [连续流形的小失真极限](derivations/geometric_bounds.md) | Gibbs 体积积分、对偶下界、显式可达核、极限交换条件 |
| [有限问题的证书与信道](derivations/finite_rate.md) | BA 更新、时分与对偶核验、权重敏感性、信息率和物理载荷 |
| [损失预算、对象代价常数、任务秩律与信道极限](derivations/structure_and_limits.md) | 五项对数损失恒等式、凸性刻画、reach 上界、混合曲率局部常数、$\tfrac r2\ln(1/D)$ 高分辨率律及度量/知识推论、链路逆定理、中断与模拟映射功率分配 |

| 已有结果 | 实现 | 观测记录 |
|---|---|---|
| 24 点后验圆和两次实高斯信道使用 | [数值复现](toy_geometry_check.py) · [率失真求解器](manifold_rd.py) | [有限界](results/finite_rd_bounds.csv) · [零率](results/toy_zero_rate.csv) · [物理信道](results/gaussian_link.csv) |
| 齐性流形的精确率失真、四种知识情形的增益、环面上任务/联合/广义三条曲线、模拟映射的信道诱导失真 | [三层率失真的数值核验](framework_checks.py) · [单元测试](tests/test_theory.py) | [全部数值](results/framework_checks.json) · [主稿图 3](../figures/framework.png) |
| 损失预算枚举、720 态等信息间隙、任务秩渐近线、容量处失真下限、知识节省 | [结构结论的数值核验](theory_checks.py) · [单元测试](tests/test_theory.py) | [全部数值](results/theory_checks.json) · [主稿图 4](../figures/gap_highres.png) |
| 两个骨干、8,192 个官方测试 ID、256 个锚点；无知识 / 共享粗类 / 共享细类三种条件 | [真实特征上的知识增益](knowledge_cifar100.py) · [单元测试](tests/test_knowledge.py) | [增益区间与嵌套界](results/knowledge_cifar100.json) · 主稿表 `tab:kreal` |
| 两个骨干、$c=2,3,4$ 类任务的后验点云，$n=4096/8192$ | [真实后验上的秩律](rank_cifar100.py) | [斜率](results/rank_cifar100.json) · [8192 点](results/rank_cifar100_8192.json) · 主稿表 `tab:rank` |
| 两个骨干、2,048 个官方测试 ID、固定索引载荷 | [SEMASIA 配对实验](semasia_cifar100.py) | [校准、代理 KL、分类准确率与跨模型比较](results/semasia_cifar100.json) |
| 有限码本、测地与弦长对照、真实索引经 AWGN/Rayleigh | [训练和测试](experiments/semasia_rd.py) · [批量运行](experiments/run_matrix.py) | [短程验证与检查点哈希](results/smoke_provenance.json) |

下列命令从本目录运行。GPU 服务器先安装适配驱动的 PyTorch，再安装研究依赖。
~~~sh
python -m pip install -r requirements.txt
python tests/test_information.py
python toy_geometry_check.py --grid 24 --output runs/circle --figure ../figures/toy_posterior_manifold.png
python tests/test_theory.py
python theory_checks.py --output runs/theory --figure ../figures/gap_highres.png
python framework_checks.py --output runs/framework --figure ../figures/framework.png
python tests/test_knowledge.py
python knowledge_cifar100.py --data data/semasia --download
python experiments/semasia_rd.py train --data data/semasia --download --device cuda --book 64 --beta 2 --checkpoint runs/rd.pt
python experiments/semasia_rd.py test --data data/semasia --download --device cuda --count 2048 --checkpoint runs/rd.pt --output runs/test.json
~~~

训练时可将代价改成 `--criterion graph_geodesic --neighbors 12` 或 `--criterion task_chord`，分别保存检查点。更换骨干用 `--model mobilenetv3_small_100.lamb_in1k`。参考资料包括[开题与 Word 的判断依据](notes/orientation_at_migration.md)、[语义率失真文献](notes/01_semantic_rate_distortion_review.md)、[机器学习几何文献](notes/02_ml_geometry_interpretability_review.md)和[理论近邻](notes/05_theory_and_closest_prior_work.md)。
