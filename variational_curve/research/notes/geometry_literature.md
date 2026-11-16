# 曲面建模的参考文献对照

[研究索引](../README.md) · [几何文献详评](TCOM_REVIEW_20261004/geometry_papers.md) · [通信与数据文献](TCOM_REVIEW_20261004/communication_data_papers.md)

## 22 份文献的取舍对照

优先级针对当前曲面通信问题，不评价论文本身的价值。“可采用”均指注明出处后作为基础、基线或组件；与当前任务、信道或误差条件不一致的定理只能借证明思路。全部文献的模型、公式位置、定理条件、实验设置及原始数据清单分别见 [几何与表征附录](<TCOM_REVIEW_20261004/geometry_papers.md>) 和 [通信、数据与拓扑附录](<TCOM_REVIEW_20261004/communication_data_papers.md>)；另一个项目的逐项证据见 [外部项目审计](<TCOM_REVIEW_20261004/external_project.md>)。

| 文献（点击本地 PDF） | 优先级 | 可借用的建模或理论 | 当前项目的用途与限制 | 原文主要实验数据 |
|---|---|---|---|---|
| [Lossy Compression](<../../../refs/theory_overlap/2023_Information_and_Inference_Lossy_Compression_of_General_Random_Variables.pdf>) | 优先：理论基础 | 参考测度 Shannon 下界、平方失真量化尺度 M^(-2/m) | 当前率下界及有限码本限制；核对绝对连续、小球正则和弦距/测地距 | 球面、Grassmannian、自相似集合 |
| [RaDOGAGA](<../../../refs/ml_geometry/2020_ICML_Rate_Distortion_Optimization_Guided_Autoencoder_for_Isometric_Embedding_in_Euclidean_Latent_Space.pdf>) | 优先：编码基线 | 熵模型＋重构＋潜变量扰动，局部 scaled-isometry 条件 | 率感知等距正则；不能把理想局部推导写成训练后全局等距 | CelebA、3D→16D 混合源、异常检测集 |
| [Low Bending Distortion](<../../../refs/ml_geometry/o_2022_arXiv_Convergent_Autoencoder_Approximation_of_Low_Bending_and_Low_Distortion_Manifold_Embeddings.pdf>) | 优先：编码组件 | 近邻 stretch 与测地 midpoint 二阶差分、Hessian 能量极限 | 低弯曲基线；原收敛要求真测地数据、光滑网络及渐近全局极小者 | sundial 半球、cow/SO(3)、Klein bottle、ellipse |
| [Implicit Latent Geodesics](<../../../refs/ml_geometry/o_2026_ICML_Geodesic_Calculus_on_Implicitly_Defined_Latent_Manifolds.pdf>) | 优先：接收组件 | id−投影定义支撑，离散测地能量与增广 Lagrangian | 有限步几何接收；去噪 posterior mean 不等于精确最近投影 | SCAPE、CMU MoCap、旋转 cow |
| [Riemannian Neural OT](<../../../refs/ml_geometry/o_2026_ICML_Riemannian_Neural_Optimal_Transport.pdf>) | 优先：受限理论 | c-transform、有限输出位置逼近障碍、局部 argmin 稳定性 | 码本维数斜率及接收误差证明模板；不是原样适用的噪声通信 converse | 球面、平环面、SO(3)、紧支撑 SE(3) 实验、古大陆 |
| [SEMASIA](<../../../refs/ml_geometry/o_2026_arXiv_SEMASIA_A_Large_Scale_Dataset_of_Semantically_Structured_Latent_Representations.pdf>) | 优先：主实验 | 同图多模型特征；whitening、LS/SVD、CCA 等强线性基线 | 异构接收端及真实比特通信；必须验证样本配对，维度不是码率 | CIFAR10/100 等八个图像集合 |
| [Metric Matching](<../../../refs/ml_geometry/o_2026_ICML_Riemannian_Metric_Matching_for_Scalable_Geometric_Modeling_of_Distributions.pdf>) | 次选：几何估计 | 去噪二阶矩回归、低秩 PSD 切算子与极限 | 切向投影/局部几何估计；切算子不是任务 g，先核查局部版本的 1/2 归一化 | S^8→R^64、MNIST、CIFAR10、CelebA、FFHQ |
| [When Scores Learn Geometry](<../../../refs/ml_geometry/2026_ICLR_When_Scores_Learn_Geometry_Rate_Separations_under_the_Manifold_Hypothesis.pdf>) | 次选：机制理论 | 小噪声下支撑几何与密度恢复的不同误差阶 | 解释 geometry-only/density-only/joint；不直接给有限步接收保证 | ellipse＋von Mises；文本到图像初步实验 |
| [Learning Eigenstructures](<../../../refs/ml_geometry/o_2026_CVPR_Learning_Eigenstructures_of_Unstructured_Data_Manifolds.pdf>) | 次选：谱基对照 | 谱 minmax、探针 PCA/QR 学习低频表示 | 低计算量平滑或压缩；零模态、探针协方差尺度需处理，本地补充证明缺失 | STL10、Imagenette、CIFAR100、Caltech256 特征 |
| Statistical Point Cloud Manifold | 后续：点云路线 | 点集 KDE、Fisher 统计度量、decoder pullback | 点云任务与信道 Fisher 建模；点云样本空间不同于物体的二维表面 | 合成形状、ShapeNetCore v2、ModelNet10/40 |
| [Robust Perceptual Geodesics](<../../../refs/ml_geometry/2026_ICLR_Counterfactual_Explanations_on_Robust_Perceptual_Geodesics.pdf>) | 次选：参考特征 | 冻结多层鲁棒特征度量及生成器 pullback | 作为 ψ 或感知对照；PSD、像素攻击鲁棒性与任务后验保证分开 | AFHQ、FFHQ/CelebA、PlantVillage |
| [Readout Representation](<../../../refs/ml_geometry/2026_ICLR_Readout_Representation_Redefining_Neural_Codes_by_Input_Recovery.pdf>) | 背景：任务接口 | 固定 readout 的 fiber、可读出语义与隐表示自由度 | 明确哪些方向任务无关；最大成功位移不是任意噪声安全半径 | ImageNet/C4 小子集及 CLIP、DINO、BERT 等 |
| [Riemannian Deep Learning](<../../../refs/ml_geometry/2026_arXiv_Riemannian_Deep_Learning.pdf>) | 背景：几何工具 | SPD、相关矩阵、旋转等已知流形上的优化模块 | 需要已知矩阵/旋转域时借组件；不是未知语义流形定理，ALEM 重根定义需审慎 | HDM05、FPHA、NTU60/120、G3D 等 |
| [Domain Learning](<../../../refs/semantic_rd/2024_arXiv_Autoencoder_Based_Domain_Learning_for_Semantic_Communication_with_Conceptual_Spaces.pdf>) | 次选：原型对照 | 监督 AE 与 EMA 类别原型/Voronoi 域 | 廉价任务原型或图册初始化；没有通信或全局几何定理 | MNIST、CelebA 两属性 |
| [Causal Abstraction / SEP](<../../../refs/ml_geometry/o_2025_ICML_Causal_Abstraction_Learning_based_on_the_Semantic_Embedding_Principle.pdf>) | 次选：线性对齐 | 右逆接口、Stiefel 子类、Gaussian covariance KL | 正交 adapter 对照；边缘协方差对齐不保证后验、因果或任务兼容 | 合成 SCM、HCP 聚合 fMRI |
| [Pragmatic Information](<../../../refs/semantic_rd/2026_arXiv_A_Mathematical_Theory_of_Pragmatic_Information.pdf>) | 背景：任务动机 | Bayes 行动等价类、价值与成本 | 直接行动传输基线；谨慎使用自定义 MI/capacity，类内奖励并非自然相等 | 离散/Gaussian/POMDP 及交通示例 |
| [VoC Pareto Analysis](<../../../refs/semantic_rd/2026_TIT_Value_of_Communication_in_Goal_Oriented_Semantic_Communications_A_Pareto_Analysis.pdf>) | 优先：评价方式 | 通信频率与远端任务成本 Pareto 曲线 | 以实际成本—任务表现评价；静态网络不能照搬 Markov/CMDP 凸性定理 | 八状态 synthetic alarm source |
| [Token Embeddings](<../../../refs/ml_geometry/o_2025_NeurIPS_Token_Embeddings_Violate_the_Manifold_Hypothesis.pdf>) | 次选：假设检验 | 多尺度邻域质量与维数斜率检验 | 流形假设的负对照；拒绝某光滑模型不等于否定全部流形 | GPT2、Llama、Mistral、Pythia token vocab |
| [TopoCode](<../../../refs/topology/2024_arXiv_Topocode_Topologically_Informed_Error_Detection_and_Correction_in_Image_Communication_Systems.pdf>) | 相邻：接收基线 | 发送 PD 签名，接收端拓扑辅助纠错 | 结构侧信息基线并计保护开销；PD 稳定性没有任务风险反向界 | MNIST、Omniglot |
| [Topological Autoencoders](<../../../refs/topology/o_2020_ICML_Topological_Autoencoders.pdf>) | 相邻：正则基线 | H0/MST pairing-index 双向距离损失 | 低成本拓扑 ablation；零训练损失不保证全局单射或全部拓扑 | Nested Spheres、MNIST、FashionMNIST、CIFAR10 |
| [Topological Flow Matching](<../../../refs/topology/o_2026_ICLR_Topological_Flow_Matching.pdf>) | 另题：图信号 | 已知复形的 Hodge 谱、热核 Gaussian bridge | 若做图域信号接收再使用；已知图结构不等于未知信源曲面 | PeMSD4、IRIS、HCP、单细胞、ocean/NOAA |
| [Topology-Preserving JSCC](<../../../refs/topology/2026_CL_Topology_Preserving_Deep_Joint_Source_Channel_Coding_for_Semantic_Communication.pdf>) | 相邻：强通信对照 | 标准 JSCC 加图像/潜表示 PH 正则、AWGN/Rayleigh | 限定“几何保护 JSCC”已有边界；本地版本没有新 RD/可靠性定理 | Omniglot、DeepGlobe road masks |

建议首先组合 **Lossy Compression＋RaDOGAGA/LBD＋Implicit Geodesics＋SEMASIA**。前者给理论基础，中间两类给编码与接收组件，最后给真实端模型差异。RNOT 的有限输出位置与稳定性证明用于限定可实现类；Metric Matching 和谱基方法先作为组件对照，确认收益后再增加复杂度。四篇拓扑工作用于界定相邻工作，避免把曲面稿的主线扩成另一篇拓扑论文。

三个已有贡献不宜换名再申报：RD 引导等距编码（RaDOGAGA）、低弯曲 AE（LBD）、几何或拓扑正则的 JSCC（Topology-Preserving JSCC）。真正需要自己完成的是它们进入固定任务评价与有限通信接口后，新的误差界、代价约束和公平实验。
