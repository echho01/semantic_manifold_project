# 几何与表征论文审阅：面向 TCOM 的可复用模块、条件与待证通信贡献

审阅日期：2026-10-04。以下页码统一为本地 PDF 页序，从封面/第一页计数；专著另标印刷页码。对照材料是本项目 `variational_geometry.tex` 与 `sections/communication_geometry.tex`。这是综合评审的技术附录，未修改论文源文件。

## 0. 覆盖范围与总判断

八篇研究论文的主体建模、理论陈述、证明中与建议有关的部分、实验和数据配置均已审阅；参考文献用于追踪来源，定性图页主要检查题注。`Learning Eigenstructures` 本地 14 页版本没有它反复引用的 supplementary，因此不能声称核验了补充证明。`Riemannian Deep Learning` 实际为 424 页博士专著：已遍历目录、全书标题/关键词，重点审阅 PDF47–54 的几何基础、PDF216–225 的 Euclidean pullback/SPD 几何、PDF235–243 的 Cholesky 几何、PDF273–275 的数据配置，以及 PDF317 的 ALEM well-definedness 讨论；没有逐页通读无关章节。

**最重要的判断：现有固定任务语义度量 g、可学发送表示 F 和 h=F*δ 的区分应保留。** 这些论文可以提供几何估计、支撑约束、统计度量、隐式优化和有限码本下界，但没有一篇证明“任意改写曲面即可降低同一任务的无限制 Shannon RDF”。当前主稿已经明示 representation invariance，这是正确边界。96 点非均匀圆周上 joint 与 density-only 基本相同，6.47% 改进主要是复现 reproduction-density/BA 效应，geometry-only 约 0.12%；不能据此把收益归因于新曲率理论。

**可以合法直接复用并引用的部分**是已有公式、假设、算法和已证明定理的适用特例；应明确写成“采用/基于某论文”。**能形成自己的贡献的部分**应是加入通信对象后重新定义且证明的结论，例如有限码本/功率/信道可辨识性/有限接收计算量下的语义风险。重新命名原论文 metric matching、c-transform、Fisher pullback、readout fiber 或投影算法，不能变成自己的原创。

## 1. Riemannian Metric Matching for Scalable Geometric Modeling of Distributions（20 页）

### 模型和可直接引用的结果

- PDF3 的 carré du champ：Γ_L(f,h)=1/2[f Lh+h Lf−L(fh)]，通过局部随机核得到 ambient 几何矩阵。PDF4 Eq.(11) 将矩阵学习化为条件去噪二阶矩回归，目标写成 `(X−Y)(X−Y)^T/(2ε)`，`Y|X ~ N(X,εI)`；Thm.3.1 的 conditional/marginal loss 只差不依赖参数的常数（证明 PDF12）。这属于条件期望与平方损失分解，可引用直接实现。
- 用 Γθ=Mθ^T Mθ，Mθ∈R^{r×D} 表示低秩 PSD 算子；PDF4 Eq.(12) 的去常数损失为 `E||MM^T||_F² − ε^{-1}E||M(X−Y)||²`。无需生成 D×D 全矩阵，是实际值得采用的工程模块。
- PDF5 Thm.4.1/Cor.4.2：支撑在光滑流形、密度光滑且 ε→0，CDC 对光滑测试函数收敛为梯度内积，ambient 坐标下对应切空间 projector；Cor.4.3 是链式法则。这里 ambient projector 是退化 cometric/切算子，不能当成某坐标中可逆的固定任务 metric g。
- PDF5 Thm.4.4（证明 PDF13–14）：r≥2m−1 的 smooth factorization 可以实现所需切算子，来自切丛的稳定嵌入/Whitney 类型构造；r=m 在全局可能受拓扑阻碍。结论是几何表达能力，不能直接读成有限样本、固定神经网、固定噪声误差保证。

### 与本项目的接口和边界

在冻结的任务特征 U 上学习切空间，可用于接收优化的 tangent projection、几何诊断、局部 intrinsic dimension/normal noise 分离。固定 g 仍来自冻结校准的 q 与 ψ。训练得到的 PSD 切算子可能奇异，不应直接替代 g^{-1} 或声称保留本稿 posterior-TV 保证。需要估计 `||Pθ−PT||`、正则化、训练噪声尺度与测试信道尺度的不同，并在自己的有限步接收风险中处理这些误差。

**本地版本一致性疑点（只针对 normalization，需平坦 sanity check）：** PDF3–4 原页已视觉核验：核 `exp(−||x−y||²/(2ε))`、回归目标 `/(2ε)` 与极限 tangent projector 的尺度需统一。均匀平坦子空间、y 在支撑上时，条件协方差为 εPT，该目标给出 PT/2；应以平坦高斯算例校准损失、生成元和极限，而非照搬未校准系数。PDF3 还同时出现 `(I+Pε)/ε` 与 `(I−Pε)/ε`，实现应依据期望生成元符号检查。此疑点不否定其低秩学习思路。

### 数据和成本

PDF6–8：S^8→R^64 的合成切投影真值最适合先复验；MNIST、CIFAR-10、CelebA、FFHQ 是高维展示。附录 PDF17–18 的 CelebA rank128 / FFHQ rank256 还需要独立 posterior-mean denoiser，并有数十至数百 epoch。FFHQ/CelebA 图像上的 Inception feature 稳定性是 proxy，没有真实 tangent/任务 metric 真值。两张 A800 条件下优先在冻结的 512/768 维特征上训练低秩网络，避免一开始复现 raw1024 图像级 Jacobian。

## 2. Geodesic Calculus on Implicitly Defined Latent Manifolds（23 页）

### 模型和可复用数学

- PDF3–4：`Z={z:ζ(z)=0}`，ζ=id−Π，其中 Π 理想上是流形近邻投影，`TzZ=ker Dζ(z)`。支撑约束与 metric 分开：这是本项目最有用的接口，可以固定任务 g，同时把接收搜索限制在可信语义支撑上。
- PDF4：离散路径能量 `EK=K Σ W(z_{k−1},z_k)`，要求 `W(z,z')=dist²(z,z')+O(dist³)`。其离散测地线收敛引用 Rumpf–Wirth (2015) Cor.4.10，依赖光滑 Z、Lipschitz metric、正确局部能量以及精确约束；不是论文为学得投影网络新证明的统一收敛。
- PDF4 augmented Lagrangian：`EK−Λ:ζ+μ||ζ||²/2`，交替优化路径、乘子与 penalty；Dζ 秩为 ambient−intrinsic dimension 时 multiplier 并不唯一。PDF5 离散 Exp 用局部 optimality residual 求解，可减少 projection 网络二阶微分。
- PDF5 Eq.(13)，证明/讨论 PDF15：平方去噪训练最优 Πσ(y)=E[Z|Z+noise=y]。均匀平坦子空间上是正交投影；光滑流形内部通常误差 O(σ²)，边界可能 O(σ)。非均匀源密度会产生切向 drift。posterior mean 不是精确 nearest projection，有限训练集也不是已证明光滑流形。

### 可以形成通信扩展的缺口

PDF9 明确未完成学习误差、有限采样、σ、路径节点 K 的联合误差分析。我们的接收器若加入 `ζθ(u)`，需自己证明：约束误差/metric 误差/有限迭代如何影响固定 `d_g²` 风险，并处理信道后 `A_H F` 的局部和全局可辨识性。不要把已有 geodesic convergence 推广为“任意训练好的 Πθ 保证最短路径/全局 MAP”。固定 g 的二次 Bayes receiver 与 MAP 也不可混称。

### 实验可行性

PDF6–8、附录15–17：SCAPE（71 pose，mesh 原 12500 vertices，下采样1250，PGA 后合成10k/100k训练点）；CMU MoCap subject86 trials1–6，约52k frames，SO(3)^30，SVAE ambient latent10；30k 张128²三维 cow 旋转图像，有 SO(3) 角度真值，latent16/PCA6。**推荐 cow/SO(3) 自生成旋转数据做受控真值实验，CMU 或 NTU 骨架做现实任务。** 对序列数据应按 subject/sequence 划分，不能照搬相邻 frame 随机80/20以避免泄漏。PDF19 和实时例子中每个路径可能耗秒至分钟；需要 amortized 或固定 K 接收器，并报告实际 latency。

代码：[LatentGeodesics](https://github.com/flrneha/LatentGeodesics)。

## 3. Riemannian Neural Optimal Transport（60 页）

### 模型和核心定理

- PDF2 假设 `(M,g)` 光滑、紧、连通、无边界，**intrinsic distance、Exp/Log 与 gradient 已可计算**；µ≪vol 时，McCann 既有定理给二次成本 `c=d_g²/2` 的 a.e. 唯一 Monge map `T(x)=Exp_x(−∇ϕ(x))`。
- PDF3 Ass.2.2 要 feature φ 连续且单射；充分稠密 distance-to-landmarks 的 injective embedding 是既有 Gromov 构造。PDF5 令 `ϕ=(fθ◦φ)^c`，通过 c-transform 直接满足 c-concavity，优化 Kantorovich semi-dual；每个样本仍需内层 manifold minimization。可引用并用作**固定源和重构分布间的 compander/初始化**，不是直接产生 bits 的码器。
- **PDF4 Thm.3.1 /证明 PDF33–34**：µ、ν≪vol，T* 推µ到ν，任何输出不超过 M 个位置的 map 有 `RMSEµ(T,T*)≥C M^{-1/m}`。因此 squared loss≥C²M^{-2/m}。这是最适合引用的有限码本维数障碍；M=2^R 给出 `D≥c 2^{-2R/m}` 的相关量化斜率（常数依赖 manifold/target density）。需要明示这是有限重构点的量化/OT逼近结果，不是所有 noisy-channel decoder 的原封不动 converse。
- PDF5 Thm.4.1：网络 prepotential 一致逼近、feature 单射、µ ac，得到 transport map µ-a.e./in-probability 收敛。不是逐点 uniform 或有限训练收敛。
- PDF6 Thm.5.1：ψ∈C^{km,1}，存在适合的 φ*∈R^{2m} 和网，使 parameter `Wε=O(ε^{-4m/[3(km+1)]})`、depth `O(ε^{-2m/[3(km+1)]})`。所谓摆脱维数灾难以随 m 增长的光滑阶数为前提，隐藏几何常数、feature 构造及内层求解成本不等于端到端 polynomial compute。
- PDF6–7 Thm.5.3（证明 H.6）：prepotential sup-error≤ε，在 full-µ-measure 集每个 x 有常数 C(x)、阈值 εx，任意对应 argmin 满足 `d(Tε(x),T*(x))≤2sqrt(ε/C(x))`。可以参考 quadratic-growth argmin stability 的证明套路，为接收近似能量建立自己的风险界；必须明确局部 gap/curvature 和无远处分支混淆，不能换成全局统一 C 后直接引用。

### 对本项目的边界

连续神经 map 没有离散输出障碍，**但数字通信的有限 bit/codebook 仍然有量化障碍**。OT 的 µ→ν density shaping 与我们的 reproduction density r 学习、BA 更新十分接近；这个部分应作为可引用编码模块/基线。真正新贡献可以是限定实际代码族（有限M、几何受限码本、功率、CSI、有限计算）后的联合设计，而不是“改曲面降低无约束 RDF”。

### 数据与算力

PDF7–9：S² 古大陆分布、球面 S^m、环面 T^m 的 uniform→wrapped normal，另附录 SO(3)、compact-support SE(3)。这些很适合和当前 circle/sphere 数值轨道合并做维数斜率检验，但地质 OT 本身不匹配语义通信任务。PDF8 球面 RNOT 单次约986s，而 RCPM约37s，较低维并不省算力；KL evaluation 还假设 test samples 处 map 局部可逆。用训练时间、内层 iterations、全局 injectivity 检查解释实验，不能只展示 KL。

代码：[RNOT](https://github.com/MLGlobalHealth/riemannian-neural-optimal-transport)（本地 PDF2 声明 MIT）。

## 4. Learning Eigenstructures of Unstructured Data Manifolds（14 页）

### 模型、现成理论和使用方式

- PDF3 Thm.3.1 针对 SPD L 的 ellipsoid `CL={f:f^TLf≤1}`：前K本征向量给 worst-case K-term squared reconstruction error `1/λ_{K+1}`，有 simple spectrum 时 uniqueness。该谱 minmax 结果已有来源，不能声称自己的新通信定理。若 L 是连通闭流形 Laplacian 有零常数模态，需要加正则或限制零均值子空间，不能直接套SPD条件。
- PDF3 Thm.3.2：ellipsoid probe 的 PCA 能恢复同一 eigenbasis。**本地版本尺度一致性疑点：**文中对均匀 ellipsoid 写 `E ff^T=L^{-1}`；通常均匀 unit ellipsoid covariance 为 `L^{-1}/(n+2)`。eigenvectors 不受共同系数影响，但 eigenvalue/error 校准必须检查探针分布与归一化。
- PDF3–4 从 `Lnorm=M^{-1/2}SM^{-1/2}` 的首本征向量恢复 diagonal mass `M=diag(q1²)`；这只给离散 mass，不给一般完整 metric、scalar curvature 或固定 task metric。
- PDF4–5 学 Φθ(P)∈R^{N×K}、QR正交化，在高斯平滑 random probes 上最小化 progressive projection loss。probe ensemble 只“类似”Laplace 的 spectral statistics，没有任意未知 manifold 的严格算子一致性。

### 与本稿的接口

可以做谱基表示/geometry smoother 的廉价近似。若固定 g、采用其真实 eigenbasis 且 µ∝dVg，`||Δg F||²` 可在谱上写成 Σλk²||Fk||²；非均匀材料源 µ=p_g dVg 时会有加权耦合，不能原样写对角式。若 probe 是任务信号族，其有限K误差只描述这一族的压缩，需要重新建立通信噪声/任务风险关系。我们应把它放在计算量对照或 ablation，不作为文章主创新。

### 数据与证明覆盖

PDF6 oracle mesh，PDF7–8 STL10/Imagenette/CIFAR100/Caltech256，DINOv2 768 / CLIP512 frozen features，1500图subset与50runs，评价 NMI/ARI/embedding。可复用 CIFAR100 frozen feature 实验；聚类更好不是通信失真或频谱节约的证据。本地无 supplementary，原文指向的部分证明未核验。

代码：[learning-eigenstructures](https://github.com/royvelich/learning-eigenstructures)。

## 5. A Statistical Manifold Framework for Point Cloud Data（25 页）

### 建模与可引用公式

- PDF2–3：每个点云是 n 个互异点的无序配置，空间维数 nD；point collision 会导致 quotient/singular 问题。固定 SPD Σ 的 KDE `p(x;X)=1/(n sqrt(detΣ)) Σ_i K(Σ^{-1/2}(x−xi))`；Gaussian 对应 Σ=σ²I。
- PDF3 Prop.2.3 configuration space 流形；Prop.2.4 在严格正定 translate kernel 条件下，点云→KDE density 单射（证明 PDF13）。**单射本身不自动证明 immersion/Fisher positive definiteness**；Gaussian 的导数可辨识性或 decoder 满秩应另外核对。
- PDF3–4 Fisher `H_{ij,kl}=∫p ∂_{Xij}logp ∂_{Xkl}logp`；Prop.2.5 的 Gaussian 展开是责任权重乘位置残差的积分（证明 PDF14）。PDF4 decoder pullback `G(z)=Jf(z)^T H(f(z)) Jf(z)`；通过随机 v 的 JVP估计 `||v^TGv−c v^Tv||²`，不显式形成巨大的 H。
- Eq.(12) `E||G−cI||F²`，c=EtrG/m 的 scaled-isometry regularizer 可以作为对照，但若无重构/满秩/scale约束，常数decoder给 G=c=0，正则也是0。可学metric自由缩放不是任务失真进步。

### 自己能做的通信部分

可把“Gaussian observation likelihood 的 Fisher pullback”接到实际 channel：观测 `Y=A_HF(u)+N`，若 Σ固定则局部 channel Fisher `K_H(u)=JF^T A_H^T Σ^{-1}A_H JF`。这条公式是标准 Fisher 链式法则；原创应是**固定 task g 下有限power/CSI/decoder预算的 trace(gK_H^{-1})风险、最优分配与 global error/outage 处理**，需重新证明，不能把 point-cloud metric结论换标签。点云 KDE Fisher 描述采样分布，可作为参考 ψ 几何；它不是类标签/任务 posterior Fisher，因此不能替代 q 自动保留task-TV定理。

### 实验与数据

PDF5–8、PDF15–17：cylinder/cone/ellipsoid/box的512points+2Dlatent；真实 ShapeNetCore.v2 AE512latent，再在 ModelNet10/40分类，2048点、unit sphere，DGCNN/FoldingNet/PointCaps，500epochs batch16。文章源噪声（bounding box 1/5/10/20%）不是链路 SNR。两A800可先复用标准点云 encoder/decoder与冻结分类q，然后新增真实bits/channel-use/power、AWGN+fading、SNR变化与端到端分类/Chamfer/任务风险。ModelNet40比直接开始ShapeNet预训练更快；单纯增加Gaussian扰动再报SVM准确率不足以成为TCOM实验。

代码：[SMF-public](https://github.com/seungyeon-k/SMF-public)。

## 6. Counterfactual Explanations on Robust Perceptual Geodesics（36 页）

### 方法和复用位置

- PDF5–6：冻结 robust backbone 多层 feature `h_k`，image-space PSD metric `GR=Σk Nk^{-1} Jh_k^T Jh_k`，generator pullback `GZ=Jgen^T GR Jgen`；path energy `1/2∫ γ'^TGZγ'`。离散能量用相邻 decoded robust features 的差分，避免显式 full Jacobian。这可复用为参考 feature ψ 或 perceptual geometry baseline，并明确引用。
- PDF6 的 source→exemplar固定端点，然后联合endpoint分类损失、re-anchoring、midpoint densification，是实际优化 heuristic。没有正式 robust geodesic/task保证或全局最优定理。PDF20 用 Softplus 代替 ReLU 改善平滑性，应重新核对固定 q/feature 的校准和语义。
- PDF8 本地 metric distance 用局部二次形式近似，只能解释小扰动；不能拿它代替任意两点真 geodesic。

### 边界与风险

Jfeature较小的方向给**较低**quadratic metric cost；鲁棒模型抑制某方向不自动意味着该方向在 metric 中代价大（PDF19相关语言需要按实际公式理解）。GR/GZ可能rank-deficient，必须有 η reference immersion 或验证fullrank。robust classifier对像素攻击有经验鲁棒性不等价于通信信道/任务posterior fidelity定理。冻结robustmetric自身可以引用；我们现有 calibratedposterior Fisher+ηψ 才有明确task-TV连接。不能把“短perceptual路径、看起来真实”当成可达语义rate降低。

### 数据与成本

PDF6–8：AFHQ动物3类、FFHQ面部属性（task classifier用CelebA）、PlantVillage病叶。PDF17–22 robust L2 ResNet50ε=3，另外VGG19task网络与robustInception/VGG16评估，避免只用同一个backbone循环打分；仍然需要标签/真实任务测试。AFHQ/PlantVillage自行训练StyleGAN2约140 H100 GPUhours；PDF29样本path分钟级，最高分辨率约7.8min。主论文应优先冻结encoder+amortizedreceiver，path只做定性解释或小规模诊断。PlantVillage是清晰边缘感知通信任务备选，需leaf_id分组并承认实验室叶片背景与真实田间有domaingap。

## 7. When Scores Learn Geometry（34 页）

### 可以直接引用的理论及条件

- PDF3–4 Ass.2.1/2.2：population distribution支撑在紧、无边界C4嵌入流形，正C1density；不是任意离散样本/实际CLIPfeature已满足的结论。
- PDF5 Thm.3.1，严格式PDF18 Thm.B.2：在固定小tube内，`log pσ(x)=−dist(x,M)²/(2σ²)+log pdata(u)−(D−m)/2 log(2πσ²)−1/2 log det Hhat(u,x)+o(1)`。领先normalgeometry为σ^{-2}，density为O(1)；真实metric与坐标density的volume因子一起变化。可引用解释geometry/support与density学习分工。
- PDF6 Thm.4.1：若已恢复分布的**score sup error**在某紧、统一rectifiably pathconnected K 为 `o(σ^{-2})`，可得到支撑趋向M；O(1)误差仍允许沿M任意错density；o(1)才恢复原density。额外假设πσ集中在K。**这不是network训练L2 loss或有限 reverse-SDE sampler 的保证。** PDF20–22证明通过pathintegral控制logdensity差，再用Laplace；条件不能省。
- PDF7 Thm.5.1 gradientfield case：`error=o(σ^β)`，β>−2，选 `max{−β,0}<α<2`，tempered Langevin `dX=σ^α sσ dt+sqrt2 dW` 的stationary law趋向uniformRiemannianvolume。PDF8 Thm.5.2非gradientcase额外要求pdata C2、唯一stationarydistribution及localWKBansatz；PDF22 Ass.B.2明确假设stationarydensity有正C2prefactor收敛，并集中tube，不是无条件证明所有neuraldiffusion有此形式。
- PDF8 Thm.6.1：加有界v及tube C1likelihoodpotential，limitingdensity∝exp(−v)dVM。接收Gaussian `v=||y−AF(u)||Σ^{-1}²/2` 在compactM有限，可构造boundedextension，但ambient全空间原v并不有界。更关键：结果**采用uniformmanifoldprior**，我们材料源µ通常非均匀，需要保留真正prior，不能去除density后宣称meanrisk必降。

### 实验与本项目

PDF9/29–30 ellipse+vonMises、transformer小模型；StableDiffusion1.5家具/汽车/建筑三类prompt仅preliminary分布/CLIP多样性。数值Langevin有10000steps、10000runs，通信接收预算应另设有限K。PDF10承认未分析累计trajectoryerror、L∞→L2转换、samplecomplexity和discretization。此论文最有用的是为我们保留**geometry-only/density-only/joint**对照提供理论解释；不要把diffusionstationarylimit用作有限时延TCOMreceiver定理。

## 8. Readout Representation（46 页）

### 概念与可用方法

- PDF4 Def.1/2：固定readoutπ下，representationfiber `Hs^π={h:π(h)=s}`。这是依赖decoder的等价类，能很好解释同一语义可有冗余表示。可引用作为固定receiver任务成功集合的来源；“fiber概念”不是我们的原创。
- PDF6 放宽为输出与原始x差<t，并用representation-size `r_x=max_{h∈Hx,t} dH(h,f(x))`。这是**最大成功位移**：并不意味着半径r内全部扰动成功，也不测量集合volume/probability。不能转写为worstcase鲁棒radius、Gaussiannoise容限、容量或固定powergain。
- π通过preimage优化或learnedreadout；强生成prior的输出可能hallucinate，需测真实label/posterior/taskdistortion。PDF6使用correlationdistance，脚注承认不是metric；PDF19feature噪声方差校准目标相关距离，不是具有统一发送power和channel-use的物理SNR。

### 我们应重新建立的量

定义固定readout的任务成功域 `S_t={y:decision(y)=t}`，度量真 `Pr(Y∈S_t|u,H)`、taskregret、posteriorloss。在power约束下，发送embedding能改变这些集合的channelmargin和分支间距，这是通信价值；自由地加入维度/缩放不应成为收益。局部近似才用K_H与trace(gK_H^{-1})，全局需remotebrancherror。用同一channeluses、平均power、decoderlatency比较readout改变的成功概率，优于把原rmax当成“几何容量”。

### 数据与成本

PDF6–9、16–20：VGG19、CLIP、DINOv2 features，BERT/GPT2/OPT；64 ImageNet图、64 C4长度≤256文本，不是大规模通信benchmark。10000AdamW steps/image/sequence，单图10noiselevels约20min V100，不能直接当接收器时延方案。比较noisy中间/末层feature很有启发，可为“intrinsiccompression vs extrinsicredundancy”设negativecontrol；推荐冻结backbone、用amortizedreadout。ImageNet分类应使用有label的官方validation，不能把无标签test当分类test。

## 9. Riemannian Deep Learning（424 页博士专著）

### 适用范围

专著主要研究**已知 manifold-valued SPD/correlation/hyperbolic/rotation data 上的networkoperators**，不是从图片学未知低维semantic surface。PDF48 Def.19、PDF54 Def.33说明SPDmetric与pullback；pullback正定需immersion，diffeomorphism产生isometry。可作基础引用与实现参考，但熟悉已知SPD空间不等于证明本稿source是紧边界less2Dsurface。

- PDF53 Thm.31 Cartan–Hadamard要求complete、connected、simplyconnected、nonpositivecurvature，才能Exp全局diffeomorphism。本稿sphere/一般compactmanifold不满足，receiver不可据此假设任意两点唯一globalLog。Fréchetmean唯一也有convexball等条件。
- PDF218–219 Lemma149 /Eq.(6.6)–(6.9)：选择diffeomorphismϕ:SPD→symmetricmatrices，pullbackflatmetric后distance=`||ϕ(P)−ϕ(Q)||F`，Exp/Log/PT有closedform。**可直接用标准统一base matrixlog或log-Cholesky**压缩协方差/动作信号，作为knowngeometrybaseline；不能据此任意flatten固定曲率taskmetric并宣称RDF降。
- PDF237–242 Thm.169/170/173/174/179：CholeskyoffdiagEuclidean、diagpower/Bures-Wasserstein的productmetric及classifier/residualoperator，可做covariance/skeletonmetric模块。poweroperation需要diagonalpositive的domain；gyroform出现 `(L^β+K^β−I)^{1/β}`，任意两个L,K不保证实正，clipping会改变原数学性质。

### ALEM 的局部适用疑点

PDF220 Prop.150 /PDF221 Thm.152把按eigenvalue分别选择不同logbase的map视为全SPDdiffeomorphism。PDF317 A.4.4.1只处理eigenpair排序，并假设eigenvectorunique，没有解决重复eigenvalue。若S=2I、a1≠a2，则 `U diag(log_a1 2,log_a2 2) U^T` 随eigespaces选择变化，因此**不同base版不宜直接宣称全SPDwell-defined/smooth**。排序不能消除此问题；统一scalarbase、standardLEM或log-Cholesky避免这个边界。记录此点是为了选稳妥组件，不是对全书作否定。

### 数据

PDF273–275：HDM05（2273seq31joints）、FPHA（1175seq45classes）、NTU60（56880seq60classes）、NTU120（114480seq120classes）、G3D；另有radar合成、EEG、CIFAR/ImageNet、AFEW、graphs等。TCOM优先**NTU skeleton模态**或CMUmotion，有明确动作task、源结构和信道传输意义；SPDcovariance是可解释knownmanifold对照。其他graph/EEG不必把整套加入主线。

代码参考：[ALEM](https://github.com/GitZH-Chen/ALEM)。

## 10. 建议优先采用的三个组件，以及真正需要证明的部分

### 第一优先：固定任务 g 的信道几何与有限功率风险

引用point-cloud/Fisherpullback与readoutfiber思想，定义实际高斯链路 `Y=A_HF(U)+N`，`K_H=JF^T A_H^T Σ^{-1}A_H JF`，固定任务g衡量失真。**待证候选：**在localidentifiability、regularity、highSNR与合适effective/unbiased/localasymptoticallyefficientdecoder条件下，局部风险首项为 `E tr(g K_H^{-1})`（若Σ=σ²Σ0，则写 `σ² E tr(g K0_H^{-1})`），据此分配power/anisotropy。不能声称任何biasedMAP有限噪声都满足pointwiseCRLB，更不能只从Fembedding推出A_HFembedding。

必须同时证明/测量：local rank、全局可分性、tube/reach内事件、低SNR远处分支混淆与fadingoutage。K奇异时inverse不存在；随机衰落接近0时E K^{-1}可能发散，应条件化nonoutage或使用完整boundedrisk，不能在积分时忽略这些channelstates。这个方向和TCOM信道/power对象最一致。

### 第二优先：有限码本的几何量化界与实际编码约束

引用RNOT Thm.3.1的M^{-2/m}量化尺度，在我们固定 `(M,g,µ)` 下证明对应有限M重构点converse/achievability，并比较**同M、同实际bits、同power**的受限codebook/quantizer。densitycompanding可作为knownBA/OTcomponent。若allowedcodebook任意且F是bijection，重参数化不会改善最优任务distortion；因此必须明确限制实际codebook/search/networkfamily、发射维数、量化器或decodecompute，才有表示设计收益。神经连续map不能取消digitalquantization。

### 第三优先：隐式支撑约束的有限步接收器及误差分解

引用MetricMatching低秩tangentestimate和ImplicitGeodesics投影/augmentedLagrangian，构建固定g上的有限Kreceiver。待证项为“noise/localstatisticalrisk + geometry/projectionestimateerror + optimization/discretizationerror + globalbranch/outageerror”；不要直接引用stationaryLangevin或exactgeodesicconvergence跳过finiteK。选择sphere/torus/SO(3)真值→CMU/NTUskeleton的两级实验，可辨别收益究竟来自constraint、densityprior还是metric/preconditioner。

## 11. 推荐数据的官方访问/许可核验（2026-10-04）

以下只核验作者/机构官方页面；未下载数据，也未声称所有binary链接都完成下载校验。许可说明是所核验页面的使用条件摘要。

| 数据 | 当前入口与条件 | 我们的用途与优先级 |
|---|---|---|
| CIFAR10/100 | [Toronto官方](https://www.cs.toronto.edu/~kriz/cifar.html)现重定向作者cave.cs.toronto.edu，列官方downloads和标准50k/10ksplit；该页要求引用technicalreport，未给统一标准license标签，不能替它填MIT/CC0 | 最快现实benchmark：冻结分类q/512或768dimfeatures，有限bits/CBR/SNR+语义task。图像集合不自动满足单一光滑manifold假设 |
| ModelNet10/40 | [Princeton官方下载](https://modelnet.cs.princeton.edu/download.html)：包含标准train/test，OFF；原CAD版权归originalauthors，数据提供仅方便academicresearch，标签freelyprovided | 推荐第二现实模态；2048points固定采样策略，classifier冻结，保持测量噪声与通信噪声分开 |
| ShapeNetCore v2 | [官方ShapeNet组织](https://huggingface.co/ShapeNet)：**旧shapenet.org/about已转到官方HF分发**，各dataset需申请access；[Terms](https://huggingface.co/spaces/ShapeNet/README/blob/main/TERMS.md)限noncommercialresearch/education | 后续大规模pointcloudpretrain。不要按旧镜像说立即免申请可下载 |
| CMU MoCap | [CMU官方](https://mocap.cs.cmu.edu/)声明数据freeforalluses，允许在products使用但禁止直接resell原/converteddata，请遵守官方acknowledgment、不要crawl | 推荐knownkinematicgeometry实验；sequence/subjectsplit，有动作/关节真值与SO(3)基线 |
| NTU RGB+D60/120 | [ROSELab官方](https://rose1.ntu.edu.sg/dataset/actionRecognition/)：60/120classes、25joints；skeleton60为5.8GB，120再+4.5GB；register/request/ReleaseAgreement，noncommercialacademic，未经许可禁止redistribution/derivingnewdataset；[作者代码](https://github.com/shahroudy/NTURGB-D) | 推荐通信动作识别主benchmark，先只skeleton避免TB级RGB/depth；followofficialcrosssubject/crossview(orcrosssetup)protocol |
| PlantVillage | [作者官方仓库](https://github.com/spMohanty/PlantVillage-Dataset)链接[作者HFcard](https://huggingface.co/datasets/mohanty/PlantVillage)：card为CC-BY-SA-3.0；54,306原图，color/grayscale/segmented是不同config，80/20保持leaf_idgroup | 农业边缘task备选；group同一叶片避免泄漏，健康/疾病posteriortask清晰，domainshift另测 |
| AFHQ | [作者StarGANv2仓库](https://github.com/clovaai/stargan-v2#animal-faces-hq-dataset-afhq)：提供downloadscript和pretrained，code/models/datasetCC-BY-NC4.0；v2resize/样本和旧版不同 | 感知geometryqualitative/小规模task。三类较简单，不能单靠动物换脸图支撑TCOM主结果 |
| FFHQ | [NVIDIA官方](https://github.com/NVlabs/ffhq-dataset)：70k图，128thumb约1.95GB；datasetmetadata/scriptCC-BY-NC-SA4.0，individualimages各自license在metadata；官方不用于facialrecognitiondevelopment | 可作非身份的perceptualreconstruction诊断，需明确attributestask标签来源；不推荐为本项目首个真实benchmark |
| HDM05 | [MPI官方入口](https://resources.mpi-inf.mpg.de/HDM05/)在本次工具访问403；未核验当前download或license | 可选motiondataset，但报告访问未核验，不能承诺立即可复现 |
| 自生成 sphere/torus/SO(3) / cow rotations | 三种intrinsicgeometry可解析，cowgenerator使用资产许可须另核对；参考ImplicitGeodesics PDF7–8/17 | **首选理论验证**：knownmetric、knownrotation/task、可控density/curvature/reach与noise。不能把generator输出的30k图说成已取得第三方publicdataset |

## 12. 必须保留的通信实验对照

同一冻结taskmetric/source/taskposterior下至少：baselinegeometry+fixeddensity；learngeometry+fixeddensity；baselinegeometry+learndensity；joint；固定g的BA/任务quantizer或标准digital/JSCCbaseline。每个报告实际bits或channeluses、平均power、SNR/fading、decoderiterations/latency，语义posteriorloss/accuracy与d_g²至少两项。geometry-only增益必须独立存在并跨seeds/资源公平比较，才能定位几何贡献。冻结robust/DINO/CLIP感知metric改善只用于proxy/ablation，不自动升格为固定tasksemanticrisk保证。

现有Gaussian-volumekernel的curvature/density高分辨率系数是**所给kernel的上界或局部展开**，不应靠这些论文包装成“对任意信道/任意最优codec的普遍曲率容量定理”。如果新的有限power/finitecodebook/finitecompute结果尚未证明，正文应写成建模与待检验猜想，不能把标准几何结果与经验数据拼接成已证明的TCOMnoveltheorem。
