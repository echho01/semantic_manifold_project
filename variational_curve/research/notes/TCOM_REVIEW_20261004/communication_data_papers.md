# 通信、压缩、数据集与拓扑论文审阅（TCOM 曲面几何导向）

审阅日期：2026-10-04。本文覆盖通信、压缩、数据与拓扑方向的 13 篇 PDF，并对照当前 `variational_geometry.tex` 与 `sections/communication_geometry.tex`。页码为本地 PDF 从 1 开始的页码；不同 arXiv 修订版可能不同。结论区分“已有方法，可引用后采用”“可适配的研究工具”“需要我们另行证明的新候选命题”，不能把引用方法换名后列为原创。

阅读范围：短论文完整阅读正文；长论文遍历完整章节/附录标题，重点阅读模型、主定理、证明条件与实验。Lossy Compression（66 页）重点为 pp.1–25、26–42 的例子及相应条件，未逐行审计全部附录证明；Pragmatic Information（152 页）重点为 pp.1–4、25–27、47–49、64–69、115–118、135–136，其余章节按目录与定理结构梳理，未逐页通读。SEMASIA、Token Embeddings、Causal Abstraction 的实验/核心定理及相关附录重点读取，未逐行审计全部证明。本文不是所有定理正确性的形式验证。

## 1. 对当前项目最直接的结论

最有用的组合是：Riegler 等的广义 Shannon 下界/量化尺度作为**理论引用基础**；RaDOGAGA 与 Low Bending Distortion（LBD）作为**几何编码的强基线与正则工具**；SEMASIA 作为**异构端模型、同源样本配对的实证资源**；Pareto Analysis 作为**通信成本—任务收益评价方式**。四篇拓扑论文主要界定相邻工作与基线，不应把本文主线从曲面几何转成持久同调。

当前项目已经有固定任务几何、表示不变性、reference SLB、Gibbs 前向核、输出密度校正及 BA 特例。不能再把“低维表示降低率”“曲率进小失真项”“几何 Gibbs 核”“任务权重度量”本身写成未有过的新理论。现有 circle 反例中，geometry-only 从 1.26498538 到 1.26343962（约 0.12%），density-only 为 1.18319194，joint 为 1.18319195；主要增益来自 reproduction density，而非几何。必须保留 density-only 与非几何同容量基线。

真正可以新增的对象应是受限可实现通信器：有限码本/图册、真实熵编码和头部开销、有限符号数与功率下的噪声信道、未知/异构接收知识，以及几何正则如何改变这些受限类的误差或资源分配。应把理论写成“在该可实现类与固定任务度量下的条件性收益”，不能写成任意可逆变换降低 unrestricted RDF。

## 2. 逐篇证据与可借鉴项

### 2.1 Lossy Compression of General Random Variables — Riegler, Koliander, Bölcskei

本地 66 页，arXiv:2111.12312v4；Information and Inference 接受稿。全结构为 regularity（§2）、RDF（§3）、quantization（§4）、二者联系（§5）、紧流形/球面/Grassmannian（§6）、self-similar sets（§7）。

- **可直接引用的建模/理论**：任意源与失真空间上选择参考测度 `mu`，定义 generalized entropy；Prop.3.1（p.15，(3.1)–(3.2)）给 `R(D) >= h_mu(X) - inf_s [sD + log sup_y int exp(-s rho(x,y)) dmu(x)]`。这正是当前 reference SLB 的基础形式，不能重新宣称原创。
- **主定理条件**：Thm.3.1（pp.16–17）要求源相对 σ-finite 参考测度绝对连续、entropy 有限；小球满足 `mu(B_{rho^(1/k)}(y,delta)) <= c delta^m`。有限半径版本还需要文中归一化/截断条件。无截断显式项为 `h_mu + (m/k)log(m/(kD)) - log(c Gamma(1+m/k)) - m/k`；并非任意云团均自动具备该小球维数。
- **可用的有限码本尺度**：Thm.4.1（p.18）以源 density 的 `L^{p/(p-1)}` 控制量化下界；Thm.4.2（p.20）以 reproduction measure 的 superregularity 与有限直径条件给 `n^{-k/m}` 上界；Cor.4.4（pp.21–22）在上下正则匹配等附加条件下给量化维数。平方距离取 k=2，可用于码本数—误差的机制验证。`log2 n` 是固定长度码本索引长度，不能直接代替真实 `H(K)` 加协议头部。
- **不能搬用**：球面章节部分采用 ambient chordal distortion，而当前项目是 intrinsic squared geodesic distortion；常数、球冠积分不可混用。m 为相应参考测度的小球维数，不能拿全球 PCA rank 当 m。不是噪声 JSCC、短包或未知接收知识的定理。
- **数据**：解析 hyperspheres、Grassmannians、Cantor/self-similar sets，无下载型真实语义数据。适合理论机制/负对照，不适合单独支撑 TCOM 的应用有效性。
- **TCOM 判断**：引用价值高；若主贡献只是该下界在流形上的改写，重合风险极高。新增有限码本、actual bit cost、task metric 与 channel mismatch 的组合需要新证明。

### 2.2 A Mathematical Theory of Pragmatic Information — Niu, Zhang

本地 152 页，arXiv:2609.10986v1，2026-09-10 预印本。章节覆盖语义映射、pragmatic entropy/MI、capacity/RD、VoI/CoI、三类 coding theorem、连续信源及联合优化。

- **值得借用的对象**：Def.2（pp.25–26，(18)–(20)）按“相同 Bayes 最优行动”构建任务等价类；Def.3 加未来价值形成动态映射。可用作固定任务下直接行动传输的 baseline，以及说明“可复用状态表示”和“一次任务最优行动”的区别。
- **可用经典论证**：Thm.19（p.65）后验 Bayes 决策的 VoI 非负，是条件优化与 Jensen 类的经典事实；没有唯一行动/无 ties 时，其等号“行动等于基准”条件须改成基准行动仍最优。Thm.21（p.66）时间共享支持 concave value envelope，不能据此断言神经编码器的有限训练可行集凸。
- **明确不能直接使用的部分**：p.64–65 以“同最优行动”为由，称类内任意代表元的全部 utility 相同；这一般不成立。例如最优行动相同而奖励幅度、次优行动 regret 不同。应改用类内条件平均 utility 或更强的奖励向量等价关系。Def.12（p.47）把 down-MI 作正部 `max(0, ...)`，正部不能保证严格凸；Thm.7（p.48）up-MI 加 coarse-graining entropy 后大于普通 MI，不能当作物理信道可传 bits 的增加。
- **连续 entropy 风险**：p.115 的 (245) 写 `h(W)-E log|fiber|`。一般 `h(W)=H(A)+E h(W|A)`；该表达式只在每类 conditional uniform 等条件下等于行动 entropy，否则通常小于 `H(A)`，不能当一般定理照搬。p.67 Thm.22 显示的 dual expression 无显式 R 依赖，不能照搬成所需 frontier duality。p.135 Thm.47 的一阶 slope 平衡若宣称充分性，需要 concavity/differentiability，而非仅 bounded utility。
- **数据**：交通灯、离散/Gaussian/POMDP 例子；未识别出适合本项目的公开真实实验集。
- **TCOM 判断**：可引用问题动机和任务 quotient，不宜用其 up/down capacity 替换普通 channel capacity。当前 repo 的固定 posterior geometry + ordinary MI + physical channel capacity 更稳妥。以这些新名词“包装”主理论风险高。

### 2.3 Value of Communication in Goal-Oriented Semantic Communications: A Pareto Analysis — Luo, Li, Pappas

本地 16 页，首页标注 IEEE TIT DOI `10.1109/TIT.2026.3720415`。通信为有限不可约 Markov source、send/skip、packet drop、远端 MAP 估计。长期 estimation cost J 与 attempt frequency F 形成 Pareto frontier。

- **建模可借鉴**：Lemma 1（p.5，(12)）的足够状态为 `(X_t,Z_{t-1},Theta_{t-1},Delta_{t-1})`：包含当前状态、上次成功源状态、age 与 consecutive error duration。说明通信价值不仅看 instantaneous mismatch，也看上下文与错误持续时间。
- **主定理条件**：Assumption 1（p.6）允许 countable state，要求有限动作、非负可无界 c、bounded resource f、至少一个 stationary policy 有不可约正再生链和有限 J、成本 sublevel 集有限等。Lemma 2 的 urgency 增长还受 transition/drop 参数约束。Thm.1（p.6）给直到足够 budget 的严格递减 convex frontier；Thm.2（p.7）给 deterministic corner policies 与分段线性/time-sharing。这些结论不自动适用于静态 feature codec 或非凸 NN 搜索。
- **算法可借鉴**：SPLIT（pp.10–12，Alg.1–2）用 Lagrange slope 求相邻 corner，Thm.3 的复杂度由 corner 数与 RVI cost 控制。mixing policy 的 coin probability 与长期时间占比不同，(38) 含 recurrent cycle length；不能把每帧随机 θ 直接当 convex coefficient。
- **数据**：p.13 synthetic 8-state normal/alarm 系统，drop 0.3，指定 transition matrix 与 exponential urgency，随机法 simulation 达 `10^7` steps，无真实数据。
- **本项目用法**：主实验报告实际 bits/use/power—固定任务 loss Pareto 曲线、达到相同 task target 的最小成本与边际收益。若发展“时序曲面 refinement/send-skip”，必须重建 Markov sufficient state 与正再生条件，是另一篇更大工作。静态论文中只借评价，不搬 CMDP frontier 定理。

### 2.4 Rate-Distortion Optimization Guided Autoencoder for Isometric Embedding in Euclidean Latent Space（RaDOGAGA）— Kato 等，ICML 2020

本地 11 页，正文完整读取。核心是 rate loss、重建 loss 与注入 latent noise 后的 decoder output loss。

- **可借模型**：p.3–4，(7) 的形式为 `-log P_z(z)+lambda1 h(D(x,xhat))+lambda2 D(xhat,xbreve)`；噪声为等方差独立 uniform latent perturbation，decoder noise sensitivity 配合 entropy model 形成 scaled orthogonal 的局部几何。这是本项目最应加入的 rate-aware isometry baseline。
- **理论条件**：p.5，(13) 给 `J^T A(x) J = [2 lambda2 sigma^2]^{-1} I`，依赖小噪声 quadratic Taylor、local distortion SPD A、decoder full rank、良好重建 `xhat≈x`、density model 接近真实 density、理想优化平衡等。形式推导先作同维 M=N，再讨论 near-manifold 高维数据。它不是任意训练结果的 global isometry theorem。
- **边界**：曲率非零的一般流形不能全局同维等距嵌入 Euclidean latent；局部近等距不等于全局 geodesic 直线、可逆性或 noise robustness。pp.5 的 density Jacobian 公式 (15)–(17) 有 codimension/volume 条件；奇异支撑不能直接当 ambient Lebesgue pdf。正文提到 Appendix B，但本地 11 页未包含该附录，不据此声称已验证全部推广条件。
- **所有实验数据**：CelebA（latent 256，BCE/SSIM 与 beta-VAE，p.6）；3-D Gaussian mixture 合成数据映至 16-D（pp.7–8）；KDDCUP99（494021×121）、Thyroid（3772×6）、Arrhythmia（452×274）、KDDCUP99-Rev（121597×121）异常检测（pp.8–9）。正常样本训练 50%，给定 anomaly 比例作 threshold，20 seeds。正文链接 UCI/ODDS。
- **TCOM 用法**：把统一扰动换成指定实际 channel 并加入固定任务 metric 是适配，但“RD optimization 推出 isometry”本身已有。异常数据偏离当前视觉曲面主线；CelebA 可作多属性任务，而需 direct posterior/label transmission baseline。

### 2.5 Convergent autoencoder approximation of low bending and low distortion manifold embeddings（LBD）— Braunsmann 等

本地 30 页，arXiv:2208.10193v2。重点读取模型、Thm.1/2 与证明、实验/限制。

- **最可直接采用的正则**：pp.6–7，(2.1)，近邻 `(x,y)` 的 geodesic distance 与 midpoint；一阶 stretch 项 `gamma(||F(y)-F(x)||/d_M(x,y))` 与二阶项 `lambda ||8[(F(x)+F(y))/2-F(mid_M(x,y))]/d_M(x,y)^2||^2`。前者最小于伸缩比 1，后者在 local scale 下近似 covariant Hessian directional energy。这与当前固定 g 的 low-bending 正则十分接近，必须明确引用。
- **主定理条件**：紧流形可含 boundary，但有 cone condition；小邻域唯一最短 geodesic 和正 injectivity radius；pair density D1–D3（p.8）归一化、按 epsilon^m 重标后双侧正界、存在极限 density。网络类 H1–H3（p.12）要求零均值 L2 closed、`C^{2,1}` 与随 epsilon 的 Lipschitz derivative bound、H2 density。Thm.1（p.12）在 H1–H2 下存在 minimizer 和 uniform H2 bound；Thm.2（p.14）在 H2–H3 且 `C_L=0` 下 Mosco convergence 到 (3.9)；Cor.1（p.16）是 global minimizers 的 subsequence weak H2 convergence。并非固定有限网络的 Adam convergence/global optimum，且 LeakyReLU 需要平滑化才满足所述正则性。
- **几何解释边界**：该 Hessian penalty 包含坐标参数化的 in-plane bending，不能直接称为 pure second fundamental form 或 Willmore energy；也不自动保证全局 injectivity/reach。已知 geodesic midpoint 是前提，真实 SEMASIA features 用 kNN 路径 midpoint 时需要离散近似误差分析，不能沿用原 Mosco theorem 无修改。
- **所有实验数据（pp.17–18）**：sundial hemisphere `S^2`（64² 灰度）；rotating cow `SO(3)`（128² PyTorch3D）；rotating colored arc on Klein bottle（128² Blender）；translation+orientation ellipse `S^1×[-1,1]^2`（64²）。这些生成器能给 ground-truth distance/midpoint，优于靠图像像素猜流形。p.25 指出像素化 binary images 形成 finite set 与理想 smooth manifold 的区别。
- **公开情况**：作者代码 [GitLab LBD](https://gitlab.com/jubrau/LBD) 可打开但本次网页未展现许可；论文给数据 [DOI 10.17879/48978454901](https://doi.org/10.17879/48978454901)，本次 DOI 页面获取失败，下载可达性/具体许可未核实。不能写“已下载/许可已确认”。
- **TCOM 潜力**：用它作为几何工具，另证明/实证“有限 quantizer 或 decoder 面对 channel noise 时”的收益。仅 AE 加 stretch+bending 后换名，原创性不足。

### 2.6 Autoencoder-Based Domain Learning — Wheeler, Natarajan

本地 6 页完整正文。把 Gärdenfors conceptual domain 建成 Euclidean quality dimensions，类别 property 用 convex prototype/Voronoi 区域。

- **可借模型**：pp.3–4，encoder+decoder+无额外 trainable classifier head 的 prototype distance→similarity→softmax。loss (6) 为 weighted reconstruction MSE、CE 与防 distance blow-up 正则；Alg.1 更新 class mean prototype EMA，首 epoch 关闭 CE。可作 supervised/prototype chart 初始化、二维示意与便宜 baseline。
- **理论边界**：没有全局等距/embedding/injectivity、通信或 calibration 定理。二维 prototype 能区分已选标签，不证明图像有 2-D intrinsic semantic manifold。Euclidean convex classes 不是任意语义任务的普遍结构。
- **所有实验数据**：MNIST 70k（latent 2）；CelebA 仅 male/smiling 两个 binary 属性（112² RGB，latent 2，p.5）。若复用应扩展任务并明示留出的属性/identity 不保证。
- **通信比较风险**：p.6 把 2 float32=64 bits 和 uncompressed pixels 301056 bits 比较得到 99% reduction；这是表示大小示例，缺少真实 entropy code、head、noise 与同任务编码基线。仅两个 binary 属性的已知任务可直接 2 bits，因此不能沿用该“通信增益”作为我们的实验。
- **TCOM 判断**：好用的任务原型 baseline；单独作为核心贡献弱，且 prototype semantic AE 已有。

### 2.7 SEMASIA: A Large-Scale Dataset of Semantically Structured Latent Representations — Pandolfo 等

本地 39 页，arXiv:2605.09485v1。主文 pp.1–10 与 alignment Appendix E（pp.30–33）重点完整读取。

- **可直接复用的数据结构**：paper p.4/Table 1 报 1697 pretrained timm vision encoders；head 前 embeddings、architecture/pretraining/resolution/dimension 等 metadata；不同 model 对同一 source image 形成配对。优先选 CIFAR-100 的 2–4 个模型，而非全量约 1700 模型下载。
- **必须借作基线的模型**：train-only centering + regularized Cholesky whitening（epsilon 1e-6）；linear LS `A=Z_B Z_A^dagger` 后 rank-k SVD 截断；whitened CCA；prototype/PPFE。p.7–8 显示低维 linear 经常胜过 PPFE/CCA；geometry方法必须先胜过这个强基线，而非只比无 alignment。
- **已覆盖的邻近贡献**：压缩维度 k、异构 latent alignment、MSE 与 downstream accuracy、metadata regression、TDA/prototype 是原论文已有。新增点需 actual bits/noisy channel、task-calibrated fixed geometry 与 mismatch 分解；不能把“模型间语义流形对齐通信”本身列新颖。
- **数学边界**：Appendix E p.32 的 PPFE 对 `k<d` 的 anchor matrix 写 full inverse square root 与 `F^T F=I`，有 rank 障碍；此时要 pseudoinverse 与 support projector，不能无条件 copy。global covariance effective rank 不等于 intrinsic tangent dimension；t-SNE 图不能证明 manifold 假设。
- **所有源数据/论文 reported splits（p.4；非自动保证当前 HF 相同）**：CIFAR-10、CIFAR-100 各 50k train/10k test；MNIST、Fashion-MNIST 各 60k/10k；CelebA 100k train/19867 val/19962 test；Oxford Flowers 7169 train/1020 test；Tiny-ImageNet 100k train/10k val；ImageNet-1k 50k val/100k test。Flowers 和 ImageNet split/labels 与常见原版设置不完全一致，必须检查实际 registry/label 可用性，不能把 test 预测当已知 GT。
- **官方验证**：[作者 GitHub](https://github.com/SPAICOM/semasia-datasets) 公共且代码 MIT；[官方 HF CIFAR-10 card](https://huggingface.co/datasets/spaicom-lab/semasia-cifar10) 公共，1699 当前 configs，Parquet 可按 model config/split 选取。card 可用 `load_dataset(repo, model_config, split=...)`，无需下载全部；本次 CIFAR-100 card 获取失败，**CIFAR-100 的下载条目与许可仍需在实际准备时核实**。没有执行数据下载。
- **配对 ID 必查**：paper 称同一 source id across models；当前 CIFAR-10 card 将 `id` 描述为 shard row index，unique per split/model，非跨模型唯一，同时声明 row order follows original split。必须检查 source order/label 或显式原图 index 后构成 join key；不能只凭相等 id join。所有 train/validation/test 先按原始图像身份划分再配对，避免同图不同 encoder 进入不同 split。
- **许可边界**：GitHub MIT 是代码许可，不能延伸为所有源图像/derived representations 的许可；HF card 本次未找到 license 字段。需保留原 dataset 的协议。官方 card “全集 size”不是单个 model subset 的下载成本。
- **TCOM 判断**：主实证首选。使用 frozen encoder/frozen posterior、模型对分层（相近预训练 vs 明显异构）、classification+retrieval/multi-task，task budget 对齐。维度 k 不能作 bit rate 的替身。

### 2.8 Token Embeddings Violate the Manifold Hypothesis — Robinson, Dey, Chiang，NeurIPS 2025

本地 37 页。主文、Algorithm 1、Thm.1/2 与 limitations 重点阅读。

- **可借方法**：pp.5–6 的 multiscale local neighbor-volume/log-radius slopes；Algorithm 1 对 sorted pairwise distance 邻居计数，比较相邻 windows 的 slope distribution，用 Holm–Bonferroni 修正多检验。适合检验本项目 smooth/low-curvature 假设是否可信，而不是默认所有 embedding 都是 smooth surface。
- **定理条件**：Thm.1（p.5）针对 smooth manifold/fiber bundle、已知 positive embedded reach 与 `r<reach` 的几何 neighborhood-volume；边界/纤维宽度引入尺度拐点。Thm.2（p.5）是 generic transformer/context output embedding 的维数条件（如 `m>2wd/ell`, `w>=m`），不能迁到视觉曲面/信道。
- **关键边界**：拒绝的是所设尺度内低曲率/高 reach 类，不是任意 finite set 的全部 manifold 模型。fail to reject 不证明流形。非均匀 density 时 neighbor counts 是概率质量，不是 volume；需 density control、多尺度、local PCA residual 与 bootstrap/held-out 校验，不应机械照搬 rejection count。
- **数据**：GPT-2、Llama 7B、Mistral 7B、Pythia 6.9B 的 token vocab embeddings；p.8 table rejection counts 为 66/33/40/54。不能由它推断 SEMASIA 图像 features 同样非流形。
- **代码/用途**：论文 p.33 指 [stratified_estimator](https://github.com/kb1dds/stratified_estimator)。建议做 assumption stress test/negative control，文本 vocab 不作当前 TCOM 主数据。

### 2.9 Causal Abstraction Learning based on the Semantic Embedding Principle — D’Acunto 等，ICML 2025

本地 37 页。主文、Thm.4.3/Prop.5.1、实验和相关附录重点阅读。

- **可借建模**：Def.4.1（p.4）SEP 要 `alpha beta = Id`，即高层可在低层稳定嵌入后取回；线性正交子类 `V∈St(ell,h)`、`V^T V=I` 是方便充分结构。对零均值 Gaussian 低/高层，(3)（p.5）用 `KL(N(0,Sigma_h)||N(0,V^T Sigma_l V))` 匹配 marginal covariance。
- **可引用条件**：Thm.4.3（p.5）SPD covariance 及 eigenvalue interlacing `lambda_i<=kappa_i<=lambda_{i+ell-h}` 给 spectral feasibility 条件；Prop.5.1（p.6）给 Stiefel 优化的解析梯度。可作为 orthogonal/linear latent adapter baseline，与 SEMASIA LS/CCA 同列。
- **不可搬成新结论**：一般 right inverse 不蕴含 Stiefel；正交是选择的子类，不是 SEP 的一般等价。匹配 observational Gaussian covariance 不识别 causal mapping，更不能保证任务 posterior、interventions、global geometry。paper p.4 明示 observational solution 可非唯一。p.6 MADMM convergence 未有保证；p.7 ManPG 是相应条件下 stationarity，不是 global truth recovery。
- **所有实验数据**：synthetic known abstractions，full prior `(12,2),(12,4),(12,6)`，30 GT×50 initializations；partial prior `(4,2)`，30×30。HCP resting-state fMRI 100 healthy adults/89 ROIs；展示一个 adult 左半球44+vermis共45→14 macro ROIs；high-level target 由已知 aggregation 生成，并非独立采集 causal ground truth（p.8）。
- **代码/TCOM 判断**：[作者 calsep](https://github.com/SPAICOM/calsep)。可借 stable linear interface，不建议把当前视觉 alignment 称“causal semantic communication”，除非增加真正 intervention/causal data 与评价。HCP 医疗数据协议与应用门槛高，非本项目优先。

### 2.10 TopoCode — Guo

本地 5 页，arXiv:2411.12825v3。作者 [正式出版清单](https://cse.unl.edu/~hguo/publications/)列 IEEE Communications Letters 2025。

- **可借模型**：发送原数据同时发送 persistence diagram（PD）作为 side signature；receiver 优化 PD Wasserstein mismatch 加 binary prior。gray image cubical filtration、H0/H1 metadata 与 birth/death pairs 构成额外包。p.3 的 (4) 可作为结构先验辅助纠错 baseline。
- **理论边界**：p.3 的 (2) 是给定 filtration 与范数条件下 persistence stability，上界从 data perturbation 到 diagram perturbation；没有反向保证，PD 非 injective，很多不同图像共 diagram。signature 正确也不能恢复所有 non-topological bit errors，paper p.5 自述高 SNR 这类局限。PD 梯度 only almost everywhere；非凸纠错目标没有全局恢复定理。
- **所有实验**：200 MNIST/Omniglot 图像，BPSK、AWGN 3 dB；LDPC 1000 bit codewords、degrees 5/20，conv rate 2/3（p.4）。示例 15 BD pairs=256 bits 相对 28²×8=6272 bits 的 4% overhead，要同时计 metadata、signature 的保护成本与同等 total channel uses。
- **公开代码**：[作者 TopoCode 项目页](https://cse.unl.edu/~hguo/research/topocode/)可见 TopoCode-PCD 新工作，未核实原 image TopoCode 代码；不能把 PCD 仓库当原 image 论文源码。
- **TCOM 判断**：说明“几何结构辅助纠错”已有 letter。适合作相邻工作/receiver 辅助 baseline；不能把 addition of PH signature 宣称新贡献，也不应据 PD 小误差声称 task loss 必小。

### 2.11 Topological Autoencoders — Moor, Horn, Rieck, Borgwardt，ICML 2020

本地 11 页，正文完整读取。

- **可借 loss**：pp.3–4，(1) 为 reconstruction + lambda topology；(2)/(3) 双向 loss，分别以 input/latent 的 PH pairing indices 比较两个 distance matrix 的被选 edges。实验主要 H0/MST；比 full H1 computation 便宜，适合作低成本拓扑正则 baseline。
- **条件与边界**：distance differentiability 与 generic unique distances 时 pairing 在局部固定，可微 almost everywhere；ties 需 symbolic perturbation。pairing-index loss 不等于 optimal diagram Wasserstein loss，零 loss 不保证全 PD 一样（paper p.4 明示），更不证明 global injectivity 或所有尺度 topology preserved。
- **理论可用范围**：Thm.1（p.4）把 subsampling diagram bottleneck error 用 Hausdorff-tail 控制；Thm.2（pp.4–5）expected Hausdorff 的近似依赖随机 distance/independence 条件。固定小 batch 不保证覆盖 rare components/holes；NN correlated point clouds 未必满足用于近似的独立性。
- **所有数据**：十个 small 100-spheres 加大 encompassing sphere，ambient 101-D 合成 Spheres；MNIST、Fashion-MNIST、CIFAR-10。latent 2-D visualization，以 reconstruction、density divergence、trustworthiness/continuity 评价，无实际bits/channel。p.8 指任务 classification 与保留全部 topology 可能冲突。
- **公开情况**：[作者 repo](https://github.com/BorgwardtLab/topological-autoencoders)公共归档（2022），BSD-3-Clause；版本依赖旧，采用时需 pin。dataset 许可不自动跟随代码 BSD。
- **TCOM 判断**：topological AE 不是新概念；它只能辅助曲面论文，不应以其替代 fixed task metric 及噪声可实现分析。

### 2.12 Topological Flow Matching — Wyrwal, Ceylan, Tong，ICLR 2026

本地 26 页。主文完整，证明附录按命题对应检查。

- **对象非常不同**：数据是定义在**已知有限 graph/simplicial complex 上的 node/edge signal**，不是 high-dimensional semantic points 位于未知流形。Hodge `L_k=B_k^T B_k+B_{k+1}B_{k+1}^T` 的 harmonic kernel 表示 signal 的全局结构。
- **可借建模**：p.5 的 (11) linear heat reference drift，典型 `-kappa L X`；spectral Gaussian prior `N(0,exp(-kappa L))`。zero eigenmodes 保留；bridge velocity 在各 spectral mode 显式计算。
- **理论条件**：Prop.2/3（pp.7–8）需要已知 finite PSD L、linear/Gaussian reference 的 transition/covariance；nonzero eigenmodes bridge 用 sinh 比值，zero mode 取直线 limit。Prop.4（p.8）把 exact cost 写成 covariance-weighted endpoint difference。一般 coupling 产生有效 TFM，只有相应 OT coupling 与所述 SBP 相关；σ=0 的 Schrödinger bridge 一般 ill-posed，paper 自述是 formal zero-noise limit，不能把 degenerate KL 当正常 theorem。
- **所有数据（pp.8–10）**：IRIS 1990–2017 earthquake（28 annual signals，576 nodes）；PeMSD4 traffic（17k signals，307 node graph lift 至340 edges）；HCP brain fMRI（约2k signals、360 nodes，1190 signals 的不同 frequency 版本）；embryoid-body differentiation single-cell（6k signals、18k vertices、PHATE embedding/5 time points）；ocean（20k edge samples，NOAA initial observations 之后 fit GP，target 合成 curl-free GP）。另有 CIFAR 图像试验，其与 ICFM FID 差异极小且有方差，不能单独强证。
- **公开情况**：[作者 repo](https://github.com/KacperWyrwal/topological-flow-matching)当前公共、MIT，包含 five signal experiments。原数据协议仍独立。
- **TCOM 判断**：若另做 remote traffic/flow signals，是合理 graph-signal receiver prior；当前视觉曲面主线无需引入 whole flow matching。不能写 graph Hodge topology 等于 latent manifold surface geometry。

### 2.13 Topology-Preserving Joint Source-Channel Coding — Erak, Alhussein, Fang, Muhaidat

本地 5 页，arXiv:2603.17126v1，2026-03-17。最接近本项目的通信风险论文。

- **已有通信模型**：CNN encoder 输出 complex z∈C^k，unit average power，`rho=k/(HWC)`，AWGN 或每图 block Rayleigh，decoder 读 real/imag（p.2）。本地 v1 没清楚交代 perfect CSI/no CSI，不能替它填入无 pilot 假设。
- **已有 loss**：p.3 (3) 为 MSE + image cubical PH(H0/H1) Wasserstein2 + batch latent pre/post channel 的 Rips PH；lambda image 1e-4、latent 1e-5，anneal。训练为 five conv PReLU DeepJSCC、Adam1e-4、batch128、随机 SNR 0/5/10/15/20，多个 rho（p.4）。推断不需实时 PH，是 training regularization。
- **理论边界**：没有新 RDF、reliability/finite-blocklength/收敛定理；batch latent topology 与真实样本 task topology 不相同。低 diagram error 不自动保 retrieval/classification/task connectivity。
- **所有数据**：Omniglot（1623 characters、50 alphabets）；DeepGlobe Road Extraction binary road masks；指标 PSNR 与 PH diagram Wasserstein，10 runs（pp.3–4）。没有直接 routing downstream test。
- **公开验证**：[官方 arXiv](https://arxiv.org/abs/2603.17126)可核实 v1；[作者项目页](https://omarsababha.github.io/NETAI/projects/semantic-task-oriented.html)把该工作列为 preprint，未找到该项 codelink。该页随后 Contrastive Privacy 项的 code link 不属于 Topology-Preserving JSCC。未确认公开官方源码，也未以二手网页作期刊接受证据；本次不能宣称 TCOM 已发表。
- **数据许可**：[官方 Omniglot](https://github.com/brendenlake/omniglot)公共、MIT，建议 alphabet holdout，作者原划分为 background30/evaluation20 alphabets。DeepGlobe [官方 resources](https://deepglobe.org/resources.html)检索到的条款为原 challenge only、后续开放许可计划；是历史页面，不能据此断言当前绝不可研究，但尚未确立当前开放研究许可，须查实际下载协议。第三方镜像 metadata 不能替代原协议。
- **TCOM 判断**：仅“manifold/topology preserving JSCC + PH regularization + noise robust”极容易和本论文重合。我们要明确不同对象：固定 task geometry、可实现 finitebit codec、density/geometry 分离、heterogeneous knowledge mismatch、可证明的几何局部噪声界。PH 损失可作 baseline，不能作为主新颖点。

## 3. 原论文所有实验数据汇总与本项目优先级

| 数据/资源 | 原论文 | 对本项目的具体用途 | 优先级/限制 |
|---|---|---|---|
| SEMASIA CIFAR-100、CIFAR-10 | SEMASIA | frozen feature、同图多 encoder、model mismatch；100 fine / 20 coarse label 可做任务尺度 | **首选主实证**；raw-image identity split，真实码率；当前 CIFAR-100 下载条目准备时核实 |
| SEMASIA CelebA | SEMASIA；RaDOGAGA/Domain 也用原 CelebA | 多属性/多任务 posterior、改变任务权重与知识兼容性 | 次选；人脸许可、任务/identity 区分，原 feature 的 train split 不是默认162770 |
| SEMASIA Flowers、Tiny-ImageNet、ImageNet1k | SEMASIA | 细粒度/层级任务和外部分布 sanity check | 次选；reported split 与 GT 标签必查，ImageNet资源大且协议单独 |
| MNIST、Fashion-MNIST、CIFAR-10 raw | Domain、TopoAE、SEMASIA | pipeline debug、TopoJSCC reproduction contrast | 不能作为全部主实证；MNIST过易 |
| Sundial hemisphere、rotating cow SO3、colored arc Klein bottle、ellipse S1×square | LBD | 真实因子、准确 geodesic/midpoint，局部弯曲与 chart 机制 | **强机制候选**；生成器/数据可达性许可尚未完全核实，像素有限化误差单列 |
| Hypersphere、Grassmannian、Cantor/selfsimilar | Lossy Compression | reference small-ball维数、有限 quantizer scaling 检验 | 理论机制，可解析，无真实 task |
| 3-D Gaussian mixture→16-D | RaDOGAGA | intrinsic/ambient维数、局部正交、density correction 负对照 | cheap synthetic，不能认作真实 smooth connected compact manifold |
| KDDCUP99/Rev、Thyroid、Arrhythmia | RaDOGAGA | anomaly score baseline | 任务不一致/低样本，不建议主实验 |
| Synthetic 8-state alarm Markov source | Pareto Analysis | 时序send/skip、urgency Pareto | 仅拓展方向，当前静态曲面论文无需并入 |
| GPT2/Llama/Mistral/Pythia token vocab | Token Embeddings | 高曲率/stratified negativecontrol | 只作假设检验，不作视觉曲面主数据 |
| Synthetic Gaussian SCM abstractions、HCP fMRI ROIs | Causal Abstraction | orthogonal adapter/known aggregation，causal对照 | 只借 baseline；无干预不写因果贡献 |
| Omniglot images+strokes | TopoCode、TopoJSCC | 原 topology baseline；可用 strokes 的结构任务 | 可公开获取；alphabet holdout，需要真正任务指标而非PD selfmetric |
| DeepGlobe Road masks | TopoJSCC | topology-sensitive road connectivity任务 | 当前许可未确认，优先寻找清楚协议再采用 |
| Nested Spheres ambient101 | TopoAE | global topology vs finitebatch覆盖负对照 | cheap synthetic，拓扑主线不替曲面 |
| PeMSD4 graph traffic、IRIS earthquake、NOAA/GP ocean、HCP fMRI、embryoid singlecell | TFM | graph-signal估计或结构先验 | 若改题才有价值；known graph与unknown manifold不能混用 |

本表只汇总这 13 篇的实验集。SmallNORB/dSprites/3dshapes 的官方下载与许可由综合评审独立核查；可作为“已知因子机制 + SEMASIA 主实证”的替代 LBD 机制集合，不需要全量使用上述清单。

## 4. 官方通信基线与新增最近邻

网页核验截至 2026-10-04，只确认公开状态和页面所示许可，没有运行/下载训练数据。各 dataset 许可不随 code license 改变。

| 基线/最近邻 | 官方来源及状态 | 应如何比较 |
|---|---|---|
| WITT（ICASSP2023） | [作者 repo](https://github.com/KeYang8/WITT)，公开，有训练代码和 pretrained model 链接；根页未显示 LICENSE，复用代码前核对许可 | AWGN/Rayleigh、同 symbol count/power/task decoder；参考 BPG+5G LDPC/Capacity |
| SwinJSCC（TCCN2024，DOI10.1109/TCCN.2024.3424842） | [官方 repo](https://github.com/semcomm/SwinJSCC)，公开、明确 official；有 rate/SNR adaptation 和 pretrained 链接，根页未显示 LICENSE | modern adaptive backbone；pin环境，README注明部分新版PyTorch inference不一致 |
| DeepJSCC-f（JSAIT，feedback） | [作者 repo](https://github.com/kurka/deepJSCC-feedback)，公开GPL-3.0 | feedback task 与当前无反馈不能无条件等预算；原2019 feedback-free DeepJSCC官方源码本次未核实 |
| Topology-Preserving JSCC | [arXiv](https://arxiv.org/abs/2603.17126)，官方代码未核实 | 同 backbone、MSE+PH两loss重现；全部 overhead/CSI与channel model声明 |
| Topological AE | [作者 repo](https://github.com/BorgwardtLab/topological-autoencoders)，BSD-3-Clause，已归档 | 便宜 H0/MST 正则 baseline，同 minibatch / latent dimension |
| Shared Semantic-Aware Codebook，TCOM2025 | [作者 SJTU repo](https://github.com/SJTU-mxtao/semantic-communication-w-codebook)，MIT，DOI10.1109/TCOMM.2024.3450877 | weighted data-semantic codebook、label codebook、variational loss、AWGN/fading；这是“共享语义码本”直接最近邻 |

最后一项是 paper 文件夹外此次发现的必要补充文献。官方 README 标 TCOM73(2):1217–1232，2025；已做 semantic-aware weighted distance codebook 与 label-based indexes，并在 WITT backbone 的 noisy/fading channel 比较 reconstruction/classification。因此“码本共享先验 + 语义距离 + channel训练”不能当本项目全新贡献。若我们采用 atlas/codebook，应说明固定任务 Riemannian geometry、曲面局部 chart、真实 entropy 与 model mismatch 如何超出该方法，且提供该基线或清楚的同预算改编。

## 5. 可真正研究的通信增量（尚未证明，不写成已有结果）

### A. 有限码本下曲面局部建模的误差—成本界（推荐理论核心候选）

固定 source law 和 task `g`；限定 codec 是有限 atlas/有限 reproduction points、局部线性/二阶 decoder 与熵编码。把总误差写成量化误差、离散/图距离误差、channel decoding error 和 decoder approximation，而不是用微分 entropy 的坐标变化当 rate gain。理论可尝试：在 lower reach、bounded Hessian、bi-Lipschitz constants、local chart radius 内，给局部线性误差的阶或 quantizer cell 条件；再与同尺寸 Euclidean/PCA codec 比较。LBD只供应正则估计工具，Riegler只供应 n^-2/m 尺度，二者均不自动证明新增界。

关键困难：bending penalty 的 L2 bound 不自动给 pointwise curvature/reach；若需 pointwise guarantee，必须加适当 Sobolev regularity/维数 embedding 或显式 pointwise bound。必须检测 injectivity 与 chart overlap；所有 chart labels、图/码本同步 overhead 计入比特。

### B. 固定任务几何与端知识不匹配的 finitebit 通信（推荐实证主线）

SEMASIA 同图 A/B 模型是 natural testbed；分别量化“source feature compression”“receiver model alignment”“noisy channel”所带误差。冻 posterior/q 与 reference geometry 防止训练把距离缩小。用 LS/CCA、Stiefel linear adapter、PPFE、RaDOGAGA/LBD 和 density-only 做 baseline。最小必要对照是“无需重构 feature 的直接任务 posterior/decision transmission”。若单任务标签传输已达到目标，曲面表示应以 reusable multitask/retrieval 或未知将来任务的额外价值证明必要性。

潜在理论：显式知道/估计 receiver error distribution 或 posterior mismatch，推到 task regret 的额外项与 actual rate-cost；但“marginal Gaussian covariance 匹配”不等于 posterior compatibility，需要单独条件。不可照搬 causal title。

### C. 跟任务与信道相关的局部资源分配（可作为算法增量）

在相同实际预算给重要且难传的 chart/patch 多 bit/power，可借 Pareto marginal slope 而非盲目给高曲率更多资源。当前局部 `sqrt(w_j c_j)` law 的 novelty 要和经典 high-rate/convex allocation 区分；新处若存在应在可估 task sensitivity、channel noise covariance 和 chart distortion 的联合设计、以及约束实际bitlength。PH重要性与curvature重要性不必等于taskimportance。

### D. 时序send/skip 或graph-signal拓展（本轮不建议主线）

Pareto CMDP、Hodge/TFM 都可另开方向。但加上时序Markov state/graphsignal会引入新的问题设定与数据，容易稀释曲面论文。当前优先把A/B实现和因果控制实验做扎实。

## 6. 建议的最小实验矩阵

1. **解析/已知因子机制**：uniform 与 nonuniform circles/spheres；已知因子的 SmallNORB/dSprites 或 LBD generator。固定 source 与 task geometry，分别观察 bending/curvature、density、finite codebook 和noise；明确boundary/quantizedpixels不满足closedsmooth假设的程度。
2. **主实证**：SEMASIA CIFAR-100 2–4模型对，sharedsource identity split；分类（fine/coarse至少两任务）+ retrieval，并将预训练相近/明显不同组区分。所有 mean/covariance、metric、posterior calibration 和码本仅train fitting。
3. **数字路径**：actualentropy encoded bytes/headers/atlas index+同 channelcode/modulation，任务loss对 total channeluses 或 latency；浮点维数仅辅助x轴。
4. **模拟路径**：固定kcomplexuses/avgpower、AWGN与blockRayleigh、SNRtrain/testshift、明确CSI/pilot成本。用同backbone比较 vanilla、LBD/RaDOGAGA regularizer、PH regularizer、ourgeometry、density-only/joint；只选择能识别目标增量的必要组合。
5. **能力控制**：同 parameter count/latentdimension/trainingbudget；直接label/posterior、LS+quantizer、PCA/orthogonal、entropy-only、geometry-only、joint；低弯曲与任务误差的关系须持出样本上报告，而非训练loss自证。

目前最有说服力的 TCOM 叙事不是“已有几何理论换语义名称”，而是“在明确 finite resource、task 与 receiver知识约束下，曲面局部结构何时改变实现代价；何时所有收益只是 density/linear alignment；用已知因子与异构模型配对实验验证边界”。
