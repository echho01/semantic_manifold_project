# 05｜理论基础与最接近交叉文献：逐篇调研

核查截止：2026-10-26。以下 15 篇用于识别已有工作、选取证明工具和区分操作模型，并非全部属于语义通信或 TIT。G11 未找到公开全文，仅依据官方摘要作有限分析。

研究切口见 [03](03_crossover_and_feasibility.md)，完整方案见 [当前论文](../../main.tex)。

## G01｜High-resolution source coding for non-difference distortion measures: the rate-distortion function

**作者：** Tamás Linder；Ram Zamir。**年份/场所：** 1999，IEEE Transactions on Information Theory。**状态：** 正式期刊论文，45(2):533–547

**背景：** 常见高分辨率结果依赖差失真；感知或间接观测问题的有效失真通常随源点位置变化，最小可达失真也可能非零。

**创新方法/理论工具：** 对非差失真在最佳重构附近作局部二次分析，给出率失真函数的高分辨率形式，并把间接源编码通过条件期望失真归约到该框架。

**实验数据或理论例证：** 理论论文；主要证据是定理、规则性条件和解析例子，不是深度学习数据集基准。已核对期刊元数据与作者全文；现已从 Ram Zamir 官方目录取得 37 页、1998 年修订作者稿，并转换为 PDF 供阅读。

**结论：** 远程源的R(D)可通过d(x,a)=E[ρ(S,a)|X=x]研究，D趋近Bayes失真下限时也已有高分辨率理论。

**局限与总结：** 需要特定光滑性、最优重构和非退化条件；不能直接当作奇异/分层语义空间的完整解。

**对拟稿的作用：** 最重要的反重叠文献：不能把Bayes失真下限、局部Hessian和高分辨率公式的组合直接宣称新贡献。

**核查范围：** 已核对作者全文中的高分辨率非差失真与间接源编码部分；已取得公开作者稿。

**版本说明：** 包中为作者官网 nd.ps.gz 转换的 37 页 1998 revised 作者稿，对应 1999 年 TIT 正式论文；不是出版社排版版。原始作者链接：https://www.eng.tau.ac.il/~zamir/papers/nd.ps.gz 。保留既有 pdf_url 仅用于记录最初检索到的旧作者端点。

**来源：** [原始记录](https://ieeexplore.ieee.org/document/749001/)；[DOI](https://doi.org/10.1109/18.749001)。

**包内全文：** [PDF](../../../refs/theory_overlap/1999_TIT_High_resolution_source_coding_for_non_difference_distortion_measures_the_rate_distortion_function.pdf)，37 页。

## G02｜Entropy and Source Coding for Integer-Dimensional Singular Random Variables

**作者：** Günther Koliander；Georg Pichler；Erwin Riegler；Franz Hlawatsch。**年份/场所：** 2016，IEEE Transactions on Information Theory。**状态：** 正式期刊论文，62(11):6124–6154；下载为arXiv v3作者版，2017-01-03更新

**背景：** 低维流形上的随机变量相对于环境空间Lebesgue测度是奇异的，普通微分熵不能直接使用。

**创新方法/理论工具：** 以Hausdorff测度建立整数维奇异变量的熵、联合熵和条件熵，研究Lipschitz变换和互信息，推出量化编码长度与Shannon下界。

**实验数据或理论例证：** 纯数学推导及奇异变量构造示例，无现实机器学习数据集。

**结论：** 为流形支持源建立合适的信息量表达；信息量应相对于实际支持维数和参考测度计算。

**局限与总结：** 主要针对整数维可求长结构；不直接提供任务语义定义，也不等价于任意深度特征都满足其正则条件。

**对拟稿的作用：** 定义流形熵与理论下界时必须引用；不能对低维支持的高维embedding直接代入普通环境空间微分熵。

**来源：** [原始记录](https://arxiv.org/abs/1505.03337)；[DOI](https://doi.org/10.1109/TIT.2016.2604248)。

**包内全文：** [PDF](../../../refs/theory_overlap/2016_TIT_Entropy_and_Source_Coding_for_Integer_Dimensional_Singular_Random_Variables.pdf)，31 页。

## G03｜Asymptotic quantization for probability measures on Riemannian manifolds

**作者：** Mikaela Iacobelli。**年份/场所：** 2016，ESAIM: Control, Optimisation and Calculus of Variations。**状态：** 正式期刊论文，22(3):770–785

**背景：** 欧氏最优量化理论不能无条件延拓到一般非紧黎曼流形，几何增长与概率尾部会影响渐近结果。

**创新方法/理论工具：** 在与曲率有关的无穷远增长条件下，证明黎曼流形概率测度的量化误差渐近，并构造说明条件尖锐性的例子。

**实验数据或理论例证：** 理论证明；流形及概率测度的解析例子，无神经网络训练数据。

**结论：** 最优有限码本量化的幂律指数与流形维数有关；曲率还通过全局增长条件影响结果是否适用。

**局限与总结：** 固定码本数的量化误差不能不加说明地等同互信息定义的Shannon R(D)；也不能据此推出一个普适标量曲率二阶修正。

**对拟稿的作用：** N^(−2/k)并非新的语义理论；可借鉴量化上界和正则条件，同时清楚区分固定率向量量化与无限块长率失真函数。

**来源：** [原始记录](https://numdam.org/articles/10.1051/cocv/2015025/)；[DOI](https://doi.org/10.1051/cocv/2015025)。

**包内全文：** [PDF](../../../refs/theory_overlap/2016_ESAIM_COCV_Asymptotic_quantization_for_probability_measures_on_Riemannian_manifolds.pdf)，16 页。

## G04｜Lossy Compression of General Random Variables

**作者：** Erwin Riegler；Günther Koliander；Helmut Bölcskei。**年份/场所：** 2023，Information and Inference。**状态：** 正式期刊论文，12(3):1759–1829；不是TIT

**背景：** 数据可取值于流形、分形和一般可测空间，需要超出欧氏密度的统一有损压缩分析。

**创新方法/理论工具：** 利用参考测度的小球概率正则性，建立一般Shannon型率失真下界和量化误差上下界；允许不同源与重构空间及非对称失真。

**实验数据或理论例证：** 超球面、Grassmann流形、自相似分形等解析案例和界的比较，不以现实语义通信训练数据为主。

**结论：** 几何可以通过小球体积和测度进入压缩极限；其框架覆盖大量已有特例。

**局限与总结：** 一般界不自动给精确最优函数，也不自动解释何种任务决定语义等价；参数依赖的统一尖锐界仍需另证。

**对拟稿的作用：** 一般流形RD和小球界已经成熟。拟稿必须提供明确操作问题及超出套用该框架的定理，例如比较流形内外重构最优曲线的差。

**来源：** [原始记录](https://academic.oup.com/imaiai/article/12/3/1759/7189915)；[DOI](https://doi.org/10.1093/imaiai/iaac035)。

**包内全文：** [PDF](../../../refs/theory_overlap/2023_Information_and_Inference_Lossy_Compression_of_General_Random_Variables.pdf)，72 页。

## G05｜On the Entropy of Rectifiable and Stratified Measures

**作者：** Juan Pablo Vigneaux。**年份/场所：** 2023，GSI。**状态：** GSI2023会议论文，LNCS14071:338–346；下载为arXiv作者版

**背景：** 单一维数的流形假设不能表达由不同维数成分混合而成的奇异分布。

**创新方法/理论工具：** 综合几何测度论与熵分解，研究可求长测度及其混合构成的分层测度，并给出链式法则与典型维数的AEP。

**实验数据或理论例证：** 理论工作，没有现实数据集实验。

**结论：** 分层分布的熵和典型维数有既有数学框架；离散连续混合只是其特例。

**局限与总结：** 不是语义通信R(D)的完整解，不能把熵/AEP直接当特定KL失真的精确率失真公式。

**对拟稿的作用：** 若采用秩变化或多维语义层作为切入点，必须与此及信息维数理论对照；平均维数领先阶通常不足以支撑首创。

**来源：** [原始记录](https://arxiv.org/abs/2306.00305)。

**包内全文：** [PDF](../../../refs/theory_overlap/2023_GSI_On_the_Entropy_of_Rectifiable_and_Stratified_Measures.pdf)，9 页。

## G06｜Curved representational Bregman divergences and their applications

**作者：** Frank Nielsen。**年份/场所：** 2025，GSI。**状态：** GSI2025有会议初版；下载为2026-03-27的arXiv v5扩充版，33页，非TIT

**背景：** 把Bregman散度限制在非仿射模型族上后，普通算术质心可能不属于允许模型族。

**创新方法/理论工具：** 定义curved representational Bregman散度，证明受限重心等于全空间重心的右Bregman投影；给出约束损失与Bregman信息的分解。

**实验数据或理论例证：** 数学推导及圆、对称化散度、复高斯和α散度球交等示例，非通信数据集基准。

**结论：** 非凸约束产生的质心投影及其附加散度已经有明确理论；概率单纯形中的非线性模型族是自然应用对象。

**局限与总结：** 该工作并未据本次阅读直接解决编码器也联合优化时流形内外重构R(D)差的尖锐几何界。

**对拟稿的作用：** 主方案最接近的数学先例。固定编码器的KL质心分解、投影以及圆上的正代价都不能单独当创新；真正贡献须是优化后的操作率失真差及非平凡上下界/渐近。

**来源：** [原始记录](https://arxiv.org/abs/2504.05654)。

**包内全文：** [PDF](../../../refs/theory_overlap/2025_GSI_Curved_representational_Bregman_divergences_and_their_applications.pdf)，33 页。

## G07｜Latent Semantic Manifolds in Large Language Models

**作者：** Mohamed A. Mabrok。**年份/场所：** 2026，arXiv。**状态：** 2026-03-17预印本，未核实正式录用

**背景：** 语言模型内部连续表征与有限词表离散输出之间有几何不匹配。

**创新方法/理论工具：** 以Fisher度量赋予潜在语义流形几何，定义词表投影及expressibility gap，提出量化失真下界与边界附近体积增长结果。

**实验数据或理论例证：** WikiText-103验证集，每层约1800 hidden vectors；GPT-2/OPT/Pythia六模型，124M至1.5B；用TWO-NN和MLE维数估计、曲率与margin诊断。

**结论：** 作者报告中间层维数扩张、末层收缩及gap线性增长，斜率0.87–1.12、R²>0.985；提出有限词表N^(−2/k)下界。

**局限与总结：** 预印本的通用下界证明需特别复核：把小球渐近直接写成全局精确欧氏常数、以及吸收体积常数的步骤不宜直接继承；经验支持不等于证明真实hidden states是光滑流形。

**对拟稿的作用：** 已经直接覆盖‘语义流形+Fisher+失真+率失真’这一宽泛组合。可引用为近期相关工作并谨慎检验，不能声称该视角从未有人研究。

**来源：** [原始记录](https://arxiv.org/abs/2603.22301)。

**包内全文：** [PDF](../../../refs/theory_overlap/2026_arXiv_Latent_Semantic_Manifolds_in_Large_Language_Models.pdf)，38 页。

## G08｜GeoIB: Geometry-Aware Information Bottleneck via Statistical-Manifold Compression

**作者：** Weiqi Wang；Zhiyi Tian；Chenhan Zhang；Shui Yu。**年份/场所：** 2026，arXiv。**状态：** 2026-02-03预印本；正文有未替换ACM模板会议信息，不能据此认定2018或顶会发表

**背景：** 深度信息瓶颈依赖互信息代理与估计器，优化可能不稳定且压缩度受代理误差影响。

**创新方法/理论工具：** 用到独立性流形的KL投影解释IB，组合Fisher–Rao局部二次项、Jacobian–Frobenius项及自然梯度，提供信息几何训练目标。

**实验数据或理论例证：** MNIST、CIFAR-10、CelebA，MLP/CNN/ResNet18；K=128、β=10^(−4)主表，比较VIB/SIB/MINE/AIB，含MINE互信息估计和反演/成员攻击。

**结论：** 主表MNIST/CIFAR-10准确率99.28%/85.54%；但CIFAR-10报告MI=1.01，高于若干对照，运行时间也更长，不能表述成所有指标普遍优越。

**局限与总结：** Fisher局部近似不等于精确Shannon压缩率；后验噪声和局部化条件关系需要复核。未证实正式同行评审。

**对拟稿的作用：** 信息几何IB并非空白；可以借鉴JVP/Hutchinson计算，不能拿其神经代理直接当理论R(D)或编码长度。

**来源：** [原始记录](https://arxiv.org/abs/2602.03906)。

**包内全文：** [PDF](../../../refs/theory_overlap/2026_arXiv_GeoIB_Geometry_Aware_Information_Bottleneck_via_Statistical_Manifold_Compression.pdf)，12 页。

## G09｜Resource-efficient Semantic Coding Schemes with Manifold-constrained Hyper-connections

**作者：** Jingwen Fu；Ming Xiao。**年份/场所：** 2026，arXiv。**状态：** 2026-08-13 arXiv v1 预印本；未核实正式发表；已下载全文并核对实验和限制。

**背景：** 语义或任务通信需要在鲁棒性、带宽、功率和网络复杂度之间权衡。

**创新方法/理论工具：** 多残差流通过双随机混合矩阵的manifold-constrained hyper-connections结合entropy bottleneck，实现显式rate-distortion/task训练。

**实验数据或理论例证：** FineWeb10B 文本训练、WikiText-103 分布外测试、AG News 四分类；6/48 层 GPT-2 风格模型，隐藏维数 288/150、上下文 1024、信道瓶颈 k=64；AWGN 5–15 dB 训练 5000 步，测试 AWGN/Rayleigh/Rician 和不完美 CSI。Table VI 总参数约 43.3M/45.7M（正文约 20M 指骨干口径，不应与总参数混同）。

**结论：** 作者报告：48 层 Rayleigh 10 dB 的 FineWeb PPL 为 639，HC 为 917，标准 residual 为 1005；48 层 AWGN 10 dB 的 WikiText-103 PPL 为 2316，对照 HC 5235、residual 3623。6 层 AG News 在 10 dB 的 mHC 准确率为 47.6%/37.5%（AWGN/Rayleigh），低于 HC 的 51.9%/49.6%；48 层 mHC 为 31.5%/31.4%。

**局限与总结：** manifold 指双随机混合矩阵约束，不是语义后验支持流形。低 SNR 下 mHC 并非一致更优；跨组件比较中 Dense 在 Rayleigh 10 dB 的 PPL 527 优于 mHC 639。报告的码率按学习熵模型的负对数概率估计；不宜等同端到端实测比特流。差分熵到离散码率依赖高分辨率近似，RD 优越性定理假设假设类包含及全局优化。训练时间相对 residual 显著增加，不能仅凭小参数开销称总计算无额外成本。

**对拟稿的作用：** 证明语义通信工程文献已使用manifold与RD术语，需解释本方案研究的是后验重构可实现集合的几何和操作极限。

**核查范围：** 全文核对：系统模型、理论结论、VI 节实验、Tables II–VII。

**版本说明：** 已阅读根代理下载的 arXiv v1 PDF；本条取代早期“摘要级证据”记录。

**来源：** [原始记录](https://arxiv.org/abs/2608.13253)。

**包内全文：** [PDF](../../../refs/theory_overlap/2026_arXiv_Resource_efficient_Semantic_Coding_Schemes_with_Manifold_constrained_Hyper_connections.pdf)，13 页。

## G10｜Clustering with Bregman Divergences

**作者：** Arindam Banerjee；Srujana Merugu；Inderjit S. Dhillon；Joydeep Ghosh。**年份/场所：** 2005，JMLR。**状态：** 正式期刊论文，6:1705–1749

**背景：** 欧氏、相对熵等不同聚类损失需要统一理解，尤其是压缩表示时的质心与信息损失。

**创新方法/理论工具：** 利用Bregman信息及其损失统一硬/软聚类，建立与规则指数族的对应，并讨论率失真权衡。

**实验数据或理论例证：** 以理论与聚类方法为主体；本次未逐项复核全文实验数据集，因此不列未经核验的基准或数值。

**结论：** 算术质心、Bregman偏差—方差分解和压缩的信息损失具有成熟框架。

**局限与总结：** 这是基础数学/算法，不是流形受限语义通信的完整编码理论。

**对拟稿的作用：** KL质心恒等式来自已知Bregman结构。拟议模型需引用并把贡献放在受限与自由最优编码曲线比较。

**来源：** [原始记录](https://jmlr.org/papers/v6/banerjee05b.html)。

**包内全文：** [PDF](../../../refs/theory_overlap/2005_JMLR_Clustering_with_Bregman_Divergences.pdf)，45 页。

## G11｜Rate-Distortion in Non-Convex Families

**作者：** Hila Ratson；Ram Zamir。**年份/场所：** 2023，ITW。**状态：** 正式会议论文，ITW2023:87–91；未找到已核验公开PDF

**背景：** 参数化Blahut算法面对非凸码字位置/分布族时，传统全局收敛结论不再自动成立。

**创新方法/理论工具：** 构造非凸族下参数化Blahut失去全局最优的例子，并结合随机退火改进搜索。

**实验数据或理论例证：** 本次只核对官方机构摘要与书目信息，未获取全文实验细节；不声称完成全文精读。

**结论：** 参数化Blahut求得的局部解不能自动当真实R(D)。

**局限与总结：** 算法研究；并非本方案的几何差值定理。

**对拟稿的作用：** 对连续非凸语义流形的计算需网格误差控制、上下界或证书；多初始值优化仅是经验工具。

**来源：** [原始记录](https://cris.iucc.ac.il/en/publications/rate-distortion-in-non-convex-families/)；[DOI](https://doi.org/10.1109/ITW55543.2023.10160233)。

**全文：** 本次未取得公开 PDF，摘要级阅读不能替代证明审查。

## G12｜Estimating the reach of a manifold via its convexity defect function

**作者：** Clément Berenfeld；John Harvey；Marc Hoffmann；Krishnan Shankar。**年份/场所：** 2022，Discrete and Computational Geometry。**状态：** 正式期刊论文，67:403–438；2021在线，2022卷期；下载为arXiv作者版

**背景：** 流形学习需要描述局部弯曲和全局自靠近的reach参数，并从点云估计它。

**创新方法/理论工具：** 连接reach与凸性缺陷函数，结合流形重构建立估计量及一致期望误差界，并比较极小极大下界。

**实验数据或理论例证：** 主要为数学证明；构造几何例子，没有语义通信或LLM基准。

**结论：** 凸组合离开流形的程度与reach相关，已有正式数学理论；C³/C⁴模型下估计率接近相应下界。

**局限与总结：** 研究对象是几何参数估计，不是率失真函数；把几何界转换到优化后的KL率差需要另外证明。

**对拟稿的作用：** 为质心离开流形的上界提供严谨先例，说明reach/convexity-defect本身也不是拟稿新贡献。

**来源：** [原始记录](https://arxiv.org/abs/2001.08006)；[DOI](https://doi.org/10.1007/s00454-021-00290-8)。

**包内全文：** [PDF](../../../refs/theory_overlap/2022_Discrete_and_Computational_Geometry_Estimating_the_reach_of_a_manifold_via_its_convexity_defect_function.pdf)，35 页。

## G13｜From Source Reconstruction to Predictive State Preservation: An Information-Theoretic Framework for AI-Native Communication

**作者：** Yi Wang；Linglong Dai。**年份/场所：** 2026，arXiv。**状态：** 2026-09-01 arXiv v1 预印本；未核实正式发表。

**背景：** AI 接收端以预测为目标，重构源样本未必是恰当保真要求；需区分传输损失与预测模型自身限制。

**创新方法/理论工具：** 以条件未来分布为 predictive state，建立任务充分性与最小性；用 proper loss 定义失真，分解通信损失、接收模型族损失与部署损失。有限字母、无记忆、匹配失真下给出源与后验状态编码的等价关系。

**实验数据或理论例证：** 理论论文；有独立 Bernoulli 变量的解析例子，把同样 1 bit 的预测损失分别归因为通信、模型族或部署。无现实数据集基准。

**结论：** log-loss 下通信项是丢失的条件互信息，模型族项是接收后验到允许预测族的 KL 近似误差；允许完整概率单纯形输出时连接 IB/log-loss 率失真。

**局限与总结：** 单字母编码结论有有限字母等限定；动态部分主要给预测状态表示性质，不是一般动态通信编码定理。

**对拟稿的作用：** 直接近邻：后验商空间、后验编码等价、family gap 都已存在。拟议工作必须推进到流形受限重构重新优化后的 RD 差值及几何界。

**核查范围：** 已核对 arXiv 元数据及 HTML 全文 IV、V、VI 节。

**来源：** [原始记录](https://arxiv.org/abs/2609.01131)。

**包内全文：** [PDF](../../../refs/theory_overlap/2026_arXiv_From_Source_Reconstruction_to_Predictive_State_Preservation_An_Information_Theoretic_Framework_for_AI_Native_Communication.pdf)，13 页。

## G14｜Source-Side Sufficiency for the Information Bottleneck: Exact Reduction and Finite-Block Equivalence

**作者：** Joss Armstrong。**年份/场所：** 2026，arXiv。**状态：** 2026-04-29 首稿；2026-08-03 arXiv v2；未核实正式发表。

**背景：** 任务无关的输入变化仍可能消耗编码率；需确定充分统计替换原源何时保留整个 IB 贸易曲线。

**创新方法/理论工具：** 在 C–Z–T 且 Z=phi(T) 条件下对源纤维内编码器作条件平均，保持任务互信息，并精确移除 I(X;T|Z)；在标准 Borel 空间证明 IB 曲线与最优解归约。

**实验数据或理论例证：** 理论证明、Gaussian IB 特例与精确有限离散例子；无真实深度模型数据集或实测通信码率。

**结论：** 有限目标字母及 log-loss 下，以 Z^n 替代 T^n 在每个块长和消息预算下都保持最优远程失真，不仅是渐近 RDF 等价。

**局限与总结：** 要求精确任务充分性；近似后验、模型失配或新增任务不自动继承等价。本文没有分析非凸可实现后验集相对自由预测的优化 RD 差值；源侧充分性归约本身加固定动作约束也通常仍成立，不能借此制造新颖性。

**对拟稿的作用：** 不能把后验充分性、IB 归约或有限块长源侧等价当创新；可实现集合限制的优化几何代价仍须单独研究。

**核查范围：** 已核对 arXiv v2 元数据与 HTML 全文归约定理、log-loss 操作形式及限制。

**来源：** [原始记录](https://arxiv.org/abs/2604.26744)。

**包内全文：** [PDF](../../../refs/theory_overlap/2026_arXiv_Source_Side_Sufficiency_for_the_Information_Bottleneck_Exact_Reduction_and_Finite_Block_Equivalence.pdf)，15 页。

## G15｜Task-Sufficient Contraction: Source Selection for Machine Information Interfaces

**作者：** Joss Armstrong。**年份/场所：** 2026，arXiv。**状态：** 2026-10-06 arXiv v1；注释说明更早版本于2026-09发于Zenodo；未核实正式发表。

**背景：** 接收机器的任务和允许动作可能在通信率、码本及求解器之前就确定；问题是在所有后续操作点上都保留任务能力的最小源描述。

**创新方法/理论工具：** 将所有允许动作产生相同regret的源状态合并；有限动作集下给出整条一步rate–regret曲线的保持。对二次损失、仿射可行动作集，用可行动作差分方向的投影刻画规范约简。

**实验数据或理论例证：** 19页理论稿，含解析构造和比较；本次未见真实数据集实验。元数据标明1幅图、1张表。

**结论：** 任务源约简可以先于具体失真目标固定；能源预算例子需要保留中心化负载，仅保留某个最优water-filling动作过粗。

**局限与总结：** 保持性依赖固定任务、动作集合和损失；主要一步有限动作及特定二次仿射结果，不等于任意动态任务、任意深度模型的通用通信定理。

**对拟稿的作用：** 最新近邻再次说明动作约束下的任务等价化与完整RD贸易曲线保持已有研究；本项目须限定为两类重构集合各自优化后的几何差值。

**核查范围：** 本代理核对原始arXiv元数据与摘要；semantic_rd代理另核PDF第3–4节。

**来源：** [原始记录](https://arxiv.org/abs/2610.08884)。

**包内全文：** [PDF](../../../refs/theory_overlap/2026_arXiv_Task_Sufficient_Contraction_Source_Selection_for_Machine_Information_Interfaces.pdf)，19 页。
