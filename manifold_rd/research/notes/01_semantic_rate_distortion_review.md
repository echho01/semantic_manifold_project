# 语义率失真：核心文献调研（截至2026-10-26）

本轮优先核查2025–2026研究，保留建立问题口径所需的2022–2024文献。这里的17篇是相关性筛选，不宣称覆盖全部已发表工作。全文阅读深度为问题、定理陈述、方法和结果章节核查，不等于已独立复现或审计每个证明。

**下载版本特别说明：S07 的正式期刊记录为 TIT 2025，但所附公开 PDF 是 NeurIPS 2021 先行版，文件名按 2021_NeurIPS 标识。2026 年 10 月的充分约简新稿 G15 放在理论交叉文献中。**

## 首先区分的五种问题

1. **潜在语义的间接RD**：编码器观察X，语义S不可直接观察。\(\bar d_s(x,\hat s)=\mathbb E[d_s(S,\hat s)\mid X=x]\)，再求\(\inf I(X;\hat S,\hat X)\)。这是经典RD的有意义任务实例，不是超越Shannon容量。
2. **语义后验分布保真**：保持\(p(S\mid X)\)，而非只保持某一分类标签；Zhao等2025已有直接定理。
3. **感知/生成分布约束**：保持输出分布看起来真实；它不保证每个样本的内容正确。TIT的RDP论文可借证明技术，但不能一概叫语义RD。
4. **战略、计算、侧信息与保密约束**：操作模型改变带来新区域或边界，2026已有多条进展。
5. **任务范围与允许重构集合**：2026-09的最新工作已研究多任务可实现预测和restricted feature的rate penalty，任何‘语义流形受限重构’必须逐项对照。

## 对拟投TIT工作的判断（本调研的分析）

‘语义表示在流形上’、‘使用Fisher/KL’、‘把同义样本取商’、‘重新运行BA’均不能单独支撑强创新。可考虑的实质问题是：在明确任务后验模型和同一操作码率下，可实现语义分布族\(\mathcal M\)相对自由概率单纯形重构引入多少代价；何时为零，何时严格正，是否能被局部曲率/二阶基本形式及后验混合偏离控制。该提问本身也需与S06、S09、S15及经典受限重构字母表RD严格区别。应从一个能证明的有限或低维指数族实例入手，再扩到流形。

下面所有‘结论’是所引论文的结果，‘局限/可借鉴’包含本调研的独立判断。

## S01 · An Indirect Rate-Distortion Characterization for Semantic Sources: General Model and the Case of Gaussian Observation

- **作者 / 年份 / 场所**：Jiakun Liu; Shuo Shao; Wenyi Zhang; H. Vincent Poor；2022；IEEE Transactions on Communications。
- **状态与版本**：正式发表，70(9):5946–5959；下载为arXiv作者稿；本次核对版本日期：2022-06-01。
- **背景**：编码器只能看到观测X，语义S是不可直接观察的内在状态；既要恢复语义，又可能要求恢复观测。
- **创新方法**：把双失真问题化为间接率失真：以条件期望失真替代未知S，求I(X;Ŝ,X̂)的约束最小化；线性语义与高斯观测下化为误差协方差凸优化，满足共同可对角化条件时得到加权反向注水。
- **实验数据**：理论及合成高斯/矩阵数值例子，无自然语言或大规模真实语义数据。
- **主要结论**：建立联合语义—观测保真限，并区分不可观测性带来的误差与可由增加码率消除的编码误差。
- **局限与总结**：线性高斯闭式不覆盖一般非线性流形；语义的物理含义由任务和失真外部指定。
- **对拟稿的作用**：任何流形远程语义模型必须退化回此类经典间接RD；换掉欧氏失真本身不足构成新定理。
- **核对范围**：已核对原始摘要及正文相关方法/结果章节，未逐行审计全部证明
- **原始来源**：[来源1](https://arxiv.org/html/2201.12477v2)；[来源2](https://collaborate.princeton.edu/en/publications/an-indirect-rate-distortion-characterization-for-semantic-sources/)。DOI：[10.1109/TCOMM.2022.3194978](https://doi.org/10.1109/TCOMM.2022.3194978)。

## S02 · The Role of Fidelity in Goal-Oriented Semantic Communication: A Rate Distortion Approach

- **作者 / 年份 / 场所**：Photios A. Stavrou; Marios Kountouris；2023；IEEE Transactions on Communications。
- **状态与版本**：正式发表，71(7):3918–3931；公开作者稿；本次核对版本日期：2023-05-09 publication record。
- **背景**：同一潜在语义与观测系统，采用不同保真准则会改变最优传输内容。
- **创新方法**：研究双失真率失真的上下界、紧致条件和参数化最优测试信道；给出语义感知Blahut–Arimoto算法；分析Hamming、擦除失真及向量高斯MSE。
- **实验数据**：有限字母表二元信源和合成联合高斯信源的解析/数值例子，无真实图像语料。
- **主要结论**：失真类型和约束激活关系可以决定系统优先传语义还是观测，语义指标具有实质设计后果。
- **局限与总结**：语义感知BA与多约束率失真已有成熟先例；本文并未自动学习语义含义。
- **对拟稿的作用**：必须解释新几何失真对应什么任务损失，并证明超出既有任意single-letter失真框架的结构结果。
- **核对范围**：已核对原始摘要及正文相关方法/结果章节，未逐行审计全部证明
- **原始来源**：[来源1](https://www.eurecom.fr/en/publication/6943)；[来源2](https://www.eurecom.fr/publication/6943/download/comsys-publi-6943.pdf)。DOI：[10.1109/TCOMM.2023.3274122](https://doi.org/10.1109/TCOMM.2023.3274122)。

## S03 · Rate-Distortion-Perception Tradeoff Based on the Conditional-Distribution Perception Measure

- **作者 / 年份 / 场所**：Sadaf Salehkalaibar; Jun Chen; Ashish Khisti; Wei Yu；2024；IEEE Transactions on Information Theory。
- **状态与版本**：正式发表，70(12):8432–8454；arXiv v2含期刊接收信息；本次核对版本日期：2025-04-25。
- **背景**：边缘分布真实性与逐样本重构不同；无公共随机性时如何给出可操作的感知率失真定理？
- **创新方法**：把感知约束定义在给定编码输出的条件分布上，推导离散和连续信源的single-letter表征；连续MSE/Wasserstein情形可用解码端加噪；高斯向量反向注水。
- **实验数据**：离散、向量高斯与高斯混合的理论/数值分析，无实际视觉benchmark。
- **主要结论**：其条件感知模型具有可实现编码刻画；零感知损失下与对应边缘感知模型一致，并有文中条件下3 dB最小失真代价。
- **局限与总结**：3 dB结论有模型和失真条件；分布真实性不能自动保证样本语义正确。
- **对拟稿的作用**：借鉴条件分布、公共随机性假设与编码定理写法；明确语义后验约束和感知分布约束不同。
- **核对范围**：已核对原始摘要及正文相关方法/结果章节，未逐行审计全部证明
- **原始来源**：[来源1](https://arxiv.org/html/2401.12207v2)；[来源2](https://experts.mcmaster.ca/scholarly-works/3478154)。DOI：[10.1109/TIT.2024.3467282](https://doi.org/10.1109/TIT.2024.3467282)。

## S04 · On the Rate-Distortion Theory for Task-Specific Semantic Communication

- **作者 / 年份 / 场所**：Jingxuan Chai; Huixiang Zhu; Yong Xiao; Guangming Shi; Ping Zhang；2025；Entropy。
- **状态与版本**：正式开放发表，27(8):775；本次核对版本日期：2025-07-23。
- **背景**：分类和生成任务对错误的容忍方式不同，单纯MSE不能同时描述它们。
- **创新方法**：在率失真中加入任务相关统计散度；对高斯源分析Wasserstein、KL及反向KL约束，并构建对应学习目标。
- **实验数据**：理论高斯源；MNIST与CIFAR-10生成实验，MNIST分类实验；对照MSE/VAE式目标。
- **主要结论**：文中生成任务里反向KL更有利于感知质量，分类任务里KL更有利于分类表现；三者不存在脱离任务的统一优劣排序。
- **局限与总结**：特定网络和数据集结果不构成所有语义场景中散度优劣的普遍定理。
- **对拟稿的作用**：不能仅以使用KL、Wasserstein或增加一个语义损失项主张创新；可将其作为任务/感知基线。
- **核对范围**：已核对原始摘要及正文相关方法/结果章节，未逐行审计全部证明
- **原始来源**：[来源1](https://pmc.ncbi.nlm.nih.gov/articles/PMC12385448/)；[来源2](https://www.mdpi.com/1099-4300/27/8/775)。DOI：[10.3390/e27080775](https://doi.org/10.3390/e27080775)。

## S05 · Extended Blahut–Arimoto Algorithm for Semantic Rate-Distortion Function

- **作者 / 年份 / 场所**：Yuxin Han; Yang Liu; Yaping Sun; Kai Niu; Nan Ma; Shuguang Cui; Ping Zhang；2025；Entropy。
- **状态与版本**：正式开放发表，27(6):651；本次核对版本日期：2025-06-18。
- **背景**：同义映射把多个语法符号归为同一语义，但一般源的语义RD缺少可计算表达。
- **创新方法**：固定同义映射下交替更新测试信道和重构分布，形成EBA；未知映射时结合模拟退火合并语法符号，寻找满足约束的大同义集合，并把语义知识库SKB视为映射实例。
- **实验数据**：高斯源N(0,1)在[-10,10]离散为201点；CUB鸟类特征实验取50类各2个样本，共100样本，对比10/20/50类SKB。属于小规模特征数值实验，无端到端语言理解评测。
- **主要结论**：论文定义下，扩大同义集合或增大已共享SKB可降低所计算的语义RD曲线；同义数为1时回到其经典对照。
- **局限与总结**：§3.2明确其语义RD可为负，因此不能直接解释成非负的物理通信码长；SKB收益须计入共享知识假设。同义映射改变保真需求，模拟退火也不自动提供全局最优证书。
- **对拟稿的作用**：作为语义等价类/商空间及数值BA已有工作的对照；单纯在流形网格上跑BA不足为主创新。
- **核对范围**：已核对下载PDF全文的方法、算法和实验章节，未独立复现证明与代码
- **原始来源**：[来源1](https://www.mdpi.com/1099-4300/27/6/651)；[来源2](https://pmc.ncbi.nlm.nih.gov/articles/PMC12192163/)。DOI：[10.3390/e27060651](https://doi.org/10.3390/e27060651)。

## S06 · Semantic Rate-Distortion Theory with Applications

- **作者 / 年份 / 场所**：Yi-Qun Zhao; Zhi-Ming Ma; Geoffrey Ye Li; Shuai Yuan; Tong Ye; Chuan Zhou；2025；arXiv。
- **状态与版本**：预印本；未从原始出版来源核实正式接收；本次核对版本日期：2025-09-12。
- **背景**：语义有歧义和多义，点估计不能充分描述理解是否保留。
- **创新方法**：约束p(S|X)与由重构Y诱导的p(S|Y)之间分歧，同时约束符号失真；给出编码表征及二元闭式，MNIST自编码器加入后验KL损失。
- **实验数据**：二元源解析与MNIST；MLP自编码器，d维L级量化，报告R=d log L这一编码率上界。
- **主要结论**：报告12 bit设定下分类准确率从纯MSE的62.41%提升至加入小权重语义约束后的90%以上；这是特定实验结果。
- **局限与总结**：使用同一预训练分类器训练/评估、近似单点标签后验及码率上界；不是对真实语义或实际熵编码极限的全面验证。
- **对拟稿的作用**：与‘后验语义流形’高度重叠：后验分布保真和其RD定理不是空白；其p(S|Y)也隐含可实现重构限制。
- **核对范围**：已核对原始摘要及正文相关方法/结果章节，未逐行审计全部证明
- **原始来源**：[来源1](https://arxiv.org/abs/2509.10061)；[来源2](https://arxiv.org/html/2509.10061v1)。DOI：[10.48550/arXiv.2509.10061](https://doi.org/10.48550/arXiv.2509.10061)。

## S07 · Universal Rate-Distortion-Perception Representations for Lossy Compression

- **作者 / 年份 / 场所**：George Zhang; Jingjing Qian; Jun Chen; Ashish Khisti；2025；IEEE Transactions on Information Theory。
- **状态与版本**：TIT正式发表，71(11):8633–8653；公开arXiv仅核到2021 NeurIPS版本，须区分版本；本次核对版本日期：2025-07-14 journal online; 2021-12-22 arXiv。
- **背景**：每个失真—感知偏好都重训编码器会造成部署开销，能否固定表示只改解码器？
- **创新方法**：定义通用RDP表示并研究可实现性、固定表示的失真—感知域和渐进细化；高斯MSE下同码率单编码器可以覆盖相应权衡。
- **实验数据**：理论高斯源；已核读会议前身实验为MNIST和SVHN，固定与可变编码器对照；期刊摘要也报告这两组数据。
- **主要结论**：通用表示在文中高斯模型可不损失对应最优权衡，在所测图像模型仅有较小代价。
- **局限与总结**：通用性依赖源模型、损失和允许解码器族；不等于一份语义表示对任意任务无损通用。
- **对拟稿的作用**：提出多任务/多几何重构的通用码时必须对比该框架；先明示decoder family再比较rate penalty。
- **核对范围**：已核对原始摘要及正文相关方法/结果章节，未逐行审计全部证明
- **原始来源**：[来源1](https://ieeexplore.ieee.org/document/11080131/)；[来源2](https://arxiv.org/html/2106.10311v3)；[来源3](https://experts.mcmaster.ca/scholarly-works/3617719)。DOI：[10.1109/TIT.2025.3584291](https://doi.org/10.1109/TIT.2025.3584291)。

## S08 · Value of Communication in Goal-Oriented Semantic Communications: A Pareto Analysis

- **作者 / 年份 / 场所**：Jiping Luo; Bowen Li; Nikolaos Pappas；2026；IEEE Transactions on Information Theory。
- **状态与版本**：原作者arXiv确认TIT录用；机构记录检索到正式DOI，未将未核实卷页写入；本次核对版本日期：2026-07-31。
- **背景**：动态系统不必耗尽通信预算，核心是每次通信对任务的边际价值。
- **创新方法**：以双目标MDP同时优化历史相关语义估计代价和通信频率；刻画凸分段线性Pareto前沿，提出SPLIT算法。
- **实验数据**：八状态Markov系统（四个正常、四个告警状态），丢包率0.3；随机策略评估使用10^7时隙。纯合成动态系统实验。
- **主要结论**：前沿角点由确定性策略给出，邻近角点混合构造其他前沿点；通信收益存在递减。
- **局限与总结**：这是动态目标导向通信的代价前沿，不是静态语义流形Shannon RD函数，不能混用码率单位。
- **对拟稿的作用**：展示TIT语义研究的合适贡献强度：可证明结构、最优性与算法；历史依赖任务是另一个独立方向。
- **核对范围**：已核对原始摘要及正文相关方法/结果章节，未逐行审计全部证明
- **原始来源**：[来源1](https://arxiv.org/html/2512.01454v2)；[来源2](https://arxiv.org/abs/2512.01454)；[来源3](https://www.diva-portal.org/smash/record.jsf?pid=diva2%3A2108259)。DOI：[10.1109/TIT.2026.3720415](https://doi.org/10.1109/TIT.2026.3720415)。

## S09 · Log-Likelihood Loss for Semantic Compression

- **作者 / 年份 / 场所**：Anuj Kumar Yadav; Dan Song; Yanina Y. Shkel; Ayfer Özgür；2026；IEEE International Symposium on Information Theory。
- **状态与版本**：ISIT2026正式会议，EPFL作者机构记录核实；下载长版arXiv；本次核对版本日期：2026-01-23。
- **背景**：接收端的输出可以是能生成原信号的抽象语义描述，而不是像素级近似。
- **创新方法**：固定条件生成分布P(X|U)，用d(x,u)=-log P(x|u)定义失真，刻画对应RD及其与log-loss、一般RD和零感知损失的联系。
- **实验数据**：二元/Hamming、高斯/平方误差等解析例子与四幅理论示意；未报告大规模真实图像/语言训练基准。
- **主要结论**：在文中条件下可把部分经典失真RD转换成log-likelihood问题，也构造达到完美感知的编码策略。
- **局限与总结**：固定似然模型的正确性和重构字母表是模型假设；负对数似然不自动等同人类语义，连续密度也有单位/基准测度问题。
- **对拟稿的作用**：生成器诱导语义流形失真和受限decoder需要直接比较此文，避免把固定生成模型损失误称全新。
- **核对范围**：已核对原始摘要及正文相关方法/结果章节，未逐行审计全部证明
- **原始来源**：[来源1](https://arxiv.org/html/2601.16461v1)；[来源2](https://www.epfl.ch/schools/ic/ipg/publications/)。DOI：[10.1109/ISIT62367.2026.11654145](https://doi.org/10.1109/ISIT62367.2026.11654145)。

## S10 · Semantic Rate Distortion and Posterior Design: Compute Constraints, Multimodality, and Strategic Inference

- **作者 / 年份 / 场所**：Emrah Akyol；2026；IEEE Journal on Selected Areas in Information Theory。
- **状态与版本**：正式发表，7:246–262（SUNY机构记录）；下载2026-02作者预印本，不是最终排版版；本次核对版本日期：2026-02-03 arXiv; journal 2026。
- **背景**：发送方和接收方可有不同语义目标，观测和计算资源也受限。
- **创新方法**：在线性高斯、二次目标和MMSE接收机下，把战略通信化为后验协方差设计和log-det信息预算，研究直接、远程、完整信息三种情形以及语义注水。
- **实验数据**：公开长版为理论推导与解析模型，无真实LLM、图像或硬件训练实验。
- **主要结论**：刻画观测模式和目标错位如何改变RD；多模态增加可观测成分；在每步信息预算模型中得到超额失真随总预算指数衰减。
- **局限与总结**：compute-as-rate是抽象假设，不能直接当成现实FLOPs或LLM推理深度的普遍定律；非高斯与近似推断尚需扩展。
- **对拟稿的作用**：后验几何、后验设计和语义waterfilling已出现；可借鉴远程误差地板及明确可观测量，不能据此宣称流形方法首创后验建模。
- **核对范围**：已核对原始摘要及正文相关方法/结果章节，未逐行审计全部证明
- **原始来源**：[来源1](https://researchconnect.suny.edu/en/publications/semantic-rate-distortion-and-posterior-design-compute-constraints/)；[来源2](https://arxiv.org/html/2602.03949v1)。DOI：[10.1109/JSAIT.2026.3699650](https://doi.org/10.1109/JSAIT.2026.3699650)。

## S11 · On the Rate-Distortion-Complexity Tradeoff for Semantic Communication

- **作者 / 年份 / 场所**：Jingxuan Chai; Yong Xiao; Guangming Shi；2026；IEEE Internet of Things Journal。
- **状态与版本**：原作者arXiv确认接收，并给出期刊DOI；本次核对版本日期：2026-05-07。
- **背景**：神经语义编解码的通信节省需要计算资源，传统RD未明确这一成本。
- **创新方法**：引入I(X;U)型复杂度约束和统计语义距离，求高斯/二元RDC闭式，并提出变分RDC学习目标。
- **实验数据**：MNIST/CIFAR-10分类和生成；对比DeepJSCC及取消复杂度项版本，测FLOPs、准确率和Inception Score。
- **主要结论**：特定低计算量设定中，复杂度约束版本有更好生成表现；论文报告CIFAR-10在10.683M FLOPs处IS改善7.27%。
- **局限与总结**：I(X;U)不是计算复杂度的恒等式；原文也指出其与FLOPs关系依数据集而变且可能饱和。v2正文为图像实验，不据摘要泛称已做视频验证。
- **对拟稿的作用**：可借鉴可测复杂度对照；主TIT方案宜先解决几何RD，避免把未验证的神经计算代理同时包装成新定理。
- **核对范围**：已核对原始摘要及正文相关方法/结果章节，未逐行审计全部证明
- **原始来源**：[来源1](https://arxiv.org/abs/2602.14481)；[来源2](https://arxiv.org/html/2602.14481v2)。DOI：[10.1109/JIOT.2026.3689652](https://doi.org/10.1109/JIOT.2026.3689652)。

## S12 · Rate–Distortion Limits for Task-Oriented Compression with Side Information

- **作者 / 年份 / 场所**：Tao Guo; Zhangyao Song; Huihui Wu; Yang Li；2026；Entropy。
- **状态与版本**：正式开放发表，28(6):593；本次核对版本日期：2026-05-26。
- **背景**：发送接收端共享背景信息时，恢复潜在语义与观测的所需码率应改变。
- **创新方法**：在S—X1—(X2,Y)及两端共享Y的模型下求I(X1,X2;X̂1,X̂2,Ŝ|Y)的三失真约束最小化；给出整数奇偶分类、高斯等特例，构造AE分类压缩验证。
- **实验数据**：合成二元/整数/高斯信源；真实MNIST、SVHN。图像中心区域定义为X1，外围为X2，原图顶部1/3作为两端免费侧信息Y；报告率R=k log2 L是固定量化码长上界。
- **主要结论**：共享Y在实验低码率区改善MSE与分类，收益随码率提高缩小；只用中心重构分类与完整重构分类呈不同质量/准确率权衡。
- **局限与总结**：理论假设perfect classifier，并允许解码器直接输出Ŝ，未强制固定任务头Ŝ=f(X̂)。共享原图片段成本未计入码率，且两端共享Y不同于仅解码端拥有Y。
- **对拟稿的作用**：固定任务头或可实现后验族对自由Ŝ带来的额外代价可成为动机；若研究上下文流形，必须固定侧信息可用性及成本，不能把本式当成任意受限解码器的答案。
- **核对范围**：已核对下载PDF全文的问题定义、定理和MNIST/SVHN实验章节，未独立复现证明与代码
- **原始来源**：[来源1](https://www.mdpi.com/1099-4300/28/6/593)；[来源2](https://pmc.ncbi.nlm.nih.gov/articles/PMC13298143/)；[来源3](https://research.nottingham.edu.cn/en/publications/ratedistortion-limits-for-task-oriented-compression-with-side-inf/)。DOI：[10.3390/e28060593](https://doi.org/10.3390/e28060593)。

## S13 · Secure Semantic Communication over Wiretap Channels: Rate-Distortion-Equivocation Tradeoff

- **作者 / 年份 / 场所**：Denis Kozlov; Mahtab Mirmohseni; Rahim Tafazolli；2026；IEEE Journal on Selected Areas in Information Theory。
- **状态与版本**：正式发表，7:120–133；arXiv列出刊发卷页及DOI；本次核对版本日期：2026-05-24。
- **背景**：语义和观测的保密需求不同；编码器是否能直接访问语义也会影响安全通信。
- **创新方法**：对记忆无关语义源/窃听信道设独立保真和保密约束，推导single-letter内外界；使用含两个私有部分的随机叠加编码。
- **实验数据**：高斯源—高斯窃听信道及Bernoulli源—二元窃听信道的数值边界，无现实网络实测。
- **主要结论**：得到rate–distortion–equivocation权衡界并比较不同语义可访问条件。
- **局限与总结**：一般问题给界而非全情形精确解；保密/多终端约束会显著提高拟稿难度。
- **对拟稿的作用**：展示最新语义RD向操作约束扩展的方向；当前几何方案不宜同时加入保密，除非可导出强新定理。
- **核对范围**：已核对原始摘要及正文相关方法/结果章节，未逐行审计全部证明
- **原始来源**：[来源1](https://arxiv.org/abs/2509.12142)；[来源2](https://arxiv.org/html/2509.12142v3)。DOI：[10.1109/JSAIT.2026.3676713](https://doi.org/10.1109/JSAIT.2026.3676713)。

## S14 · Semantic Rate-Distortion for Bounded Multi-Agent Communication: Capacity-Derived Semantic Spaces and the Communication Cost of Alignment

- **作者 / 年份 / 场所**：Anthony T. Nixon；2026；arXiv。
- **状态与版本**：预印本，未核实同行评审发表；本次核对版本日期：2026-04-10。
- **背景**：能力不同的智能体可能形成不同语义字母表，而非仅用不同码率压缩同一字母表。
- **创新方法**：用能力限制诱导的商POMDP表示语义，讨论语义失配阈值、单向无记忆Wyner–Ziv化约和对齐遍历界。
- **实验数据**：八个POMDP环境，包括RockSample(4,4)；作者报告结构化策略在特定基准相比计数界有最高19倍差异。
- **主要结论**：提出语义字母表/能力差异本身也产生通信成本；完整双向推广仍区分为猜想或额外条件。
- **局限与总结**：理论依商结构定义与假设；尚未核实同行评审，不能把作者宣称的所有定理视为已被独立确认。
- **对拟稿的作用**：‘语义商空间’和‘对齐成本’已有直接先例；其商POMDP非一般黎曼统计流形，需逐假设比较。
- **核对范围**：已核对原始摘要及正文相关方法/结果章节，未逐行审计全部证明
- **原始来源**：[来源1](https://arxiv.org/abs/2604.09521)；[来源2](https://arxiv.org/html/2604.09521v1)。DOI：[10.48550/arXiv.2604.09521](https://doi.org/10.48550/arXiv.2604.09521)。

## S15 · On Task Scope and Information Retention in Source Coding

- **作者 / 年份 / 场所**：Alireza Furutanpey; Kerstin Bunte；2026；arXiv。
- **状态与版本**：预印本；2026-09-29最新直接相关工作；本次核对版本日期：2026-09-29。
- **背景**：面向机器的编码不一定比面向人的编码便宜，真正决定所需信息的是任务集合、损失、风险容忍度、可用观测与允许解码。
- **创新方法**：扩展到有限多任务和可联合实现预测元组；刻画source/feature编码等价条件，以最优测试信道KL给出限制观测造成的码率惩罚界。
- **实验数据**：有限源与学习block code；ImageNet、CIFAR-100、Taskonomy及鸟类图像核查，关注实际发送消息与任务保持。
- **主要结论**：原任务上保持最优预测，不保证增加任务后仍与原信号编码同率；分别可实现的输出不保证联合可实现。
- **局限与总结**：未做一般黎曼流形曲率/高分辨率展开；当前仍是预印本。
- **对拟稿的作用**：‘可实现语义重构集合导致RD代价’最近的直接对照。应比较其有限预测集合与连续流形约束，避免只重命名允许的reproduction alphabet。
- **核对范围**：已核对原始摘要及正文相关方法/结果章节，未逐行审计全部证明
- **原始来源**：[来源1](https://arxiv.org/html/2609.37575v1)。DOI：[10.48550/arXiv.2609.37575](https://doi.org/10.48550/arXiv.2609.37575)。

## S16 · Fundamental Limitation of Semantic Communications: Neural Estimation for Rate-Distortion

- **作者 / 年份 / 场所**：Dongxu Li; Jianhao Huang; Chuan Huang; Xiaoqi Qin; Han Zhang; Ping Zhang；2023；Journal of Communications and Information Networks。
- **状态与版本**：2023期刊工作；公开arXiv于2024发布，正文确认2023-11-15接收；本次核对版本日期：2024-01-02。
- **背景**：高维实际源的分布未知，解析语义RD或离散BA难以直接计算。
- **创新方法**：提出生成网络驱动的神经语义RD估计NESRD；已知分布下使用通用BA；语义是观测确定函数时采用生成器—分类器级联。
- **实验数据**：合成高斯观测3维、语义2维；50000个观测、每个1000个条件语义样本；真实实验为MNIST和SVHN，与超先验压缩模型比较。
- **主要结论**：数值估计与可计算高斯基准吻合，并给图像语义RD参考；说明BA复杂度随高维字母表迅速增长。
- **局限与总结**：神经估计受有限样本、优化和函数族限制；原文指出可估计信息量受O(log样本数)限制，不能直接当成严格全局下界或实际codec。
- **对拟稿的作用**：几何RD实验应将神经估计、离散BA与实用编码率分开呈现；可借其高斯核对与已知分布基准。
- **核对范围**：已核对原始摘要及正文相关方法/结果章节，未逐行审计全部证明
- **原始来源**：[来源1](https://arxiv.org/html/2401.01176v1)；[来源2](https://arxiv.org/abs/2401.01176)。DOI：[10.23919/JCIN.2023.10387242](https://doi.org/10.23919/JCIN.2023.10387242)。

## S17 · UniTAC: Universal Task-Aware Compression via Weighted Distortion Measures

- **作者 / 年份 / 场所**：Homa Esfahanizadeh; Matin Mortaheb; Adeel Mahmood; Jinfeng Du; Harish Viswanathan；2026；arXiv。
- **状态与版本**：预印本；未从原始出版来源核实正式接收；本次核对版本日期：2026-09-11。
- **背景**：任务变化时逐任务重训codec成本高，希望一个可重构图像的编码器按当前任务分配码率。
- **创新方法**：以J_fᵀJ_f刻画局部任务二次失真，讨论对角权重的task-consistency；用Integrated Gradients生成重要性图，条件化ViT编解码器，运行时换图适配任务。
- **实验数据**：AffectNet训练codec、CelebA测试两个属性分类；嘴部任务在报告0.034 bpp时准确率91.4%，通用codec为76.9%。
- **主要结论**：单模型可接近逐任务codec；它直接连接解释性归因与任务加权率失真。
- **局限与总结**：报告bpp为熵模型估计且省略16×16权重图传输。线性精确一致性需要交叉误差项抵消；非线性黑图基线IG只是代理，同一分类器用于归因和评估。
- **对拟稿的作用**：拉回度量及IG驱动语义码率分配已有先例；可借对角/完整任务度量消融，主创新应落在可实现后验族的严格RD差价。
- **核对范围**：已核对arXiv v2的理论假设、方法及实验章节，未独立复现
- **原始来源**：[来源1](https://arxiv.org/abs/2608.16696)；[来源2](https://arxiv.org/html/2608.16696v2)。DOI：[10.48550/arXiv.2608.16696](https://doi.org/10.48550/arXiv.2608.16696)。


## 包内 PDF 快速链接

- **S01**：[An Indirect Rate-Distortion Characterization for Semantic Sources: General Model and the Case of Gaussian Observation](../../../refs/semantic_rd/2022_TCOM_An_Indirect_Rate_Distortion_Characterization_for_Semantic_Sources_General_Model_and_the_Case_of_Gaussian_Observation.pdf)
- **S02**：[The Role of Fidelity in Goal-Oriented Semantic Communication: A Rate Distortion Approach](../../../refs/semantic_rd/2023_TCOM_The_Role_of_Fidelity_in_Goal_Oriented_Semantic_Communication_A_Rate_Distortion_Approach.pdf)
- **S03**：[Rate-Distortion-Perception Tradeoff Based on the Conditional-Distribution Perception Measure](../../../refs/semantic_rd/2024_TIT_Rate_Distortion_Perception_Tradeoff_Based_on_the_Conditional_Distribution_Perception_Measure.pdf)
- **S04**：[On the Rate-Distortion Theory for Task-Specific Semantic Communication](../../../refs/semantic_rd/2025_Entropy_On_the_Rate_Distortion_Theory_for_Task_Specific_Semantic_Communication.pdf)
- **S05**：[Extended Blahut–Arimoto Algorithm for Semantic Rate-Distortion Function](../../../refs/semantic_rd/2025_Entropy_Extended_BlahutArimoto_Algorithm_for_Semantic_Rate_Distortion_Function.pdf)
- **S06**：[Semantic Rate-Distortion Theory with Applications](../../../refs/semantic_rd/2025_arXiv_Semantic_Rate_Distortion_Theory_with_Applications.pdf)
- **S07**：[Universal Rate-Distortion-Perception Representations for Lossy Compression](../../../refs/semantic_rd/o_2021_NeurIPS_Universal_Rate_Distortion_Perception_Representations_for_Lossy_Compression.pdf)
- **S08**：[Value of Communication in Goal-Oriented Semantic Communications: A Pareto Analysis](../../../refs/semantic_rd/2026_TIT_Value_of_Communication_in_Goal_Oriented_Semantic_Communications_A_Pareto_Analysis.pdf)
- **S09**：[Log-Likelihood Loss for Semantic Compression](../../../refs/semantic_rd/2026_ISIT_Log_Likelihood_Loss_for_Semantic_Compression.pdf)
- **S10**：[Semantic Rate Distortion and Posterior Design: Compute Constraints, Multimodality, and Strategic Inference](../../../refs/semantic_rd/2026_JSAIT_Semantic_Rate_Distortion_and_Posterior_Design_Compute_Constraints_Multimodality_and_Strategic_Inference.pdf)
- **S11**：[On the Rate-Distortion-Complexity Tradeoff for Semantic Communication](../../../refs/semantic_rd/2026_IoTJ_On_the_Rate_Distortion_Complexity_Tradeoff_for_Semantic_Communication.pdf)
- **S12**：[Rate–Distortion Limits for Task-Oriented Compression with Side Information](../../../refs/semantic_rd/2026_Entropy_RateDistortion_Limits_for_Task_Oriented_Compression_with_Side_Information.pdf)
- **S13**：[Secure Semantic Communication over Wiretap Channels: Rate-Distortion-Equivocation Tradeoff](../../../refs/semantic_rd/2026_JSAIT_Secure_Semantic_Communication_over_Wiretap_Channels_Rate_Distortion_Equivocation_Tradeoff.pdf)
- **S14**：[Semantic Rate-Distortion for Bounded Multi-Agent Communication: Capacity-Derived Semantic Spaces and the Communication Cost of Alignment](../../../refs/semantic_rd/o_2026_arXiv_Semantic_Rate_Distortion_for_Bounded_Multi_Agent_Communication_Capacity_Derived_Semantic_Spaces_and_the_Communication_Cost_of_Alignment.pdf)
- **S15**：[On Task Scope and Information Retention in Source Coding](../../../refs/semantic_rd/o_2026_arXiv_On_Task_Scope_and_Information_Retention_in_Source_Coding.pdf)
- **S16**：[Fundamental Limitation of Semantic Communications: Neural Estimation for Rate-Distortion](../../../refs/semantic_rd/2023_JCIN_Fundamental_Limitation_of_Semantic_Communications_Neural_Estimation_for_Rate_Distortion.pdf)
- **S17**：[UniTAC: Universal Task-Aware Compression via Weighted Distortion Measures](../../../refs/semantic_rd/2026_arXiv_UniTAC_Universal_Task_Aware_Compression_via_Weighted_Distortion_Measures.pdf)
