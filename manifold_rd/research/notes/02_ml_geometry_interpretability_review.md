# AI/ML会议流形、黎曼几何与可解释性文献调研

核查截止：**2026年10月26日**。本清单精选18篇正式主会论文，覆盖ICML、NeurIPS、ICLR，以及统计机器学习强相关会议AISTATS和计算机视觉主会CVPR/ECCV。它不是穷尽目录；所选论文分别承担度量、图册、维数诊断、概念解释和对照验证的作用。年份按正式会议年份填写，不把arXiv首次提交年或workshop年份冒充主会发表年。正文中的结果是原论文已报告的结果；末节的语义通信方案是本次调研提出的借鉴方向，未开展实验。

## 1. 阅读路线和主题映射

| 优先级 | 文献 | 借鉴重点 | 对拟议TIT工作的作用 |
|---|---|---|---|
| 第一批 | M04、M15 | 概率模型Fisher拉回、KL/Bregman、任务后验几何 | 确定失真与几何来源，识别已有理论 |
| 第一批 | M14、M09 | 多图册、局部维数及切/法空间 | 实现可实现重构并验证流形假设 |
| 第一批 | M11、M02 | 多维周期概念、概念干预 | 构造可控语义与验证语义保真 |
| 第二批 | M18、M13、M12 | 切向约束、密度测地线、离流形惩罚 | 区分支持约束和数据密度效应 |
| 第二批 | M08、M10、M16、M17 | 反事实、层级、组合泛化、相关性反例 | 防止对语义几何作过强解释 |
| 背景 | M01、M06、M03、M05、M07 | 潜空间几何、局部基、率代理、概念库、白盒网络 | 建立有意义的算法与解释基线 |

## 2. 逐篇调研

### M01. Latent Space Oddity: on the Curvature of Deep Generative Models

**ICLR 2018｜正式发表：主会论文｜生成潜空间黎曼几何**

作者：Georgios Arvanitidis；Lars Kai Hansen；Søren Hauberg。

**背景与问题：** 生成模型潜变量常用欧氏距离，但非线性解码和潜空间重参数化会使直线插值穿过低支持区域，欧氏聚类也未必对应观测变化。

**创新方法：** 从随机解码器的均值和标准差诱导期望拉回度量 G(z)=J_mu(z)^T J_mu(z)+J_sigma(z)^T J_sigma(z)，并用能反映离开训练支持后不确定性增大的方差模型构造测地线、黎曼聚类与采样。

**数据与实验设置：** MNIST三个三类子集{0,1,2}、{3,4,7}、{5,6,9}；每类随机取1,000张图，多次划分；另有视频帧卷积VAE、数字插值与潜空间随机游走。

**已报告结果与结论：** MNIST {0,1,2} 的聚类F分数由欧氏基线77.57±0.87%提高到黎曼方法94.28±1.14%；另两组也改善。测地线避开高解码不确定性区域，但效果依赖不确定性模型。

**局限与边界：** 解码空间中的几何相似性不自动等于任务语义相似性；较可信的方差外推是关键。论文没有定义语义通信的操作率失真函数。

**可借鉴到语义流形率失真的方法：** 作为欧氏距离与解码器拉回距离的历史基线；检验重构原型是否离开可实现支持。本文的度量本身不能作为新论文的理论创新。

**版本说明：** 正式发表年份为2018；arXiv首次提交2017，下载的公开作者版本可能带2021修订日期。

来源：[正式出版记录](https://openreview.net/forum?id=SJzRZ-WCZ)；[公开全文/作者稿](https://arxiv.org/html/1710.11379v3)；[PDF下载](https://arxiv.org/pdf/1710.11379)。

### M02. Concept Bottleneck Models

**ICML 2020｜正式发表：主会论文｜概念可解释性**

作者：Pang Wei Koh；Thao Nguyen；Yew Siang Tang；Stephen Mussmann；Emma Pierson；Been Kim；Percy Liang。

**背景与问题：** 预测网络难以让人理解、纠正中间决策。若语义概念可观测，希望通过明确概念变量连接输入与任务输出。

**创新方法：** 采用x→c→y概念瓶颈，比较独立、顺序、联合训练；测试时直接纠正概念并重新计算预测，检验概念层对决策的作用。

**数据与实验设置：** OAI膝关节X光36,369个数据点，10个临床概念和关节炎严重程度；CUB 11,788张图、200种鸟，使用112个处理后的二元属性。

**已报告结果与结论：** OAI顺序/联合模型RMSE约0.418，标准模型约0.441；CUB联合瓶颈错误率约0.199，标准模型约0.175，说明准确率代价依任务而异。概念纠正能降低任务误差。

**局限与边界：** 依赖概念标签与概念完整性；概念监督不保证概念间因果独立。CUB使用类级属性处理，不能等同于每个样本的精确语义。

**可借鉴到语义流形率失真的方法：** 提供可核验的语义变量和干预协议：传输后同时评估概念重构、下游任务以及人工纠错效果。可用CUB和合成可控属性构造语义原型；OAI并非最方便的开放起步数据。

**版本说明：** 使用ICML2020正式版本。

来源：[正式出版记录](https://proceedings.mlr.press/v119/koh20a.html)；[公开全文/作者稿](https://proceedings.mlr.press/v119/koh20a/koh20a.pdf)；[PDF下载](https://proceedings.mlr.press/v119/koh20a/koh20a.pdf)。

### M03. Learning Diverse and Discriminative Representations via the Principle of Maximal Coding Rate Reduction

**NeurIPS 2020｜正式发表：主会论文｜编码率代理与表示几何**

作者：Yaodong Yu；Kwan Ho Ryan Chan；Chong You；Chaobing Song；Yi Ma。

**背景与问题：** 交叉熵分类可能导致类内特征塌缩，单纯判别准确率不能描述可解释的低维结构及类内多样性。

**创新方法：** 提出最大编码率缩减MCR²：使用log-det近似编码率，最大化总体编码率减去按类加权的类内编码率。在假设下鼓励不同类别占据近正交低维子空间，同时保留类内变化。

**数据与实验设置：** CIFAR-10监督表示学习使用ResNet-18、128维表示；另有CIFAR-10/100、STL-10自监督和聚类实验，比较交叉熵及不同标签噪声。

**已报告结果与结论：** 在论文CIFAR-10、标签噪声比例0.5的配置中，MCR²准确率84.30%，交叉熵60.37%；展示更清晰的子空间结构。自监督直接目标会出现饱和，作者还使用调整后的CTRL策略。

**局限与边界：** log-det是特定分布/量化近似下的编码率代理，并非已经证明的语义信源操作码率；真实语义未必属于类条件线性子空间。

**可借鉴到语义流形率失真的方法：** 可以作为表示谱、类内多样性、子空间混叠的诊断和训练基线；必须另外报告真实熵编码比特，不能把MCR²目标直接当作R(D)。

**版本说明：** NeurIPS2020主会论文。

来源：[正式出版记录](https://papers.nips.cc/paper_files/paper/2020/hash/6ad4174eba19ecb5fed17411a34ff5e6-Abstract.html)；[公开全文/作者稿](https://papers.neurips.cc/paper_files/paper/2020/file/6ad4174eba19ecb5fed17411a34ff5e6-Paper.pdf)；[PDF下载](https://papers.neurips.cc/paper_files/paper/2020/file/6ad4174eba19ecb5fed17411a34ff5e6-Paper.pdf)。

### M04. Pulling back information geometry

**AISTATS 2022｜正式发表：主会论文｜Fisher拉回度量**

作者：Georgios Arvanitidis；Miguel González-Duque；Alison Pouplin；Dimitris Kalatzis；Søren Hauberg。

**背景与问题：** 欧氏观测误差无法统一处理离散、方向及其他概率型解码器；需要把解码分布可辨识性转成潜空间几何。

**创新方法：** 将输出统计模型的Fisher信息通过解码参数映射h(z)拉回，G(z)=J_h^T I(h(z))J_h；利用局部KL的二阶导近似度量，支持多种分布而不局限于高斯VAE。

**数据与实验设置：** MNIST数字1的四类VAE对比；含Normal、Bernoulli、Beta、Dirichlet、Exponential等解码模型的合成圆；CMU动作捕捉69_06序列，用26段骨骼的方向分布建模。

**已报告结果与结论：** 信息几何路径能够随模型分布解释插值，在动作数据中产生更合理的周期运动；对不确定性的正则化显著影响结果。数值度量近似在数据边界可能失稳。

**局限与边界：** Fisher度量仅刻画给定概率模型的局部可辨识性；学习到的输出分布失配会改变几何，数值近似也不自动正定。没有给出语义重构约束的率失真定理。

**可借鉴到语义流形率失真的方法：** 最直接的度量来源：把h换为任务后验q(y|z)，用任务KL而非像素差异定义局部成本，并单独检验退化秩和模型校准。新颖性必须落在可实现重构约束产生的额外代价。

**版本说明：** AISTATS2022正式主会；AISTATS是统计机器学习强相关会议，未标为NeurIPS/ICML/ICLR。

来源：[正式出版记录](https://proceedings.mlr.press/v151/arvanitidis22b.html)；[公开全文/作者稿](https://proceedings.mlr.press/v151/arvanitidis22b/arvanitidis22b.pdf)；[PDF下载](https://proceedings.mlr.press/v151/arvanitidis22b/arvanitidis22b.pdf)。

### M05. Post-hoc Concept Bottleneck Models

**ICLR 2023｜正式发表：主会论文｜事后概念可解释性**

作者：Mert Yuksekgonul；Maggie Wang；James Zou。

**背景与问题：** 从头训练概念瓶颈需要大量概念标注，而且难以利用现成高性能模型。

**创新方法：** 在冻结表示上用概念激活向量或CLIP文本构造概念库，训练稀疏概念预测器；PCBM-h再加入残差分支恢复准确率，支持对概念权重进行全局编辑。

**数据与实验设置：** CIFAR-10/100、COCO-Stuff、CUB以及HAM10000、ISIC皮肤图像；主要结果跨10个随机种子，分类与医学任务分别使用准确率、mAP、AUROC等指标。

**已报告结果与结论：** CIFAR-10原模型准确率0.888，纯PCBM为0.777±0.003，带残差PCBM-h为0.871±0.001；COCO原模型mAP0.770，PCBM-h为0.768。结果显示保真与可解释性的实际权衡。

**局限与边界：** 概念库遗漏会降低性能；残差分支恢复预测力但使全部决策不能再由概念解释。概念投影也不等于证明存在光滑语义流形。

**可借鉴到语义流形率失真的方法：** 低成本构造可解释码字/原型与反事实测试；若采用残差传输，应把概念比特和残差比特分别计数，同时报告未解释的任务收益。

**版本说明：** ICLR2023主会spotlight；arXiv首次发表于2022，不能按预印本年命名会议。

来源：[正式出版记录](https://openreview.net/forum?id=nA5AZ8CEyow)；[公开全文/作者稿](https://arxiv.org/html/2205.15480v2)；[PDF下载](https://arxiv.org/pdf/2205.15480)。

### M06. Understanding the Latent Space of Diffusion Models through the Lens of Riemannian Geometry

**NeurIPS 2023｜正式发表：主会论文｜扩散局部切空间**

作者：Yong-Hyun Park；Mingi Kwon；Jaewoong Choi；Junghyo Jo；Youngjung Uh。

**背景与问题：** 扩散模型的逐时刻潜空间没有像VAE那样直接给定的低维语义坐标，单一全局方向难以描述图像编辑。

**创新方法：** 用U-Net瓶颈映射的Jacobian构造拉回几何，并以SVD获取局部潜空间基；分析不同时间、样本和文本条件的切空间变化，用局部基及平行移动实现编辑。

**数据与实验设置：** CelebA-HQ、Flowers、FFHQ、ImageNet、LSUN等预训练扩散模型以及Stable Diffusion；条件几何比较还使用50条MS-COCO文本描述。

**已报告结果与结论：** 局部方向随扩散阶段呈现由粗到细的变化；复杂数据对应更强的切空间差异，文本条件的影响也随时间变化。论文展示编辑效果及切空间相似性分析。

**局限与边界：** 语义方向存在属性纠缠，例如胡须和性别；局部基不是唯一因果因素。不同模型和时刻间的几何可比性需要验证，未分析编码极限。

**可借鉴到语义流形率失真的方法：** 借鉴局部Jacobian谱、有效秩和Grassmann切空间距离，检验压缩前后及噪声前后语义子空间是否稳定；作为几何诊断而非语义定义。

**版本说明：** NeurIPS2023主会正式PDF。

来源：[正式出版记录](https://proceedings.neurips.cc/paper_files/paper/2023/hash/4bfcebedf7a2967c410b64670f27f904-Abstract.html)；[公开全文/作者稿](https://proceedings.neurips.cc/paper_files/paper/2023/file/4bfcebedf7a2967c410b64670f27f904-Paper-Conference.pdf)；[PDF下载](https://proceedings.neurips.cc/paper_files/paper/2023/file/4bfcebedf7a2967c410b64670f27f904-Paper-Conference.pdf)。

### M07. White-Box Transformers via Sparse Rate Reduction

**NeurIPS 2023｜正式发表：主会论文｜白盒表示与编码代理**

作者：Yaodong Yu；Sam Buchanan；Druv Pai；Tianzhe Chu；Ziyang Wu；Shengbang Tong；Benjamin D. Haeffele；Yi Ma。

**背景与问题：** Transformer编码器的模块通常靠经验设计，缺乏统一且可直接检查的优化目标解释。

**创新方法：** 把稀疏率缩减目标的迭代优化展开为CRATE网络；多头子空间注意力对应压缩步骤，前馈模块对应类似ISTA的稀疏化步骤。

**数据与实验设置：** ImageNet-1K及ReaL、CIFAR-10/100、Flowers102和Pets等分类/迁移实验，同时跟踪网络内部的率缩减和稀疏性。

**已报告结果与结论：** 相近约22M参数下，论文表1的CRATE-B在ImageNet-1K为70.8%，ViT-S为72.4%；具有可解释优化结构的网络达到有竞争力但并非全面领先的性能。

**局限与边界：** 架构可解释不等于人类概念或因果机制可解释；率缩减仍是代理目标，不能替代真实通信比特与操作RD。

**可借鉴到语义流形率失真的方法：** 适合作为具有可检查矩阵谱与压缩目标的编码器对照；不建议将复杂新网络作为TIT论文主创新，核心应保持数学问题清晰。

**版本说明：** 此处为NeurIPS2023版本；不是标题带“Compression Is All There Is?”的JMLR扩展稿。

来源：[正式出版记录](https://papers.neurips.cc/paper_files/paper/2023/hash/1e118ba9ee76c20df728b42a35fb4704-Abstract-Conference.html)；[公开全文/作者稿](https://proceedings.neurips.cc/paper_files/paper/2023/file/1e118ba9ee76c20df728b42a35fb4704-Paper-Conference.pdf)；[PDF下载](https://proceedings.neurips.cc/paper_files/paper/2023/file/1e118ba9ee76c20df728b42a35fb4704-Paper-Conference.pdf)。

### M08. The Linear Representation Hypothesis and the Geometry of Large Language Models

**ICML 2024｜正式发表：主会论文｜概念几何与反事实**

作者：Kiho Park；Yo Joong Choe；Victor Veitch。

**背景与问题：** 线性表示假说混合了“方向表示”“线性探针”“加向量干预”三种不同含义，欧氏余弦也未必具有语义依据。

**创新方法：** 以反事实词对形式化输入/输出空间的概念表示，证明在假设下探针与干预的联系；定义尊重概念可分离性的因果内积，并用词向量协方差逆矩阵构造一个选择。

**数据与实验设置：** LLaMA-2 7B、32K词表；27类概念，包括BATS中的22类、4组语言转换及词频概念；另有法语/西语Wikipedia上下文和构造的干预提示。

**已报告结果与结论：** 多数反事实词对方向在适当内积下较随机词对更集中，可分离概念更接近正交；thing→part是线性假说的反例之一。探针与控制实验支持理论联系。

**局限与边界：** 结论依赖反事实概念、可分离性与线性结构假设；内积未由条件唯一确定，不能推广成“所有语义均为一维方向”。

**可借鉴到语义流形率失真的方法：** 设计语义保持/改变的反事实对及非目标概念保持指标；检查选用欧氏距离是否合理，为任务敏感几何提供对照。

**版本说明：** ICML2024正式版；arXiv首次为2023。

来源：[正式出版记录](https://proceedings.mlr.press/v235/park24c.html)；[公开全文/作者稿](https://arxiv.org/html/2311.03658v2)；[PDF下载](https://raw.githubusercontent.com/mlresearch/v235/main/assets/park24c/park24c.pdf)。

### M09. Diffusion Models Encode the Intrinsic Dimension of Data Manifolds

**ICML 2024｜正式发表：主会论文｜内在维数与法空间**

作者：Jan Pawel Stanczuk；Georgios Batzolis；Teo Deveney；Carola-Bibiane Schönlieb。

**背景与问题：** 高维数据虽可能集中在流形上，但真实有效维数难以估计；传统邻域方法在较高内在维数时容易失准。

**创新方法：** 证明低噪声score主要指向数据流形的法空间；在某个样本附近多次加噪并收集score向量，用SVD谱间隙估计法空间维数，进而估计局部内在维数。

**数据与实验设置：** 100维空间中的10/50维球面、spaghetti曲线、不同维数球面并集；10/20/100维的方块与高斯斑点合成图像流形；MNIST及不同维度自编码器验证。

**已报告结果与结论：** 一维spaghetti曲线估计为1，PPCA为98；100维高斯斑点图像流形估计为98，ID-NF为56.3。MNIST整体给出152，但真实维数未知，这个数是模型与尺度相关估计而非真值。

**局限与边界：** 小噪声理论要求流形支持及较准确score；有限噪声、曲率、数据非均匀性和训练误差都会影响谱间隙。原始图像内在维数不等于给定任务的语义有效维数。

**可借鉴到语义流形率失真的方法：** 把ID-diff作为独立的几何诊断；与任务后验Fisher秩、可控真实维数和经验率失真斜率交叉核对，避免由可视化直接指定k。

**版本说明：** 已核对最终主会PDF第5—9页及表1/2；arXiv早期标题为“Your diffusion model secretly knows the dimension of the data manifold”。

来源：[正式出版记录](https://proceedings.mlr.press/v235/stanczuk24a.html)；[公开全文/作者稿](https://raw.githubusercontent.com/mlresearch/v235/main/assets/stanczuk24a/stanczuk24a.pdf)；[PDF下载](https://raw.githubusercontent.com/mlresearch/v235/main/assets/stanczuk24a/stanczuk24a.pdf)。

### M10. The Geometry of Categorical and Hierarchical Concepts in Large Language Models

**ICLR 2025｜正式发表：主会论文｜多类别与层级几何**

作者：Kiho Park；Yo Joong Choe；Yibo Jiang；Victor Veitch。

**背景与问题：** 二元概念的单方向描述不能完整表达动物种类、语言等多类别概念以及它们的层次关系。

**创新方法：** 在概念几何框架中把多类别概念表示为多面体，并描述层级概念间的子空间/正交关系；几何比较使用与因果概念结构相容的内积。

**数据与实验设置：** 在Gemma和LLaMA-3模型上检验由WordNet构造的900余个层级概念，研究类别间对比以及上位/下位概念关系。

**已报告结果与结论：** 观察到与理论相符的多面体和层级几何结构，为超越单方向的概念分析提供依据；结果针对指定模型与概念体系。

**局限与边界：** 离散类别多面体、层次和光滑流形是不同对象；不能由这些实验直接假定所有语义空间都是固定维数的单个黎曼流形。

**可借鉴到语义流形率失真的方法：** 在实验中加入分层任务、离散分支与连续因子；必要时使用流形并集或分层空间，并把类别/图册索引比特计入码率。

**版本说明：** ICLR2025主会oral；不能将早期ICML2024 workshop版本当成ICML主会。

来源：[正式出版记录](https://proceedings.iclr.cc/paper_files/paper/2025/hash/be7430d22a4dae8516894e32f2fcc6db-Abstract-Conference.html)；[公开全文/作者稿](https://proceedings.iclr.cc/paper_files/paper/2025/file/be7430d22a4dae8516894e32f2fcc6db-Paper-Conference.pdf)；[PDF下载](https://proceedings.iclr.cc/paper_files/paper/2025/file/be7430d22a4dae8516894e32f2fcc6db-Paper-Conference.pdf)。

### M11. Not All Language Model Features Are One-Dimensionally Linear

**ICLR 2025｜正式发表：主会论文｜多维概念与因果验证**

作者：Joshua Engels；Eric J. Michaud；Isaac Liao；Wes Gurnee；Max Tegmark。

**背景与问题：** 稀疏自编码器等解释技术常把每个概念压到单一线性方向，但周期性变量可能天然需要多维结构。

**创新方法：** 提出不可约多维特征的定义，通过聚类稀疏自编码器特征发现圆形表示，并用模型内干预检验这些方向对周期任务的因果作用。

**数据与实验设置：** GPT-2与Mistral-7B的特征分析；星期/月以及相关模运算任务；Mistral-7B和Llama-3 8B的内部干预与电路验证。

**已报告结果与结论：** 发现星期/月的二维圆形特征，并证明相关表示参与任务计算；结论支持某些语义需要多维几何，而非仅仅二维投影看起来像圆。

**局限与边界：** 研究的是具体周期概念和模型，不能据此推断一般语言语义均为圆形、光滑或低曲率流形；稀疏特征恢复本身也有估计误差。

**可借鉴到语义流形率失真的方法：** 最合适的可控验证例之一：以圆/环面上的周期语义构造可实现后验族，比较流形内原型和自由后验均值，并用干预确认语义保持。

**版本说明：** ICLR2025正式标题含“One-Dimensionally”；早期预印本标题较短。

来源：[正式出版记录](https://proceedings.iclr.cc/paper_files/paper/2025/hash/d3221cdb27e49d9c1cd35ad254feccfe-Abstract-Conference.html)；[公开全文/作者稿](https://proceedings.iclr.cc/paper_files/paper/2025/file/d3221cdb27e49d9c1cd35ad254feccfe-Paper-Conference.pdf)；[PDF下载](https://proceedings.iclr.cc/paper_files/paper/2025/file/d3221cdb27e49d9c1cd35ad254feccfe-Paper-Conference.pdf)。

### M12. Geometry-Aware Generative Autoencoders for Warped Riemannian Metric Learning and Generative Modeling on Data Manifolds

**AISTATS 2025｜正式发表：主会论文｜学习支持与测地线**

作者：Xingzhi Sun；Danqi Liao；Kincaid MacDonald；Yanlei Zhang；Guillaume Huguet；Guy Wolf；Ian Adelstein；Tim G. J. Rudner；Smita Krishnaswamy。

**背景与问题：** 低维可视化能暴露数据几何，但往往无法在保持流形支持的同时生成新点、计算可靠路径或模拟群体变化。

**创新方法：** GAGA结合几何感知自编码器与扭曲拉回度量；离流形样本辅助训练用于惩罚偏离支持，进一步进行体积引导采样、测地线计算和群体插值。

**数据与实验设置：** 椭球、半球、鞍面、环面等合成流形及噪声/高维嵌入；27天Embryoid Body单细胞数据；CITE-seq与Multiome留一时间点预测，在50/100个PCA维度评价。

**已报告结果与结论：** 公开全文表4中，CITE-seq的50维1-Wasserstein距离为23.29±0.83，所列最强比较方法为34.16±0.04。合成测地线多数情形改善，但鞍面上Dijkstra优于GAGA。

**局限与边界：** 流形支持和几何由训练与正则化估计，密度不均衡与噪声会造成偏差；图像/任务语义并非本文主对象，实验优势不意味着普遍优于图算法。

**可借鉴到语义流形率失真的方法：** 借鉴离流形负样本、支持一致性评价及扭曲度量基线；与任务KL几何分别对照，区分“真实可生成”与“对任务相近”。

**版本说明：** AISTATS2025主会；作者名单采用PMLR正式书目，公开arXiv v4作者列表有差异，引用时以所用版本为准。

来源：[正式出版记录](https://proceedings.mlr.press/v258/sun25c.html)；[公开全文/作者稿](https://arxiv.org/html/2410.12779v4)；[PDF下载](https://raw.githubusercontent.com/mlresearch/v258/main/assets/sun25c/sun25c.pdf)。

### M13. Probability Density Geodesics in Image Diffusion Latent Space

**CVPR 2025｜正式发表：主会论文｜密度测地线对照**

作者：Qingtao Yu；Jaskirat Singh；Zhaoyuan Yang；Peter Henry Tu；Jing Zhang；Hongdong Li；Richard Hartley；Dylan Campbell。

**背景与问题：** 扩散图像插值中的欧氏或球面直线可能生成不自然的中间帧，需要借助模型分布约束变化路径。

**创新方法：** 构造依赖概率密度的共形度量，利用score写出测地线初值/边值问题；在预训练扩散模型上计算偏向高密度区域的路径，无需再训练生成模型。

**数据与实验设置：** Stable Diffusion 2.1 base；40个每段100帧、单因素变化的CLEVR视频；MorphBench 90对、AFHQ和CelebA-HQ各50对，以及20对网络图像；每对生成17帧序列。

**已报告结果与结论：** 简单可控视频轨迹与所定义测地线较一致，多项插值连贯性/自然度指标有竞争力；论文同时报告某些质量指标并不占优，不能概括为全部提升。

**局限与边界：** 高密度路径不等于真实支持的最短路径，也不必保留指定语义；远距离端点、局部最优及数值计算开销仍是问题。

**可借鉴到语义流形率失真的方法：** 重要对照：分别比较欧氏、密度、切空间及任务后验几何，验证收益到底来自偏爱常见样本，还是来自可实现性约束。

**版本说明：** CVPR2025主会；以CVF正式书目命名，公开作者稿作为可访问备份。

来源：[正式出版记录](https://openaccess.thecvf.com/content/CVPR2025/html/Yu_Probability_Density_Geodesics_in_Image_Diffusion_Latent_Space_CVPR_2025_paper.html)；[公开全文/作者稿](https://arxiv.org/html/2504.06675v2)；[PDF下载](https://openaccess.thecvf.com/content/CVPR2025/papers/Yu_Probability_Density_Geodesics_in_Image_Diffusion_Latent_Space_CVPR_2025_paper.pdf)。

### M14. Atlas-based Manifold Representations for Interpretable Riemannian Machine Learning

**AISTATS 2026｜正式发表：主会论文｜多图册与可解释几何**

作者：Ryan Allen Robinett；Sophia Madejski；Kyle Ruark；Samantha J. Riesenfeld；Lorenzo Orecchia。

**背景与问题：** 单个全局低维坐标可能无法保持非平凡拓扑，而常用可视化嵌入也不能直接支持黎曼运算。

**创新方法：** 用显式可微图册存储局部坐标、图间转换与局部几何；Atlas-Learn从点云拟合局部二次图，支持近似测地更新、平行移动及受限优化。

**数据与实验设置：** Grassmann流形上的Fréchet均值；具有Klein瓶结构的高对比度图像块及凸/凹边界分类；人类造血单细胞数据，30个PCA环境维度中学习5维、30图册模型，积分309条初始轨迹。

**已报告结果与结论：** 在给定案例保持拓扑和可解释路径；单细胞Atlas受限流比直接环境空间积分更接近观测数据。图像块分类用64个图册中心构造可解释边界。

**局限与边界：** 尚不支持维数随位置变化的广义流形；图册中心自适应、转换误差及近似测地优化的严格误差控制仍有不足。

**可借鉴到语义流形率失真的方法：** 最直接的实现途径：图册索引加局部坐标量化、图册内原型重构、切空间一致性检验；图册索引和转换误差必须计入通信性能。

**版本说明：** AISTATS2026正式主会，PMLR300:3961–3969；arXiv首次为2025。

来源：[正式出版记录](https://proceedings.mlr.press/v300/robinett26a.html)；[公开全文/作者稿](https://arxiv.org/html/2510.17772v1)；[PDF下载](https://raw.githubusercontent.com/mlresearch/v300/main/assets/robinett26a/robinett26a.pdf)。

### M15. The Information Geometry of Softmax: Probing and Steering

**ICML 2026｜正式发表：主会论文｜语义信息几何**

作者：Kiho Park；Todd Nief；Yo Joong Choe；Victor Veitch。

**背景与问题：** 用全局欧氏方向控制softmax模型容易改变非目标概念，需要将概率分布变化与表示空间几何严格联系。

**创新方法：** 将softmax写成指数族，连接KL与对数配分函数的Bregman散度、Fisher信息与输出向量协方差，并在均值/对偶坐标中进行干预；特定可分解条件下控制非目标分布变化。

**数据与实验设置：** Gemma-3 4B与C4上下文，测试性别、时态、英法语言等概念；MetaCLIP-2配合合成颜色/形状和COCO图像，评估目标概念概率、非目标KL、概率质量与排序变化。

**已报告结果与结论：** 对偶坐标干预在目标改变与非目标分布保持之间表现较好；建立了表示、探针、概率分布和信息几何之间的直接理论联系，是语义几何方向的重要既有研究。

**局限与边界：** 依赖概念线性/可分解假设，主要研究输出softmax；集中分布下协方差可能病态，需要数值正则。e-/m-联络测地线不能直接称作Levi-Civita最短测地线。

**可借鉴到语义流形率失真的方法：** 首要新颖性对照。使用任务后验KL和非目标损失诊断；研究真正的可实现后验集合与自由概率重构之间差距。均值坐标可平均并不保证概率mixture仍属于非饱和模型族。

**版本说明：** ICML2026正式主会，PMLR306:96977–97000；不是仅arXiv预印本。

来源：[正式出版记录](https://proceedings.mlr.press/v306/park26p.html)；[公开全文/作者稿](https://arxiv.org/html/2602.15293v2)；[PDF下载](https://raw.githubusercontent.com/mlresearch/v306/main/assets/park26p/park26p.pdf)。

### M16. Compositional Generalization Requires Linear, Orthogonal Representations in Vision Embedding Models

**ICML 2026｜正式发表：主会论文｜组合语义表示**

作者：Arnas Uselis；Andrea Dittadi；Seong Joon Oh。

**背景与问题：** 视觉基础模型在未见过的属性组合上泛化有限，需要说明什么表示几何有助于组合泛化。

**创新方法：** 通过可分解、可迁移、稳定等形式化要求，在指定学习设定下推导线性正交因子表示的必要性与维数约束；用条件平均及投影解释率测量模型对该结构的符合程度。

**数据与实验设置：** PUG-Animal、dSprites、MPI3D、ImageNet-AO；比较CLIP/OpenCLIP/MetaCLIP、SigLIP和DINO系列视觉表示，使用部分属性组合训练、独立组合测试。

**已报告结果与结论：** 训练良好的视觉表示具有一定因子化与正交结构，但仍存在较大未解释残差；理论必要条件是在给定定义和优化条件下成立，不是所有模型的无条件规律。

**局限与边界：** 理论要求较强，真实任务存在相关因素和交互作用；线性探针与有限数据的评估不能穷尽组合推理能力。

**可借鉴到语义流形率失真的方法：** 为语义流形实验提供有真实生成因子的公开数据与组合分布外测试；把因子表示残差和未见组合失真作为额外检验。

**版本说明：** ICML2026正式主会，PMLR306:123040–123107。

来源：[正式出版记录](https://proceedings.mlr.press/v306/uselis26a.html)；[公开全文/作者稿](https://arxiv.org/pdf/2602.24264)；[PDF下载](https://raw.githubusercontent.com/mlresearch/v306/main/assets/uselis26a/uselis26a.pdf)。

### M17. From Data Statistics to Feature Geometry: How Correlations Shape Superposition

**ICLR 2026｜正式发表：主会论文｜几何可解释性的反证检查**

作者：Lucas Prieto；Edward Stevinson；Melih Barsbey；Tolga Birdal；Pedro A. M. Mediano。

**背景与问题：** 概念在表示中形成簇或圆形，既可能反映计算结构，也可能仅由数据共现与训练约束造成。

**创新方法：** 提出Bag-of-Words Superposition受控实验；已知二元词特征的相关性在低秩瓶颈中产生建设性干涉，区分线性/非线性叠加以及presence-coding/value-coding。

**数据与实验设置：** WikiText-103构建10,000词词表、每样本20条记录，训练1,621,198例、验证180,133例；OpenWebText复验；另有循环相关合成特征、模运算和相对位置等受控任务。

**已报告结果与结论：** 紧瓶颈或权重衰减下，相关特征可形成语义簇和周期结构，非线性自编码器也可能利用线性叠加；几何图形本身不足以确定其语义或因果解释。

**局限与边界：** 主要是受控自编码器和小型任务分析，不能直接推广为所有LLM内部机制的充分描述。

**可借鉴到语义流形率失真的方法：** 必须纳入负对照：打乱概念共现、改变采样密度、控制瓶颈维数，同时观察流形几何和率失真曲线；避免将数据相关性误解释为语义结构。

**版本说明：** ICLR2026主会正式版本。

来源：[正式出版记录](https://proceedings.iclr.cc/paper_files/paper/2026/hash/85cae4ea32697bbdc2a48e7ca0d15988-Abstract-Conference.html)；[公开全文/作者稿](https://arxiv.org/html/2603.09972v1)；[PDF下载](https://proceedings.iclr.cc/paper_files/paper/2026/file/85cae4ea32697bbdc2a48e7ca0d15988-Paper-Conference.pdf)。

### M18. Be Tangential to Manifold: Discovering Riemannian Metric for Diffusion Models

**ECCV 2026｜正式发表：主会论文｜切空间与支持约束**

作者：Shinnosuke Saito；Takashi Matsubara。

**背景与问题：** 密度引导插值易偏向高概率区域，而扩散模型没有显式低维潜变量，支持流形及切方向难以直接提取。

**创新方法：** 从score的Jacobian谱分离切向/法向，构造无需再训练的各向异性黎曼度量；用它进行测地插值以及限制classifier-free guidance偏离流形。

**数据与实验设置：** 合成C形数据100,000点；Stable Diffusion 2.1；MorphBench、AFHQ、CelebA-HQ插值；DAVIS、HumanPexels、RealEstate10K视频中帧；30,000条COCO2014验证集文本进行生成评价。

**已报告结果与结论：** 作者报告插值感知自然度改善，并在保持文本图像对齐的同时提升生成质量；核心贡献是区分向高密度移动与沿数据流形切向移动。

**局限与边界：** 真实score与谱分离均由模型近似；噪声尺度、谱间隙、计算成本及局部支持误差影响效果；模型支持几何仍不自动等同于任务后验几何。

**可借鉴到语义流形率失真的方法：** 最新且相关的支持约束基线：比较切向保真、离流形偏移、任务KL和原型可实现性，以区别“语义重构合法”与“样本更常见”。

**版本说明：** ECCV2026主会，Springer正式章于2026-09-12上线；下载使用作者开放版本。

来源：[正式出版记录](https://link.springer.com/chapter/10.1007/978-3-032-37016-7_4)；[公开全文/作者稿](https://arxiv.org/html/2510.05509v3)；[PDF下载](https://arxiv.org/pdf/2510.05509)。

## 3. 对“可实现语义重构约束的率失真代价”的方法建议

### 3.1 应当明确区分的三个对象

- **数据支持流形**描述哪些输入或生成样本是合法的。M01、M12、M13、M14、M18主要研究这类对象。
- **表示几何**描述神经网络中的局部方向、概念关系、维数与拓扑。M06、M08—M11、M16、M17提供诊断；它未必与输入支持流形同构。
- **任务后验的可实现集合**描述接收端允许输出的语义对象对应哪些预测分布。若解码器必须输出原型、合法对象或指定生成族，可定义该集合为 \(\mathcal M=\{q(\cdot\mid \hat x):\hat x\in\widehat{\mathcal X}_{\rm valid}\}\)。这是拟议信息论问题真正需要的重构集合。

把这三者统称“语义流形”会使创新与结论不清楚。尤其是原始图像的内在维数，不能直接替代任务后验有效维数；类别层级也可能需要分层空间或流形并集。

### 3.2 一个能连接几何与操作问题的精确定义

以下恒等式是本次调研的数学整理，不是声称上述18篇已解决了拟议通信问题。令 \(q_X=q(\cdot\mid X)\)，编码器输出 \(U\)，自由重构的后验均值为 \(\bar q_U=\mathbb E[q_X\mid U]\)。在相关KL有限时，对任意只依赖 \(U\) 的重构分布 \(r_U\)，有

\[
\mathbb E D_{\rm KL}(q_X\|r_U)
=\mathbb E D_{\rm KL}(q_X\|\bar q_U)
+\mathbb E D_{\rm KL}(\bar q_U\|r_U).
\]

因此，**固定编码器**时，强制 \(r_U\in\mathcal M\) 的额外损失就是对 \(\bar q_U\) 的可实现性投影代价

\[
\mathbb E\inf_{r\in\mathcal M}D_{\rm KL}(\bar q_U\|r).
\]

这个分解给出可检验对象，但尚不是对两个最优率失真函数之差的闭式解：加入约束后最优编码器也会改变。后续定理需要处理该联合优化，以及何时额外代价为零、何时严格为正、局部曲率或统计模型结构怎样影响其量级。

不要把“对偶均值参数能够取平均”与“概率混合分布仍属于模型族”混为一谈。非饱和指数族的均值参数平均可以合法，却不一定使原概率mixture属于同一族。若可实现集合本身对mixture闭合，上述固定编码器缺口可消失；集合不闭合也不代表每一种源和每个编码器都一定产生严格缺口。

### 3.3 可以直接迁移的方法

| 研究环节 | 可借鉴论文 | 具体做法 | 必须另行验证 |
|---|---|---|---|
| 任务失真 | M04、M15 | 从任务后验KL拉回度量，检查Fisher谱与校准 | 退化秩、概率零点、模型失配以及KL方向 |
| 可实现集合 | M14、M12 | 图册拟合或显式生成族，图内量化并输出合法原型 | 图册覆盖、转换误差、样本支持误差 |
| 维数诊断 | M09、M06 | score法空间估计、Jacobian谱与切空间距离 | 尺度依赖和有限样本偏差；数据维数≠任务维数 |
| 反事实评估 | M02、M05、M08、M11 | 改变一个概念、保持其余概念，进行模型内或输入级干预 | 相关不等于因果；残差分支不能隐藏在解释之外 |
| 支持与密度对照 | M13、M18、M12 | 同时比较密度路径、切向路径和任务后验路径 | 更常见不等于更合法，更合法不等于任务更接近 |
| 泛化与负对照 | M16、M17 | 留出属性组合，打乱共现，改变数据采样密度 | 几何收益是否只是训练分布捷径 |
| 真实码率 | M03、M07 | 保留log-det代理为训练对照，独立做熵编码 | 代理编码率、互信息和实际比特不可互相替代 |

### 3.4 推荐的验证层次

1. **解析模型。** 在可明确写出后验族的圆、环面、曲线或低维指数子族上，比较自由后验重构与流形内重构；加入mixture闭合集合作为零缺口对照。此处优先采用可精确计算的离散概率任务，避免先把误差归咎于深度模型。
2. **可控生成因素。** 借鉴M11的周期语义和M16的dSprites/MPI3D等因素数据，区分任务相关因素与无关因素；报告任务KL、概念重构、可实现性和实际比特。
3. **可解释真实数据。** 用CUB概念或预训练表示上的PCBM构造可解释原型；比较自由概率输出、约束原型输出、欧氏量化、Fisher/图册量化和随机原型。若有残差通道，单独计入比特并测量其贡献。
4. **稳定性。** 改变任务集合、采样密度、概率校准、图册数量和噪声；检查观察到的代价是否仍来自合法重构约束。补充非目标概念KL与反事实一致性，避免只报告分类准确率。

## 4. 新颖性风险和投稿定位

Fisher拉回、KL局部二次展开、潜空间测地线、内在维数估计以及圆形概念均已有明确研究。单纯把这些方法用于一个语义通信网络，不足以支撑“首次从流形解释语义”的表述，也很难单凭实验成为TIT的主要贡献。

更清晰的切入点是：定义接收端语义对象必须可实现的操作约束，证明它相对自由预测后验重构的代价，并刻画零缺口条件、严格正缺口例子以及可计算界。几何算法用于验证理论，并量化学习几何和真实对象支持之间的误差。本文调研只能确认相关ML方法和近邻，不能单独证明信息论问题没有被其他传统约束重构、量化或率失真文献覆盖；最终新颖性仍需要结合信息论文献检索。

**还应知晓的两篇主会近邻（未计入上述18篇精选包）：**

- [The Spacetime of Diffusion Models: An Information Geometry Perspective，ICLR2026](https://proceedings.iclr.cc/paper_files/paper/2026/hash/ad2b3b0963114adddbd2c97d6e74abb9-Abstract-Conference.html)。正式摘要指出，确定性概率流ODE的朴素拉回会使测地线解码为环境空间直线；作者引入带噪声时间的随机后验族，构造Fisher几何和Diffusion Edit Distance。它提醒我们必须核查所用“流形度量”是否真正包含数据几何。
- [Score-based Pullback Riemannian Geometry: Extracting the Data Manifold Geometry using Anisotropic Flows，ICML2025](https://proceedings.mlr.press/v267/diepeveen25a.html)。这是score、生成流及数据驱动黎曼几何的直接近邻，也说明“以score学习流形几何”本身已不是空白。此处只列出近邻定位，未把摘要核查当作与18篇等深度的逐篇阅读。

## 5. 数据和版本核查说明

所有入选论文均有正式会议或出版机构记录及公开全文来源。定量结果是各论文特定实验配置中的数值，不构成跨论文公平排名；它们也不是本次新实验的结果。下载采用正式会议PDF优先、公开作者稿备用的策略。作者稿与最终版的作者名单、标题、页码或修订日期若有差异，需同时保留正式出版身份与实际下载版本的说明。

主会与预印本状态未确认或互相矛盾的条目没有混入精选18篇。例如“Origins of Representation Manifolds”检索到作者页面与审稿页面的不同状态，因此没有仅凭二手页面把它归为某一主会正式论文。


## 包内 PDF 快速链接

- **M01**：[Latent Space Oddity: on the Curvature of Deep Generative Models](../../../refs/ml_geometry/o_2018_ICLR_Latent_Space_Oddity_on_the_Curvature_of_Deep_Generative_Models.pdf)
- **M02**：[Concept Bottleneck Models](../../../refs/ml_geometry/o_2020_ICML_Concept_Bottleneck_Models.pdf)
- **M03**：[Learning Diverse and Discriminative Representations via the Principle of Maximal Coding Rate Reduction](../../../refs/ml_geometry/o_2020_NeurIPS_Learning_Diverse_and_Discriminative_Representations_via_the_Principle_of_Maximal_Coding_Rate_Reduction.pdf)
- **M04**：[Pulling back information geometry](../../../refs/ml_geometry/o_2022_AISTATS_Pulling_back_information_geometry.pdf)
- **M05**：[Post-hoc Concept Bottleneck Models](../../../refs/ml_geometry/o_2023_ICLR_Post_hoc_Concept_Bottleneck_Models.pdf)
- **M06**：[Understanding the Latent Space of Diffusion Models through the Lens of Riemannian Geometry](../../../refs/ml_geometry/o_2023_NeurIPS_Understanding_the_Latent_Space_of_Diffusion_Models_through_the_Lens_of_Riemannian_Geometry.pdf)
- **M07**：[White-Box Transformers via Sparse Rate Reduction](../../../refs/ml_geometry/o_2023_NeurIPS_White_Box_Transformers_via_Sparse_Rate_Reduction.pdf)
- **M08**：[The Linear Representation Hypothesis and the Geometry of Large Language Models](../../../refs/ml_geometry/o_2024_ICML_The_Linear_Representation_Hypothesis_and_the_Geometry_of_Large_Language_Models.pdf)
- **M09**：[Diffusion Models Encode the Intrinsic Dimension of Data Manifolds](../../../refs/ml_geometry/o_2024_ICML_Diffusion_Models_Encode_the_Intrinsic_Dimension_of_Data_Manifolds.pdf)
- **M10**：[The Geometry of Categorical and Hierarchical Concepts in Large Language Models](../../../refs/ml_geometry/o_2025_ICLR_The_Geometry_of_Categorical_and_Hierarchical_Concepts_in_Large_Language_Models.pdf)
- **M11**：[Not All Language Model Features Are One-Dimensionally Linear](../../../refs/ml_geometry/o_2025_ICLR_Not_All_Language_Model_Features_Are_One_Dimensionally_Linear.pdf)
- **M12**：[Geometry-Aware Generative Autoencoders for Warped Riemannian Metric Learning and Generative Modeling on Data Manifolds](../../../refs/ml_geometry/o_2025_AISTATS_Geometry_Aware_Generative_Autoencoders_for_Warped_Riemannian_Metric_Learning_and_Generative_Modeling_on_Data_Manifolds.pdf)
- **M13**：[Probability Density Geodesics in Image Diffusion Latent Space](../../../refs/ml_geometry/o_2025_CVPR_Probability_Density_Geodesics_in_Image_Diffusion_Latent_Space.pdf)
- **M14**：[Atlas-based Manifold Representations for Interpretable Riemannian Machine Learning](../../../refs/ml_geometry/o_2026_AISTATS_Atlas_based_Manifold_Representations_for_Interpretable_Riemannian_Machine_Learning.pdf)
- **M15**：[The Information Geometry of Softmax: Probing and Steering](../../../refs/ml_geometry/o_2026_ICML_The_Information_Geometry_of_Softmax_Probing_and_Steering.pdf)
- **M16**：[Compositional Generalization Requires Linear, Orthogonal Representations in Vision Embedding Models](../../../refs/ml_geometry/o_2026_ICML_Compositional_Generalization_Requires_Linear_Orthogonal_Representations_in_Vision_Embedding_Models.pdf)
- **M17**：[From Data Statistics to Feature Geometry: How Correlations Shape Superposition](../../../refs/ml_geometry/o_2026_ICLR_From_Data_Statistics_to_Feature_Geometry_How_Correlations_Shape_Superposition.pdf)
- **M18**：[Be Tangential to Manifold: Discovering Riemannian Metric for Diffusion Models](../../../refs/ml_geometry/2026_ECCV_Be_Tangential_to_Manifold_Discovering_Riemannian_Metric_for_Diffusion_Models.pdf)
- **M19（补录）**：[A Geometric Lens on Physics-Aligned Data Compression](../../../refs/ml_geometry/2026_ICML_A_Geometric_Lens_on_Physics_Aligned_Data_Compression.pdf)
- **M20（补录）**：[Geometric Rate–Distortion Invariance for Domain Generalization](../../../refs/ml_geometry/2026_ICML_Geometric_Rate_Distortion_Invariance_for_Domain_Generalization.pdf)

## M19 · 补录：A Geometric Lens on Physics-Aligned Data Compression

**2026-10-31 补录｜ICML 2026 正式主会论文｜PMLR 306:122906–122926**

- **作者：** Aleix Segui Ugalde；Wesley Armour（按出版页写法；PDF 第一作者署名为 Aleix Seguí）。
- **背景与方法：** 研究固定码率下物理可观测量保真与常规重构失真的冲突。以熵模型、物理可观测量和信号失真诱导的潜空间局部度量刻画各向异性误差分配，提出局部切空间率失真分析与主特征子空间重叠的对齐诊断。
- **实验范围：** 计算流体动力学、宇宙学模拟及大脑皮层体积电子显微镜测量。
- **课题关系与边界：** 是任务/物理量保真、潜空间几何和压缩的直接近邻。其局部压缩噪声配置问题不应直接等同于全局操作性语义率失真或可实现后验流形重构的代价；后续新颖性分析需明确区别。
- **核查深度：** 已核对正式出版页、PDF 摘要与引言；未逐项审计证明、实验数值或复现代码。作为原18篇精选之外的补录，不宣称完成同等深度调研。
- **代码状态：** 正式 PDF 和作者发布页未找到该文配套代码链接，暂不加 `o_`。
- **来源：** [正式出版页](https://proceedings.mlr.press/v306/ugalde26a.html)；[作者发布页](https://aleixsegui.dev/research)。
- **包内全文：** [PDF](../../../refs/ml_geometry/2026_ICML_A_Geometric_Lens_on_Physics_Aligned_Data_Compression.pdf)，21 页。

## M20 · 补录：Geometric Rate–Distortion Invariance for Domain Generalization

**2026-10-31 补录｜ICML 2026 正式主会论文｜PMLR 306:77298–77343**

- **作者：** Tong Liu；Sen Liang；Shuo Bai。
- **背景与方法：** 面向未见目标域的域泛化，RDI 将类条件表示建模为 Grassmann 流形上的低维子空间，联合优化跨域子空间对齐的几何失真与谱—体积复杂度的容量正则化率项，以减少分布偏移下的判别几何退化。
- **理论与实验：** 正式摘要及引言报告了有界偏移下的有限样本稳定性保证，并在 DomainBed 上与域泛化基线比较、消融对齐和复杂度控制；未在此补录中逐项核实数值或审计证明。
- **课题关系与边界：** 是“流形几何＋率失真式目标”的直接近邻。该文的率项用于表示复杂度/容量正则化，不应未经操作性编码论证就等同于实际通信码率；域泛化子空间对齐也不同于接收端可实现后验重构约束。
- **核查深度：** 已核对正式出版页、PDF 摘要与引言，并检查全篇代码链接。作为原18篇精选之外的补录，不宣称完成同等深度调研或代码复现。
- **代码状态：** 正式 PDF 与本次检索结果未找到作者配套代码；OpenReview 页面要求浏览器验证，未进一步访问。因此暂不带 `o_`，不等于已证明作者从未开源。
- **来源：** [正式出版页](https://proceedings.mlr.press/v306/liu26cc.html)；[OpenReview](https://openreview.net/forum?id=96ma1Gyrnj)。
- **包内全文：** [PDF](../../../refs/ml_geometry/2026_ICML_Geometric_Rate_Distortion_Invariance_for_Domain_Generalization.pdf)，46 页。


