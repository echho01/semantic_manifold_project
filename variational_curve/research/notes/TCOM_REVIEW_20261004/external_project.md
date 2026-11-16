# /Users/gaofan.86/Documents/semantic-manifold 只读审计与 TCOM 借鉴建议

日期：2026-10-04。审计范围为该目录全部 13 个研究文本/脚本（common、两篇 main.tex 与证明笔记、3 个 synthetic 脚本、STATUS、code/README、related-work、references.bib），并对照当前工作区曲面正文、通信扩展、GPU 计划和已保存结果。外部编译日志/PDF属于生成物，未重新编译；外部 refs 的两份论文与当前 paper 完全相同（SHA-256 已比对），由本次文献评审统一处理。本报告没有改外部项目或当前仓库，没有安装依赖。

## 1. 总体判断

**该外部项目可作“问题组织与机制实验清单”，不能作当前论文的成熟理论底座或独立验证。当前仓库在固定任务度量、编码计费、可达核与理论边界上明显更严谨，应保留当前框架，逐项吸收外部项目中经过修正的内容。**

外部曲率稿是骨架：`/Users/gaofan.86/Documents/semantic-manifold/paper1-curvature/main.tex:78` 起各主定理环境仅有注释，正式陈述和证明尚未填入；`:54`–`:56` 摘要却声称 CLIP/BERT 已验证，实际 `/Users/gaofan.86/Documents/semantic-manifold/code/README.md:22` 明确真实嵌入实验待做。两篇稿 `main.tex:24` 写 companion paper “submitted”，该文件自身不足以证明投稿状态，整合时不得直接沿用。STATUS 中的“✅证明闭合/novelty 成立”应视作作者工作笔记，必须独立检验。

当前关键对照是：

- `/Users/gaofan.86/Seafile/博四/语义通信/语义流形曲面拓扑/variational-semantic-manifold/variational_geometry.tex:25` 固定任务度量；`:57` 表示不变性；`:81` 率上下界；`:109` 指定可达核的展开；`:182` 真实码率计费；`:198` 似然匹配接收器。
- `/Users/gaofan.86/Seafile/博四/语义通信/语义流形曲面拓扑/variational-semantic-manifold/sections/communication_geometry.tex:32` 保留 density-only 强负对照；`:38` 明示真实数据实验待做。
- `/Users/gaofan.86/Seafile/博四/语义通信/语义流形曲面拓扑/variational-semantic-manifold/notes/RESEARCH_AND_GPU_PLAN.md:57`–`:89` 已有比外部 CLIP/BERT 泛化计划更完整的实际 bit、信道使用、功率、CSI、种子、停止条件。

## 2. 借鉴分级表

| 外部对象 | 可借鉴的部分 | 判断与修正 | 对当前 TCOM 稿的价值 |
|---|---|---|---|
| common/framework.tex 的 reach、维数、曲率、源密度 | 建模假设清单、内蕴/外显失真区分 | 可引用经典定义；真实嵌入是否满足必须测；诱导欧氏度量不能自动代表任务语义 | 为局部图册/接收器补充有效域和可检验假设 |
| T1.1 曲率修正 RD | 球面/平环面可解源对照 | 原目标公式不可靠；fixed-rate 与 Shannon 混用，二阶系数未证明且相互矛盾 | 保留当前 Gibbs-kernel 与 exact homogeneous specialization；不要替换为该“普适展开” |
| T1.2 弦弧关系 | reach 控制局部失真换算 | 目标结论可由更好的现有定理修复，当前证明不成立 | 写成引用引理，用于任务 metric/编码 chord 的误差预算 |
| T1.3 码点密度变分 | 固定码本大小下 p^{d/(d+2)} 分配 | 首项是经典 Zador/Gersho；熵约束版本错；曲率系数依赖 λ 时漏导数 | 作为量化/密度基线；原创需来自限制性几何族、通信资源和可计算近似差 |
| T1.4 流形增益 | 小噪声切向投影、条件误差/跨分支错误分开 | d/D 是后处理去噪维数比，非无条件通信 SNR 增益；cliff 只能条件化定义 | 补机制图，真实主图仍需任务错误—SNR、bit 和功率 |
| T1.5 task 与 metric mismatch | 任务 Jacobian、退化方向、metric 等价 | max rank 不足以保证规则商流形；失配需相对谱界/固定解释 | 当前 Fisher+η metric 更稳健；可补尺度/度量失配消融 |
| Paper 2 思路 | 结构任务、PH、feature margin、噪声证书 | 主要定理多处不成立，不应作为曲面稿的附加理论 | 仅取反例/结构实验构思，严格使用当前任务参考证书 |
| common/notation、proofs、synthetic、STATUS 分工 | 一个结论对应假设、证明、脚本、原始输出 | 可直接借项目组织，取消未证结论的“通过”标记 | 给当前实验增加“proved / computed / pending / falsified”证据台账 |

## 3. 曲率 RD：不能回迁原公式

外部 `/Users/gaofan.86/Documents/semantic-manifold/paper1-curvature/proofs/T1-proof-notes.md:28`–`:32` 用了

\[
R(D)=\frac d2\log_2(1/D)+h_g(p)+c_d+\frac{C_d\bar S_p}{\ln2}D+o(D),
\qquad C_d=\frac{d}{12(d+2)},
\]

同时称 fixed-rate operational 版本与 Shannon RDF “首项+常数项一致”。这不是正确的理论区分。

### 3.1 固定大小码本与 Shannon RDF 的常数不同

设 N 个码点、平方测地失真，令 Q_d 为欧氏单位立方体最优量化的平方误差常数。在紧致规则流形，既有 Zador 首项为

\[
D_N\sim Q_d N^{-2/d}\left(\int p^{d/(d+2)}\,dV_g\right)^{(d+2)/d}.
\]

因此 R=log2 N 对应

\[
R_{\rm FR}(D)=\frac d2\log_2\frac{Q_d}{D}
+\frac{d+2}{2}\log_2\int p^{d/(d+2)}dV_g+o(1).
\]

密度常数是 Rényi 型项，不是 h_g(p)。源为均匀分布时它简化为 log2 Vol + (d/2)log2 Q_d。Shannon 小失真主项则是

\[
R_g(D)=h_g(p)-\frac d2\log\frac{2\pi eD}{d}+o(1)
\]

（此式为 nats；再除 ln2）。两者连均匀源一般也存在量化形状/熵约束常数差；例如均匀一维圆周 Q₁=1/12，固定码本量化常数比 Shannon 主项多 ½ln(2πe/12)≈0.17649 nats≈0.25461 bits。首项斜率相同不足以推出相同常数，更不能把一个码本实验当 Shannon 二阶项证明。

外部笔记 `:46`、`:124` 所称 Graf–Luschgy Thm.6.2 “流形量化系数”也不精确：该经典结果为欧氏量化首项；流形推广应查 Kloeckner/Iacobelli 等原文。已核实 [Iacobelli, Theorem 1.4 / Corollary 1.5](https://arxiv.org/pdf/1412.4026)：给流形 Zador 首项和最优码点分布，不给该笔记所声称的普适线性曲率 RD 系数。

### 3.2 原稿球面系数自相矛盾

原 `:30` 若 d=2、S=2/ρ²，所得系数为 D/(12ρ² ln2)。但 `:51` 又写 D/(24ρ² ln2)，自差 2 倍；相对于当前均匀球面 **Shannon RDF** 的真系数 D/(6ρ² ln2)，分别差 2 倍和 4 倍。这些并非可用“符号/单位不同”消解的同一公式。

当前 `/Users/gaofan.86/Seafile/博四/语义通信/语义流形曲面拓扑/variational-semantic-manifold/variational_geometry.tex:112` 的 **指定 forward kernel** 展开为

\[
I_h(s_D)=B_g(D)+\Gamma_F+\frac{aD}{6m}(\bar S_h+3\mathcal J_h)+O(D^2),
\]

并在 `:115` 明确齐性均匀源精确 RDF 的曲率系数 S/(6m)。单位 S²：S=2、m=2，因此为 1/6 nats/D。非均匀源的 Fisher density 项只属于所选 forward kernel，正文 `:140` 正确地没有把它宣称为 unrestricted RDF 真二阶系数。应保留这一限定。

本次外部项目审计用现有 scipy 对球面精确径向积分再做一次轻量核查（nats）：s=100、1000、10000 时 `(R-B)/D` 分别为 0.1667778042、0.1666777780、0.1666677778，趋于 1/6。本次评审另已重新执行当前 `validate_implementation.py`，球面系数 0.1666703704，有限核梯度误差 1.73e-11，所有断言通过。该验证与上述独立积分是两项不同的数值检查。

### 3.3 球体积不能直接决定最优量化二阶系数

外部 `:38`–`:47` 从 geodesic ball 的体积展开直接推精确编码上下界。ε-net 的 covering 常数、Voronoi cell 形状、密度导数、曲率导致的 cell 变形和全局几何缺陷均可能贡献二阶项；G–L 首项不填补这一缺口。一般 C¹ 密度也不足以支撑写成 O(r⁴) 的带 Laplacian 展开。正确用途是容量/覆盖的阶界或明确给定 cell 模型的启发式；不是“只差一般 p 常数、常曲率已闭合”。

## 4. 密度分配与失真准则：可用但必须修正

### 4.1 T1.3 的 λ 定义要固定

外部 `/Users/gaofan.86/Documents/semantic-manifold/paper1-curvature/proofs/T1-proof-notes.md:58` 将 λ 定义为每单位体积码点数量；`:63` 固定约束 ∫λ=N。若 m(x) 为给定、正且不依赖 λ，则

\[
\min_{\lambda>0}\int p(x)m(x)\lambda^{-2/d}dV,
\quad\int\lambda dV=N
\quad\Rightarrow\quad
\lambda^*(x)=N\frac{(p m)^{d/(d+2)}}{\int(p m)^{d/(d+2)}dV}.
\]

这是严格凸可解的经典首阶分配，值得用于 baseline。原 `:58` 的 mκ 本身含 λ^{-2/d}，则不能在变分时把 mκ 当常数：D=∫p[Qλ^{-2/d}+aSλ^{-4/d}+…] 时 derivative 还有第二项；`:64` 的封闭比例因此只可作首阶/固定 m 模型。

原 `:66` 声称 entropy-constrained λ∝p 错。该定义下 H(K)≈h(p)+∫p logλ，变分

\[
-\frac2d pm\lambda^{-2/d-1}+\mu p/\lambda=0
\]

给 λ∝m^{d/2}；m 常数时码点密度近似 **均匀**，cell 大小近似相同。λ∝p 指等概率 cell 之类不同目标，不等价于熵约束平方误差最优。若定义 λ 为归一化分布或 cell probability，又必须重新写约束。

将其叫“manifold reverse water-filling”只能作为类比，不能当新的编码定理；经典 Gaussian reverse water-filling 有独立 spectral mode 与零分配阈值，连续码点密度 Euler–Lagrange 比例本身没有该结构。

### 4.2 T1.2 的证明错误可引用既有 reach 定理修复

外部 `:10`–`:15` 想证 d_E≤d_g≤d_E+2d_E³/(3τ²)。其引用 NSW 的真实界是

\[
d_g\le\tau[1-\sqrt{1-2d_E/\tau}],\quad d_E\le\tau/2.
\]

但右端 Taylor 为 d_E+d_E²/(2τ)+d_E³/(2τ²)+…，有二次项，不能直接推出三次余项。已对照 [NSW 原文 Proposition 6.3](https://math.uchicago.edu/~shmuel/NSW1.pdf)，第16页。可以改引 [Boissonnat–Lieutier–Wintraecken, SoCG2018 Theorem 1 / Lemma 5](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SoCG.2018.10)：正 reach 下

\[
d_E\le d_g\le2\tau\arcsin(d_E/(2\tau))
=d_E+d_E^3/(24\tau^2)+O(d_E^5/\tau^4),\quad d_E<2\tau.
\]

这一修复为既有定理的引用，不能当原创。原 `:19` 还漏了尺度：D(1+O(D/τ²)) 对 RDF 的变化应为 O(D/τ²)，因为 δD=O(D²/τ²) 而 |R'(D)|≈d/(2D)。仅靠平均平方误差 D→0 也不能推出 E e⁴=O(D²)，需要尾部/最大 cell 半径或具体 kernel 控制；原证据链未给。

## 5. 信道增益、outage 与有限块长

### 5.1 小噪声切向投影是可靠机制，收益口径不能扩大

外部 `/Users/gaofan.86/Documents/semantic-manifold/paper1-curvature/proofs/T1-proof-notes.md:77`–`:84` 若 N∼N(0,σ_tot² I_D/D)，并在有效 tubular neighborhood 中投影，则一阶

\[
\pi_M(s+N)=s+P_{T_sM}N+O(\|N\|^2/\tau),
\quad E\|P_TN\|^2=d\sigma_{\rm tot}^2/D.
\]

可以做解析机制图。条件为真实输入流形、投影/解码器可知、局部噪声、协方差约定和语义 metric 不变。各向异性一般式应是 tr(G_task J_rec Σ J_recᵀ)，不只是 d/D；当前 `variational_geometry.tex:213` 的 local sensitivity allocation 更贴通信实现。

但 10log10(D/d) 不是通用通信 SNR、容量或端到端节能增益。增加完全无用的 ambient 坐标并给它们加噪，可以任意增大 raw MSE/投影 MSE 比值；没有改变任务所需信息量。必须同样 source、channel uses、实际 power、CSI、量化和 decoding budget 比较。

外部 `:83` 声称普适正二阶修正 `1+σ_ax²||II||²_HS/d` 也不应直接用。法向 shape operator 的逆只给部分项，还需切向扰动、内蕴距离和投影的高阶贡献。一个可检验反例推导为：S^m(r) 上各向同性 ambient noise、σ_ax²=t，测地平方误差

\[
E d_g^2=m t+\frac{m(5-2m)}{3r^2}t^2+O(t^3).
\]

该式由 `r² arctan²(||N_T||/(r+N_n))` 的局部展开得到（远端尾部指数小）；m=2 的相对二阶系数为 t/(3r²)，原式为 t/r²；m≥3 甚至可为负。此处属于本审计解析反例，没有将它当发表定理或数值复现结论。保守建议只留一阶式及阶余项，二阶另证。

### 5.2 t14 已重跑，能复现表格但不是论文性能验证

已在外部现有 `.venv` 执行 `/Users/gaofan.86/Documents/semantic-manifold/code/synthetic/t14_outage_threshold.py`。结果与 STATUS 表基本一致：σ_tot=.05 时 gain=15.06dB，σ_tot=1 时 gain=14.96dB、P_out=.0465。源为两条同心圆、D_amb=32、Δ=.6，30 个坐标仅噪声；`script:52`–`:56` 的 geoMSE **只在正确分支条件下计算**，没有计入跨分支概念错误。

该脚本 `:4` 把 reach 设为 ρ，但两圆并的 reach 实为 min(ρ,Δ/2)=.3。投影 `:38`–`:40` 只依赖前两维；出平面噪声不影响最近分支，因此总噪声范数 bound 可以极其保守。第一阶径向判错尺度由 σ_ax，而非二维噪声总幅度 sqrt(2)σ_ax 决定。

此例精确错误率可用 Rice 半径分布：

\[
P_e=\tfrac12[1-F_{\rm Rice}(1.3;1,\sigma_{ax})
+F_{\rm Rice}(1.3;1.6,\sigma_{ax})].
\]

本次轻量 scipy 积分/分布计算得 σ_tot=.3/.5/.7/1 的 Pe 分别约 7.87e-9、.00035197、.00784912、.04600768，和 Monte Carlo 是一个模型内的重复核对。σ>0 时跨支概率一般正且光滑；有限样本报 0 不等于概率严格 0，没有真实数学“不连续跳变”。

外部 `T1 notes:87` 将离正确分支距离超过 reach 或离开 tube 当作必然投影跳支，也不成立：外部 tube 可仍有唯一、正确最近点。reach 是足够安全域，不是每个源点实际 task decision margin，更不能无条件宣称 noise above reach ⇒ semantic cliff。

### 5.3 原端到端“三项相加”是模型，尚非等式定理

外部 `T1 notes:91` 写 D_e2e=D_quant+(d/D)σ²(1+…)+d_cross P_out。量化与噪声的 cross term、curved receiver 的 bias、条件正确分支误差、wrong branch 内连续误差都需处理。合理写法为

\[
E d_{sem}=(1-P_e)E[d_g^2\mid correct]+P_eE[d_{cross}\mid wrong],
\]

先定义 wrong branch cost，再对第一项作条件上界或局部 approximation。如果设跨分支 +∞，任何 Pe>0 都导致平均失真∞，就不能再使用有限 RD 和该相加式；要使用固定 finite penalty 或独立 excess-distortion/concept-error constraint。

### 5.4 外部项目没有有限块长理论，TCOM 可补但需区分问题

common `:52`–`:55` 只有实 AWGN 与 per-block sigma scaling，无 power budget、source block n、channel block k、codebook reliability、CSI/导频、BLER、dispersion。故不能把上述 geometric branch error 叫 finite-block channel outage。

可以将当前严格率夹逼接到标准准静态信道：L_g(D)≤R_g(D)≤I_h(s)，v=n/k，C(H) 与 R 均用同单位，则容量 outage 满足

\[
\Pr[C(H)<vL_g(D)]\le\Pr[C(H)<vR_g(D)]
\le\Pr[C(H)<vI_h(s)].
\]

对单位平均 Rayleigh gain、固定 P/N0=γ、complex SISO C=log2(1+γ|H|²)，中间式为

\[
P_{out}^{\infty}=1-\exp[-(2^{vR_g(D)}-1)/\gamma].
\]

这只计渐近容量中断；有限块长 excess-distortion 应接 [Kostina–Verdú lossy JSCC](https://arxiv.org/abs/1209.1317) 的 joint source/channel dispersion，典型 regular memoryless 设置的近似为

\[
kC-nR(D)\simeq\sqrt{kV_c+nV_s(D)}Q^{-1}(\epsilon).
\]

必须验证 d-tilted information、有限 dispersion、分离/联合编码假设；不能用首项 R(D)+curvature term 就宣称 finite block reliability。若把 conditional-on-H 正态近似再平均，须明确这是工程近似，不能无说明当整个 quasi-static fading channel 的新定理。CLIP/feature reconstruction errors 是否真为 isotropic Gaussian 也需 residual covariance 校准。

## 6. 任务 metric、知识失配：保留当前校准模型

外部 `/Users/gaofan.86/Documents/semantic-manifold/common/framework.tex:96` 取 d_f=max rank df；`T1-proof-notes.md:99`–`:104` 由 f* h 的半度量推 smooth quotient N 和其 curvature RDF。这需要 **constant rank、规则 fiber、适当 quotient、输出推前密度正则**；max rank 单独不够，奇异 rank strata 的正概率质量甚至会改变 RD dimension 的系数。O'Neill curvature identities 需要 Riemannian submersion，不是任意任务 smooth map 即可调用。f* h 位于原 manifold 且可退化，不能把它当推前源 N 上的非退化 volume measure 使用。

分类 f:M→finite labels 更不能一般光滑且非恒定地定义在 connected M 上：离散值连续映射在每个连通分量内常值，但真实 class decision boundary 可以穿过 connected feature manifold。将 β0 当 class alphabet 必须另立 task 条件。

当前 `variational_geometry.tex:28`–`:41` 的 frozen posterior Fisher metric+η positive regularizer，给 TV/regret/IB loss 联系，任务尺度固定；这是比外部“rank 越小越省率”更可操作的基础。要优先精读 task calibration 的文献支撑及是否真满足 sufficiency，保留直接 label/posterior coding baseline。

metric mismatch 外部 `T1 notes:106`–`:108` 只用绝对 coordinate operator norm ε_g，缺少参考 metric 和 uniform ellipticity；也没有证明无限码率能克服 task label/decoder mapping 失配。可用的严格限定是

\[
(1-\epsilon)g_T\preceq g_R\preceq(1+\epsilon)g_T
\Rightarrow (1-\epsilon)d_T^2\le d_R^2\le(1+\epsilon)d_T^2.
\]

固定相同解释、相同点集/码本、仅 metric 近邻失配下可推 multiplicative distortion，随 cell error 趋零。不能扩成“所有 KB metric mismatch 无底板”。当前 `communication_geometry.tex:18`–`:24` 的 conditional posterior TV/KL risk bound 已正确指出 marginal matching 不够；现有 circle label rotation counterexample 更有说服力。

## 7. 拓扑稿仅作为反例清单，不合并其未证定理

这部分已读取所有源，但当前用户聚焦曲面 TCOM，建议不把它扩成额外论文主线。

1. `/Users/gaofan.86/Documents/semantic-manifold/paper2-topology/proofs/T2-proof-notes.md:11`–`:17` 称 bottleneck 对同伦等价空间为0并 canonical unique，错误。PH 保的是 filtered topology/metric scale；不同半径圆同伦等价，offset H1 bar death 不同，d_B>0。`:66` 同伦等价 KB ⇒ 底板0 同样错。取0失真或某些缩放的稳定 pseudometric也满足不少弱公理，不能直接推出唯一 bottleneck metric。
2. `:20` 的概念无差错率下界 H(comp(S)) 正确（指定 source task）；主稿 `/Users/gaofan.86/Documents/semantic-manifold/paper2-topology/main.tex:39` 却写 R≥log2 β0，只有 equiprobable labels 才对。一般 Fano 应为 R≥H(B)−h2(Pe)−Pe log2(Nc−1)。
3. `T2 notes:27`–`:31` 的 Gaussian max norm σsqrt(2log n) 漏 ambient dimension 与 failure probability；正确典型界 σ[√D+√(2log(n/α))]。Vietoris–Rips 以边长为 filtration 时由 pairwise perturbation 得常数2ε，而不是任意 convention 都1ε；offset filtration 才对应 Hausdorff 1-Lipschitz。可引用 [Chazal–de Silva–Oudot](https://arxiv.org/abs/1207.3885)，不要不分对象直接引 scalar function stability。
4. `:36` 称 γ_min≤2reach 恒成立。framework 任意 sublevel task filtration 可乘尺度常数而不改变 manifold reach，立即反例；即便限定距离 filtration 也不能从“颈部直觉”证明 universal relation。bar gap 是过滤尺度，reach 是外在几何量，必须统一对象和单位。所谓 topology universally more fragile 不宜用。
5. `:43`–`:51` 用 N=2^R 点云覆盖阈值 d log(V^{1/d}/τ)，然后 converse 给 V/τ^d bits，**阶不匹配**，不能叫 order optimal。更关键一条 source symbol 的码本大小不意味着一段接收样本实际覆盖所有 N 个码点；需 source sample count、density lower bound、sampling confidence。NSW noisy Proposition7.1 实际 requires tube and r-density，原文 noise 模型是 bounded tube，不是高斯无界 raw distribution 可直接套。笔记“KB 已知可0 bit”反例提醒可借鉴，须分清传样本 vs 同步未知 manifold。
6. `:61`–`:64` 的 mismatch floor 只有在强制输出diagram距离固定 receiver KB≤ε(R)才由 triangle 成立；任意可通信 decoder 可以发送/重建 transmitter point cloud，何来固定KB diagram底板？所谓子结构/形变没有自动给所需 closeness。若 output admissible set A_R 已固定，严格 floor 是 inf_{D∈A_R}d_B(D,D_T)，而非必等 d_B(D_R,D_T)。diagram metric也不能独自保证 labels/task意思相同。
7. `:79`–`:82` 任务 f:M→R 时直接编码 f(S) 即可；Reeb edge index 仅当任务需要辨认不同 level-set 连通分支才收费。homology class在ker f_*也不意味着可以无误差丢掉某些可测 task values。不能仅函子性推“safe semantic feature deletion”。
8. `/Users/gaofan.86/Documents/semantic-manifold/code/synthetic/t23_dumbbell_cliff.py:41` 实际是两条分离圆，没有连接颈部，不是 connected dumbbell；`:63` 每个噪声级重采 clean points，混合 sampling与channel effect；未计算宣称的 bottleneck distance；significance threshold=.15 是选择条件而非γ_min 本身。false-positive bars、merge 和 measured task decision必须分开。未重跑该脚本，源审计不能当运行结论。

## 8. 实验与证据审计

### 8.1 t11 不能证明曲率二阶项

`/Users/gaofan.86/Documents/semantic-manifold/code/synthetic/t11_sphere_rd.py:18`–`:21` 比球面 Fibonacci 与平环面 square lattice（不同 cell class）；`:79`–`:90` 用两参数 a+b/N 拟合。常数 offset 可由不同 cell shape解释，但 N^{-1} 项也可能来自 lattice defects、Fibonacci离最优差、全局 unavoidable defects，不能只归scalar curvature。`:59` 实际环面码本 N_eff=round(sqrt N)^2，不是所有点都同N。回归没有标准误/种子/置信区间，更无均方残差或 alternative fits。

`:93`–`:100` 的 PASS 是固定打印的文字，没有 assert、容差判据或 exit failure。STATUS的“通过”只是作者宣告。M=200k，N=16384 时预期 b/N≈3.4e-5，相对 Monte Carlo noise 可远大于该项；chunk=20000 的 cos matrix 超2.6GB，甚至可能同时分配额外矩阵，不能称轻量充分核验。本审计没有运行该大脚本。

替代：沿用当前 exact spherical radial RDF（明确是 analytic source）、再加同一source/fixed task cost、同 lattice class、优化码本、多个 N/ρ/seeds、精确数值积分/CI 和 residual模型对照。若研究真实 quantizer 二阶項，要明确为该 quantizer、不要当 exact RDF。

### 8.2 当前结果不能伪装 GPU 性能，也不能视外部项目为 independent confirmation

当前 `/Users/gaofan.86/Seafile/博四/语义通信/语义流形曲面拓扑/variational-semantic-manifold/results/communication_final/summary.json:3` 标 CPU executed/GPU pending，`:5` 为 test-channel mutual information；`:9`–`:18` 的 geometry-only .12% vs density解释全部joint收益；`:56` 指定 prototype synthetic、完美ZF和fixed gain，无real image/random fading。这些边界应该完整保留。

external refs 两份与current论文文件完全相同：

- Riemannian Deep Learning SHA-256：201bff85529c0a2351aea2d973e2d499a1c24c30c122557476a6a491d380196c。
- Riemannian Metric Matching SHA-256：ba0b7698a11ead65adb3357261b59b0168c315b7cb97c5b6bcabaf782cee2f44。

它们不是额外两条独立文献。两个项目的 sphere/投影/PH toy scenes 对同一基础建模的重复不是独立 real-data validation；本次 t14+Rice 只是相同统计模型的 MC与解析一致性，不能升格 novel communication result。

`/Users/gaofan.86/Documents/semantic-manifold/refs/related-work.md:7`–`:8` 和 `:37`–`:42` 的 novelty确认、关键词零结果没有完备检索或论证；arXiv conjunction zero hits不能证明不存在related geometry RD/PH quantization/learned compression。该项目自己 `:44`–`:48` 也列 IEEE Xplore、G–L、Mabrok精读未完。TCOM 应对近邻逐条判断对象、distortion、rate、channel、assumptions、guarantee，不能把“无某关键词”当originality。

## 9. 建议怎么参考：保留当前主线的最短路线

1. **建模只取 reach/局部近似有效域。** 保留fixed posterior task metric g、learned nondegenerate encoding h、material source measure、true rate accounting。补 residual/noise covariance、injectivity/margin是否满足的可测指标；不假设CLIP空间天然是compact smooth connected manifold。
2. **理论只新增通信接口。** 已有率夹逼→fixed H throughput→Rayleigh asymptotic outage上下界；receiver Jacobian与task metric→local sensitivity/power；reach/margin→条件projection safety。经典引理清楚引用，新的贡献放在这些接口的受限可计算结果，而非改名包装Zador/Gersho/PH stability。
3. **机制实验增量有限。** circle/sphere projection可作为receiver机制sanity，必须计同source、power、channel uses、错branch代价；补Rice analytic control和CI。sphere exact RDF保持当前版本，不能搬external coefficient或“validated”摘要。
4. **真实主实验沿用当前GPU计划。** SEMASIA小特征子集、classification+declared retrieval/direct posterior baseline、digital actual-bit轨道和analog fixed-symbol轨道、density/isometry/bending/joint分解、AWGN+Rayleigh+CSI residual。这里外部项目没有可复用dataset loader、trained models或real result。当前 geometry gain若仍被density-only解释，不应“包装”成independent manifold performance gain。
5. **工程组织可直接吸收。** common模型/记号、proof dependency、每定理一个script、显式assumption checklist值得保留，但要新增结果哈希、运行命令、raw CSV、参数、seed、CI、actual-failure exits，以及proved/computed/pending/falsified台账。当前代码/计划已有这些大部分，整合时避免形成另一套互相矛盾的draft。

可向用户直接概括：外部项目提供很好的研究问题清单，却没有可直接复制的成熟新定理；最有用的迁移是“局部投影机制与reach有效域、受限码点密度基线、通信outage接口、proof-script组织”。曲率普适RD、维数比SNR增益、拓扑KB普适底板等强主张不应移入当前TCOM稿。

## 10. 本次独立核验的一手来源

- [Iacobelli: Asymptotic quantization for probability measures on Riemannian manifolds](https://arxiv.org/pdf/1412.4026)：Thm1.4/Cor1.5，流形量化首項与point density；不是普适second-order curvature RDF。
- [NSW: Finding the Homology of Submanifolds with High Confidence from Random Samples](https://math.uchicago.edu/~shmuel/NSW1.pdf)：Prop6.3真实弦弧界与Prop7.1 noisy tube/density条件。
- [Boissonnat–Lieutier–Wintraecken: The Reach, Metric Distortion, Geodesic Convexity and the Variation of Tangent Spaces](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SoCG.2018.10)：Thm1/Lemma5，sharp reach/chord-geodesic bound。
- [Chazal–de Silva–Oudot: Persistence stability for geometric complexes](https://arxiv.org/abs/1207.3885)：几何filtered complex/GH stability，应区分Rips/Čech/offset尺度。
- [Kostina–Verdú: Lossy joint source-channel coding in the finite blocklength regime](https://arxiv.org/abs/1209.1317)：finite-blocklength excess-distortion含source/channel dispersion。
- [Riegler–Koliander–Bölcskei: Lossy compression of general random variables，作者ETH页面](https://www.mins.ee.ethz.ch/pubs/p/lossy-2021)：reference-measure/small-ball general compression已有理论，应作为基础引用。

以上网络核查只使用作者/出版社/arXiv原始来源。数学修正和反例是根据源公式的独立推导，不将网页摘要视作新证明。当前GPU计划与结果表只据保存内容和主审计实际核验，不声称本次训练过模型。

轻量核查脚本：`/private/tmp/semantic-manifold-review-20261004/recheck_other_project.py`；实际输出：`/private/tmp/semantic-manifold-review-20261004/other_project_recheck.txt`。用外部现有环境执行，未安装依赖，未改项目文件。t14 原脚本的 stdout 已在本审计工具记录中逐行核对。
