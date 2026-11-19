# 语义流形建模、率失真与变分学习

三个目录对应开题报告第一项研究的三篇论文。先读每项工作的研究索引可以找到定理假设、完整证明、对应实验和已测结果。

本目录使用一个 Git 仓库统一维护，三篇论文保持独立目录。提交范围、复现材料与多机同步规则见 [Git 维护约定](GIT_WORKFLOW.md)；本地文献、检查点和数据缓存的初始校验清单见 [LOCAL_ARTIFACTS.sha256](LOCAL_ARTIFACTS.sha256)。

| 工作 | 论文形式 | 主稿 | 研究过程 |
|---|---|---|---|
| 语义流形率失真 | IEEE Transactions（拟投 IEEE Transactions on Communications，现 16 页） | [LaTeX](manifold_rd/main.tex) · [PDF](manifold_rd/build/main.pdf) | [证明、真实特征与信道实验](manifold_rd/research/README.md) · [主线：三层语义率失真](manifold_rd/research/derivations/framework.md) · [接收机字母表与信道极限](manifold_rd/research/derivations/structure_and_limits.md) · [对 Word 与 PPT 的批判性对照](manifold_rd/research/notes/06_framework_critique.md) |
| 变分编码曲面 | IEEE Letter（拟投 IEEE Communications Letters，现 5 页） | [LaTeX](variational_curve/main.tex) · [PDF](variational_curve/build/main.pdf) | [曲面变分、码字概率与四组对照](variational_curve/research/README.md) · [主线：曲面变分泛函（等距亏损、维数、面积–reach–功率）](variational_curve/research/derivations/variational_surface.md) · [有限字母表：超额恒等式与交互秩](variational_curve/research/derivations/interaction_rank.md) |
| 变分拓扑 | IEEE Conference（拟投 IEEE GLOBECOM 2027，现 6 页） | [LaTeX](variational_topology/main.tex) · [PDF](variational_topology/build/main.pdf) | [连通询问、知识字典与通信实验](variational_topology/research/README.md) · [主线：关系泛函的一阶变分、临界负载与变分分配](variational_topology/research/derivations/relation_variation.md) · [对称轨道率失真与衰落逆定理](variational_topology/research/derivations/orbit_rdf_and_fading.md) |

每个目录的 `main.tex` 包含该文从摘要到结论的全部正文，`ref.bib` 保存书目，`figures/` 保存主稿实际引用的图。在对应目录运行 `python build.py` 即可重建 `build/main.pdf`。编译使用 pdfLaTeX 和 BibTeX，各目录随附 IEEEtran 文档类与书目样式。研究中的推导按数学问题维护在 `research/derivations/*.md`，训练、测试与测量记录由 `research/README.md` 索引。实验目录中的图、检查点和逐样本数据用于复现论文里的数值。

[开题报告](../高凡-开题报告-1114.pptx)第 19–20 页提出固定任务下的失真和资源边界，第 21–22 页关注编码映射和码字概率，第 24–25 页将结构误判与通信资源联系起来。[老师的 Word](../语义通信的语义流形理论.docx)提供分析思路，任务定义、局部度量和接收扰动的取舍见[开题材料对照](manifold_rd/research/notes/orientation_at_migration.md)。文献库按主题存放在 [refs](refs/)。

SEMASIA 提供真实图像的预训练特征。工作一和工作二固定各自的任务头、合法重构和评价成本，在同一官方划分内按样本 ID 对齐。每个骨干对应自己的特征分布。实验既可直接分析有限点云，也可检验把它近似为连续流形所需的假设。[数据论文](https://arxiv.org/abs/2605.09485) · [CIFAR-100 特征](https://huggingface.co/datasets/spaicom-lab/semasia-cifar100)

投稿目标与页数限制写在各目录的 `research/ccfa.yaml`，投稿前请对照期刊与会议的最新作者须知核对版式、页数和截稿日期。新增定理的数值核验只依赖 NumPy/SciPy：`manifold_rd/research/tests/test_theory.py`、`variational_curve/research/tests/test_variational.py` 与 `test_interaction.py`、`variational_topology/research/tests/test_variation.py` 与 `test_orbit.py`。
