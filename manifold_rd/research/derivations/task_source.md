# 任务、共享知识与局部度量

从观测变量推到条件信源，区分观察损失、任务超额对数损失和对象位移。体积熵与小失真的双边推导见同目录的连续几何笔记。

[研究索引](../README.md) · [符号和引用所用书目](../../ref.bib)

## An operational source and its geometry

Let $X$ be the observation, $C=(K,t)$ the context shared before encoding, and $Y_t$ the unknown task variable. The reference $U=E_0(X,C)$, the source law $\mu_c$, the allowed output set $\mathcal A_c$, and the evaluation cost $\delta_c$ are fixed before fitting a communication system. For an i.i.d. source--context pair, additive distortion and context known at both ends, the conditional single-letter characterization in bits is 

$$
R_\delta(D\mid C)=\inf_{W:\,\mathbb{E}\delta_C(U,A)\le D} I_2(U;A\mid C).
$$

This is an asymptotic source-coding statement. For indirect observation the effective cost for an action is $\mathbb{E}[\ell(Y_t,a)\mid X,C]$ [`liu2022indirect`]; replacing $X$ by $U$ can lose information. Under Bayes log loss the irrecoverable increase is $H_2(Y_t\mid U,C)-H_2(Y_t\mid X,C)=I_2(Y_t;X\mid U,C)$. The task identifier is shared; the sample's true label is not side information.

Three specializations correspond to the proposal. (i) On a smooth stratum with a specified metric $g_0$, evaluate squared geodesic distance to a manifold-valued reconstruction. (ii) When $C$ changes the source, use $(\mu_c,g_c,\mathcal A_c)$. Conditioning on $c$ produces a family of feasible kernels, hence 

$$
R(D\mid C)=\inf_{\sum p_cD_c\le D}\sum p_c R_c(D_c).
$$

A feasible global kernel gives the lower inequality by conditioning; choosing each conditional kernel within $\epsilon$ of its own optimum gives the reverse. Supported interior points have the same slope, not necessarily the same $D_c$. (iii) For a task posterior $q(u,c)$ use directed $D_{\mathrm{KL}}(q(u,c)\Vert q(a,c))$ or allow a free simplex action $a$. These are distinct reproduction alphabets.

Let $q_i>0$ be twice differentiable. Expanding $\ln q_i(u+v)$ and using $\sum_iq_i=1$ gives 

$$
D_{\mathrm{KL}}(q(u)\Vert q(u+v))=\tfrac12 v^\top J_q^\top\operatorname{diag}(1/q_i)J_qv+o(\|v\|^2)
$$

in nats. Null directions leave the posterior unchanged to first order. A quotient requires regularity; otherwise keep the semimetric or fix $g=J_q^\top\operatorname{diag}(1/q_i)J_q+\epsilon g_0$ with $\epsilon>0$. The pullback is a local construction [`arvanitidis2022pullinginformationgeometry`]; the KL is generally asymmetric and is not a global squared geodesic. Point-to-support distance is zero even for two far-apart on-manifold states.

For a smooth positive density $p$ on a compact regular $m$-manifold without boundary, squared geodesic loss and on-manifold reproduction, the small-distortion benchmark in nats has leading form 

$$
h_g(p)-\tfrac m2\ln(2\pi eD/m)+o(1),\qquad h_g(p)=-\int p\ln p\,dV_g .
$$

Local Gaussian cells have covariance $(D/m)I$; the rate is source entropy relative to volume minus the cell entropy. Use rigorous bounds for irregular measures instead of assigning a dimension by inspection [`linder1999high`, `riegler2023lossy`]. The [geometric bounds](geometric_bounds.md) give matching leading terms. A finite dataset saturates at finite source entropy; its last few slope estimates cannot establish a continuum theorem.
