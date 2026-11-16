# 位置询问的通信资源下界

接收机得到噪声信号，评估器按相同位置询问收发两侧。对离线询问条件化，才能把二元错误熵转成信道互信息下界。

[研究索引](../README.md) · [符号和引用所用书目](../../ref.bib)

### A converse for preserving identifiable relations

Let $L=(a,b,t)$ specify a query drawn before observing channel noise. Given $X$ and shared context $C$, the query sampler and the physical channel use independent randomness. The offline evaluator knows $L$ when it applies the same readout to $X$ and $\hat X$. Let $B=A(f_X,L)$ and $\hat B=A(f_{\hat X},L)$, with error probability $e$.

**Proposition.** If the decoder output depends only on the received signal $V$ (including receiver CSI) and $C$, then 

$$
[H_2(B\mid L,C)-h_2(e)]_+ \le I_2(X;V\mid C).
$$

For $d$ named queries on the same image, replace the entropy by $H_2(B_1,\ldots,B_d\mid L_{1:d},C)$ and subtract $\sum_j h_2(e_j)$. This is a vector bound for one transmission, not $d$ independent copies of the channel.

**Proof.** Given $L,C$, predicting the binary bit from $\hat B$ leaves conditional entropy at most $h_2(e)$, since $B$ is determined by $\hat B$ and the error indicator. Conditional data processing bounds $I(B;\hat B\mid L,C)$ by $I(X;V\mid L,C)$. Conditional independence of $L$ and $V$ given $X,C$ implies $I(X;V\mid L,C)=I(X;V\mid C)-I(L;V\mid C)\le I(X;V\mid C)$. For the vector, condition jointly and upper bound the residual entropy by the sum of individual error-indicator entropies. Normalize by the number of source samples only after defining what one transmission represents.

For complex AWGN with $k$ uses of unit mean energy this is at most $k\log_2(1+\mathrm{SNR}_{\rm lin})$ by data processing and the Gaussian capacity bound. In a fading link the appropriate bound must use the declared CSI, power policy and outage or averaging assumption. A pooled $H(B)$ is *not* the conditional entropy in the proposition. For an empirical converse, preselect a small finite, source-independent dictionary of aligned endpoints and thresholds, estimate the joint answer law on held-out images and report its sampling uncertainty. The hybrid foreground sampler is useful for training but can itself reveal information about $X$; the theorem explicitly conditions on the query. This converse connects the branch's finite dictionary reasoning to the new image model without assuming a source-event taxonomy.
