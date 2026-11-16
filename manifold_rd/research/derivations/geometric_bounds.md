# 连续流形的 Gibbs 界与小失真极限

测度为黎曼体积，损失为平方测地距离。由条件熵界和显式测试信道分别推出两侧。

[研究索引](../README.md) · [符号和引用所用书目](../../ref.bib)

### Geometry-only finite-distortion bounds and the leading slope

For the squared geodesic cost, define $Z_s(a)=\int_\mathcal{M}e^{-s d_g(u,a)^2}dV_g(u)$ with $s>0$. Comparing the conditional density of $U|A=a$ to $e^{-s d_g(u,a)^2}/Z_s(a)$ yields 

$$
h_g(U|A)\le sD+\ln\sup_a Z_s(a),\quad
R_{\mathrm{nat}}(D)\ge \sup_{s>0}\{h_g(U)-sD-\ln\sup_a Z_s(a)\}.
$$

This holds when the differential entropies exist; a singular conditional law with infinite information cannot improve a finite optimum. Conversely send according to the explicit continuous kernel $W_s(da|u)=e^{-s d_g^2(u,a)}dV_g(a)/Z_s(u)$. It gives 

$$
D_s=-\mathbb{E}\,\partial_s\ln Z_s(U),\qquad
R_{\mathrm{nat}}(D_s)\le h_g(A_s)-\mathbb{E}\ln Z_s(U)-sD_s .
$$

Both follow by substituting the kernel in its normalized expectation and conditional entropy; neither requires identifying the optimal kernel at arbitrary $D$. For smooth positive density on a smooth compact boundaryless manifold with positive injectivity radius, normal coordinates $v$ have $d_g(u,a)^2=\|v\|^2$ locally and $dV=(1+O(\|v\|^2))dv$ uniformly. Contributions beyond a fixed normal neighborhood are exponentially small. Gaussian integration gives $Z_s=(\pi/s)^{m/2}(1+O(1/s))$ and, integrating the second moment independently, $D_s=m/(2s)+O(1/s^2)$. The symmetric kernel and uniform expansion make the output density an approximate identity applied to $p$ (including the factor $1/Z_s$); it converges uniformly to $p$. Compactness and positivity then yield $h_g(A_s)\to h_g(p)$. Taking $s=m/(2D)$ in the lower bound and inverting the continuous decreasing $D_s$ for the upper bound recovers $h_g(p)-\tfrac m2\ln(2\pi eD/m)+o(1)$. Boundaries, singular strata, nonpositive densities or another loss need separate control. In particular, the finite empirical graph below only approximates distances; it has no literal continuum high-resolution limit at fixed sample count.
