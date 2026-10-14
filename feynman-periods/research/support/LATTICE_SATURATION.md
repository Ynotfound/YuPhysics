# 边界格的本原性与指定纤维的完整超越格

本稿补上 Picard 证书之后缺少的整数饱和论证。采用已有 Hessian 格定理，并逐项核对本模型的适用条件，得到：目标 \(m=(1,2,4,8),s=-64\) 的二十条边界曲线生成整个几何 Néron–Severi 格，完整超越格为
\[
 T(X)\simeq U\oplus U(2)\oplus A_2(-2),\qquad \det T(X)=48.
\]
两个旧 Picard 样本也有同一结论。这没有认证候选 Q6 所标记的六个循环；它先确定独立于该候选矩阵的几何目标格。

## 一、从 banana 方程到 Sylvester Hessian

设 \(\mu_i=m_i^2\ne0\)、\(s\ne0\)，并记
\[
 F=M U_G-s\alpha_1\alpha_2\alpha_3\alpha_4,
 \quad M=\sum_i\mu_i\alpha_i,
 \quad U_G=\sum_i\prod_{j\ne i}\alpha_j.
\]
取射影坐标
\[
 y_i=\mu_i\alpha_i\ (i\le4),\qquad y_5=-M,
 \qquad H=\{\textstyle\sum_{i=1}^5y_i=0\}\simeq\mathbb P^3.
\]
这是 \(\mathbb P^3_\alpha\) 到 H 的线性同构。令
\[
 \lambda=(\mu_1,\mu_2,\mu_3,\mu_4,s),\qquad
 Q_\lambda=\sum_i\lambda_i\prod_{j\ne i}y_j.
\]
直接得到多项式恒等式
\[
 Q_\lambda(\mu_1\alpha_1,\ldots,\mu_4\alpha_4,-M)
 =-\Bigl(\prod_{i=1}^4\mu_i\Bigr)F. \tag{1}
\]

考虑 H 上的三次曲面
\[
 S_\lambda:\quad \sum_{i=1}^5\frac{y_i^3}{\lambda_i}=0.
\]
它的 Sylvester 三次系数是 **1/λ_i**，不是 λ_i。这一点避免了将原文三次系数与本项目四次系数混用。
用 y₁,…,y₄ 为坐标，Hessian 为
\(6\operatorname{diag}(y_i/\lambda_i)+6(y_5/\lambda_5)\mathbf1\mathbf1^t\)。
对角阵加秩一阵的行列式恒等式给出
\[
 \det\operatorname{Hess}(S_\lambda)
 =\frac{6^4}{\prod_{i=1}^5\lambda_i}\,Q_\lambda. \tag{2}
\]
因此本项目四次曲面确为该三次曲面的 Hessian。

三次曲面奇异当且仅当存在符号使 \(\sum_i\epsilon_i\sqrt{\lambda_i}=0\)：Lagrange 乘子条件为 \(y_i^2/\lambda_i\) 全相等。对正质量和实 s<0，前四项实、第五项非零纯虚，故不可能相消。等价地，坏参数为 s=0 或 \(s=(\sum_{i=1}^4\epsilon_i m_i)^2\)，以及某个 μ_i=0。三个已经认证 Picard 数的负参数样本均满足光滑三次曲面条件。

## 二、二十条曲线和整数基的逐项对应

用 1,…,5 为文献标记（原文 0,…,4 整体加一）。令 N_ab 为直线 \(y_a=y_b=0\) 的严格变换，N_abc 为节点 \(y_a=y_b=y_c=0\) 上的例外曲线。

若 \(\{i,j,k,l\}=\{1,2,3,4\}\)，则原稿曲线与 Hessian 曲线的对应为

| 本项目曲线 | Hessian 标记 | 几何理由 |
|---|---|---|
| V_i | N_jkl | 原坐标顶点仅 α_i 非零，故 y_j=y_k=y_l=0 |
| E_ij | N_kl5 | 边上节点满足 M=0，且 α_k=α_l=0 |
| C_ij | N_kl | 原坐标边 α_k=α_l=0 |
| F_i | N_i5 | 原边界线 α_i=M=0 |

两组各十条 (−2) 曲线内部互不相交；N_ab 与 N_cde 相交一次，当且仅当 \(\{a,b\}\subset\{c,d,e\}\)。这逐项复现已经保存的二十阶 Gram 矩阵。

本项目选取的十六元基删去 E12、F2、C12、F3。它们恰对应 N345、N25、N34、N35，也就是原文删去的 N234、N14、N23、N24 整体加一。因而不仅抽象交图相同，所用整数基也与文献一致。全部其他曲线的整数坐标此前已在 boundary_projector.json 给出。

## 三、文献输入和沿光滑族保持本原性的证明

使用 [Dardanelli–van Geemen，§§1.5–1.6，第2–3页](https://arxiv.org/pdf/math/0409322) 的结论：一般光滑 Sylvester 三次曲面的 Hessian 解消，其 Néron–Severi 格由上述二十条曲线生成；秩为十六，判别式绝对值48，正交补为 \(U\oplus U(2)\oplus A_2(-2)\)。这是已有格定理的应用，不声称重新证明其分类。原 PDF 公式和基标记已直接查看原页，快照保存在 tmp/pdfs/dardanelli-van-geemen-p2.png、p3.png。

为了避免把“一般纤维”未经证明地推广到指定纤维，再给出如下变形论论证。

取复参数空间 B⊂P⁴，为所有 λ_i 非零且 S_λ 光滑的开集。B 是不可约复代数簇的非空 Zariski 开集，故解析路径连通。Hessian 族有十个固定节点 p_ij=[e_i−e_j]。沿十个互不相交的节点截面吹起，得到光滑射影 K3 族 \(\pi:\mathcal X\to B\)；局部非退化二次切锥保证吹起后各纤维光滑。十条例外曲线和十条边的严格变换构成带标记的相对除子。

因此二十个整上同调类为 \(R^2\pi_*\mathbb Z\) 的平坦截面，给出秩十六子局部系统 L。沿任何路径的整平行运输同时识别
\[
 (H^2(X_b,\mathbb Z),L_b)
 \simeq(H^2(X_{b'},\mathbb Z),L_{b'}).
\]
商群 H²/L 的挠部分及本原闭包指数因而局部常数。

在文献定理给出的一般纤维，L=NS。对复 K3，NS 在 H²(Z) 中本原：若 n v 是 (1,1) 整类，则 v 也是 (1,1) 整类，由 Lefschetz (1,1) 定理 v∈NS。所以该一般纤维的 H²/L 无挠。路径连通性推出 **B 的每个纤维中 L 都本原**，无须假设每个纤维 Picard 数相同。

特别地，若某个 b 满足 ρ(X_b)=16，则 L_b⊂NS(X_b) 同秩，而 L_b 在 H²(Z) 本原，故
\[
 \boxed{NS(X_b)=L_b.} \tag{3}
\]
由于本原嵌入及其正交补沿整平行运输保持，所有纤维的边界正交格都同构于文献的一般正交格。在 ρ=16 的纤维，它就是完整 T(X_b)。在 ρ>16 的特殊纤维，它仍是秩六边界正交格，但不能称为完整超越格。

本论证补充的是整数嵌入的信息；单独使用秩十六和 det L=−48 无法代替它。

## 四、三个指定纤维的结论

[PICARD_CERTIFICATE](PICARD_CERTIFICATE.md) 已分别证明
\[
 (m,s)=((1,2,4,8),-1),\quad((1,3,7,15),-2),
 \quad((1,2,4,8),-64)
\]
的几何 Picard 数为十六。由 (3)，它们均有 NS=L，故 NS 的 Gram 行列式为 −48，超越格可以取
\[
 Q_T=\begin{pmatrix}0&1\\1&0\end{pmatrix}
 \oplus\begin{pmatrix}0&2\\2&0\end{pmatrix}
 \oplus\begin{pmatrix}-4&2\\2&-4\end{pmatrix}. \tag{4}
\]
其签名为 (2,4)，行列式48，Smith 不变量为 (1,1,2,2,2,6)。因此判别群为 \((\mathbb Z/2)^3\oplus\mathbb Z/6\)。此处的格同构为无标记整格同构；尚未构造将原六个数值循环送到式 (4) 基的几何变换。

## 五、对 Q6 的精确影响

目标纤维的真实 det T=48 已不再是条件。Q6_GEOMETRIC_CERTIFICATE现已证明六个规定几何循环属于T且Gram为Q6（det192），故它们生成的格L₆满足
\[
 [T:L_6]^2=192/48=4,\qquad[T:L_6]=2. \tag{5}
\]
这也证明它们不是完整整基。先前有限代数计算给出的十一种偶指数二补全中，只有 w=(e₁+e₂+e₃+e₄)/2 的类被四个候选质量矩阵保持。因此若再识别这四个矩阵为真实几何运输，该唯一性就确定真实补全。

实际循环配对与十五对迹、六三迹现已认证；质量运输标记和具体补全类仍未认证。绝对周期矩阵另需几何初值与 de Rham 归一化。不能因为候选补全与 (4) 抽象同构，就反推全部标记已经正确。

有限代数输入核对见 [hessian_lattice.json](data/hessian_lattice.json)。它记录方程恒等式、二十曲线与相同基的对应、Smith 形及文献 PDF 哈希。本原性由第三节的文献输入和族论证承担，脚本不自行认证这一几何命题。

最新MARKED_LATTICE已认证实际T0及全部八节点循环；全射线单值同时保持三种偶指数二补全，故仍不能据此选出w。六循环的绝对周期尺度另由ABSOLUTE_COMPARISON完成，不依赖先选出补全类。
