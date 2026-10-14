# 指定质量射线的纯几何单值群

本稿基于 Q6_GEOMETRIC_CERTIFICATE 已认证的实际六循环 Gram，不再以候选矩阵为前提。设 V=T_(s=−64)⊗Q，Q=Q₆，Γ 为 m=(1,2,4,8) 射线在光滑参数域上的纯 K3 几何单值群。

**定理。**
\[
 \boxed{\ \overline{\Gamma}^{\,\mathrm{Zar}}=O(Q)\ },
 \qquad
 (\overline{\Gamma}^{\,\mathrm{Zar}})^0=SO(Q).
                                                               \tag{1}
\]
等式为有理代数群等式；复化后为六维完整正交群。它不声称 Γ 等于整格全部等距群，也不计算该算术群中的指数；它不识别固定纤维或相对混合对象的绝对动机群。

## 1. 可精确验证的李代数见证

在六循环有理基中，六个实际反射为
\[
 S_i=I+e_i e_i^tQ,\qquad S_i^2=I,\quad
 S_i^tQS_i=Q,\quad\det S_i=-1.
\]
S₆=T₀S₁T₀⁻¹ 也属于 Γ。取
\[
 U=S_1S_2,\quad A=U-I,\quad X=A-\tfrac12A^2.
\]
精确乘法验证
\[
 A^3=0,\quad \operatorname{rank}A=2,\quad
 \operatorname{rank}A^2=1,\quad
 X^3=0,\quad e^X=I+X+\tfrac12X^2=U.                    \tag{2}
\]

令 G 为 Γ 的 Zariski 闭包。对任意在 G 上消失的正则函数 f，多项式
\[
 t\longmapsto f(I+tX+\tfrac12t^2X^2)
\]
在所有整数 t 上为零，因为该矩阵等于 U^t。故它恒为零，于是整条一参数子群 exp(tX) 属于 G，特别 X∈Lie(G)。此处直接证明所需闭包事实，不把数值“接近幂零”当作幂零。

Lie(G) 对所有 Ad(S_i) 稳定。按列表从左至右递归作用 Y←S_iYS_i，以空字代表 X，取以下十五个字：
\[
 \begin{gathered}
 (),\ (3),\ (4),\ (5),\ (6),\
 (3,1),\ (3,4),\ (3,6),\\
 (4,1),\ (4,5),\ (4,6),\
 (5,1),\ (5,6),\ (6,1),\ (3,1,5).
 \end{gathered}
\]
每个所得 Y 满足 YᵗQ+QY=0。以 QY 的十五个严格上三角元素按字典顺序为坐标，将这些向量作为列，得到整数15×15矩阵 C，满足
\[
 \det C=210119944214597861376\ne0.                       \tag{3}
\]
完整矩阵、字和逐项结果已保存在 geometric_monodromy.json。

## 2. 从李代数维数到群等式

因为 Q 非退化，Y↦QY 把 so(Q) 识别为反对称矩阵空间，维数为15。式(3)证明 Lie(G)=so(Q)。在特征零中代数群光滑，故 dim G⁰=15。SO(Q) 连通，而 G⁰⊂SO(Q) 为同维闭连通子群，所以 G⁰=SO(Q)。Γ 含任意一个 det=−1 的 S_i，遂 G 同时含正交群的两个分支，得到式(1)。

所用正交群连通分支事实可参见 [Milne, Algebraic Groups, §2.43](https://www.jmilne.org/math/Books/iAG2022.pdf#page=64)：非退化二次型在特征不为二时，SO 为连通群，O 的分支由行列式±1区分。本稿的十五维见证是实际 Q₆ 的精确计算，不由该一般文献代替。

## 3. 绝对不可约性的直接证明

也可不使用群分类，直接证明这六个反射在 V⊗C 上绝对不可约。设非零子空间 W 对所有 S_i 不变。若每个 (S_i−I)W=0，则 W 与所有 δ_i 正交；六循环张成非退化空间，迫使 W=0，矛盾。因此某个 δ_i∈W。若 q_ij≠0，则 (S_j−I)δ_i=q_ijδ_j，故 δ_j∈W。Q 的非零配对图连通，最终六个 δ_j 均属于 W，因此 W=V⊗C。

这给出指定射线纯六维局部系统的不可约性，不等价于任意固定代数纤维下层 regulator 非零。

## 4. 证书和下一步

生成器使用精确符号矩阵；独立验证器只用 Fraction 有理运算、直接矩阵乘法及高斯消元，重算式(2)、全部共轭和式(3)，且重复一列的负对照给行列式零。

- data/banana_pure_ray/geometric_monodromy.json，SHA256=e5cb8de592258eefc1d37d9cd716db28b37bf291542343691edd805004e9ad06。
- data/banana_pure_ray/geometric_monodromy_replay.json。
- code/verify_geometric_monodromy.py。

该群等式使此前“Zariski闭包至少有两个分支”的结论得到加强。绝对动机群、混合扩张的幺幂根、指定纤维比较条目、质量方向的整格标记及绝对周期初值仍是独立任务。尤其不可把参数族的完整几何单值群直接当作 s=−64 单个动机的群。
