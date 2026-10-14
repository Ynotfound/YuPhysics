# 只由消失特征向量认证节点反射

本稿给出供 Q6 后续闭路认证使用的理论简化，以及不依赖浮点尾项的显式误差引理。**后续实际计算已认证全部八个节点的 Frobenius 向量、尾界和局部反射**，详见 [NODE_REFLECTION_CERTIFICATES](NODE_REFLECTION_CERTIFICATES.md)；本稿保留引理的推导。五个目标节点的支路、共同基点反射及十二项迹现已在 [BASED_NODE_TRANSPORT](BASED_NODE_TRANSPORT.md) 认证。开放干路见 [HORIZONTAL_TRANSPORT](HORIZONTAL_TRANSPORT.md)。

## 一、把节点 residue 对角化

在常数 D 基中写连接为 A(s)。对任一八个节点 a，已精确证明

\[
A(s)=\frac{R_a}{s-a}+A_{\rm reg}(s),\qquad
R_a^2=-\tfrac12R_a,\quad\operatorname{rank}R_a=1.
\]

取有理可逆矩阵 Q，第一列为 −1/2 特征向量，其余五列为 ker Rₐ 的基，便有

\[
Q^{-1}R_aQ=\Lambda=\operatorname{diag}(-1/2,0,0,0,0).
\]

取非零高斯有理 h，令 s=a+ht，所用闭单位圆盘不能含其他极点。在 Q 基中寻找 y(t)=t⁻¹ᐟ²g(t)。则

\[
t g'=D_0g+tH(t)g,\quad
D_0=\Lambda+\tfrac12I=\operatorname{diag}(0,1/2,1/2,1/2,1/2),
\quad H(t)=hQ^{-1}A_{\rm reg}(a+ht)Q. \tag{1}
\]

固定 g(0)=e₁。对每个正整数 n，nI−D₀ 可逆，因而形式解析解唯一。下面的残差引理还直接给出其收敛和端点包围，不需要把形式级数的存在误称为收敛证书。

## 二、加权 Frobenius 残差引理

设 H=N/q 在闭单位圆盘的邻域全纯，|q|≥ℓ>0。给出高斯有理向量多项式 P，P(0)=e₁，满足

\[
q(t)tP'(t)-q(t)D_0P(t)-tN(t)P(t)=t^mS(t),\qquad m\ge1. \tag{2}
\]

这里整除 tᵐ 须精确检查，不能靠浮点小量判定。普通精确 Frobenius 递推到次数m−1可产生这种见证。

以矩阵无穷范数和向量无穷范数分别取有理上界

\[
M\ge\sup\|H\|,\qquad
\delta\ge\sup\|S/q\|,
\qquad \kappa:=\frac{M}{m+1/2}<1.
\]

**引理。** 式 (1) 有唯一指定初值的解析解，且

\[
\boxed{\|g(t)-P(t)\|_\infty
\le |t|^m\frac{\delta}{(m-1/2)(1-\kappa)}}
\quad(|t|\le1). \tag{3}
\]

证明：令 g=P+tᵐu。将式 (2) 代入 (1)，得到

\[
t u'=(D_0-mI)u+tH(t)u-S(t)/q(t).
\]

其解析解等价于下列不动点方程：

\[
u(t)=\int_0^1 x^{mI-D_0-I}
\left[txH(tx)u(tx)-S(tx)/q(tx)\right]dx. \tag{4}
\]

因 D₀ 已对角化，对0<x≤1有 ∥x^{mI−D₀−I}∥∞=x^{m−3/2}。故右侧的线性部分范数至多 M∫₀¹x^{m−1/2}dx=M/(m+1/2)=κ，非齐次部分范数至多 δ/(m−1/2)。在闭圆盘上连续、内部全纯的向量函数空间中，它是收缩映射；不动点全纯，范数至多 δ/((m−1/2)(1−κ))。微分式 (4) 还原原方程。正整数 n 的 nI−D₀ 可逆给指定初值解析解的唯一性。于是得到式 (3)。

全部界可由清分母多项式的系数绝对值和、其他极点的有理距离下界取得。因子 h 已纳入 H，不得遗漏。选择 m 或 h 直到严格满足 κ<1；如果 κ≥1，只能判定本估计不足，不能宣布解析解不存在。

此界没有普通 Gronwall 的3^ceil(M)放大，代价是精确的低阶残差消失和已证明的节点指数。舍入后的近似通常破坏 tᵐ 整除性，必须先修正低阶系数或另给低阶缺陷估计；不能直接把普通步的量化多项式塞入式 (3)。

## 三、一个特征向量决定整个反射

设 t₀≠0 为局部种子，v=Qg(t₀)。因比例 t₀⁻¹ᐟ²不影响特征直线，v 张成局部单值的 −1 特征线。其他五个指数是0，与−1/2无整数差，因此局部单值可对角化，另有五维 +1 特征空间。

令 K=K_D(a+ht₀) 为该种子的精确逆配对矩阵；它非退化，并由实际连接保持。不同单值特征值的两个空间正交，所以 +1 特征空间为 v⊥，而 d=vᵗKv≠0。后一非零性可直接由非退化性推出：若 v 本身也与 v 正交，则 v 与全部六维空间正交，矛盾。因此局部单值为

\[
\boxed{T_{a,\rm loc}=I-2\frac{vv^tK}{v^tKv}.} \tag{5}
\]

更换 v 的非零倍数或把 K 乘以任何非零常数均不改变式 (5)。这说明计算反射不必先知道消失球的绝对周期幅度，也不必先校准极化的整体几何常数。不过指定整循环与绝对周期的认证仍需要这些额外比较。

## 四、把向量包围变成反射包围

设 ∥v−V∥∞≤ε，V 高斯有理，K 精确。取 k≥∥K∥∞、u≥∥V∥∞，并精确计算 d_c=VᵗKV。令

\[
\Delta=6k\epsilon(2u+\epsilon),\quad
d_0=\max(|\Re d_c|,|\Im d_c|),\quad d_-=d_0-\Delta.
\]

三角不等式给 |d−d_c|≤Δ，也给 ∥vvᵗK−VVᵗK∥∞≤Δ。如果 d₋>0，就得到包含真正反射的全有理误差界

\[
\boxed{
\left\|T_{a,\rm loc}-\left(I-2\frac{VV^tK}{d_c}\right)\right\|_\infty
\le2\left(\frac\Delta{d_-}
 +\frac{\|VV^tK\|_\infty\Delta}{d_0d_-}\right).} \tag{6}
\]

式 (6) 来自对分子和倒数分别作差。理论上 d≠0 不保证给定区间已经分离零；若 d₋≤0，必须提高向量精度或调整局部基，不能静默除法。

## 五、与共同基点及迹判据连接

设 Hₐ 为从高基点 c 到该局部种子的严格运输。共同基点的反射为 Hₐ⁻¹T_{a,loc}Hₐ，消失特征线可用 vₐ=Hₐ⁻¹v 表示。所有 vₐ 的配对都在同一个 K_D(c) 下计算。置 dₐ=vₐᵗK_D(c)vₐ、bᵢⱼ=vᵢᵗK_D(c)vⱼ，有

\[
\operatorname{tr}((T_i-I)(T_j-I))
=\frac{4b_{ij}^2}{d_id_j},
\]
\[
\operatorname{tr}((T_i-I)(T_j-I)(T_k-I))
=-\frac{8b_{ij}b_{jk}b_{ki}}{d_id_jd_k}. \tag{7}
\]

这是秩一算子的直接乘法恒等式，故不依赖每条消失线的比例。在几何消失球平方均为−2的归一化下，它们还原为旧稿的 qᵢⱼ² 与 qᵢⱼqⱼₖqₖᵢ。

零点单值 T₀ 保持同一配对，故第六反射 T₀T₁T₀⁻¹ 的特征线为 T₀v₁，且平方等于 d₁。因此 Q6 只需五条节点消失线以及零点作用的认证，再用式 (7) 输出真实十五个对迹和六个三迹区间。

八个实际节点在各自 a+i/4 的向量和反射包围已保存并独立重放，统一局部反射误差≤2^(-160)，详见 NODE_REFLECTION_CERTIFICATES。进一步的 BASED_NODE_TRANSPORT 已把五个目标节点的消失向量运输到共同基点，认证十个对迹、两个三迹及五阶Gram子块。T₀及涉及第六循环的九项迹仍未认证。零点有不同的共振结构，不能套用 D₀=diag(0,1/2,…,1/2) 的引理。完整群等号、固定代数纤维 regulator 非零及绝对比较矩阵也不由上述秩一恒等式推出。

`code/check_frobenius_reflection_lemma.py` 已在独立的耦合模型 g=(1,t/(1−t/4),0,…,0) 中核对该引理：m=4、κ=20/81，端点真实误差1/48被界297/6832包含。该模型确实使用了+1/2特征空间。它还精确核对反射的正交性、平方为单位、二重与三重迹的符号，以及一个复向量扰动的反射包围。结果在 [frobenius_reflection_lemma_check.json](data/frobenius_reflection_lemma_check.json)；这些模型检查不声称是 banana 节点证书。
