# 指定活动对象的十三维幺幂群与实际 coaction 结构

本稿固定 m=(1,2,4,8)、s=−64 的秩八活动对象 E。LOWER_REGULATOR 已证明其下层扩张在实 Hodge 实现中非零；上层质量 Kummer 扩张亦非零。由这两项和纯核简单性，可以识别 Hodge 群和 Nori 动机群的全部幺幂部分，而不必先识别固定 K3 纤维的纯动机群。

## 1. 记号及不能省略的输入

按权重由低至高排列
\[
 A=T(X)(3),\qquad B=\mathbb Q(1),\qquad C=\mathbb Q(0),
 \quad \dim(A,B,C)=(6,1,1).
\]
权重依次为 −4、−2、0。令 N=W_{-2}E、K=E/A。两项实际扩张为
\[
 0\to A\to N\to B\to0,\qquad 0\to B\to K\to C\to0.       \tag{1}
\]
第一项非零由 LOWER_REGULATOR 的实分裂障碍认证。第二项非零由实际平方质量点给出：其活动单位是2的一个非零有理幂，非零实对数不属于 2πiℚ。A 的有理 Hodge 结构因 Picard 数十六而简单；B、C 简单。

记 G_H(E) 为有理混合 Hodge 结构 E 的 Mumford–Tate 群，
\[
 R_H=G_H(A\oplus B),\qquad R_N=G_N(A\oplus B).
\]
这里 G_N 是既定 Nori 范畴和 Betti 纤维函子的实际群。R_H、R_N 均为各自的实际纯分级群，并未被替换为射线几何单值群 O(Q6)。

选择一个有理权重分裂。保持过滤且在分级上恒等的最大群为
\[
 U_{\max}=\left\{
 \begin{pmatrix}I_6&x&z\\0&1&y\\0&0&1\end{pmatrix}:
 x,z\in\mathbb Q^6,\ y\in\mathbb Q
 \right\}.                                             \tag{2}
\]
式 (2) 按任意 ℚ 代数解释为群概形；不是只取其有理点。维数十三，中心为六维 z 块，是二步幺幂群。

## 2. Hodge 群的相邻块全部出现

先考虑一般的两步扩张 0→A₁→V→A₂→0，其中纯分级半单。其 Hodge 群到分级群的核是 Hom(A₂,A₁) 中的向量子群。若此 Hom 表示简单且扩张非零，则核必为整个 Hom：

核若为零，V 的群为约化群，其特征零表示半单，扩张随之分裂，矛盾；非零子表示由简单性只能是全部。这个论证针对实际群，不由任意形式矩阵推断群等号。

应用于式 (1)，下层给出全部 Hom(B,A)≅T(X)(2)，维数六；上层给出全部 Hom(C,B)≅ℚ(1)，维数一。

令 U_H 为 G_H(E) 在分级上的核。分级是可极化纯 Hodge 结构的直和，群 R_H 约化，所以 U_H 正是幺幂根基。Tannakian 商到 N、K 的群是满射；特征零下满射把幺幂根基满射到商的幺幂根基。因此 Lie(U_H) 在两个相邻块上的投影分别满射。

其相邻块联合像是
\[
 S=\operatorname{Hom}(B,A)\oplus\operatorname{Hom}(C,B)
\]
的 R_H 子表示。S 的两个直和因子是维数六和一的非同构简单表示。R_H 约化，故子表示半单；分别满射迫使联合像等于整个 S。只凭两个投影非零一般并不足够，这里“非同构简单”排除了两块之间的对角关系。

## 3. 交换子迫使全部中心出现

写幂零矩阵
\[
 n(x,y,z)=\begin{pmatrix}0&x&z\\0&0&y\\0&0&0\end{pmatrix}.
\]
直接乘法给
\[
 [n(x,y,z),n(x',y',z')]
       =n(0,0,xy'-x'y).                                 \tag{3}
\]
因第二节联合像为全部 S，可分别提升 (x,0)、(0,1) 到 Lie(U_H)；提升中的中心项不影响式 (3)。其交换子遍历 Hom(C,A) 的全部六个方向。故 Lie(U_H) 同时有七维相邻商和六维中心核，等于 Lie(U_max)。特征零幺幂群连通且由 Lie 代数确定，因此
\[
 \boxed{U_H=U_{\max},\qquad\dim U_H=13.}                  \tag{4}
\]

这是混合 Hodge 群的结论，与 GEOMETRIC_MONODROMY_GROUP 的纯参数单值闭包不同。一般三步扩张的框架可对照 [Eskandari，Theorem A](https://arxiv.org/html/2407.01379#S1)；上面的简单性、相邻块及交换子论证给出本对象的具体应用。没有把该文的“很一般”结论用于固定纤维。

## 4. Nori 群的幺幂部分也完全确定

Hodge 实现是忠实正合张量函子，并与所用 Betti 纤维函子相容。因此实际矩阵群满足
\[
 G_H(E)\subseteq G_N(E)\subseteq P,
\]
其中 P 为权重抛物群。在两侧取“作用于分级恒等”的核，得到
\[
 U_{\max}=U_H
 \subseteq\ker(G_N(E)\to R_N)\subseteq U_{\max}.
\]
于是中间核等于 U_max。

这里也可确认它就是 G_N(E) 的幺幂根基，而非只给出某个幺幂子群：A 的 Hodge 实现简单，忠实正合性使 A 在 Nori 范畴中也简单；B 同样简单。因此 A⊕B 是半单且在其生成群上为忠实表示，R_N 约化。结论为
\[
 \boxed{R_uG_N(E)=R_uG_H(E)=U_{\max}.}                   \tag{5}
\]
这一步不假设所有纯 Nori 动机均半单，也不假设一般 Hodge 实现全忠实。

进一步，由到实际分级群的满射及完整的核，在选定权重分裂后有
\[
 \boxed{
 G_H(E)=P^{-1}(R_H),\qquad
 G_N(E)=P^{-1}(R_N).
 }                                                     \tag{6}
\]
P→GL(A)×G_m 为对角块投影，顶层 C 的作用恒等。式 (6) 意味着对固定的实际分级群，没有额外削小 x、y、z 块的混合群关系。

R_H 与 R_N 的等号、固定纯 K3 的全部 Hodge/动机张量及其具体群仍待独立识别。式 (6) 是相对于这些实际分级群的完整群描述，不能省略 R_N 而直接写为某个预设正交群与 U 的半直积。

## 5. 明确的坐标 Hopf 结构

在上述分裂下，实际群元素写为
\[
 g=\begin{pmatrix}H&x&z\\0&t&y\\0&0&1\end{pmatrix},
 \qquad (H,t)\in R_N .
\]
由实际群等号 (6)，x、y、z 是完整幺幂方向。其坐标余乘法为
\[
 \begin{aligned}
 \Delta H_{ij}&=\sum_kH_{ik}\otimes H_{kj},&
 \Delta t&=t\otimes t,\\
 \Delta x_i&=\sum_jH_{ij}\otimes x_j+x_i\otimes t,&
 \Delta y&=t\otimes y+y\otimes1,\\
 \Delta z_i&=\sum_jH_{ij}\otimes z_j+x_i\otimes y+z_i\otimes1 .
 \end{aligned}                                         \tag{7}
\]
H、t 之间仍须满足 R_N 的真实纯关系；本稿没有删去这些关系。

若将普适比较写作 dR→Betti 的矩阵，并以右复合 de Rham 群规定右 coaction，可用兼容权重的基写
\[
 \mathcal P=\begin{pmatrix}P_A&p&f\\0&\tau&q\\0&0&1\end{pmatrix},
 \qquad \Delta\mathcal P=\mathcal P\otimes g_{\rm dR}.
\]
从而
\[
 \begin{aligned}
 \Delta P_A&=P_A\otimes H,&\Delta\tau&=\tau\otimes t,\\
 \Delta p&=P_A\otimes x+p\otimes t,&
 \Delta q&=\tau\otimes y+q\otimes1,\\
 \Delta f&=P_A\otimes z+p\otimes y+f\otimes1 .
 \end{aligned}                                         \tag{8}
\]
另一左右作用约定会转置或调整因子顺序，使用时必须先固定约定。P_A 属于 T(X)(3) 的相应比较，不能未经 Tate 扭转与对偶变换就等同于裸二形式的 P6。

式 (7)–(8) 现是实际秩八对象的结构公式，而不仅是允许的形式块群。但它们尚未求值所有比较条目：当前只认证了纯六循环 P6 以及五个球上的 β 周期；第六混合周期、顶层相对条目及完整基转换仍需计算。也不据幺幂维数推断复周期值的代数独立性。
