# 相对费曼动机中纯 K3 核的识别

> 2026-10-14专项状态：旧假设11.2的逐项依赖与实际物理双标架等式现见HYPOTHESIS_11_2_AUDIT.tex（PDF同名）。本稿后部的部分“待证”措辞记录撰写时状态；固定下层非零、P8及固定群的后续证据分别在LOWER_REGULATOR、RELATIVE_TOP_COMPARISON与FIXED_PURE_GROUP，不由本稿或新补篇重新数值认证。新补篇未将结论推广到所有参数，也尚未回写统一全文。

在正质量欧氏区域及已固定的二十曲线模型下，秩十五物理相对 Nori 动机的最大一动机商之核，可由一个明确的 K3 代数对应实现。本稿给出对应和进入相对对象的实际态射。结论在有理 Nori 动机范畴内成立；它不计算扩张周期，也不证明候选块群等于整个动机群。

## 对象及主命题

沿用 [几何比较稿](GEOMETRIC_COMPARISON.md) 的记号。设 \(k\subset\mathbb R\) 为数域，质量平方 \(\mu_i>0\)、\(s<0\) 属于 \(k\)。记 \(X_s\) 为四线 banana 四次曲面的最小光滑 K3 解消，\(U_s\subset X_s\) 为其环面部分，\(D=X_s\setminus U_s\) 为二十条边界曲线。

取光滑射影三维模型 \(P\)，严格变换 \(X\) 与边界 \(Y\) 组成正常交叉除子，\(Z=X\cap Y\)。这里 \(X\) 允许是 \(X_s\) 的进一步边界吹起，且 \(P\setminus Y=T\)、\(X\setminus Z=U_s\)。定义

\[
 \mathsf E=H^{\mathrm N}_3(P\setminus X,Y\setminus Z),\qquad
 \mathsf K=(\mathsf E)_{\leq1},\qquad
 \mathsf A=\ker(\mathsf E\to\mathsf K).
\]

上一稿已经证明 \(R_H(\mathsf A)=A_s(3)\)，其中

\[
 A_s=H^2(X_s)/\langle[D_i]\rangle,\qquad
 \dim A_s=6.
\]

**主命题。** 存在下文显式构造的 Chow 投影子 \(\pi_{\mathrm{res}}\)，使其 Nori 实现 \(\mathsf H_{\mathrm{res}}\) 为 \(H^2_{\mathrm N}(X_s)\) 的秩六纯直和因子，并有实际 Nori 动机同构

\[
 \boxed{\mathsf H_{\mathrm{res}}(3)\xrightarrow{\sim}\mathsf A.}
\]

由此得到实际短正合列

\[
 0\longrightarrow\mathsf H_{\mathrm{res}}(3)
 \longrightarrow\mathsf E\longrightarrow\mathsf K\longrightarrow0.
\]

这里 \(\mathsf H_{\mathrm{res}}\) 指边界除子跨度的正交因子。如果 \(\rho(X_s)>16\)，它仍有秩六，但包含额外代数方向。仅在 \(\rho(X_s)=16\) 时，可将它识别为完整超越 K3 因子的 Nori 实现。

本批 [Picard 证书](PICARD_CERTIFICATE.md) 已对 (m,s)=((1,2,4,8),−1)、((1,3,7,15),−2) 及目标基点 ((1,2,4,8),−64) 证明 rho=16。因此在这些指定纤维上，纯核确为完整 T(X)(3)。后续 [LATTICE_SATURATION](LATTICE_SATURATION.md) 另证边界格本原性，得到完整 T=U⊕U(2)⊕A2(−2)，det48；六循环的几何标记仍未认证。

## 相对对偶性的动机版本

以下一般工具给出所需的动机箭头。Ivorra–Morel 的 Proposition 2.11 将基域上的 perverse Nori 范畴与通常 Nori 范畴识别；Tubach 的 Theorem 4.5 和 Corollary 4.8 提供与 Betti、Hodge 实现相容的六运算。采用的版本分别为 [arXiv 1901.02096v3](https://arxiv.org/pdf/1901.02096v3) 与 [arXiv 2309.11999v4](https://arxiv.org/html/2309.11999v4)。只使用这些无条件结论，不假定 Voevodsky 动机上存在猜想的 motivic t-structure。

**引理一。** 设 \(P\) 为光滑射影 n 维簇，\(X,Y\) 无共同不可约分量，\(X\cup Y\) 正常交叉。则在 Nori 动机范畴中有与实现相容的自然同构

\[
 H^i_{\mathrm N}(P\setminus Y,X\setminus Z)(n)
 \simeq H^{\mathrm N}_{2n-i}(P\setminus X,Y\setminus Z).
\]

**证明。** 记 \(W=P\setminus(X\cup Y)\)，并记开浸入

\[
 j:P\setminus Y\hookrightarrow P,\quad
 k:P\setminus X\hookrightarrow P,\quad
 v:W\hookrightarrow P\setminus Y,\quad
 u:W\hookrightarrow P\setminus X.
\]

在 Nori 导出范畴内，由开基变换和伴随构造交换态射

\[
 j_!v_*\mathbf1_W\longrightarrow k_*u_!\mathbf1_W.
\]

其 Betti 实现由正常交叉的局部乘积模型检验为同构。具体地，在 \(X\cap Y\) 附近，把局部坐标分成控制 X 与 Y 的两组，两边分别由同一个外张量积计算；离开交集时交换是直接的。这个拓扑论证见 [Bloch–Kerr–Vanhove，Lemma 6.1.1 至 6.1.4](https://arxiv.org/html/1406.2664v3#S6.SS1)。Nori 的 Betti 实现正合且忠实，其导出函子保守，因此已经构造的交换态射本身为同构。

令 \(a:P\to\operatorname{Spec}k\)。物理相对上同调复形为 \(a_*k_*u_!\mathbf1_W\)，交换 X、Y 后的复形为 \(a_*j_*v_!\mathbf1_W\)。相对上同调与这些复形的识别来自局部化三角及闭基变换；相应交换在同一个局部乘积模型中检验。利用 \(a_!=a_*\)、Verdier 对偶交换星号与惊叹号，以及 \(\mathbb D\mathbf1_W=\mathbf1_W(n)[2n]\)，得到

\[
 \mathbb D(a_*k_*u_!\mathbf1_W)
 \simeq (a_*j_*v_!\mathbf1_W)(n)[2n].
\]

取次数即为所需同构。所有箭头来自六运算、伴随和局部化；没有从一个任意的 Hodge 同构反推其动机来源。证毕。

应用 n=i=3，得到实际比较

\[
 d:\mathsf E\xrightarrow{\sim}H^3_{\mathrm N}(T,U_s)(3).
\]

本引理处理完整相对复形。不能直接套用只针对 very good pair 的对偶引理，因为本模型还有非零的 \(H^2(T,U_s)\)。这里使用的是基域上通常 Nori 对象的比较和相对六运算，不要求证明所有不同构造的相对导出实现函子全局等价。

## 边界对应与 K3 投影子

按顺序记全部边界曲线为

\[
 (V_1,\ldots,V_4,E_{12},\ldots,E_{34},
 C_{12},\ldots,C_{34},F_1,\ldots,F_4),
\]

其中二元下标取字典序。从中删去 \(E_{12},F_2,C_{12},F_3\)，得到十六条曲线 \(B_1,\ldots,B_{16}\)。令

\[
 G=(B_i\cdot B_j)_{i,j=1}^{16},\qquad \det G=-48.
\]

边界的四个主除子关系及这个非退化子式共同证明：\(B_i\) 的类是全部边界除子有理跨度的一组基。四个关系来自单位 \(a/d,b/d,c/d,M/d\)，其中 \(M=\sum\mu_i\alpha_i\)。因此这里的秩十六同时有实际除子关系与交叉矩阵下界支持。

定义代数对应

\[
 \pi_D=\sum_{i,j=1}^{16}(G^{-1})_{ij}\,[B_i\times B_j]
 \in\operatorname{CH}^2(X_s\times X_s)_{\mathbb Q}.
\]

采用对应作用 \(\Gamma_*(x)=p_{2*}(p_1^*x\cdot\Gamma)\)。其作用是

\[
 (\pi_D)_*x=\sum_{i,j}(G^{-1})_{ij}(x\cdot B_i)B_j.
\]

由乘积对应的复合公式，\(\pi_D^2=\pi_D\) 在 Chow 群中成立：复合系数即 \(G^{-1}GG^{-1}=G^{-1}\)。它在每个边界类上为恒等，且对称于交换两个因子。

选取一个次数为一的有理零圈 \(o\in\operatorname{CH}^2(X_s)_{\mathbb Q}\)。例如任一闭点除以其次数即可，无需额外假定有 k 有理点。令

\[
 \pi_0=o\times X_s,\quad \pi_4=X_s\times o,\quad
 \pi_2=\Delta_{X_s}-\pi_0-\pi_4,\quad
 \pi_{\mathrm{res}}=\pi_2-\pi_D.
\]

乘积对应的维数计算给出 \(\pi_0,\pi_4,\pi_D\) 两两正交，故 \(\pi_{\mathrm{res}}\) 也是 Chow 投影子。K3 的一、三次上同调为零，所以其 Nori 实现仅位于二次，定义

\[
 \mathsf H_{\mathrm{res}}
 =\operatorname{im}\big((1-\pi_D)_*:H^2_{\mathrm N}(X_s)
 \to H^2_{\mathrm N}(X_s)\big).
\]

该对象纯权重二，Hodge 型为 (2,0)、(1,1)、(0,2)，重数分别为 1、4、1。其 Betti 实现是边界跨度的正交补，并通过商映射同构于 \(A_s\)。

这是一个具体 Chow 对应在 Nori 范畴中的实现。后面证明的是它与混合对象之核的 Nori 同构；不声称在纯 Chow 范畴中定义或分类整个混合扩张。

## 进入物理对象的实际态射

令 \(r\) 为开浸入 \(U_s\hookrightarrow X_s\) 的限制，\(\delta\) 为相对对 \((T,U_s)\) 的连接态射。定义

\[
 \iota:\mathsf H_{\mathrm{res}}(3)
 \hookrightarrow H^2_{\mathrm N}(X_s)(3)
 \xrightarrow{r}H^2_{\mathrm N}(U_s)(3)
 \xrightarrow{\delta}H^3_{\mathrm N}(T,U_s)(3)
 \xrightarrow{d^{-1}}\mathsf E.
\]

这是代数限制、相对长正合列和引理一的复合，因而已是动机态射。

**引理二。** \(\iota\) 为单射，其 Hodge 像恰为 \(W_{-4}E_s\)。

**证明。** 边界 Gysin 序列给出

\[
 \ker(H^2(X_s)\xrightarrow r H^2(U_s))
 =\langle[D_i]\rangle.
\]

因此 r 在边界正交补上单射，其像为开曲面的权重二部分。另一方面，相对长正合列给出

\[
 \ker\delta=\operatorname{im}(H^2(T)\to H^2(U_s)).
\]

\(H^2(T)=\mathbb Q(-2)^3\) 为纯权重四，所以其像也是纯权重四。它与 r 的纯权重二像交为零。这里使用 MHS 态射对权过滤的严格性，而不是假定开曲面已经有理分裂。于是 \(\delta r\) 在六维正交补上单射。扭转 (3) 后，其像是 E 中六维的权重负四部分。Betti 实现忠实且正合，故 \(\iota\) 在 Nori 范畴中也单射。证毕。

令 \(q:\mathsf E\to\mathsf K\) 为最大一动机商。由于 \(R_H(\mathsf K)\) 权重只为 -2、0，纯权重 -4 的像在其中为零。因此 \(R_H(q\iota)=0\)，由忠实性有 \(q\iota=0\)。这给出唯一实际态射

\[
 \widetilde\iota:\mathsf H_{\mathrm{res}}(3)\longrightarrow\mathsf A.
\]

它的 Betti 实现是两个六维空间间的同构，故其核与余核的 Betti 实现皆零。正合忠实实现使核与余核本身为零，\(\widetilde\iota\) 为同构。主命题证毕。

## 实际环境扩张的 motivic Yoneda 消失

这个构造还把环境三步对象的 Yoneda 消失落实到 Nori 范畴。将实际一动机写成 \([L_s\to T_s]\)，记其环面子对象为 \(\mathsf B_s=[0\to T_s]\)，格商为 \(\mathsf L_s=[L_s\to0]\)。于是

\[
 0\to\mathsf B_s\to\mathsf K\to\mathsf L_s\to0.
\]

定义 \(\mathsf N=\mathsf E\times_{\mathsf K}\mathsf B_s\)，得到实际短正合列

\[
 0\to\mathsf H_{\mathrm{res}}(3)\to\mathsf N\to\mathsf B_s\to0,
 \qquad 0\to\mathsf N\to\mathsf E\to\mathsf L_s\to0.
\]

若第一列的类为 \(\eta_{\mathrm{env}}\)，而 \(\mathsf K\) 的类为 \(\xi_{\mathrm{env}}\)，则

\[
 \eta_{\mathrm{env}}\circ\xi_{\mathrm{env}}=0
 \quad\text{于 }\operatorname{Ext}^2_{\mathrm{Nori}}
 (\mathsf L_s,\mathsf H_{\mathrm{res}}(3)).
\]

证明是扩张长正合列：\(\xi_{\mathrm{env}}\) 已由 \(\mathsf E\) 提升到 \(\operatorname{Ext}^1(\mathsf L_s,\mathsf N)\)，故连接障碍为零。这不宣称整个 Ext² 群为零。后续 [KUMMER_BOUNDARY.md](KUMMER_BOUNDARY.md) 计算实际质量点并构造活动拉回，同一论证再给出指定活动类 eta_kappa、xi_kappa 的 Yoneda 消失。

## 比较解决了什么

现在可以在 Nori 范畴中把纯层写成指定的 K3 边界正交因子，并把它与物理相对对象的嵌入一并给出。这个识别同时保留 Betti、Hodge 和 de Rham 实现的自然比较。

它没有给出比较矩阵的数值条目或特殊函数表示，也没有判断所选点的 Picard 数。质量比 Kummer 比较由独立的 [边界稿](KUMMER_BOUNDARY.md) 计算；不能用本稿的纯层识别代替那里的单位杯积与扩张比较。

对后续论文，假设 11.2 应拆成实际动机态射、边界点比较和待认证的具体周期数据，逐项指向对应证明。活动最小性与完整群等号不能由本命题推出。

## 对应系数的精确核对

`data/boundary_projector.json` 保存十六曲线基、G、G 的有理逆、全部二十曲线的基坐标，以及每个基关系作为四个主除子关系之组合的系数。`code/verify_boundary_projector.py` 重新构造交叉表并验证这些恒等式。

文件中的二十维矩阵作用于边界曲线的生成元空间。其四维核是除子关系空间，不是六维 K3 核。六维来自实际 \(H^2(X_s)\) 的维数二十二减去除子跨度十六。有限代数核对不承担六运算、动机对偶性或 Picard 数的证明。
