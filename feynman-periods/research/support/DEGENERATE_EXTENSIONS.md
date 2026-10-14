# 下层扩张退化、完整纯核与等质量 banana

本稿继续纯理论主线，补齐 THEORY_MANUSCRIPT 第7节明确留出的 r=0 情形，并把 Picard 跳跃时的纯核处理落到实际代数对应。最后应用于等质量四传播子 banana：借助已知一般 Picard 数十九及相对 Cartan 公式，得到很一般欧氏参数的秩四最小对象、深度一和七维实际群。本批没有新增验证程序、周期运输或数值积分。

## 1. 假设与“扩张为零”的含义

在有理中性 Tannakian 范畴 C 中，考虑实际三步对象
\[
 0\subset A\subset N\subset E,\qquad
 N/A=B\otimes V,\quad E/N=\mathbf1,
 \qquad V=\mathbb Q^\kappa.                              \tag{1.1}
\]
A 为简单纯对象，维数 n>1、权重−4；B=Q(1)的权重为−2。分级群约化，权重余特征和相容的严格 Hodge 实现均存在。记
\[
 \eta\in\operatorname{Ext}^1_C(B\otimes V,A),\qquad
 \xi\in\operatorname{Ext}^1_C(\mathbf1,B\otimes V).
\]
当 κ>0 时，ξ 的 κ 个坐标在 Ext¹_C(1,B) 中有理独立。

本稿先假定 \(\eta=0\) **于所讨论的实际范畴 C 中**。这等价于前稿下层 D_A 秩 r=0。若只知道 Hodge regulator 消失，可以在 Hodge 范畴应用本稿，不能未经证明把结论套给 Nori 群。有限精度的数值零更不满足这里的假设。

物理框架为顶层非零的 \(\sigma\in E_B\) 和
\(\omega\in F^3E_{\rm dR}^{\vee}\)，且 \(\omega|A\ne0\)。最小对象始终指原有理范畴中的实际子商；其上同调形式周期对象是本文同调对象的对偶。

## 2. r=0时，独立的顶到底扩张是规范对象

由于 Hom_C(B,A)=0，分裂的下层序列有唯一截面。记其像为
\(\widetilde B=B\otimes V\subset N\)。因此
\[
               N=A\oplus\widetilde B                    \tag{2.1}
\]
是由实际态射确定的分解，不需选择任意向量空间分裂。

定义
\[
 Z=E/\widetilde B,\qquad K=E/A,
 \quad\zeta=[0\to A\to Z\to\mathbf1\to0].              \tag{2.2}
\]
商映射给出规范同构
\[
                    E\simeq Z\times_{\mathbf1}K.        \tag{2.3}
\]
确实，两者到顶层的核都识别为 A⊕B⊗V，态射在核与商上为恒等。等价地，Ext¹ 对有限直和的可加性把 E 的顶层扩张写为(ζ,ξ)。

这明确了前稿“还需研究顶到底扩张”中的对象：ζ在η=0时是规范的，不是依赖某个非唯一低层分裂的待定常数。

若 κ=0，则 E=Z，本节仍成立，不存在任何 Kummer 中间层。

## 3. 两支实际群与最小对象

设 ε=1 当 ζ≠0，ε=0 当 ζ=0。若 κ>0，记 R=G(A⊕B)；若 κ=0，记 R=G(A)。不要在没有中间层时人为把一个不在 E 生成范畴中的 Tate 根加入分级群。

**定理3.1。** 在(1.1)的假设及η=0下，
\[
                U(E)\simeq\mathbb V(A)^\varepsilon
                              \times\mathbb V(B\otimes V),
 \qquad \dim U(E)=\varepsilon n+\kappa.                 \tag{3.1}
\]
κ=0时删去第二因子。U(E)为交换向量群。

**证明。** 非分裂二步扩张 Z 的幺幂像是 A 的非零 R 子表示；若该像为零，Z 的表示降到约化群，半单性使 ζ=0。因此 ζ≠0时，A的简单性迫使该像为整个A；ζ=0时该像为零。上层独立性同样使 K 的幺幂像为全部 B⊗V。A与B的权重不同，两个像没有同构简单因子；联合像不能是非平凡对角子表示，故为直和。式(2.3)给忠实表示
\[
 u(z,y)=\begin{pmatrix}I_A&0&z\\0&I_V&y\\0&0&1\end{pmatrix}.
                                                               \tag{3.2}
\]
这些矩阵相乘只把 z、y 相加，没有 xy' 项。证毕。

### 3.1 ζ非零：最小秩n+1，深度一

全部 z、y 平移把任意顶层非零 σ 生成整个 E。若子对象 J 被ω消去，则 J∩A=0；若它在顶层非零，生成性将迫使 J=E，矛盾。所以 J⊂N，并因(2.1)及 Hom(B,A)=0而包含于 \(\widetilde B\)。反之，\(F^3(B\otimes V)^\vee=0\)，ω确实消去它。故
\[
      E_{\min}=Z,
      \qquad \operatorname{rank}E_{\min}=n+1,
      \qquad U(E_{\min})=\mathbb V(A).                  \tag{3.3}
\]
在 \(Z^\vee\) 上，幺幂作用为
\((\alpha,\gamma)\mapsto(\alpha,\gamma-\alpha z)\)。因为α≠0，作用非平凡；任意两次约化作用为零。因此物理形式周期深度恰为一，而非仅给出上界。

与权重环面交换后，End(E)的元素是diag(a,h,λ)。与全部 z 交换给a=λI_A；有中间层时，与全部 y 交换再给h=λI_V。于是
\[
                    \operatorname{End}(E)=\mathbb Q,
 \qquad \operatorname{End}(Z)=\mathbb Q.                \tag{3.4}
\]
因此η=0并不强迫整个E分解；独立的ζ足以保持不可分解性。

### 3.2 ζ为零：物理周期由纯核承载

Hom(1,A)=0使 Z 的截面唯一，故(2.3)给
\[
                         E\simeq A\oplus K.             \tag{3.5}
\]
\(F^3K^\vee=0\)，物理读出只作用于 A。写σ=σ_A+σ_K。若σ_A≠0，A简单使σ_A生成A，且ω|A≠0排除其非零不可见子对象。因此
\[
                 E_{\min}=A,\qquad
                 \operatorname{rank}E_{\min}=n,
                 \operatorname{depth}I^{\mathfrak m}=0. \tag{3.6}
\]
欧氏积分I>0保证σ_A≠0。若讨论任意抽象框架且σ_A=0，形式系数就是零，最小对象为零；不能省略这个条件。

跨 A 与 K 的 Hom 均为零，且活动 K 的自同态只有标量，故
\[
                      \operatorname{End}(E)=D_A\times\mathbb Q.
                                                               \tag{3.7}
\]
即使κ>0且环境对象仍有κ维幺幂群，物理周期也可能完全看不见它。深度零表示纯对象承载，不表示周期是代数数。

### 3.3 三种深度与最小秩

合并 THEORY_MANUSCRIPT 第7节的r>0结果，在纯核简单、上层活动及非零物理框架范围内：

| 扩张情形 | 最小同调对象的秩 | 物理幺幂深度 | 原E的自同态 |
|---|---:|---:|---|
| r>0，d=dim_Q ker L_eta | n+κ+1−d | 2 | Q |
| r=0，ζ≠0 | n+1 | 1 | Q |
| r=0，ζ=0，σ_A≠0 | n | 0 | D_A×Q |

这给前稿一般秩公式中δ分支一个实际范畴解释。深度不是圈数，单纯看E是否不可分解也分不出第一、二行。

## 4. 删光Tate边界时，群核可能含有限部分

第3节的群维数不足以确定 G(E)→G(E_min) 的全部核。以下假定κ>0，并定义
\[
 R_A=G(A),\qquad
 T_0=\ker\{R=G(A\oplus B)\longrightarrow R_A\}.        \tag{4.1}
\]
T₀在A上平凡、在B上按标量作用。由于R对A⊕B忠实，它嵌入Gm且在R中中心。

**命题4.1。** 在第3节的两个非零物理分支中，都有
\[
 1\to\mathbb V(B\otimes V)\to
 \ker\{G(E)\to G(E_{\min})\}\to T_0\to1.              \tag{4.2}
\]
选Levi分裂后，核写成 \(\mathbb V(B\otimes V)\rtimes T_0\)，T₀按B的角色作用。

**证明。** ζ≠0时，G(E)的矩阵由(H,t)∈R和自由z∈A、y∈B⊗V组成，而G(Z)只保留H、z。其核要求H=I、z=0，尚有t∈T₀、任意y。ζ=0时从开始就没有z，映到G(A)给完全相同的核。矩阵乘法为(y,t)(y',t')=(y+ty',tt')。证毕。

在banana的 \(A=T(X)(3)\) 上，非退化极化取值于 \(B^4=\mathbb Q(4)\)。因此T₀的标量t满足t⁴=1，即
\[
                              T_0\subseteq\mu_4.        \tag{4.3}
\]
这里是有限群概形，不能只按有理点数计数。κ=0时B原本就不必属于E生成范畴，此时不额外加(4.2)的有限核。

一个明确的表示论例子是前稿的分级形式R=SO(q)×Gm、\(\dim A=6\)，作用A=t²g、B=t。若t²g=I，则g=t⁻²I；正交性给t⁴=1，六维的行列式条件不再增加限制。故T₀=μ₄。若R=O(q)×Gm，同样如此。对三维SO(q)情形，行列式另给t⁶=1，与t⁴=1合起来只有μ₂。

这些是所述分级表示的核计算，不宣称已找到具有η=0的新的指定非等质量纤维。它说明第7节r>0商定理为什么要使用κ−d>0：有剩余B时不会遗失Tate根，删掉全部B时就不能跳过(4.1)。

## 5. Picard跳跃时的实际超越投影

此前边界投影仅剔除了固定十六维除子空间。因此当ρ(X)>16时，秩六的H_res含额外代数方向，不能把整个H_res当成简单纯核。本节给出不依赖新数值的处理。

令X/k为当前光滑K3。取有限Galois扩域L/k，使几何Néron–Severi群的一组基D₁,…,D_ρ都由L上的除子给出。其交叉矩阵G_NS非退化。设
\[
 \pi_{\rm NS,L}=\sum_{i,j}(G_{\rm NS}^{-1})_{ij}
                        [D_i\times D_j].                \tag{5.1}
\]
组合律给 \(\pi_{\rm NS,L}^2=\pi_{\rm NS,L}\)。K3的Pic⁰为零，因此Galois对NS基的变化也给对应除子类的变化；(5.1)是同一非退化空间的正交投影，故在有理Chow群中Galois不变。

取corestriction并除以[L:k]，得到k上的对应π_NS，其限制为π_NS,L。限制在有理Chow群中单射，因为corestriction与限制复合为[L:k]；因此投影子的恒等式也在k上成立。这一步解决下降，不能只说“所有Hodge(1,1)类都是代数的”便省略实际对应。

用相容的次数零、四投影定义π₂。令π_D为已有边界投影，边界空间非退化且包含在NS中。于是
\[
 \pi_{\rm alg}=\pi_{\rm NS}-\pi_D,
 \qquad \pi_T=\pi_2-\pi_{\rm NS},\qquad
 \pi_{\rm res}=\pi_T+\pi_{\rm alg}.                     \tag{5.2}
\]
这些投影两两正交，等式由嵌套非退化空间的正交投影及相同的交叉矩阵组合律验证。Nori实现给实际直和
\[
 A_{\rm res}=H_{\rm res}(3)=A_T\oplus C_{\rm alg},
 \quad A_T=T(X)(3),\quad\dim A_T=22-\rho(X),\quad
 \dim C_{\rm alg}=\rho(X)-16.                            \tag{5.3}
\]
C_alg是Artin表示扭转Q(2)，一般不在k上分裂成平凡Tate直和。

物理ω属于F³E∨，而C_alg∨的Hodge型全为(2,2)，所以ω消去C_alg。它通过实际商
\[
                         E_T=E/C_{\rm alg}              \tag{5.4}
\]
下降。链取商像后矩阵系数不变，故(5.4)与原E具有相同的最小物理对象。其纯核现在为实际简单对象A_T：超越K3 Hodge结构简单，忠实Hodge实现遂排除Nori中的真子对象。

因此前述群与最小性判据可用于E_T，只需把n改成22−ρ，并把η推到A_T；原始质量κ不变。这里得到的是点上的构造与条件判据，没有自动计算ρ、找到全部除子或求出各投影regulator。原κ、ρ16、指定−64的结论均未被改动。

## 6. 等质量族：不靠数值证明很一般的顶层非零

在正质量区域，κ=0等价于四个平方质量相等：正实单位根只有1。通过共同尺度归一化，取μ₁=…=μ₄=1、参数t=s/m²<0。定义函数为
\[
       f=(1+a+b+c)(1+a^{-1}+b^{-1}+c^{-1}),\qquad f=t.   \tag{6.1}
\]
物理参数积分只是一般尺度下的非零代数倍，尺度不改变本节的最小对象秩。

本节唯一新增的族几何外部输入是：这个等质量族的一般Picard数为十九，并有十九维由除子生成的局部子系统。采用BKV §3.2陈述的Verrill族结果；其正文与(1.2)给出的正坐标方程相符，§3.2的负坐标通过a=−x等代换得到。[1] 本稿未独立重证该Picard数，也不把它当作所有特殊纤维的逐点结论。

因此在去掉额外Picard跳跃点后，(5.3)有dimC_alg=3、dimA_T=3。持久的除子类给平坦的代数子系统；在局部标记图表或所需有限覆盖上，(5.4)成为实际相对变分的商，其序列为
\[
                 0\to A_T\to Z\to\mathbf1\to0.         \tag{6.2}
\]
无需把旧秩八比较矩阵补零或删行来构造它。

**定理6.1。** 对很一般实欧氏参数t<0，(6.2)在Hodge实现中非分裂，因而在实际Nori范畴中非分裂。

**证明。** 取相对形式e=(dloga∧dlogb∧dlogc,0)，并按顶层生成元归一化。扭转后，F⁰A_T=0，因此e的像是顶层1的唯一Hodge提升e_F。已推导的相对Cartan公式在归一化参数下为
\[
                       \nabla_{\partial_t}e_F
                               =\delta\omega_t,         \tag{6.3}
\]
其中ω_t是非零全纯K3二形式。额外代数类均为(1,1)，投到T后ω_t仍非零；δ把纯核单射到相对对象。因此(6.3)在光滑局部图表中非零，e_F不可能在任何非空开集成为平坦截面。

在单连通图表固定有理Betti标记。具有顶层1的有理平坦提升q只有可数多个。某纤维Hodge分裂当且仅当某个q在该点等于e_F。对固定q，q−e_F是全纯向量函数，且非恒零；一复变量中其共同零点离散。对可数个q和图表取并，再并入额外Picard跳跃集合，仍是可数例外集。故其外的实t<0非分裂。

如果实际Nori扩张分裂，其Hodge实现必分裂；逆否命题得到同一点的Nori非零。并未使用“Hodge零推出Nori零”。证毕。

该证明没有把某个点的导数非零直接当作该点扩张非零；导数排除开集分裂，可数平坦提升才给很一般纤维结论。任何指定代数t仍须单独判断。

结合第3节及物理正链，得到
\[
           \boxed{\operatorname{rank}M(I^{\mathfrak m}_{\rm equal})=4,
                  \qquad\operatorname{depth}I^{\mathfrak m}_{\rm equal}=1}
                                                               \tag{6.4}
\]
对同一很一般欧氏集合成立。这里从原秩七活动对象先商掉三个代数方向，而不是错误套用秩六简单核。

其Hodge过滤维数为(\(\dim F^0,\dim F^{-1},\dim F^{-2},\dim F^{-3}\))=(1,2,3,4)：纯核的三个型分别为(−1,−3)、(−2,−2)、(−3,−1)，顶层为(0,0)。

## 7. 等质量最小对象的实际群

当ρ=19时，T的有理秩为三。Zarhin定理给End_HdgT为全实域或CM域。[2] 其次数整除三，所以CM情形不可能。若全实域次数为三，T在每个实嵌入的特征空间复维一；但H²⁰与H⁰²在同一实嵌入下的特征值相同，且是两个独立方向，矛盾。因此
\[
                         \operatorname{End}_{\rm Hdg}T=\mathbb Q.
                                                               \tag{7.1}
\]
相应权重零Hodge群为SO(q)，在A_T=T(3)上连同权重标量得到GO(q)。这是三维的特殊便利：若HᵗqH=λq，令
\[
                  c=\det H/\lambda,\qquad g=c^{-1}H.
\]
由(detH)²=λ³可得c²=λ、gᵗqg=q、detg=1，故
\[
                         GO(q)\simeq\mathbb G_m\times SO(q)
                                                               \tag{7.2}
\]
为连通的四维群。没有在这个没有B中间层的对象上额外指定一个Tate平方根。

杯积给实际极化A_T⊗A_T→Q(4)；作为该张量的像，Q(4)已属于A_T生成范畴。故实际Nori纯群包含Hodge下界GO(q)，又包含于相同的相似群上界，二者相等。这里不需要识别T(1)单独的可能Artin行列式角色。

对第6节的非分裂点，Hodge幺幂核已是整个A_T，Nori核也只能介于它和同一上界A_T之间。于是选Levi分裂后
\[
 \boxed{G_H(Z)=G_N(Z)=
 \left\{\begin{pmatrix}H&z\\0&1\end{pmatrix}:
                 H\in GO(q),\ z\in\mathbb Q^3\right\},
             \quad\dim G=7.}                            \tag{7.3}
\]
群按任意系数代数解释；等式是在相容Betti中性化中的表示群。三维幺幂根基与四维约化群均已由理论确定，两群连通。此结论限定同一很一般等质量欧氏纤维；它不计算其绝对比较矩阵或特定t值。

若C_Z的块为(P_A,f;0,1)，则实际结构coaction为
\[
                         \Delta f=P_A\otimes z+f\otimes1.
                                                               \tag{7.4}
\]
没有中间质量对数块，也没有前稿的p⊗y项。按物理对偶读出可写四项基展开；这不是新的特殊函数求值公式。

## 8. 对主稿及后续工作的影响

1. r=0与κ=0的群、最小性和深度现在有完整条件性分类；不再仅留下“另外分析ζ”的空缺。
2. 删光中间层时必须同时检查Tate根是否仍由纯核生成；有限核T₀不由群维数检验排除。
3. 任意光滑纤维上，额外代数方向可由下降后的实际NS投影处理；其是否在k上分裂为Tate不是必要条件。
4. 等质量很一般欧氏纤维提供实际秩四/深度一/七维群应用。已有一般Picard十九归于文献，非分裂和最小性通过本文的实际相对框架推导，不宣称该经典族研究的首创性。
5. 指定非等质量m=(1,2,4,8)、s=−64仍为原来的秩八、深度二与29维群；它的证书及证明没有修改。

尚未给出的是所有特殊代数参数的ζ非零判定、完整NS的有效求取、一般图观察到几何态射的桥梁以及整体论文剩余附篇。很一般等质量结论不自动覆盖某个CM纤维。用户要求的纯理论优先继续执行。

## 文献与证明责任

[1] Bloch–Kerr–Vanhove, *A Feynman integral via higher normal functions*, arXiv:1406.2664v3，(1.2)、§3.2以及§6.2：
https://arxiv.org/html/1406.2664v3#S3.SS2
https://arxiv.org/html/1406.2664v3#S6.SS2

本批核对作者HTML：等质量方程、十九维除子子系统、一般ρ19及Hodge唯一提升约定。Verrill原论文未在本批直接读取，故明确把族几何作为BKV所陈述的已知输入；NS下降、r=0群分类与本文最小化推论给出各自论证。

[2] Charles, *On the Picard number of K3 surfaces over number fields*, §2.2（Zarhin结果）：
https://www.imo.universite-paris-saclay.fr/~francois.charles/Picard.pdf

[3] Brown, *Notes on motivic periods*, §2.4–2.5及§3.8：
https://arxiv.org/html/1512.06410#S2.SS5

本批原始HTML核对了coradical的C₀不变量定义与深度零/一的范围；具体(3.2)作用和最小商在文中直接证明。Eskandari v3的本次网页读取超时，未把失败读取列为新的原文核对；前稿已读一般框架不作为本批新定理的替代证明。

内部依赖：THEORY_MANUSCRIPT、REGULATOR_RANK_AND_GROUP、VISIBILITY_AND_CORADICAL、MINIMALITY、K3_KERNEL、RELATIVE_TOP_COMPARISON、MIXED_PARAMETER_CONNECTION。第6节Cartan恒等式沿用实际相对复形证明；全部零/非零及群结论按所用范畴分别陈述。
