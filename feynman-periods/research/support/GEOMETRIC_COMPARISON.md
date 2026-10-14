# 四线 banana 的物理相对动机与边界一动机

> 2026-10-14专项状态：旧假设11.2的逐项依赖与实际物理双标架等式现见HYPOTHESIS_11_2_AUDIT.tex（PDF同名）。本稿后部的部分“待证”措辞记录撰写时状态；固定下层非零、P8及固定群的后续证据分别在LOWER_REGULATOR、RELATIVE_TOP_COMPARISON与FIXED_PURE_GROUP，不由本稿或新补篇重新数值认证。新补篇未将结论推广到所有参数，也尚未回写统一全文。

本稿把相对 Hodge 对象落实为有效 Nori 同调动机，并证明它有一个实际的一动机商。结论适用于正代数质量、非奇异欧氏参数及 HODGE.md 中的二十曲线模型。[纯核比较稿](K3_KERNEL.md) 进一步识别边界正交 K3 因子；[边界 Kummer 比较稿](KUMMER_BOUNDARY.md) 用实际留数与同伦纤维计算质量点。完整群识别仍待证明。下述结果是已有一般定理在本模型上的推导，不宣称一般定理本身或这些结论的文献新颖性。

## 模型与约定

取数域 \(k\) 的一个实嵌入，\(\mu_i=m_i^2>0\) 和 \(s<0\) 均在 \(k\) 中。令

\[
 U_G=abc+abd+acd+bcd,\qquad
 F_s=(\mu_1a+\mu_2b+\mu_3c+\mu_4d)U_G-sabcd.
\]

记 T=P^3\{abcd=0}、U_s={F_s=0}∩T。本文的积分是去掉约定性动量测度因子的参数积分

\[
 I(s)=\int_{\sigma}\omega_s,\qquad
 \omega_s=\Omega/F_s,
 \quad \sigma=\{[a:b:c:d]:a,b,c,d\geq0\}.
\]

定向选择使积分为正。采用 wt Q(r)=-2r，定义

\[
 E_s=H^3(T,U_s;\mathbb Q)(3),\qquad
 A_s=H^2(X_s;\mathbb Q)/\langle[D_i]\rangle,
\]

其中 X_s 是最小光滑 K3 解消，D_i 是二十个边界分量。[Hodge 推导](HODGE.md) 给出 rank E_s=15、W_{-4}E_s=A_s(3)，以及 E_s/A_s(3) 的权重层 Q(1)^8、Q(0)。这些秩依赖指定的边界模型；不要求特殊纤维 Picard 数恰为十六。

## 正 Schwinger 链的显式提升

**命题一。** 在上述参数范围，正单纯形在一个定义于 k 的光滑射影模型 P 中提升为闭的相对链，避开曲面的严格变换 X，边界落在坐标边界的总变换 Y 中。可以取 X∪Y 为正常交叉除子，P\Y=T、X\(X∩Y)=U_s。

**证明。** 对任一排列 p=(p_0,p_1,p_2,p_3)，在相应 Hepp 扇区取

\[
 \alpha_{p_0}=1,\quad \alpha_{p_1}=t_1,\quad
 \alpha_{p_2}=t_1t_2,\quad \alpha_{p_3}=t_1t_2t_3,
 \qquad 0\leq t_i\leq1.
\]

令

\[
 q=t_1t_2t_3,\quad V=1+t_3+t_2t_3+q,
 \quad M_p=\mu_{p_0}+\mu_{p_1}t_1+\mu_{p_2}t_1t_2+\mu_{p_3}q.
\]

直接代入得

\[
 F_s=t_1^2t_2 D_p(s,t),\quad
 D_p=M_pV-sq,\quad
 \left|\frac{\partial(\alpha_{p_1},\alpha_{p_2},\alpha_{p_3})}
 {\partial(t_1,t_2,t_3)}\right|=t_1^2t_2.
\]

因此在一致定向下

\[
 \omega_s=\frac{dt_1dt_2dt_3}{D_p},\qquad
 D_p\geq\mu_{p_0}>0\quad\text{于整个闭立方体}.
\]

吹起四个坐标顶点及六条坐标线后，这二十四个单项式坐标图覆盖提升链。t_i=1 的面是扇区之间的接缝，在链之和中抵消；t_i=0 的外边界落在 Y。严格变换由 D_p=0 给出，故与链不交。

几何模型也可先吹起四个顶点和六个坐标边上的节点，再吹起坐标线的严格变换。这六个额外节点满足相应的 M=0；正质量使它们远离正单纯形。因此该模型在提升链附近与上述扇区模型一致。

这里的节点条件可以直接检验。坐标顶点处的二次项正比于 bc+bd+cd，非退化。在边 ab 上的额外节点附近，以 x=M、c、d 为局部坐标，二次项正比于 x(c+d)-scd，因 s≠0 而非退化。边界其他点光滑。在 T 内，奇点条件迫使 alpha_i=epsilon_i C/m_i，从而 s=(sum epsilon_i m_i)^2，epsilon_i=±1；s<0 排除了这些点。因而只有上述十个普通二重点。

曲面节点解消后为光滑曲面。继续在边界内作相对解消，且保持已有正常交叉开集不变，即可使 X∪Y 正常交叉；这些后续中心与提升链不交。所有中心都在坐标边界上，故保持 T 和 U_s 不变。

形式也合法：沿光滑中心的吹起，若中心余维为 c、曲面重数为 m，拉回三形式沿新例外除子的阶数为 c-1-m。节点处 c=3、m=2；坐标线一般点处 c=2、m=1。两者均为零。曲面光滑之后的中心有 m≤1、c≥2，不产生额外极点。于是 omega_s 在 P\X 正则，并给出相对 de Rham 类。其在 Y 的限制为零，可在正常交叉边界的相对 de Rham 复形中取零边界分量。证毕。

这个证明用到了本模型的普通二重点几何及边界解消，不把任意费曼图的发散积分包括进来。

**解析推论。** 以上公式同时给出

\[
 I(s)=\sum_{p\in S_4}\int_{[0,1]^3}\frac{dt}{D_p(s,t)},
 \qquad 0<I(s)\leq6\sum_{i=1}^4\mu_i^{-1}.
\]

固定 s_0<0，若 |z|<mu_min=min_i mu_i，则

\[
 I(s_0+z)=\sum_{n=0}^{N}c_nz^n+R_N(z),\quad
 c_n=\sum_p\int_{[0,1]^3}\frac{q^n\,dt}{D_p(s_0,t)^{n+1}},
\]

\[
 |R_N(z)|\leq
 6\sum_{i=1}^4\frac{(|z|/\mu_i)^{N+1}}{\mu_i-|z|}.
\]

这是几何级数余项在闭立方体上的一致估计，不使用数值拟合。对 J(tau)=I(-tau)、tau>0，还有 (-1)^n J^{(n)}(tau)>0。该结果可检查欧氏物理分支及局部积分逼近，但不认证六维周期矩阵或全局单值运输。

## 有效相对动机与实际物理框架

令 Z=X∩Y，构造有理有效 Nori 同调动机

\[
 \mathsf E_s=H^{\mathrm N}_3(P\setminus X,Y\setminus Z;\mathbb Q).
\]

这里使用实际代数对的同调，因此有效性来自构造。不能仅凭 Hodge 型为负或写下一个 Tate twist 就断言某个动机有效。

正常交叉条件下的相对 Verdier 对偶性给出

\[
 R_{\mathrm{Hdg}}(\mathsf E_s)
 \simeq H^3(P\setminus Y,X\setminus Z;\mathbb Q)(3)
 =E_s.
\]

所用对偶形式见 [Bloch–Kerr–Vanhove，Lemma 6.1.4](https://arxiv.org/html/1406.2664v3#S6.SS1)。[纯核比较稿的引理一](K3_KERNEL.md) 用相对六运算构造了实际 Nori 对偶同构；其同构性通过已构造交换态射的 Betti 实现检验，不是从任意 Hodge 同构反推。

提升链给出 \(\sigma_B\in R_B(\mathsf E_s)\)。全部构造可在实嵌入下进行，链上的点被复共轭逐点固定，故几何复共轭 \(c_\infty\) 满足

\[
 c_\infty(\sigma_B)=\sigma_B.
\]

选择顶层框架使 sigma_B 映到 1∈Gr_0^W E_s。此归一化具有拓扑含义：正环面链与三个单位圆的乘积恰有一个横截交点，选择相容定向后交数为一。Tate 因子通过对偶性包含在识别中，不能再把 Q(3) 的共轭符号额外乘一次。

令 \(\mathsf M_s=\mathsf E_s^\vee\)，在允许 Tate 逆的 Nori 动机范畴中取对偶。参数积分是实际框架的周期

\[
 I^{\mathfrak m}(s)=[\mathsf M_s,\sigma_B,[\omega_s]]^{\mathfrak m},
 \qquad\operatorname{per} I^{\mathfrak m}(s)=I(s).
\]

由 F^0N_s=0，顶层的 Hodge 提升 e_F 唯一。对偶物理形式属于 F^3，故它与 e_F 的配对为零，得到 I(s)=<omega_s,sigma_B-e_F>。HODGE.md 定义 nu_E=[e_F-e_B]；当 e_B=sigma_B 时，这里配对的是 -nu_E。明确这个符号后，regulator 比较才有一致的归一化。

## 最大一动机商的构造

**所用一般定理。** 有理有效 Nori 动机到有效一动机的包含有左伴随 (−)_≤1，且与最大一动机型 Hodge 商相容。这是 [Ayoub–Barbieri-Viale，Nori 1-motives，arXiv 1206.5923v1，Theorem 8.12](https://arxiv.org/pdf/1206.5923) 的内容。这里固定版本，因为其他版本的定理编号不同。

**命题二。** 在上述几何和边界模型下，定义

\[
 \mathsf K_s=(\mathsf E_s)_{\leq1},\qquad
 \mathsf A_s=\ker(\mathsf E_s\longrightarrow\mathsf K_s).
\]

则存在实际的 Nori 动机短正合列

\[
 \boxed{0\longrightarrow\mathsf A_s\longrightarrow\mathsf E_s
 \longrightarrow\mathsf K_s\longrightarrow0},
 \qquad(\operatorname{rank}\mathsf A_s,\operatorname{rank}\mathsf E_s,
 \operatorname{rank}\mathsf K_s)=(6,15,9),
\]

并且

\[
 R_{\mathrm{Hdg}}(\mathsf A_s)=A_s(3),\qquad
 R_{\mathrm{Hdg}}(\mathsf K_s)=E_s/A_s(3).
\]

**证明。** 一动机型 Hodge 结构只含权重 -2、-1、0。严格性使 E_s 到任一此类对象的态射在纯权重 -4 子对象 A_s(3) 上为零。因此所有这些态射都通过 E_s/A_s(3) 因子化。该商只有型 (-1,-1)、(0,0)，自身就是一动机型，故它恰为最大一动机型 Hodge 商。应用上述左伴随的实现相容性，其伴随单位在 Hodge 实现上为满射。Hodge 实现忠实且正合，因而原态射也是满射；其核的实现与秩随即确定。证毕。

这个论证避免了“任意 Hodge 子对象都可提升为动机子对象”的无依据假设：实际态射先由一般定理构造，再计算其核。

\(\mathsf K_s\) 可由一个一动机 \([L_s\to T_s]\) 表示，其几何格秩为一、环面维数为八，无阿贝尔簇部分。单个 Hodge 实现不能排除算术 Artin 扭曲。[边界比较稿第一节](KUMMER_BOUNDARY.md) 通过定义在 k 上的交点、两次 Gysin 连接及实际 H^3_N(T)(3) 证明两层在 k 上分别为 Q(1)^8 和 Q(0)，从而补足分裂与下降。

\(\mathsf A_s\) 的 Hodge 实现为纯 K3 商 \(A_s(3)\)。[K3_KERNEL.md](K3_KERNEL.md) 已通过边界 Chow 投影子、开限制、相对连接和动机对偶构造实际态射，证明它同构于边界正交 K3 因子的 Nori 实现扭转 (3)。这不计算其特殊函数比较条目；完整超越因子的识别仍需 Picard 数十六。

## 实际质量比 Kummer 比较

[KUMMER_BOUNDARY.md](KUMMER_BOUNDARY.md) 在指定基域、余字符框架及同源意义下推导了

\[
 \mathsf K_s\simeq
 \left[\mathbb Z\longrightarrow\mathbb G_m^8\right]\otimes\mathbb Q,
 \quad
 1\longmapsto
 (\mu_1/\mu_4,\mu_2/\mu_4,\mu_3/\mu_4,1,1,1,1,1).
\]

顶层生成元由三个坐标单位的有序杯积固定；边界留数整数基为 P，一动机同调余字符基为 −P。这个整体负号由单单位相对周期校准。三十个交点中十二个旗标非零，旗标矩阵在基 P 中给出三个质量比及五个零单位类。

只计算 tame residue 不足以验证同构。比较稿第三、四节构造单位杯积到 RΓ_N(T)(3)[3] 的实际态射，再用同伦纤维引理证明符号扩张向物理相对扩张的推出关系。与局部化相容的实现把八个边界推出计算成 Kummer 单位类。这个构造补出了 [BKV 式 6.2.9](https://arxiv.org/html/1406.2664v3#S6.SS2) 在本模型中的尺度归一化；它没有固定纯 K3 方向的全部扩张周期。

质量比生成的最小子环面给出实际活动子一动机。对正实质量，令 kappa 为三个平方质量比在 k^×⊗Q 中的秩，则活动环面维数为 kappa，活动相对对象的秩为 7+kappa。对 m=(1,2,4,8)，得到秩八的实际实现上界；尚未证明其物理最小性。

## 物理链的活动因子化

**命题三。** 在正质量欧氏区域，取上述实际 Kummer 点生成的子环面及其在 E_s 中的拉回 E_act。则物理 Betti 链 sigma_B 属于 R_B(E_act)。

**证明。** Kummer 比较稿给出 E_s/E_act=Q(1)^(8-kappa)，且为 k 上分裂的实际商。几何复共轭在该商上为 -1，在 sigma_B 上为 +1。设 sigma_bar 是它在商中的像，则 sigma_bar=-sigma_bar；有理系数下 sigma_bar=0。证毕。

此处同时使用真实链的实性、相对对偶和实际分裂 Tate 商。抽象 Hodge 型 (-1,-1) 本身不能替代最后一项。

这个因子化结论不证明最小性。最小性仍需证明指定形式和链确实看见完整六维纯核与全部活动环面方向，并说明是在泛参数变分还是某个固定代数纤维上工作。

## 实际 motivic coaction 与待定群的区别

对已构造的 \(\mathsf M_s\)，在其生成的张量范畴中，取任意 de Rham 基 \(e_1,\ldots,e_{15}\)，规范余作用已经有意义：

\[
 \Delta I^{\mathfrak m}(s)=
 \sum_{j=1}^{15}
 [\mathsf M_s,\sigma_B,e_j]^{\mathfrak m}
 \otimes
 [\mathsf M_s,e_j^\vee,[\omega_s]]^{\mathrm{dr}}.
\]

这是 [Brown，Notes on motivic periods，式 2.3](https://arxiv.org/html/1512.06410#S2.SS2) 的一般基无关公式在实际对象上的应用。对一般数域需相应使用 Betti 与 de Rham 纤维函子的基域约定；上述公式在原实例 k=Q 时可直接按通常记号理解。

计算全部特殊函数值并非定义此余作用的前提。当前仍未得到的是：指定 K3 与质量单位基下的完整比较矩阵、所有关系、候选块群恰等于实际动机群的证明，以及由此导出的显式最小余作用。实际余作用的存在与这些显式识别问题应在论文中分开陈述。

## 假设 11.2 的当前拆分

| 内容 | 当前结论 | 下一项可检验目标 |
|---|---|---|
| 相对动机存在与物理链 | 命题一及有效同调构造，欧氏正质量范围 | 推广到所需解析延拓分支时跟踪链 |
| 边界一动机存在 | 命题二，实际商秩九 | 与后续显式群描述合并 |
| 指定纯 K3 核 | 已在 Nori 范畴识别边界正交因子，见 K3_KERNEL | 特殊纤维 Picard 数、框架与比较矩阵条目 |
| 质量 Kummer 比较 | KUMMER_BOUNDARY 的实际留数、纤维及局部化推导 | 统一主论文的杯积、余字符与周期符号 |
| 物理链活动因子化 | 命题三，实际分裂实商已识别 | 独立的最小性下界 |
| 规范 motivic coaction | 对实际对象可定义 | 指定群、基、关系及显式比较条目 |
| 全局 K3 认证 | 独立计算支线 | 完整闭路、21 个迹及绝对周期归一化 |

Q6 的指数二补全在张量 Q 后不改变对象，无法替代这些扩张比较。二十一项迹即使全部严格认证，也只识别相应单值和 Gram 数据；它们不单独确定 de Rham 归一化后的绝对周期矩阵。

## 有限代数核对

`code/verify_sector_geometry.py` 从原始 F_s 对二十四个排列分别计算代换、Jacobian 与分母正性多项式，结果记录在 `data/sector_geometry_checks.json`。这些精确计算验证扇区恒等式；对偶性、动机左伴随与下降条件由上面的数学论证承担。
