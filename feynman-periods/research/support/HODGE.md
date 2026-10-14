# 四线 banana 的相对 Hodge 结构与扩张问题

当前状态更新（2026-10-13 07:05）：LOWER_REGULATOR已独立认证指定m=(1,2,4,8)、s=−64的下层实Hodge非零和最小秩八；MIXED_UNIPOTENT_GROUP给出活动对象的完整13维幺幂核。下文较早“固定非零未证”的文字仅记录本稿原范围；完整混合比较和固定纯群仍待完成。

本稿从具体几何推出环境相对对象的权重层及 Hodge 维数，再说明识别物理周期最小对象所需的比较。本稿中的范畴论论证先在有理 mixed Hodge structures 中进行；[后续几何比较](GEOMETRIC_COMPARISON.md) 已另行构造有效 Nori 同调动机及其实际最大一动机商，并明确哪些比较尚未提升。

主要项目来源为 [相对粘合稿第 3 至 9 页](E:/document/feynman_period/banana_relative_motivic_gluing.pdf) 和 [综合稿第 11 节](E:/document/feynman_period/feynman_periods_expanded.tex)。前者对比较提出了具体论证，后者把它保留为假设 11.2。本稿逐项分析论证，不将两种证据状态混为一谈。

## 几何对象和适用范围

取非零质量平方 mu_i=m_i^2，四线 banana 的方程为

\[
 U_G=abc+abd+acd+bcd,\quad
 F_s=(\mu_1a+\mu_2b+\mu_3c+\mu_4d)U_G-sabcd.
\]

排除 s=0 和所有带符号质量和的平方。记 X_s 为四次曲面的最小光滑解消，T=(G_m)^3，U_s={F_s=0}∩T。U_s 是光滑仿射曲面。记 D=X_s\U_s。

令 y_i=mu_i alpha_i、y_5=-sum(mu_i alpha_i)。在 sum y_i=0 的超平面上，lambda=(mu_1,mu_2,mu_3,mu_4,s) 满足

\[
 \sum_{i=1}^5\lambda_i\prod_{j\ne i}y_j
 =-\left(\prod_{i=1}^4\mu_i\right)F_s.
\]

这将模型识别为 Sylvester Hessian。所用的十节点解消及二十条曲线配置，与 [Dardanelli 和 van Geemen 第 1.5 至 1.6 节](https://arxiv.org/pdf/math/0409322) 一致。其一般 Picard 秩十六层的超越格是 U⊕U(2)⊕A2(-2)；特定参数点是否在此层，必须单独确认。

本节之后的 Hodge 推导使用该光滑 K3 解消和正常交叉边界模型。它不要求每个特殊纤维的 Picard 数都等于十六。

## 边界的四个精确秩

边界由四条 V_i、六条 E_ij、六条 C_ij 和四条 F_i 构成，都是有理的 (-2) 曲线。非零异分量配对为

\[
 C_{ij}\cdot V_i=C_{ij}\cdot V_j=C_{ij}\cdot E_{ij}=1,
 \qquad E_{ij}\cdot F_k=E_{ij}\cdot F_l=1,
\]

其中 {k,l} 是 {i,j} 的补集。因此对偶图有 20 个顶点、30 条边，连通，且

\[
 \operatorname{rank}\partial=19,\qquad b_1(\Gamma_D)=30-20+1=11.
\]

交叉矩阵 Q_D 的秩为 16。四个单位 a/d、b/d、c/d、M/d，其中 M=sum mu_i alpha_i，给出四个独立的主除子关系。这里 M/d 在 U_s 上可逆，因为 s≠0 且 abcd≠0。将它们的 valuation 排成 V，有 Q_D V=0、rank V=4。

由此，二十个除子类的实际跨度至多为 16；交叉矩阵秩又给出至少 16，故其跨度恰为 16。这个结论比仅报告 Gram 矩阵秩多了“所有关系确实来自除子类”的一步。

本轮从相交表和 valuation 重新构造矩阵，以精确有理数核对了这些秩。还核对了整除子基的行列式 -48；这是指定曲线格的结论，不单独证明特殊纤维没有额外代数类。

## 开曲面的 Hodge 结构

采用 wt Q(r)=-2r。正常交叉边界的权谱序列、K3 的 H^1=H^3=0，以及边界分量的 H^1=0，给出

\[
 H^1(U_s)\simeq\mathbb Q(-1)^4,
\]

\[
 0\longrightarrow A_s\longrightarrow H^2(U_s)
 \longrightarrow H_1(\Gamma_D,\mathbb Q)(-2)\longrightarrow0,
 \qquad A_s=H^2(X_s)/\langle[D_i]\rangle.
\]

具体地，weight 2 部分是除子 Gysin 映射的余核；weight 3 部分因各 D_i 有理而消失；weight 4 部分是交点至分量的 incidence 映射之核。故

\[
 \dim A_s=22-16=6,\quad h(A_s)=(1,4,1),\quad \dim H^2(U_s)=6+11=17.
\]

由于 U_s 光滑仿射且复维数为二，H^j(U_s)=0 对 j>2 成立。Euler 数交叉核对为 chi(U_s)=24-(40-30)=14=1-4+17。

只有在 rho(X_s)=16 时 A_s 才等于完整超越部分。若 rho 增加，A_s 仍然是秩六的纯 Hodge 商，但包含额外的代数方向。

## 环面限制映射为何恰好去掉三个方向

H^2(T)=Q(-2)^3 由三个坐标 dlog 的二次外积生成。令 v_i^p 表示第 p 个坐标单位在 D_i 上的 valuation，双留数矩阵为

\[
 Z_{(i,j),(p,q)}=v_i^p v_j^q-v_i^q v_j^p.
\]

在一致的交点定向下，partial Z=0；选择 (V_1,C_12)、(V_1,C_13)、(V_2,C_23) 三行得到单位矩阵。因此 H^2(T)→H^2(U_s) 为单射，在 weight 4 部分秩为三。本轮已独立核对这个单位子式。

定义 N_s=coker(H^2(T)→H^2(U_s))(3)，则

\[
 0\longrightarrow A_s(3)\longrightarrow N_s
 \longrightarrow\mathbb Q(1)^8\longrightarrow0.
\]

这里的八来自 11-3，不能直接替换为质量单位的秩 kappa。

## 环境相对对象的完整分层

定义

\[
 E_s=H^3(T,U_s)(3).
\]

相对复形为 Cone(RGamma(T)→RGamma(U_s))[-1]。长正合列、H^3(U_s)=0 和 H^3(T)=Q(-3) 给出

\[
 0\longrightarrow N_s\longrightarrow E_s\longrightarrow\mathbb Q(0)\longrightarrow0.
\]

所以 E_s 的权重过滤由实际对象构造出来：

\[
 W_{-4}E_s=A_s(3),\quad W_{-2}E_s=N_s,\quad W_0E_s=E_s,
\]

\[
 \operatorname{Gr}_{-4}^WE_s=A_s(3),\quad
 \operatorname{Gr}_{-2}^WE_s=\mathbb Q(1)^8,\quad
 \operatorname{Gr}_0^WE_s=\mathbb Q(0).
\]

维数依次为 6、8、1，rank E_s=15。A_s 原本权重为 2，twist (3) 后为 -4；无需把未扭转的 T_K3 直接写成负四权层。

各层的 Hodge 型及重数如下。

| 权重层 | Hodge 型和重数 |
|---|---|
| -4 | (-1,-3):1，(-2,-2):4，(-3,-1):1 |
| -2 | (-1,-1):8 |
| 0 | (0,0):1 |

由 Hodge 过滤的严格性得到

| p | 1 | 0 | -1 | -2 | -3 |
|---|---|---|---|---|---|
| dim F^p E_s | 0 | 1 | 10 | 14 | 15 |

这些维数不指定 F^p 在 E_B⊗C 中的具体位置；实际位置携带扩张周期。也不声称存在有理的 Hodge 直和分裂。

另一个容易遗漏的结果是 H^2(T,U_s)=Q(-1)：H^1(T)=Q(-1)^3 注入 H^1(U_s)=Q(-1)^4，且 H^2(T)→H^2(U_s) 单射。因而整个相对复形不集中于三次；本文的秩十五结论专指 H^3(T,U_s)(3)。

## 从八个边界方向到 kappa 个活动方向

令 A=A_s(3)，B=Q(1)^8，K=E_s/A。实际对象给出 0→B→K→Q(0)→0。[KUMMER_BOUNDARY.md](KUMMER_BOUNDARY.md) 通过实际留数、单位杯积和同伦纤维推出，识别 K 的质量点。在该稿的留数基 P 及 nu=e_F-e_B 约定中，其坐标为

\[
 (\log(\mu_1/\mu_4),\log(\mu_2/\mu_4),\log(\mu_3/\mu_4),0,0,0,0,0)
 \pmod{(2\pi i\mathbb Q)^8}.
\]

该稿第三、四节提供实际比较映射；tame residue 表负责计算该映射的单位值。一动机同调余字符基为 −P，故其点为三个质量比及五个 1；不能忽略这个整体符号。

在上述正实代数质量及几何条件下，令 kappa 为三个质量平方比在 k^×⊗Q 中的秩。由实对数，质量的乘法关系恰对应这些扩张坐标的有理线性关系。

取 B_act⊂B 为扩张类的最小有理 Tate 支撑。投影到 B/B_act 后的扩张分裂，且 Hom_MHS(Q(0),B/B_act)=0，所以分裂唯一。由此得到 K_act⊂K，并定义

\[
 E_{\rm act}=E_s\times_K K_{\rm act}.
\]

它是实际子对象，满足

\[
 \operatorname{Gr}^W E_{\rm act}
 =A_s(3),\ \mathbb Q(1)^\kappa,\ \mathbb Q(0),
 \qquad\operatorname{rank}E_{\rm act}=7+\kappa.
\]

相应 Hodge 过滤维数是 (dim F^0,dim F^-1,dim F^-2,dim F^-3)=(1,kappa+2,kappa+6,kappa+7)。对 (1,2,4,8)，mu=(1,4,16,64)，kappa=1，得到秩八的实际活动实现；这不是秩十五环境对象的另一种计数，也不是最小性证明。

Kummer 比较通过 k 有理交点的实际 Gysin 映射排除 Artin 扭曲，给出 E_s/E_act=Q(1)^(8-kappa) 的分裂实商。几何复共轭在商上作用为 -1，在真实正链上作用为 +1，故物理链因子化成立。这个结论同时使用 [几何比较稿](GEOMETRIC_COMPARISON.md) 的实际链与对偶，以及边界稿的下降；仅有抽象 Hodge 型 (-1,-1) 不足以推出它。

即使链因子化成立，最小性也需要额外下界：完整六维纯核、活动边界是否被指定框架看见，以及所用常数域。对泛参数变分成立的最小性，不自动转移到每个特殊纤维的数值周期。

## 剩余问题可以精确地表为扩张纤维

这一节给出当前最有用的理论拆分。固定已构造的 N=N_s 及序列 0→A→N→B→0。由于 A、N、B 的 Hodge 型第一指标均为负，F^0A=F^0N=F^0B=0。

对一个扩张 0→N→E→Q(0)→0，取 1 的有理提升 e_B。严格性使 F^0E→C 成为同构，因而有唯一 Hodge 提升 e_F。差

\[
 \nu_E=[e_F-e_B]\in N_\mathbb C/N_\mathbb Q
\]

正是该 Hodge 扩张类。换有理提升只改变 N_Q 中的元素。反过来，任取 N_C 中的代表，在 N⊕Q 上把 F^0 定义为 C(e_B+nu)，并由 F^pN 与这条直线生成其余 p≤0 的过滤，即构造具有给定 graded pieces 的扩张。因此这个描述不仅是计数。

投影 nu_E 到 B_C/B_Q 得到边界 Kummer 类。即使这个影像已经完全算出，nu_E 在 N_C/N_Q 中仍未确定。固定该影像且存在一个提升时，所有提升构成

\[
 \operatorname{Ext}^1_{\rm MHS}(\mathbb Q(0),A)
 \simeq A_\mathbb C/A_\mathbb Q
\]

的挠子；这里使用 Hom(Q(0),B)=0 和扩张长正合列。这是固定 N、固定边界类之后的陈述；N 自身的扩张也需几何识别。

这个剩余对象是六维纯层对应的 Hodge 扩张自由度，并非“六个已证明独立的数值常数”。在变分中，几何来源、连接和 Griffiths 横截性还会约束其变化。它指出实际目标：从相对积分或 regulator 明确构造 nu_E，验证其边界投影、连接和物理初始分支。仅有 Hodge 数、边界对数或一个 Q6 不能代替此计算。

由实际三步对象得到的两类 eta∈Ext^1(B,A)、xi∈Ext^1(Q(0),B) 满足 eta∘xi=0，因为 E 就给出了 xi 向 N 的提升。这一形式结果不确定提升 nu_E 的具体值，也不自动完成 motivic 比较。

## 与实际动机比较的衔接

[GEOMETRIC_COMPARISON.md](GEOMETRIC_COMPARISON.md) 已从物理对的同调构造有效 Nori 动机 mathsf E_s，并通过最大一动机商定理证明实际商 mathsf K_s 的 Hodge 实现就是 E_s/A_s(3)。该结论不需要先猜质量对数，也不假定任意 Hodge 子对象都来自动机。

[K3_KERNEL.md](K3_KERNEL.md) 已构造边界正交 K3 因子到纯核的实际 Nori 同构。其关键是：限制到开曲面只消去边界除子；相对连接的核来自纯权重四的 H^2(T)，不与纯权重二的剩余 K3 因子相交。

实际一动机的字符格、下降和质量点已由 KUMMER_BOUNDARY 给出证明稿。PICARD_CERTIFICATE 证明目标 s=−64 的 rho=16；MINIMALITY 给出链生成性和最小商判据。REGULATOR_COMPARISON 已推导下层活动扩张的 Chow regulator 比较及三个 Tate 直和因子；DE_RHAM_MINIMALITY 给出显式 F² 提升、质量连接等式和很一般参数最小性。固定 s=−64 的 eta_act 非零、nu_E 的纯核提升求值、全部特殊函数及完整群识别仍未完成。

Q6 的指数二补全属于积分结构支线。它在张量 Q 后是同构，不能消除上述有理扩张问题，也不应成为继续这一理论主线的前置阻塞。
