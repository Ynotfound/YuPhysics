# 物理框架的生成性与最小性判据

本稿把“秩八活动实现”与“秩八最小形式周期对象”的差距化为一个具体扩张类。以下命题在前述实际 Nori 构造及质量比较下成立；它们不证明数值周期的代数独立性。后续 REGULATOR_COMPARISON.md 已给出下层扩张与高阶 Chow regulator 的实际比较，DE_RHAM_MINIMALITY.md 给出很一般参数的非零性证明稿；固定 s=−64 的非零性现由 LOWER_REGULATOR.md 的实分裂障碍另证，故该点最小秩八。

## 设置

固定正质量欧氏纤维，先假定 \(\rho(X)=16\)。由 [纯核比较](K3_KERNEL.md) 和 [质量比较](KUMMER_BOUNDARY.md)，有实际对象
\[
 0\longrightarrow A=T(X)(3)\longrightarrow E=E_{\rm act}
   \longrightarrow K\longrightarrow0,\qquad
 0\longrightarrow B=\mathbb Q(1)^\kappa\longrightarrow K
   \longrightarrow\mathbb Q(0)\longrightarrow0.
\]
记 \(N=W_{-2}E\)，其扩张类为
\(\eta\in\operatorname{Ext}^1_{\mathrm N}(B,A)\)。
K 是单位点生成的最小活动一动机。物理框架是
\[
 \sigma\in E_B,\qquad \omega\in F^3(E^\vee_{\rm dR}),\qquad
 I=\langle\omega,\sigma\rangle>0.
\]
采用先前的链归一化，σ 在顶层 Q(0) 的像为 1。ω 在 A 上的限制非零：相对残数给出 K3 的非零全纯二形式，且 \(F^3K^\vee=0\)。

A 的 Hodge 实现为简单对象，故 A 在 Nori 范畴也简单，因为 Hodge 实现忠实且正合。超越 K3 Hodge 结构的简单性可见 [Huybrechts，Chapter 3，Lemma 2.7 和 Lemma 3.1](https://www.math.uni-bonn.de/people/huybrech/K3Global.pdf)：一个真子 Hodge 结构若包含 (2,0) 方向便违反超越部分的最小性；若不包含，则为纯 (1,1) 型，违反其与 Néron–Severi 群的正交性。

[PICARD_CERTIFICATE](PICARD_CERTIFICATE.md) 已证明本项目的三个样本包括 \(m=(1,2,4,8),s=-64\) 满足此 Picard 条件。对于 Picard 数跳跃的其他纤维，本稿不直接套用秩六简单核；应先构造完整 Néron–Severi 投影子、剔除额外代数方向并说明下降。

## 一、K 的顶层不能由真子对象承载

**引理。** 若 \(K'\subset K\) 满射到 Q(0)，则 \(K'=K\)。

**证明。** 令 \(B'=K'\cap B\)。推出到 B/B' 后，K' 给出顶层的截面，所以 K 的 Kummer 类在 B/B' 中为零。但 B 是该类的最小 Tate 支撑，故 B'=B，从而 K'=K。

这里也可直接用质量单位的有理乘法独立性核对：任何非零商 B→Q(1) 对应非平凡的单位组合，其对数不落入 \(2\pi i\mathbb Q\)。正实代数单位在该群中的核仅为 1。无须从一般 Hodge 扩张的消失推断 Nori 扩张的消失。证毕。

## 二、物理 Betti 向量生成整个活动对象

**命题一。** 包含 σ 的最小 Nori 子对象就是 E。

**证明。** 设 R⊂E 包含 σ。R 在 K 中的像满射到顶层，因而由引理等于 K。因 A 简单，R∩A 为零或 A。

若为 A，则 R 同时包含核并满射到商，故 R=E。
若为零，则 R≅K。态射在 Hodge 过滤上严格，ω 限制到 R 属于 \(F^3R^\vee=F^3K^\vee=0\)，所以 \(I=\langle\omega,\sigma\rangle=0\)，与欧氏积分的正性矛盾。故只能 R=E。证毕。

这排除了由较小子对象承载物理链的可能性，尚未排除 ω 杀死某个子对象后得到较小商。正性不能代替后者的扩张检验。

## 三、最小商由一个扩张映射控制

定义有理线性映射
\[
 L_\eta:\operatorname{Hom}_{\mathrm N}(\mathbb Q(1),B)
 \longrightarrow \operatorname{Ext}^1_{\mathrm N}(\mathbb Q(1),A),
 \qquad v\longmapsto v^*\eta.
\]
左侧与 Q^κ 同构。其“单射”指这些扩张类在实际 Nori Ext 群中有理线性独立。

**命题二。** 若 κ>0，则 E 是该框架的最小形式周期对象，当且仅当 \(L_\eta\) 单射。

**证明。** 由命题一，只需判断是否存在非零子对象 J⊂E 被 ω 杀死。这是形式矩阵系数的最小子商构造，见 [Brown，§2.4，Lemma 2.4](https://arxiv.org/html/1512.06410#S2.SS4)；对通常以 E∨ 表示的动机周期取对偶不改变秩。

因 A 简单且 ω|A≠0，J∩A=0，故 J 注入 K。
如果其像在 B 内，记为 B'，则 J 给出 η 沿 B' 的分裂。非零 B' 含一个 Q(1) 直和因子，故 \(L_\eta\) 有非零核。
如果其像不在 B 内，引理给出像为 K，于是 J≅K 为 E→K 的截面，特别地 η=0。这在 κ>0、\(L_\eta\) 单射时不可能。因此单射蕴含 J=0。

反之，如果 \(v^*\eta=0\) 且 v≠0，则相应 Q(1) 可提升为 N⊂E 中的子对象。提升唯一，因为 \(\operatorname{Hom}(\mathbb Q(1),A)=0\)，可由纯权重 −2、−4 的不同及忠实 Hodge 实现看出。ω 限制到该 Tate 子对象为零，因为 \(F^3\mathbb Q(-1)=0\)。因此存在严格较小的承载商，E 不是最小对象。证毕。

非零 Hodge 扩张足以证明相应 Nori 扩张非零；反向使用必须另有理由。特别地，某个 Hodge regulator 求值为零不能自动证明动机扩张分裂。

## 四、一般的秩公式与 κ=0 分支

令 \(d=\dim_{\mathbb Q}\ker L_\eta\)，并记对应 Tate 子对象为 B₀⊂B。η 沿 B₀ 分裂，唯一的提升给出 \(\widetilde B_0\subset E\)。上述证明同时说明：所有被 ω 杀死且像包含在 B 内的子对象，恰被这个秩 d 的对象包含。

令 \(\bar E=E/\widetilde B_0\)，其余边界秩为 κ−d。若 κ−d>0，诱导的扩张映射单射，\(\bar E\) 已最小。

若 κ−d=0，则
\[
 0\to A\to\bar E\to\mathbb Q(0)\to0
\]
还有一个顶层扩张类 \(\zeta\in\operatorname{Ext}^1_{\mathrm N}(\mathbb Q(0),A)\)。
若 ζ≠0，无可再商去的非零子对象；若 ζ=0，唯一 Q(0) 截面被 ω 杀死，可以再降秩一。综上
\[
 \boxed{\operatorname{rank}M(I^{\mathfrak m})
 =\operatorname{rank}A+\kappa+1-d-\delta,}
\]
其中 δ=1 当且仅当 κ−d=0 且 ζ=0，否则 δ=0。κ=0 一开始就走这个顶层分支；不能把空映射“单射”误用成满秩结论。

该公式的秩针对固定代数纤维、固定 Nori 范畴中的形式周期最小对象。它不等同于参数微分方程的最小阶数，也不等同于求值后复数间的线性独立数目。

## 五、目标样本的精确剩余问题

对 \(m=(1,2,4,8),s=-64\)，已证明 r=rank A=6，且 κ=1。故
\[
 \operatorname{rank}M(I^{\mathfrak m})=8
 \quad\Longleftrightarrow\quad
 \eta_{\rm act}\ne0
 \ \text{in }\operatorname{Ext}^1_{\mathrm N}(\mathbb Q(1),T(X)(3)).
\]
如果该类为零，最小秩降为 7 或 6，还要看 ζ。LOWER_REGULATOR现已用实际混合周期认证该类非零，所以本目标纤维的秩八最小性成立。质量 log 2 非零只证明上层 Kummer 扩张非零，不证明下层 \(\eta_{\rm act}\ne0\)。

利用 Tate 扭转，可把目标群写成
\[
 \operatorname{Ext}^1_{\mathrm N}(\mathbb Q(0),T(X)(2)).
\]
下文给出已核对除子相消的高阶 Chow 循环。REGULATOR_COMPARISON.md 随后用支持上同调与同伦纤维完成其与 \(\eta_{\rm act}\) 的比较；DE_RHAM_MINIMALITY.md 推进了很一般参数非零性。固定代数基点的非零性由LOWER_REGULATOR独立完成。

## 六、质量射线的一个显式高阶 Chow 循环

采用 KUMMER_BOUNDARY 的整数留数基 P，活动方向由
\(c=3P_1+2P_2+P_3\) 张成，因为质量对数列为 \(-2\log 2\,c\)。
对 F_i，令 j<k<l 为其余下标，设
\[
 u_i=-\frac{\mu_j\alpha_j}{\mu_l\alpha_l}.
\]
其零点、1 点、极点分别为 \(F_i\cap E_{kl}\)、\(F_i\cap E_{jl}\)、\(F_i\cap E_{jk}\)。
在 E_ij 的例外二次曲线上，令 k<l 为互补下标，取
\(v_{ij}=\alpha_k/\alpha_l\) 的严格变换限制。它在 F_k 处为零、在 F_l 处为极点。在 C_ij 交点两坐标同时消失，约去公共零点后其值为 −1，故那里无除子贡献。

定义
\[
\begin{aligned}
 f_{F_1}&=(1-u_1)^2/u_1,&
 f_{F_2}&=u_2/(1-u_2)^3,\\
 f_{F_3}&=(1-u_3)^3/u_3^2,&
 f_{F_4}&=u_4/(1-u_4)^2,\\
 (f_{E_{12}},f_{E_{13}},f_{E_{14}},f_{E_{23}},f_{E_{24}},f_{E_{34}})
 &=(v_{12},v_{13}^{-2},v_{14}^3,v_{23},v_{24}^{-2},v_{34}).
\end{aligned}
\]
其他边界分量取常数 1。按照每条边 E→F 的方向，E 上的除子系数为 c_e，F 上为 −c_e。十二个非零旗标逐点抵消，因此
\[
 Z_c=\sum_D(D,f_D),\qquad
 \sum_D\operatorname{div}_D(f_D)=0.
\]
这给出光滑曲面上高阶 Chow 群 \(\operatorname{CH}^2(X,1)_{\mathbb Q}\) 的一个循环代表；本节的除子检验本身不证明非零，LOWER_REGULATOR现另证固定纤维的不可分解类非零。这里采用曲线函数的除子复形描述，参见 [Sato，§2.1](https://arxiv.org/html/2408.09102#S2.SS1)。其 §2.2–2.4 给出 regulator 的链积分实现，但该文的非消失定理针对另一类 K3，不能直接用于本模型。

检验见 [除子核对数据](data/minimality_cycle.json) 和 [验证程序](code/verify_minimality_cycle.py)。[REGULATOR_COMPARISON](REGULATOR_COMPARISON.md) 已构造并证明
\[
 (\pi_{\rm res})_*\,r_{\rm H}(Z_c)
 \quad\text{与}\quad
 R_H(c^*\eta)(-1)
\]
在 \(\operatorname{Ext}^1_{\rm MHS}(\mathbb Q(0),T(X)(2))\) 中相等，并给出实际 Nori 层面的等式。符号采用本稿 P 框架和归一化留数。除子相消本身不意味着 regulator 非零；很一般参数与固定基点的非消失仍须分别处理。
