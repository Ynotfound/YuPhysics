# 纯 K3 层的六维框架与连接重建公式

本稿从已有质量导数和四维参数族的局部 Torelli 出发，构造未扭转纯核 K 的六维 de Rham 框架，并把全部 Gauss–Manin 连接系数归约到一个四阶有理配对矩阵 Y 及其导数。这里给出重建公式。后续 BANANA_PAIRING 已恢复文献 Frobenius 归一化下的十个有理函数，PURE_RAY_SYSTEM 已恢复指定质量射线的实际纯系统；绝对几何 trace 尺度、原标量六阶算子及绝对周期仍未完成。

## 一、纯 F¹ 的五元基

沿用 [DE_RHAM_MINIMALITY](DE_RHAM_MINIMALITY.md) 的记号：K 是权重二的秩六边界正交核，物理纯核是 K(3)。取 μ₄ 为固定非零尺度，在四维参数图表上令
\[
 D_1=\mu_1\partial_{\mu_1},\quad
 D_2=\mu_2\partial_{\mu_2},\quad
 D_3=\mu_3\partial_{\mu_3},\quad
 D_4=\partial_s,
 \qquad v_a=\nabla_{D_a}\omega.
\]
四个向量场彼此交换。取非零全纯二形式 ω 的原有归一化。十节点族的 Kodaira–Spencer 单射和局部 Torelli 给出 v₁,…,v₄ 在 Gr¹_F K 中的基；该商秩为四。因此
\[
 \boxed{F^1K=\langle\omega,v_1,v_2,v_3,v_4\rangle.} \tag{1}
\]
由已有恒等式 \(\nabla_s\beta_i=\mu_i\nabla_{\mu_i}\omega\)，前三项正是下层质量提升的纯导数。这给出在 m=(1,2,4,8)、s=−64 同样有效的纯层 F¹ 框架；它本身不证明该纤维扩张非零。

## 二、利用极化补足第六个方向

用曲面杯积和 de Rham trace 固定未扭转配对 \(\langle\ ,\ \rangle\)。它平坦、对称、非退化。Hodge 正交性给出
\[
 \langle\omega,\omega\rangle=\langle\omega,v_a\rangle=0,
 \qquad F^1K=\omega^\perp.
\]
定义对称矩阵
\[
 Y_{ab}=\langle v_a,v_b\rangle. \tag{2}
\]
F¹/F² 上的诱导配对非退化，而 v_a 是该商的一组基，故 det Y≠0。所有条目由代数 de Rham 上同调及代数 trace 定义，属于参数函数域；无需先知道绝对周期。

由平坦性，
\[
 \langle\omega,\nabla_{D_a}v_b\rangle=-Y_{ab}.
\]
记 Y 的逆为 \((Y^{ab})\)，则
\[
 \tau=-\frac14\sum_{a,b=1}^4Y^{ab}\nabla_{D_a}v_b,
 \qquad\langle\omega,\tau\rangle=1. \tag{3}
\]
因此 {ω,v₁,…,v₄,τ} 已是 K 的六元基。为了简化连接，先令
\[
 b_a=\langle v_a,\tau\rangle,\qquad
 \widetilde w=\tau-\sum_{a,b}v_aY^{ab}b_b,
 \qquad w=\widetilde w-\tfrac12\langle\widetilde w,\widetilde w\rangle\omega.
\]
于是
\[
 \langle\omega,w\rangle=1,\quad
 \langle v_a,w\rangle=0,\quad\langle w,w\rangle=0.
\]
在框架 \(\mathcal B=(\omega,v_1,v_2,v_3,v_4,w)\) 中，极化矩阵为
\[
 Q_{\mathcal B}=\begin{pmatrix}0&0&1\\0&Y&0\\1&0&0\end{pmatrix}. \tag{4}
\]
这是一组有理 de Rham 基；它不是整 Betti 基，不能与 LATTICE_SATURATION 中的整格矩阵混同。

## 三、全部连接系数由 Y 决定

因 D_a 交换、Gauss–Manin 平坦，\(\nabla_{D_a}v_b=\nabla_{D_b}v_a\)。将其在上述基中展开：
\[
 \begin{aligned}
 \nabla_{D_a}\omega&=v_a,\\
 \nabla_{D_a}v_b&=\sum_c\Gamma^c_{ab}v_c+C_{ab}\omega-Y_{ab}w,\\
 \nabla_{D_a}w&=-\sum_{b,c}Y^{bc}C_{ac}v_b.
 \end{aligned} \tag{5}
\]
第二式 w 的系数由与 ω 配对确定；第三式由分别与 ω、w、v_b 配对确定。因此 C 对称，Γ 对下标 a,b 对称。

对 Y 求导并使用配对平坦性，得到唯一的系数
\[
 \Gamma^c_{ab}=\frac12\sum_dY^{cd}
 (D_aY_{bd}+D_bY_{ad}-D_dY_{ab}). \tag{6}
\]
这正是复代数度量 Y 在这组交换向量场中的 Levi–Civita 公式；这里不要求 Y 是正定实度量。

C 也无需引入独立周期函数。定义曲率约定
\[
 R^d{}_{cab}=D_a\Gamma^d_{bc}-D_b\Gamma^d_{ac}
 +\sum_e(\Gamma^e_{bc}\Gamma^d_{ae}-\Gamma^e_{ac}\Gamma^d_{be}),
 \quad\operatorname{Ric}_{cb}=\sum_aR^a{}_{cab},
 \quad R=\sum_{b,c}Y^{bc}\operatorname{Ric}_{cb}.
\]
式 (5) 的平坦性在 v_d 系数上给出
\[
 R^d{}_{cab}=-C_{bc}\delta_a^d+C_{ac}\delta_b^d
 -Y_{bc}\sum_eY^{de}C_{ae}
 +Y_{ac}\sum_eY^{de}C_{be}.
\]
对 d=a 缩并（维数四）得到
\(\operatorname{Ric}_{cb}=-2C_{bc}-Y_{bc}\operatorname{tr}_Y C\)，再次缩并得到 \(R=-6\operatorname{tr}_Y C\)。故
\[
 \boxed{C_{ab}=-\frac12\left(\operatorname{Ric}_{ab}
                              -\frac R6Y_{ab}\right).} \tag{7}
\]
这是在已知来自本 K3 族的 Y 上成立的重建式。任意给定四阶对称矩阵并不自动产生平坦连接；其曲率还必须满足上式未缩并的关系以及余下平坦条件。

于是 (2)、(6)、(7) 和 (5) 给出全部纯层连接。从原形式和杯积计算出十个 Y_ab 后，可用纯有理运算和求导恢复四个六阶一阶系统矩阵。

## 四、固定质量射线及明确的剩余计算

固定 μ 时，s 方向只需取式 (5) 的 a=4。三个活动质量提升的导数分别为 v₁、v₂、v₃，第四个由
\[
 \nabla_s\beta_4=-\omega-sv_4-v_1-v_2-v_3
\]
给出。这是由尺度 Euler 恒等式得到的实际耦合列。

BANANA_PAIRING 现已给出十个配对的实际有理式；其归一化来自公开 Frobenius 数据，绝对 trace 常数仍待几何校准。PURE_RAY_SYSTEM 从文献 IBP 系统独立提取 m=(1,2,4,8) 的纯六维系统，并核对全射线极化恒等式。可继续用 (5)–(7) 作另一路重建，检查全四参数平坦性；如需标量算子，再选循环向量并比较表观因子。

这条路线无需把缺失的旧六阶算子作为理论前置输入，但不减少后续严格数值运输的要求。已恢复的配对与一阶系统是精确有理数据；原始几何 trace 常数及严格数值运输仍有独立验收要求。

[verify_polarized_connection.py](code/verify_polarized_connection.py) 在一个可直接微分的非恒定零锥模型上，独立比较原始六维平坦导数与 (5)–(7)，核对曲率缩并的符号和系数；该核对不替代 banana 的杯积计算。

## 实际射线核对更新

POLARIZED_TRANSPORT 已利用逆配对构造同一个唯一w，给出实际有理规范和不含表观极点的s方向连接。全部16个Gamma^c_(4b)已经从四变量逆度量的偏导独立重算并与连接精确一致；不能把这个沿射线核对称为全四参数平坦性证明。

## 后续定标更新（2026-10-13）

ABSOLUTE_COMPARISON已在指定质量射线、裸Poincare留数约定下完成Morse定标：P_D^t K_D P_D=8 pi² Q6；标准Tate trace下几何代数配对为公开配对的−2倍。六循环绝对P6在s=−64的明确有理B框架误差≤2^(-85)，八节点尺度一致。MARKED_LATTICE已完成全部八节点支路、实际T0标记及三种补全歧义的证明。本稿较早“尚无绝对尺度/周期”的记录不再是这些项目的当前状态。完整四参数平坦性、具体饱和半循环及固定regulator非零仍独立待证。
