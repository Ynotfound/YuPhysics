# 九个上方接入点的严格运输与共同干路消去

本稿接续 [HIGH_TRUNK_TRANSPORT](HIGH_TRUNK_TRANSPORT.md)。固定 b=−64、c=−64+64i 和常数规范 D=diag(−64,−64,−64,−64,4096,−1/64)。旧竖直证书、旧代码及其哈希均保留。新的计算从已经认证的 c 出发，沿 Im s=64 到225+64i，并在每个 a+64i，a∈{0,1,9,25,49,81,121,169,225}，保存接入点。

## 一、实际完成的73个新步骤

新计划 [horizontal_trunk/plan.json](data/banana_pure_ray/horizontal_trunk/plan.json) 采用步长至多4，在九个接入点准确停下，共73步。每个运输圆盘对全部有限极点的距离下界至少60。计划先于运输保存；旧45圆盘计划未被修改。

每步仍采用64阶、240位 dyadic 多项式，由高斯整数残差认证。改动的是累计半径的向上取整：

\[
\operatorname{up}_{240}(x)=2^{-240}\lceil2^{240}x\rceil,
\qquad x\le\operatorname{up}_{240}(x)<x+2^{-240}. \tag{1}
\]

这避免了每步向整二幂取整可能损失接近一位精度。矩阵中心的舍入误差仍逐步精确计算。浮点数只提出多项式；全部接受判据使用整数及有理数。

73份多项式见证、每步误差及输入哈希见 [transport.json](data/banana_pure_ray/horizontal_trunk/transport.json) 引用的分段文件。另一个入口逐项重算源矩阵、标量卷积残差和传播不等式，全部通过：[replay.json](data/banana_pure_ray/horizontal_trunk/replay.json)。没有重新生成旧竖直八步。

## 二、两个基点与一次性组合

令 S 为从 b 到 c 的真实竖直运输，F_a 为从 c 到 a+64i 的真实水平运输。我们同时保存：

1. 从 c 的精确单位初值出发的 F_a 包围；
2. 将原 b 初值逐步传播所得的包围。

第二种包围严格，但因反复使用各步范数，远端的估计变得很松。例如225+64i的逐步半径约9.31×10⁻⁴；这反映估计方式的放大，不是已知真实误差。

利用精确恒等式 U_{b,a}=F_a S，在每个接入点将共同种子只组合一次。若

\[
\|F_a-\widehat F_a\|\le f_a,\qquad
\|S-\widehat S\|\le e_S,
\]

则对舍入后的中心 C_a，

\[
\boxed{\|F_aS-C_a\|\le
(\|\widehat F_a\|+f_a)e_S+f_a\|\widehat S\|+\rho_a.} \tag{2}
\]

这里 ρₐ 为精确中心乘法的舍入误差。旧种子误差 e_S=2⁻¹⁰⁷ 没有被删除；只避免把同一初始误差经多个宽松范数重复放大。最终225+64i的这一半径约1.65×10⁻³⁰。小数仅供展示，证书保存完整有理数。

九个接入点以及两种基点下的18个逆矩阵均用 Neumann 残差认证，并独立重放。所有接入点统一可用的矩阵无穷范数界为：

| 基点及方向 | 九个接入点的统一误差上界 |
|---|---:|
| c 到 a+64i | 2⁻¹¹⁹ |
| a+64i 返回 c | 2⁻¹⁰¹ |
| b 到 a+64i，采用式 (2) | 2⁻⁹⁸ |
| a+64i 返回 b | 2⁻⁸⁰ |

各点实际界通常更好，详见 [landmarks.json](data/banana_pure_ray/horizontal_trunk/landmarks.json)。[landmarks_replay.json](data/banana_pure_ray/horizontal_trunk/landmarks_replay.json) 记录一次性组合、所有逆矩阵和上述统一界的独立精确复核。

## 三、二十一项迹可以在高基点计算

对从 c 出发的局部闭路 λₐ，记矩阵为 Mₐ。将它移回 b，按列向量先后作用的约定，其矩阵为

\[
T_a=S^{-1}M_aS. \tag{3}
\]

因此对任意含逆的矩阵字 W，以及 Rₐ=Tₐ−I，

\[
\operatorname{tr}W(T_a)=\operatorname{tr}W(M_a),\qquad
\operatorname{tr}(R_iR_jR_k)
=\operatorname{tr}((M_i-I)(M_j-I)(M_k-I)). \tag{4}
\]

证明只需在乘积中相消相邻 S S⁻¹，再用迹的共轭不变性。第六反射 T₆=T₀T₁T₀⁻¹ 也属于同一矩阵字，故这一化简同时适用于全部十五个对迹与六个三迹。

实际路径比较由 [POLARIZED_TRANSPORT](POLARIZED_TRANSPORT.md) 的上半平面同伦承担：共同竖直干路、水平段以及到各节点的上侧下降路径，与原稿的几何闭路类相容。必须继续保持逆时针及第六反射的次序。

式 (4) 的价值是可以在**写入区间运算之前**，精确消掉共同干路的共轭。因此后续迹计算可使用从 c 精确单位初值出发的高精度证书，而不必让旧竖直种子的误差反复进入。该化简不确定绝对周期、整 Betti 标记或 de Rham 绝对尺度；这些任务仍使用原基点及几何初值。

## 四、下一步和严格边界

本稿认证完整水平开放路径及九个上方接入点。后续 NODE_REFLECTION_CERTIFICATES 已认证八个局部种子；BASED_NODE_TRANSPORT 已完成其中五个节点的共同基点反射，以及十个对迹、两个三迹和对应五阶Gram子块。零点作用、涉及第六循环的剩余九项迹及绝对比较仍待认证。

节点无需必然运输完整六维局部框架：[FROBENIUS_REFLECTION](FROBENIUS_REFLECTION.md) 给出只认证消失特征向量、再用极化恢复反射的严格路线及加权残差界。该路线现已执行，五个节点闭路算子通过上升向量运输和水平逆运输组合得到。零点不是节点反射，仍须使用共振消去后的幂零系统或实际闭路运输。

可重放入口：

- `code/build_horizontal_plan.py`：独立新计划，旧计划保持不变。
- `code/certify_horizontal_trunk.py`：两种归一化、每步精确残差、原子检查点。
- `code/verify_horizontal_trunk.py`：源矩阵与全部73步的独立复核。
- `code/certify_horizontal_landmarks.py`：一次性组合与18个 Neumann 求逆。
- `code/verify_horizontal_landmarks.py`：不做数值求逆的有理证书复核。

固定s=−64的 regulator 非零、完整四参数比较、绝对周期矩阵、有界一般图算法及物理基准均未因这批计算而被替代。全文需继续区分开放运输、共同基点单值、整数格标记和绝对比较。

## 后续完成状态（2026-10-13）

本稿记录相应先行批次。后续ZERO_MONODROMY已完成实际T0及第六反射；Q6_GEOMETRIC_CERTIFICATE已认证全部21迹、完整Q6及实际指数二；GEOMETRIC_MONODROMY_GROUP已证明指定质量射线纯几何单值的Zariski闭包=O(Q6)。本稿中“尚未完成”所指先行批次不再代表这些项目的当前状态。指定补全类、绝对周期与固定纤维regulator非零仍待证；旧证书及其哈希未改变。
