# Morse 留数定标与六循环绝对比较矩阵

当前状态更新（2026-10-13 07:05）：LOWER_REGULATOR已独立认证指定m=(1,2,4,8)、s=−64的下层实Hodge非零和最小秩八；MIXED_UNIPOTENT_GROUP给出活动对象的完整13维幺幂核。下文较早“固定非零未证”的文字仅记录本稿原范围；完整混合比较和固定纯群仍待完成。

本稿在明确的裸 Poincare 留数约定下，固定 MARKED_LATTICE 中尚未知的共同尺度。它给出 m=(1,2,4,8)、s=−64 的六循环绝对周期矩阵；这六个整循环张成指数二子格，同时构成完整的有理 Betti 基。特定指数二补全半循环、旧标量 theta 基及固定纤维的 regulator 非零仍未认证。

## 1. 形式和框架约定

在环面图表上令
\[
 f(\alpha)=\left(\sum_{i=1}^4m_i^2\alpha_i\right)
           \left(\sum_{i=1}^4\alpha_i^{-1}\right),\qquad
 X_s:\ f=s.
\]
纯全纯二形式取
\[
 \omega_s=\operatorname{Res}_{f=s}
 \frac{d\log(\alpha_1/\alpha_4)\wedge
       d\log(\alpha_2/\alpha_4)\wedge
       d\log(\alpha_3/\alpha_4)}{f-s}.                   \tag{1}
\]
这里的二形式不包含 (2πi) 的除因子。与项目参数三形式 Ω/F 的整体符号选择相差时，必须同时调整最初循环取向；本文用下面的局部球取向直接固定它。

实际周期向量采用已保存的有理框架：先取 PURE_RAY_SYSTEM 的完整二阶导数 H 框架，再乘 polarized_system.json 中 gauge_P，得到 B 框架，最后 D=diag(−64,−64,−64,−64,4096,−1/64) 得平衡框架。gauge_P 第一行精确为 (1,0,0,0,0,0)，因此 D 框架第一项就是 −64∫ω_s。公开数据中采用规范化周期而产生的常数倍不改变这个微分系统；这里通过式(1)给实际几何周期选择绝对尺度。

框架中的第六项仍按公开 Frobenius 配对构造，是明确的有理 de Rham 形式。本稿不暗中把它换成另一种 trace 归一化的对偶形式。

## 2. 局部球的留数积分

在一个环面节点，取符号质量 r_i=ε_i m_i、S=Σr_i≠0，设 α_i/α_4=(ε_i/m_i)/(ε_4/m_4)·x_i，x₄=1。于是
\[
 f(x)=\left(\sum r_ix_i\right)\left(\sum r_i/x_i\right),
 \quad a=S^2,\quad x=(1,1,1).
\]
环面体积形式仍为 dlog x₁∧dlog x₂∧dlog x₃，在节点的一阶密度为1。三个局部坐标下的 Hessian 为
\[
 H=2(S\operatorname{diag}(r_1,r_2,r_3)-rr^t),\qquad
 \det H=8S^2r_1r_2r_3r_4.                               \tag{2}
\]

由复解析 Morse 引理，存在坐标 u 使 f−a=u₁²+u₂²+u₃²。节点处体积形式变换的 Jacobian 平方为8/det H。对 t=s−a，实标准球取通常边界取向，有
\[
 \int_{u_1^2+u_2^2+u_3^2=t}
 \operatorname{Res}\frac{du_1\wedge du_2\wedge du_3}
 {u_1^2+u_2^2+u_3^2-t}=2\pi\sqrt t.                    \tag{3}
\]
可直接验证：半径 r=√t 时，取法向量场 \(\sum u_i\partial_{u_i}/(2r^2)\)，其收缩体积形式为球面面积形式除以2r，积分等于4πr²/(2r)。对复 t 和复线性坐标连续延拓定义对应球；平方根选择决定取向。

高阶 Jacobian 项不改变首项。因此在一致的平方根和取向下，
\[
 \int_{\delta_a(s)}\omega_s
 =\frac{2\pi}{S\sqrt{r_1r_2r_3r_4}}\sqrt{s-a}
   +o(\sqrt{s-a}).                                    \tag{4}
\]
这是实际几何留数计算，不是由单值或近整数迹推测出来的周期幅度。

## 3. 节点1确定唯一共同尺度

在 a=1，取 r=(−1,−2,−4,8)，S=1，det H=−512，选择 √(∏r_i)=8i。局部种子为 h=i/4，τ=(s−1)/h；选择 √h=e^(iπ/4)/2。于是实际 D 框架第一项的首项为
\[
 -64\,\frac{2\pi\sqrt h}{8i}\,\tau^{1/2}
 =8i\pi\sqrt i\,\tau^{1/2}.                             \tag{5}
\]

既有精确 Frobenius 解是
\[
 v_D(\tau)=\tau^{-1/2}Q\sum_{n\ge0}P_n\tau^n.
\]
从已认证系数直接得到
\[
 (QP_0)_1=0,\qquad (QP_1)_1=i/128.                      \tag{6}
\]
节点反射的负特征空间一维，所以几何消失周期向量和这个解只差常数。比较式(5)、(6)给
\[
 \boxed{\rho=1024\pi\sqrt i
             =512\sqrt2\,\pi(1+i).}                    \tag{7}
\]
MARKED_LATTICE 构造的所有列与该第一列共享同一幅度，故
\[
 P_D(c)=\rho\,\widetilde P(c).                          \tag{8}
\]
式(7)选择的是明确局部球取向；其他五列沿已经认证的 Gram 取向约定确定。无需另外猜测六个独立周期初值。

## 4. 配对尺度及八节点交叉检查

已保存的 K_D 是在公开 Frobenius 归一化下的逆配对，满足 K_D'+A_DᵗK_D+K_DA_D=0。它在节点1正则。令 v₀=QP₀、v₁=QP₁，则
\[
 v_0^tK_D(1)v_0=0,\qquad
 2v_1^tK_D(1)v_0+h\,v_0^tK_D'(1)v_0=i/65536.            \tag{9}
\]
平坦性使第二式恰为整个局部解的常数范数，而非有限截断近似。因此由 δ₁²=−2 及 \(\rho^2=1048576i\pi^2\)，
\[
 \boxed{P_D^tK_DP_D=8\pi^2Q_6.}                         \tag{10}
\]
六反射绝对不可约保证平坦不变量配对的比例唯一；也可直接用 MARKED_LATTICE 中共同幅度框架的已知 Gram 推出。

其余七个节点独立计算同样的局部首项及式(9)，全部给 ν/π²=8。八个有理恒等式和种子哈希在 all_node_scale_crosscheck.json，独立逐项重放见 all_node_scale_replay.json；没有重算任何 Frobenius 级数或运输。

若代数 de Rham trace 采用标准 Tate 比较 \(\operatorname{per}(\mathbb Q(-2))=(2\pi i)^2\)，则该射线上的几何代数杯积矩阵与公开配对满足
\[
 C_{\rm dR}=-2K_D^{-1}.                                \tag{11}
\]
这来自 Riemann 双线性关系 \(P_DQ_6^{-1}P_D^t=(2\pi i)^2C_{\rm dR}\) 与式(10)。因此一般 PURE_CONNECTION 里的几何归一化各向同性对偶为 −1/2 倍公开框架的第六项；不要不经换基直接把两套第六行混为一谈。完整四参数连接的独立核对仍是另一任务。

## 5. 严格矩阵包围

π 使用已认证 Machin 区间；√2 用320位 dyadic 端点 l<u，精确验证 l²<2<u²。按式(7)进行复标量包围，将尺度、相对框架、量化和乘法误差一起传播，再使用既有共同竖直逆运输返回 s=−64，最后乘 D⁻¹。

| 矩阵 | 基点与框架 | 矩阵无穷范数误差上界 |
|---|---|---|
| P_D(c) | −64+64i，D平衡框架 | 2⁻¹⁰⁵ |
| P_B(−64) | −64，明确的有理B框架 | 2⁻⁸⁵ |

36个条目的实部和虚部均由同一范数界包围；中心为320位高斯 dyadic 数。数据保存在 data/banana_pure_ray/marked_lattice/absolute_periods.json，SHA256=910108d46e952b8cfc67d5bb0bc8c3b413f026f4c799faf6f8d5be4ba9803e68。独立验证器重算 Hessian、首项、精确局部范数、π/√2界、逆运输残差和最终矩阵误差，结果 absolute_periods_replay.json；入口 code/verify_absolute_periods.py。

这是指定六循环有理 Betti 基的绝对比较认证。尚未转换到旧稿未恢复的标量 theta 框架，也未选择指数二饱和格的具体半循环。它不证明固定s=−64的下层regulator非零或完整相对混合比较矩阵。
