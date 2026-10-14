# banana 的十个有理配对：公开数据恢复与独立核对

本批得到十个配对的实际有理表达式，不再停留在以未知 Y 表示连接。采用 [Duhr–Maggio–Porkert–Semper–Stawinski](https://arxiv.org/abs/2507.23061) 的 Frobenius 归一化，记配对为 \(\widehat Y\)。它与本项目绝对几何 trace 及整 Betti 标记的尺度比较仍保留；不能将这个归一化问题隐去。

## 一、原始数据和框架

数据来自 Sara Maggio 的 [Bonn 补充数据，2025，V1](https://doi.org/10.60507/FK2/J4OPZY)，许可证 CC BY 4.0。原始文件保存在 external/bonn-banana-2025，保留作者文件，派生数据另存。七份文件下载后均逐字节核对服务端 MD5，并保存 SHA256。

使用的关键输入是 ZCYInverse.m。它是 Mathematica 表达式文本，程序仅将其解析为数学表达式，没有执行其中代码。CBanana4M.mx 是平台相关二进制 dump，本批未加载，也没有让它成为计算依赖。四份 Diffeqm*.m 已用于恢复 m=(1,2,4,8) 射线的完整纯六维系统，见 PURE_RAY_SYSTEM.md。

文献变量为 \(x_i=m_i^2/p^2=\mu_i/s\)。第11–12个 PDF 原页已查看，正文式 (3.13)、(3.15)、(3.16) 固定 Frobenius 配对及微分框架；HTML 的相应编号是 (30)、(32)、(33)。两种编号不同，引用时勿据文本编号猜页码。相关框架的前五项为全纯周期及四个一阶质量导数。

设 Z 是此六元框架的配对矩阵。由 Hodge 正交性，Z 有块形
\[
 Z=\begin{pmatrix}0&0&z\\0&Y_x&b\\z&b^t&c\end{pmatrix},
 \qquad (Z^{-1})_{\{2,3,4,5\},\{2,3,4,5\}}=Y_x^{-1}. \tag{1}
\]
最后一项取二阶导数或其一半不改变这个中间块。因此不需要反演整个长表达式，也不会把最后一项的因子二带入十个配对。

## 二、紧凑的完整有理表达式

令 D_x=diag(x₁,x₂,x₃,x₄)。从原数据提取并精确化简得到
\[
 \boxed{Y_x^{-1}=\tfrac13D_x C(x)D_x,\qquad
        Y_x=3D_x^{-1}C(x)^{-1}D_x^{-1}.} \tag{2}
\]
对角项由余下三个指标组成：
\[
 C_{ii}=2\left[(x_i-1)^2+
 3\left(\sum_{j\ne i}x_j^2-2\!\sum_{j<k;\ j,k\ne i}x_jx_k\right)\right]. \tag{3}
\]
若 {i,j,k,l}={1,2,3,4}，非对角项为
\[
 C_{ij}=3(x_i+x_j-x_k-x_l)^2-4x_ix_j-12x_kx_l
        -2(x_i+x_j)+6(x_k+x_l)-1. \tag{4}
\]
式 (2)–(4) 已确定全部十个有理函数。data/banana_metric.json 还保存精确伴随式的分子、分母；无需巨大展开即可重放。

定义消去平方根后的多项式
\[
 \Delta(x)=\prod_{\epsilon\in\{\pm1\}^4}
       (1+\epsilon_1\sqrt{x_1}+\cdots+\epsilon_4\sqrt{x_4}).
\]
用独立十六因子展开核对，所有根变量的奇次幂相消，并得到
\[
 \det C=-27\Delta,\qquad
 \det Y_x=-\frac3{(x_1x_2x_3x_4)^2\Delta}. \tag{5}
\]
这是精确多项式恒等式，并非浮点拟合。

## 三、转换到本项目的质量与 s 框架

本项目形式满足 \(\omega_{\mu,s}=s^{-1}\omega_{x,1}\)。固定 μ₄，令 D₁,…,D₄ 为 PURE_CONNECTION 中的三个对数质量导数和 ∂s。参数变换矩阵为
\[
 J_{ai}=D_a x_i=
 \begin{cases}x_i\delta_{ai},&a\le3,\\-x_i/s,&a=4.\end{cases}
\]
对 s 求导出现的额外 ω 项与所有 F¹ 方向配对为零，故
\[
 \boxed{\widehat Y=s^{-2}JY_xJ^t.} \tag{6}
\]
也可给出更便于计算的逆矩阵：
\[
 \widehat Y^{-1}=\frac{s^2}{3}L^tC(\mu/s)L,
 \qquad
 L=\begin{pmatrix}1&0&0&0\\0&1&0&0\\0&0&1&0\\-1&-1&-1&-s\end{pmatrix}. \tag{7}
\]
右端在 μ 和 s 中是多项式矩阵。
由 det J=−x₁x₂x₃x₄/s 得到
\[
 \det\widehat Y=-\frac3{s^{10}\Delta(\mu/s)}.
\]
在 m=(1,2,4,8) 上，
\[
 P_8(s)=\prod_{t\in\{1,9,25,49,81,121,169,225\}}(s-t),
 \qquad\det\widehat Y=-\frac3{s^2P_8(s)}. \tag{8}
\]

特别地，在指定 s=−64，逆矩阵非常短：
\[
 \boxed{\widehat Y(-64)^{-1}=\frac13
 \begin{pmatrix}
 48738&-2469&-8517&-2035392\\
 -2469&46182&-9669&-1838784\\
 -8517&-9669&33078&-1052352\\
 -2035392&-1838784&-1052352&136798208
 \end{pmatrix}.} \tag{9}
\]
其配对行列式为 \(-3/353067798239177216000\)。这认证了文献归一化下的有理耦合值，尚不认证任何绝对周期值。

## 四、不依赖原矩阵的局部交叉核对

由本项目 Laurent 方程的常数项展开，全纯解为
\(F=\sum_{n\in\mathbb N^4}(|n|!/\prod n_i!)^2x^n\)。
通过带指数的 Gamma 系数微分，四个对数解的正则部分为
\[
 G_i=\sum_n2\left(\frac{|n|!}{\prod n_j!}\right)^2
        (H_{|n|}-H_{n_i})x^n,
 \qquad \tau_i=\log x_i+G_i/F.
\]
记 θ_a=x_a∂_(x_a)，S 的对角为零、非对角为−1。Frobenius 极化关系给出
\[
 (Y_\theta)_{ab}=F^2\sum_{i,j}S_{ij}(\theta_a\tau_i)(\theta_b\tau_j).
\]
脚本从阶乘与调和数独立生成四变量级数，核对
\[
 (C/3)Y_\theta=I_4+O(|x|^5)
\]
的全部十六个条目，共涉及七十个单项式的精确有理系数。没有把 ZCYInverse 的系数用作级数输入。

这个有限阶交叉核对有能力发现参数、转置、尺度和符号错误，但不单独证明无限阶恒等式。完整有理式的来源仍是已注明出处的公开矩阵，配合精确抽取、求逆与判别式核对。没有把有限匹配升级为新的全局证明。

## 五、归一化及后续连接工作的边界

文献的无 2πi 对数 Frobenius 框架不是本项目已经标记的整 Betti 基。其 Σ 也不能直接当成行列式48的整格矩阵。因此本批明确不固定相对于本项目绝对几何 trace 的常数尺度；有关字段 absolute_geometric_cup_scale_certified 保持 false。

这一尺度问题不阻塞纯微分系统：把 Y 乘以非零常数不改变 Levi–Civita 系数和 PURE_CONNECTION 的 C_ab，只改变补足方向 w 的常数尺度。可在明确的文献归一化下使用完整重建式。这条直接路径现已完成：从四份公开微分矩阵取 epsilon=0、去掉 tadpole 并核对纯六维闭合，再沿 x_i=μ_i/s 拉回。PURE_RAY_SYSTEM.md 给出实际矩阵、全射线配对相容性和有限极点的完整根包围。

仍未完成的是整循环标记、绝对比较常数、严格全局闭路以及 s=−64 的 regulator 非零性。式 (9) 的非退化性不能代替最后一项。

可重放入口为 [verify_banana_metric.py](code/verify_banana_metric.py)，结果为 [banana_metric.json](data/banana_metric.json)。

补充归一化审计：原 ZCYInverse.m 对应最后一项为半二阶导数的基；转换到完整二阶导数的 T=diag(1,1,1,1,1,2) 已由单点偏导相容性及全射线逆配对恒等式共同核对。不能将原完整六阶矩阵直接当作 W 基的配对。中间四阶块和本文 (2)–(9) 不受此因子影响。

## 后续定标更新（2026-10-13）

ABSOLUTE_COMPARISON已在指定质量射线、裸Poincare留数约定下完成Morse定标：P_D^t K_D P_D=8 pi² Q6；标准Tate trace下几何代数配对为公开配对的−2倍。六循环绝对P6在s=−64的明确有理B框架误差≤2^(-85)，八节点尺度一致。MARKED_LATTICE已完成全部八节点支路、实际T0标记及三种补全歧义的证明。本稿较早“尚无绝对尺度/周期”的记录不再是这些项目的当前状态。完整四参数平坦性、具体饱和半循环及固定regulator非零仍独立待证。
