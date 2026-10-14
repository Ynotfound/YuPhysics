# 显式 de Rham 提升及很一般参数的最小性

当前状态更新（2026-10-13 07:05）：LOWER_REGULATOR已独立认证指定m=(1,2,4,8)、s=−64的下层实Hodge非零和最小秩八；MIXED_UNIPOTENT_GROUP给出活动对象的完整13维幺幂核。下文较早“固定非零未证”的文字仅记录本稿原范围；完整混合比较和固定纯群仍待完成。

本稿给出下层质量扩张的显式对数二形式、一个 Gauss–Manin 恒等式，并据此证明固定质量射线上很一般参数的非消失。结合已有 Picard 证书和最小性判据，质量 (1,2,4,8)、(1,3,7,15) 的活动实现对很一般欧氏参数分别达到秩八、十。

“很一般”在这里允许排除一个可数集。证明不认定指定代数基点 s=−64 位于该集之外，故该固定纤维的最小性仍待认证。

## 一、十二维对数形式的具体框架

**扭转约定（审阅修正）。** 本稿形式和连接计算均在未扭转的上同调中进行，记 \(K=\mathsf H_{\rm res}\)，其权重为2、F²K 秩1、F¹K 秩5。前稿的物理纯核是 A=K(3)，regulator 的纯核是 A₂=K(2)。以下 F²K、F¹K 不能写成 F²A、F¹A；后者的过滤已经平移。扩张类作 Hodge 比较时才统一扭转 (2)。

沿用光滑十节点解消 X、二十曲线 D，以及
\[
 M=\sum_{i=1}^4\mu_i\alpha_i,\quad
 N=\sum_{i=1}^4\alpha_i^{-1},\quad MN=s\text{ on }U.
\]
在 d=1 图表设 x=a/d、y=b/d、z=c/d，并令
\[
 f=xyz(MN-s),\qquad
 \omega=\operatorname{Res}_{f=0}\frac{dx\wedge dy\wedge dz}{f}
        =\frac{dx\wedge dy}{f_z}.
\]
其在 X 上延拓为非零全纯二形式；固定上述留数次序。

每个 \(\alpha_i/\alpha_j\) 在 D 的各分量上阶数仅为 −1、0、1，且在 U 上可逆。因此它乘以 ω 是全局对数二形式。定义
\[
 q_{ij}^{\pm}=
 \mu_i\frac{\alpha_i}{\alpha_j}
 \ \pm\ \mu_j\frac{\alpha_j}{\alpha_i},\qquad i<j.
\]
它们满足唯一的线性关系
\[
 \sum_{i<j}q_{ij}^{+}=s-\sum_i\mu_i. \tag{1}
\]
于是如下十二个形式构成 \(F^2H^2_{\rm dR}(U)\) 的一组基：
\[
 \omega,\qquad
 \{q_{ij}^{-}\omega: i<j\},\qquad
 \{q_{ij}^{+}\omega: i<j,\ (i,j)\ne(3,4)\}. \tag{2}
\]

**证明。** 对数 Hodge–de Rham 谱序列退化给出
\(H^0(X,\Omega_X^2(\log D))\simeq F^2H^2_{\rm dR}(U)\)，后者由已计算的分层维数为 1+11=12。
十三个函数 1 和 \(\mu_i\alpha_i/\alpha_j\) 仅有 (1) 所对应的关系：U 的 Laurent 方程 MN−s 不可约；任何支撑包含在同一个 A₃ 根多面体中的倍式只能是常数倍，因为 Newton 多面体的 Minkowski 和不能真包含于原来的有界多面体。亦即不能有另一个具有这些单项式的独立线性关系。
将 (2) 的系数列再加上该关系列，得到十三阶矩阵，行列式为 64。因此十二个函数在 U 上独立，恰构成全部对数二形式。证毕。

## 二、六个单位杯积与全部留数的核对

四个单位为 x、y、z、h=M/d。在 d=1 图表，直接使用 df=0 计算得到
\[
\begin{aligned}
 d\log x\wedge d\log y&=(-q_{13}^{-}-q_{23}^{-}+q_{34}^{-})\omega,\\
 d\log x\wedge d\log z&=(q_{12}^{-}-q_{23}^{-}-q_{24}^{-})\omega,\\
 d\log y\wedge d\log z&=(q_{12}^{-}+q_{13}^{-}+q_{14}^{-})\omega,\\
 d\log x\wedge d\log h&=-q_{23}^{-}\omega,\\
 d\log y\wedge d\log h&=q_{13}^{-}\omega,\\
 d\log z\wedge d\log h&=-q_{12}^{-}\omega.
\end{aligned} \tag{3}
\]
这些是 U 上形式的等式；不是仅在某个数值点验证。系数矩阵可逆，所以六个反对称形式正好张成六个已知单位杯积方向。

其留数由 [REGULATOR_COMPARISON](REGULATOR_COMPARISON.md) 的 Z、W 给出。对称形式的留数也能无需数值积分而确定：在任一双点，
\(\alpha_i/\alpha_j\) 和其逆不可能同时沿两条分量具有极点。如果前者沿两支的阶数均为 −1，则 q⁺ 和 q⁻ 的双留数相同；如果逆有双极点，则互为相反数；两者都没有双极点则双留数均为零。

由此恢复所有三十行。对称形式留数秩为五，六列之和为零；加上六个反对称列后秩为十一。原始单位 valuation、杯积的六个形式恒等式、全部留数以及基矩阵行列式均已作精确核对，见 [log_form_frame.json](data/log_form_frame.json)。

## 三、质量方向的显式 Hodge 提升

令 F_mu 为 KUMMER_BOUNDARY 中的 30×4 质量旗标矩阵。定义
\[
 \boxed{\beta_i=-\mu_i\alpha_i N\,\omega
       =-\frac{s\mu_i\alpha_i}{M}\,\omega
       \quad\text{on }U.} \tag{4}
\]
它们属于 (2) 的线性跨度，因而是全局对数二形式。将
\(-\mu_i\alpha_iN=-\mu_i-\sum_{j\ne i}\mu_i\alpha_i/\alpha_j\)
代入前节的全部留数，得到精确等式
\[
 \rho(\beta_i)=(F_\mu)_{\bullet i},\qquad
 \sum_i\beta_i=-s\omega. \tag{5}
\]
所以对任何 \(r=(r_1,\dots,r_4)\in\mathbb Q^4\)，
\(\beta_r=\sum r_i\beta_i\) 是图循环 \(c_r=F_\mu r\) 的 F² 提升。
常数向量 r 给出纯全纯形式，不改变图循环。

特别地，对 m=(1,2,4,8)、r=(3,2,1,0)，得到前稿的活动方向 c=3P₁+2P₂+P₃，且在 s=−64 有
\[
 \boxed{\beta_c=
 64\,\frac{3a+8b+16c}{a+4b+16c+64d}\,\omega.} \tag{6}
\]
这给出实际 Hodge 过滤中的具体提升，但没有给出它与有理 Betti 提升的全部比较积分。两种提升之差，模掉 \(F^2A_{\mathbb C}\) 与适当扭转的有理格，才是下层扩张的 regulator。

## 四、Gauss–Manin 的质量导数恒等式

把 μ_i 与 s 暂时都作为参数。因为 f_s=−xyz，而
\(f_{\mu_i}=\alpha_i xyzN\)，对留数形式求导可得
\[
 \boxed{\nabla_s\beta_i=\mu_i\nabla_{\mu_i}\omega.} \tag{7}
\]
证明可在取留数前完成：
\[
 \partial_s\left(\frac{-\mu_i\alpha_i N}{f}\right)
 =-\frac{\mu_i\alpha_i xyzN}{f^2}
 =\mu_i\,\partial_{\mu_i}\left(\frac1f\right).
\]
取留数和限制到 U 与 Gauss–Manin 连接相容，所以 (7) 是上同调恒等式。因 (5) 的图循环不随参数变化，左侧已经落在未扭转纯核 K 中。作为独立一致性检查，求和得到
\[
 \sum_i\mu_i\nabla_{\mu_i}\omega=-\omega-s\nabla_s\omega,
\]
正是所有五个系数共同缩放时 ω 的齐次次数 −1。

以下非消失证明只使用 (7)、Hodge 型和参数族的局部 Torelli 性，不需要先恢复六阶 Picard–Fuchs 算子的巨大系数。

## 五、四维参数族的周期映射微分单射

使用已固定的 Sylvester 坐标
\[
 y_i=\mu_i\alpha_i\ (i\le4),\qquad y_5=-M,\qquad
 \sum_{i=1}^5y_i=0,\qquad
 Q_\lambda=\sum_{i=1}^5\lambda_i\prod_{j\ne i}y_j,
\]
其中 \(\lambda=(\mu_1,\mu_2,\mu_3,\mu_4,s)\)，只计共同缩放后的 P⁴ 参数。在只有十个普通节点的开集上，这些节点恒为
\[
 p_{ij}=[e_i-e_j]\quad(1\le i<j\le5).
\]
吹起十个固定节点给出光滑 K3 族。

**引理。** 此四维族的 Kodaira–Spencer 映射单射。

**证明。** 若某个一阶参数方向使抽象 K3 变形平凡，则因 H¹(O_X)=0，四次模型的拉回超平面线丛也唯一地平凡化。它的四个截面恢复到原节点四次曲面的映射，所以该一阶四次曲面变形必由一个无穷小 PGL₄ 变换产生。这里截面的基变换合法：节点为有理双点，解消的 O_X 高阶直接像为零；四次超曲面的标准正合列给出 h⁰(O_X(1))=4、h¹(O_X(1))=0，故一阶底变换下截面仍为自由秩四。

每条节点例外曲线为 (−2) 曲线，其法丛为 O_(P¹)(−2)，没有一阶移动。因此在上述平凡化下，这十条曲线及其收缩点都保持不动。无穷小 PGL₄ 变换须逐点固定十个 p_ij。

用 p_i5 (i=1,…,4) 作射影基点，它们迫使线性变换为对角阵；固定 p_12、p_13、p_14 又迫使四个对角元相同。因此在 PGL₄ 中变换为零。相应线性条件矩阵的秩为十五，已经精确核对。

最后，五个 \(\prod_{j\ne i}y_j\) 在超平面 \(\sum y_i=0\) 上线性独立：取 y_i=0 且其余坐标非零，便能单独读出第 i 个系数。因此一阶参数改变只能是 Q 的整体缩放，在 P⁴ 中为零。证毕。

由 K3 局部 Torelli，以上映射与周期映射微分通过全纯二形式的收缩相联系，故周期映射微分也单射。所用标准结果为 [Huybrechts，Chapter 6，Proposition 2.4 和 2.8](https://www.math.uni-bonn.de/people/huybrech/K3Global.pdf)。

这个论证使用带标记参数族，不要求粗模空间没有有限自同构。各例外曲线的 H⁰(O(−2))=0 排除的是一阶移动，有限置换不产生无穷小稳定子。α 与 y 的坐标变换虽依赖质量，却只改变四次曲面的射影平凡化；其 Kodaira–Spencer 类不变。全纯形式的相应标量归一化至多给一阶导数添加 F²K 方向，不影响以下模 F² 的线性无关性。

## 六、非恒定质量方向的无穷小障碍非零

现在固定质量，仅令 s 变化。记
\[
 v_s=(0,0,0,0,1),\qquad
 v_r=(r_1\mu_1,r_2\mu_2,r_3\mu_3,r_4\mu_4,0)
\]
为五个系数空间中的向量。若 r 非常数向量，则 v_r 不在
\(\langle v_s,\lambda\rangle\) 中，因为等式
\(v_r=a\lambda+bv_s\) 的前四个坐标强迫所有 r_i=a。

所以在 \(\operatorname{Gr}_F^1 K\) 中，
\(\nabla_{v_r}\omega\) 与 \(\nabla_s\omega\) 的像线性无关。结合 (7)：
\[
 \delta_r:=[\nabla_s\beta_r]
 \quad\text{in}\quad
 \frac{F^1K}{F^2K+\mathcal O\cdot\nabla_s\omega}
 \quad\text{不为零}. \tag{8}
\]
这里分母中的连接项按其在 F¹ 中的像理解；\(\omega\) 张成 F²K。十六个边界除子类为平坦截面，故这些纯导数仍与边界正交。

若图循环 c_r 的扩张在某个开圆盘上存在平坦 Hodge 截面 q，则必有
\[
 q=\beta_r+g(s)\omega
\]
（扭转 (2) 后顶层为 Q(0)）。对 s 求导并模去 F²K，得到
\([\nabla_s\beta_r]+g(s)[\nabla_s\omega]=0\)，与 (8) 矛盾。
因此任何非常数有理 r 所对应的下层扩张，都不可能作为 Hodge 变分在一个非空开集上分裂。

这证明一个变分非消失结果；它仍不排除个别纤维的扩张类恰为零。

## 七、很一般固定纤维与形式周期最小性

在一个单连通参数圆盘上，固定 c_r 的有理 Betti 提升只有可数多个。对每个平坦有理提升 q，“q 位于 F² 提升空间”是一个解析条件；若在整个圆盘上成立，则给出前节排除的平坦 Hodge 截面。因此其成立集合是离散集。

取可数圆盘覆盖、可数有理提升，再对所有非零有理质量方向 r mod Q(1,1,1,1) 取并，得到一个可数例外集。在其外，
\[
 B_{\rm mass}\longrightarrow
 \operatorname{Ext}^1_{\rm MHS}(\mathbb Q(0),K_s(2))
\]
为单射。特别地，在其中任何实际活动子空间上亦单射。Hodge 非零蕴含 Nori 扩张非零；这里没有使用逆方向。

质量 (1,2,4,8) 和 (1,3,7,15) 的射线已各有一个 rho=16 样本，由 [Picard 证书](PICARD_CERTIFICATE.md) 确立泛纤维秩十六。额外有理 (1,1) 类的出现也在一个可数例外集内。结合 [MINIMALITY](MINIMALITY.md) 的链生成性和最小商判据，得到：

**定理。** 在上述两条固定质量射线上，对很一般的欧氏参数 s<0，物理形式矩阵系数的最小 Hodge 对象秩分别为
\[
 \boxed{8\quad\text{和}\quad10.}
\]
其活动单位秩分别为 κ=1 和 κ=3。将该纤维定义在包含质量与 s 的实嵌入特征零基域上，既定实际 Nori 实现也有同样最小秩：Hodge 下界与秩 7+κ 的实际实现上界相合。所有局部化和比较构造均在这种基域上成立，不要求 s 为代数数。

这不是“所有代数参数成立”的定理。可数例外集可能包含给定的代数点；本证明不判定 s=−64 的 regulator 非零。它也不是关于 C(s) 上标量微分方程阶数的定理，后者会把质量对数视为复常数，须另作分析。

## 验证责任和后续工作

[verify_log_forms.py](code/verify_log_forms.py) 核对六个杯积形式恒等式、三十处留数、十二维框架、质量提升、导数恒等式及十节点无穷小射影稳定子的秩。程序输出明确保留 fixed_basepoint_regulator_nonvanishing_certified=false。局部 Torelli、支持比较及可数例外集论证由本稿和前稿承担。

下一项真正的固定纤维问题是：用式 (6) 与明确 Betti 提升计算扩张类，并排除它在 \(F^2K_{\mathbb C}+K_{\mathbb Q}(2)\) 中消失。单靠 (8) 在 s=−64 处非零，不足以完成这一步。这里按既定比较同构把扭转有理格放入未扭转的复空间。Q6 的全局闭路、循环标记和绝对比较仍独立保留；边界格本原性现由 [LATTICE_SATURATION](LATTICE_SATURATION.md) 另证。
