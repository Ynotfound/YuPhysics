# 边界一动机：Tate 框架、全部旗标与实际 Kummer 比较

> 2026-10-14专项状态：旧假设11.2的逐项依赖与实际物理双标架等式现见HYPOTHESIS_11_2_AUDIT.tex（PDF同名）。本稿后部的部分“待证”措辞记录撰写时状态；固定下层非零、P8及固定群的后续证据分别在LOWER_REGULATOR、RELATIVE_TOP_COMPARISON与FIXED_PURE_GROUP，不由本稿或新补篇重新数值认证。新补篇未将结论推广到所有参数，也尚未回写统一全文。

本稿在已经固定的二十曲线模型和有理 Nori 范畴中，计算实际秩九边界一动机。证明分为三步：先通过几何留数识别分级对象，再计算符号的三十个旗标，最后用实际导出态射比较扩张。有限矩阵计算只核对第二步及基变换；动机比较由第三步承担。

## 范围与结论

设 k 是实数域中的数域，mu_i>0、s<0 属于 k。沿用 [几何比较稿](GEOMETRIC_COMPARISON.md) 与 [纯核稿](K3_KERNEL.md)：

\[
 T=\mathbb G_m^3,\quad U_s\subset X_s,\quad
 \mathsf E=H^3_{\mathrm N}(T,U_s)(3),\quad
 \mathsf A=\mathsf H_{\mathrm{res}}(3),\quad
 \mathsf K=\mathsf E/\mathsf A.
\]

物理相对同调对象通过已经构造的 Nori 对偶同构识别为这里的 E。各秩为 15、6、9；此处不要求 Picard 数恰为十六。

**主命题。** 下文指定环面余字符框架及顶层生成元后，在有理一动机范畴中有

\[
 \boxed{\mathsf K\simeq
 [\mathbb Z\longrightarrow\mathbb G_m^8]\otimes\mathbb Q,\qquad
 1\longmapsto(\mu_1/\mu_4,\mu_2/\mu_4,\mu_3/\mu_4,1,1,1,1,1).}
\]

框架中的整体符号必须保留：下面的旗标矩阵使用“先沿有序边的第一个分量留数”；由此得到的八维留数基记为 P。一动机同调的余字符框架为 **−P**，理由见最后的单单位校准。若直接用 P 作为余字符框架，三个质量比须同时取逆；不能一边换基、一边保留原来的点坐标。

这识别了边界商及其实际扩张。它不计算纯 K3 比较矩阵、不确定整个动机群，也不证明物理对象的最小性。

## 一、实际留数映射排除算术扭曲

将二十条曲线按 V1,…,V4，E12,…,E34，C12,…,C34，F1,…,F4 排序；两位下标按字典序。每个交点对应从序号小的分量指向序号大的分量的边。设

\[
 \partial:\mathbb Q^{30}\to\mathbb Q^{20}
\]

为入点 +1、出点 −1 的关联矩阵。

这些曲线及交点均定义在 k 上。原四个顶点和六个边节点是 k 有理点；各严格变换与例外曲线也定义在 k 上。交点由两条定义在 k 上的曲线横截相交且交数为一给出，因此交点的剩余域为 k。这里的结论来自明确几何，不从单个复 Hodge 实现推断下降。

对 D_i，令 X_i=X_s\setminus\bigcup_{j\ne i}D_j，D_i^\circ=D_i\cap X_i。先对光滑除子 D_i^\circ⊂X_i 用 Gysin 连接，再对曲线 D_i 上删去的交点用连接，得到实际 Nori 态射

\[
 H^2_{\mathrm N}(U_s)\longrightarrow
 H^1_{\mathrm N}(D_i^\circ)(-1)\longrightarrow
 \bigoplus_{p\in D_i\setminus D_i^\circ}\mathbb Q(-2).
\]

这些箭头由局部化与纯性构造。Nori 实现与六运算、Hodge 实现的相容性采用 [Tubach，Theorem 4.5 与 Corollary 4.8，v4](https://arxiv.org/html/2309.11999v4)。基域上的通常 Nori 范畴识别沿用纯核稿的约定。

在每条有向边选择第一个分量的旗标，合成

\[
 \rho:H^2_{\mathrm N}(U_s)\longrightarrow\mathbb Q(-2)^{30}.
\]

约定 Res(dlog t∧eta)=eta；交换两次留数的顺序产生负号。每条完备有理曲线上的总留数为零，所以 ∂ρ=0。严格性及 HODGE.md 的分层给出 Betti 层面的正合列

\[
 0\to A_s\to H^2(U_s)\xrightarrow{\rho}
 \ker(\partial)(-2)\to0.
\]

这里的映射已经在 Nori 范畴构造；Betti 实现忠实正合，故这也是实际 Nori 正合列，其核为已经识别的 H_res。无需提升任意 Hodge 子对象。

令 v_i 是 a/d、b/d、c/d 在 D_i 上的阶数。环境环面二次上同调的三个生成元在交点 e=(i,j) 的留数为

\[
 Z_{e,(p,q)}=v_i^pv_j^q-v_i^qv_j^p,\qquad 1\le p<q\le3.
\]

取 R 为抽取 (V1,C12)、(V1,C13)、(V2,C23) 三行的矩阵。精确地，RZ=I3，∂Z=0。因此

\[
 C:=\ker\partial\cap\ker R,\qquad \dim C=8,
 \qquad \operatorname{pr}=I-ZR
\]

将 ker∂/im Z 明确识别为 C。扭转 (3) 并对相对长正合列取商，得到实际短正合列

\[
 0\to C\otimes\mathbb Q(1)\to\mathsf K\to\mathbb Q(0)\to0.
\]

顶层的 Q(0) 是实际 H^3_N(T)(3)，其生成元由有序坐标单位的杯积固定。环面层是八个实际分裂 Tate 对象。结合此前最大一动机商的构造，K 的格和环面均已在 k 上分裂；Artin 扭曲不再是本模型的待证项。

## 二、三十个旗标的独立计算

取单位符号

\[
 \xi=\{a/d,b/d,c/d\}\in H^3_M(U_s,\mathbb Q(3)),
 \qquad \lambda_i=[\mu_i]\in k^\times\otimes\mathbb Q.
\]

这里的符号是三个全局单位类的有序杯积。tame 映射用 ∂{pi,u2,…}= {u2,…} 的约定。

在面分量 F_i 上，令 j<k<l 为其他三个下标，epsilon_i=(−1)^(i−1)，并置

\[
 u=-\frac{\mu_j\alpha_j}{\mu_l\alpha_l},\quad
 1-u=-\frac{\mu_k\alpha_k}{\mu_l\alpha_l},\quad
 A=\mu_l/\mu_j,\quad B=\mu_l/\mu_k.
\]

模去由 −1 产生的挠元，有

\[
 \partial_{F_i}\xi
 =\epsilon_i\big(\{A,B\}+\{A,1-u\}+\{u,B\}\big).
\]

Steinberg 关系消去 {u,1−u}。在 u=0、1、∞ 的再次留数分别为

\[
 \epsilon_i(\lambda_l-\lambda_k),\quad
 \epsilon_i(\lambda_j-\lambda_l),\quad
 \epsilon_i(\lambda_k-\lambda_j).
\]

这三个点依次是与 E_kl、E_jl、E_jk 的交点。

其余零项也需要检查。在 C_ij 的一般点，两个消失坐标之比为 −1，第一次 tame 符号有一个挠元因子，故有理系数下为零。在 V1 的例外圆锥上，若 y=b/d、z=c/d，则 yz+y+z=0；相应的 {y,z} 由 z=−y/(1+y) 和 Steinberg 关系归零。其他顶点经坐标置换相同。因此所有 V–C 和 E–C 旗标为零。节点 E_ij 处的另外两个旗标由交换次序的负号确定；这也给出独立的相容性检验。

按前述边定向，E 位于 F 之前，故表中数值为 **负的面优先留数**。全部十二个非零行如下，未列出的十八行均为零：

| 有向边 | 二次 tame 值 |
|---|---|
| E12 → F3 | lambda1 − lambda2 |
| E12 → F4 | −lambda1 + lambda2 |
| E13 → F2 | −lambda1 + lambda3 |
| E13 → F4 | lambda1 − lambda3 |
| E14 → F2 | lambda1 − lambda4 |
| E14 → F3 | −lambda1 + lambda4 |
| E23 → F1 | lambda2 − lambda3 |
| E23 → F4 | −lambda2 + lambda3 |
| E24 → F1 | −lambda2 + lambda4 |
| E24 → F3 | lambda2 − lambda4 |
| E34 → F1 | lambda3 − lambda4 |
| E34 → F2 | −lambda3 + lambda4 |

把此表记为 F_mu lambda，则

\[
 \partial F_\mu=0,\quad RF_\mu=0,\quad
 \operatorname{rank}F_\mu=3,\quad F_\mu(1,1,1,1)^t=0.
\]

几何上，非零支撑是四个 F 顶点与六个 E 中点组成的细分四面体图，它的循环秩为三。不能把该秩与八维相对 Tate 层或质量单位的实际乘法秩混为一谈。

取 F_mu 的前三列作为 P 的前三列，再按随附程序给出的确定顺序补成 C 的八维基。P 为整数矩阵，以下八行的子式行列式为 1：

\[
 (V2,C12),(V3,C13),(V3,C23),(V4,C14),(V4,C24),
 (E12,F3),(E13,F2),(E14,F2).
\]

因而这也是 C∩Z30 的整数基：任一整数循环的基坐标由这个幺模子式计算，仍为整数。由于 RZ=I3，环境三方向的商同样没有由本次基变换引入的有限指数问题。这个边界结论与 Q6 的指数二补全无关。

在基 P 中

\[
 F_\mu\lambda
 =P(\lambda_1-\lambda_4,\lambda_2-\lambda_4,
       \lambda_3-\lambda_4,0,0,0,0,0)^t.
\]

## 三、从单位符号到实际相对扩张

仅计算上一表不能识别 K。本节补上实际态射。

**同伦纤维引理。** 在一个阿贝尔范畴的导出范畴内，设 V→W，H^0(W)=0，H^0(V)=Q；取 h:Q→V 诱导 H^0 上的恒等。令

\[
 F=\operatorname{Fib}(V\to W),\qquad
 F_h=\operatorname{Fib}(Q\xrightarrow{h}V\to W).
\]

则自然态射 F_h→F 在零次上同调诱导交换图

\[
 \begin{array}{ccccccccc}
 0&\to&H^{-1}(W)&\to&H^0(F_h)&\to&Q&\to&0\\
 &&\downarrow&&\downarrow&&\Vert\\
 0&\to&\operatorname{coker}(H^{-1}(V)\to H^{-1}(W))
 &\to&H^0(F)&\to&Q&\to&0.
 \end{array}
\]

下行是上行沿左侧商映射的推出。

**证明。** 对两个纤维三角取长正合列。H^(−1)(Q)=0 和 H^0(W)=0 给出上行短正合列；下行直接由 V→W 的长正合列得到。纤维构造的自然性给出图。沿左侧商映射推出上行后，所得短正合列到下行在两端均为恒等，因此中间态射是同构。此证明不使用任何 Ext² 群整体消失。

现在在基域上的 Nori 导出范畴中令

\[
 V=R\Gamma_{\mathrm N}(T)(3)[3],\quad
 W=R\Gamma_{\mathrm N}(U_s)(3)[3],\quad Q=\mathbb Q(0).
\]

取 h 为 **三个坐标单位 a/d、b/d、c/d 的有序 motivic 杯积** 经 Nori 实现所得的态射。这个类并非根据 Hodge 实现猜测出来：单位类位于 H^1_M(Gm,Q(1))，乘积位于 H^3_M(T,Q(3))，再使用与六运算相容的实现。

将 H^3_N(T)(3) 的生成元固定为该杯积的 Betti 类，h 在 H^0 上就是恒等。等价地，以各单位圆的正定向乘积归一化，三个 dlog/(2pi i) 类的积分乘积为一。这里没有额外 3! 或未定的非零有理尺度。

U_s 为光滑仿射曲面，H^3(U_s)=0，因此 H^0(W)=0。由忠实正合的 Betti 实现，这个消失也成立于 Nori 范畴。引理给出实际符号扩张 M_xi 及图

\[
 \begin{array}{ccccccccc}
 0&\to&H^2_{\mathrm N}(U_s)(3)&\to&\mathsf M_\xi&\to&\mathbb Q(0)&\to&0\\
 &&\downarrow&&\downarrow&&\Vert\\
 0&\to&\mathsf N&\to&\mathsf E&\to&\mathbb Q(0)&\to&0,
 \end{array}
\]

其中 N=coker(H^2_N(T)→H^2_N(U_s))(3)。上、下对象的秩分别是 18 和 15，核差为三个环境环面方向。

这给出了 [BKV，式 6.2.9](https://arxiv.org/html/1406.2664v3#S6.SS2) 所指比较在本模型中的具体归一化。我们使用纤维引理重新构造箭头，并未把该文保留的尺度直接当作已固定。物理积分的相对周期歧义属于提升和配对问题，不是此处已固定扩张类的额外自由参数。

## 四、局部化相容性计算推出后的点

前面的两次几何连接在复形层面给出

\[
 W\longrightarrow\mathbb Q(1)^{30}[1].
\]

与 h 的限制复合，即两个 motivic 局部化连接作用于单位符号。其值是对应的二重 tame 单位。这里使用的是乘积单位类、纯性与局部化的相容性；不是根据数值对数识别一个未知类。Milnor 边界作为 motivic 局部化边界的说明可见 [Geisser，Motivic Cohomology, K-Theory and Topological Cyclic Homology，§1](https://www2.rikkyo.ac.jp/web/geisser/Handbook.pdf)。该文采用将 uniformizer 放在末位的约定；本稿将它移到首位时加入相应的分次交换符号。

令 S 是数据文件给出的 P 的左逆，SP=I8。复合 S(I−ZR)ρ 在 H^(−1)(W) 上就是此前的边界商映射到 Q(1)^8。对同伦纤维应用这个复形态射，得到 M_xi 的推出扩张。由上一节的交换图和推出的结合律，它正是 K；其八个扩张类是单单位类的实现

\[
 \big(\kappa_{\mathrm N}(\mu_1/\mu_4),
 \kappa_{\mathrm N}(\mu_2/\mu_4),
 \kappa_{\mathrm N}(\mu_3/\mu_4),0,0,0,0,0\big).
\]

此处 κ_N(v) 明确定义为单位类 [v]:Q(0)→Q(1)[1] 的 Nori 实现，并用其同伦纤维的 H^0 表示扩张。下述校准说明如何把这个约定翻译为一动机的点。

**单单位校准。** 在有序点对 {1,v}⊂Gm 中，用第二点减第一点识别 coker(H^0(Gm)→H^0({1,v}))。固定相对连接的符号，使边界函数沿路径 1→v 的配对为 f(v)−f(1)；在相对 de Rham 复形中可写 d(alpha,beta)=(d alpha,−i*alpha−d beta)，配对为路径积分加边界取值。相对长正合列及上述纤维引理将单位 v 的限制减去单位 1 的限制识别为 H^1_N(Gm,{1,v})(1) 的扩张。取从 1 到 v 的路径与正定向圆周，dlog t 的积分分别为 log v、2pi i。按 HODGE.md 的 ν=e_F−e_B 约定，该扩张的坐标为

\[
 \nu(\kappa_{\mathrm N}(v))=+\log v
 \quad\pmod{2\pi i\mathbb Q}.
\]

另一方面，一动机 [Z→Gm,1↦v] 的同调实现由 exp z=v^n 的对 (n,z) 描述；有理提升为 (1,log v)，F^0 提升的 Lie 坐标为零，所以同调余字符的通常正圆周基下 ν=−log v。故两者以 **负的** Tate 识别相对应。这个符号是统一的，不依赖 v 或质量。

因此，在余字符框架 −P 下，K 的实际单位点就是主命题所写的三个质量比及五个 1。可用一动机 Hodge 实现验证最后这个符号：k^×⊗Q 到 C/(2pi i Q) 的对数映射单射，因为落入 2pi i Q 意味着某次幂为 1。这里只用该单射核对已知分裂一动机的单位点；没有把 Hodge 全忠实性外推到 K3 混合扩张。v=1 的分裂情形通过零单位类处理。

于是主命题成立。符号的原始基 P、同调余字符基 −P、ν=e_F−e_B，以及物理配对中的 −ν 已分别固定，不能混用。

## 五、物理活动因子化的推论与边界

令

\[
 \kappa=\dim_{\mathbb Q}\langle[\mu_1/\mu_4],
 [\mu_2/\mu_4],[\mu_3/\mu_4]\rangle\subset k^\times\otimes\mathbb Q.
\]

取实际点生成的最小连通子环面 T_mu⊂Gm8，令 K_act=[Z→T_mu]⊗Q，并取 E_act=E×_K K_act。正质量排除点在实嵌入下非平凡的根单位值；同源范畴中有限挠元也不影响构造。则

\[
 \operatorname{rank}T_\mu=\kappa,\qquad
 \operatorname{rank}\mathsf E_{\mathrm{act}}=7+\kappa,\qquad
 \mathsf E/\mathsf E_{\mathrm{act}}\simeq\mathbb Q(1)^{8-\kappa}.
\]

这些商在 k 上分裂且与实结构相容。物理正链 sigma_B 被几何复共轭固定，而其商上的共轭为 −1；故 sigma_B 在商中为零。由此，原欧氏积分的动机周期通过 E_act 的对偶实现。

对 m=(1,2,4,8)，平方质量比为 (1/64,1/16,1/4)，单位秩 κ=1，得到秩八的实际活动实现。这是**实现上界**，不是最小性或八个周期独立性的证明。纯核在特殊纤维上还可能有不可见的代数方向。

[PICARD_CERTIFICATE](PICARD_CERTIFICATE.md) 已在 s=−64 排除额外代数方向。[MINIMALITY](MINIMALITY.md) 随后证明物理链生成整个活动对象，并将秩八最小性等价地化为下层扩张 eta_act 非零。该非零性仍未证明；本稿的 log 2 非零只控制上层 Kummer 类。

实际三步对象 E_act 同时实现相邻扩张，因此相应 eta_kappa∘xi_kappa=0 在 Nori Ext² 中成立。理由是 xi_kappa 经 E_act 提升到 W_(−2)E_act 的扩张，与纯核稿对环境对象的论证相同；不声称整个 Ext² 群为零。

## 可重放核对及尚未完成的工作

`code/verify_boundary_kummer.py` 重建关联矩阵、环境留数 Z、全部三十个旗标、整数基 P 和左逆，检查 ∂F_mu=0、RF_mu=0、rank F_mu=3、共同质量缩放不变性及八维坐标。故意反转一个非零旗标会破坏全局相容性。结果保存为 `data/boundary_kummer_flags.json`。

数据文件中的 actual_extension_comparison_certified=false 是有意的：有限线性代数程序不承担本稿第三、四节的动机证明。理论推导与计算证据的责任应保持明确。

后续 REGULATOR_COMPARISON 与 DE_RHAM_MINIMALITY 已给出下层支持比较、三个 Tate 直和因子、显式质量提升和两条射线上很一般参数的最小性证明稿。仍需将各稿合并进统一正文，证明指定 s=−64 的非零性，识别纯 K3 与扩张群的实际关系，并计算指定比较条目。Q6 的全局运输、绝对周期归一化和几何标记仍是独立未完成任务。
