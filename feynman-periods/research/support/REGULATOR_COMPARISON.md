# 下层扩张的支持上同调比较与三个 Tate 直和因子

本稿完成 [MINIMALITY](MINIMALITY.md) 中候选边界循环与下层扩张的比较：在既定的实际 Nori 实现中，两者通过支持上同调的局部化三角联系。另由四个全局单位的二次杯积，给出三个可显式分裂的 Tate 方向。后续 [DE_RHAM_MINIMALITY](DE_RHAM_MINIMALITY.md) 给出很一般参数非零性的证明稿；固定 s=−64 的非零性现由LOWER_REGULATOR独立认证，结合MINIMALITY得到最小秩八；该新证据不能省去。

## 一、对象、方向与结论

沿用二十曲线模型 \(U=X\setminus D\)。令
\[
 A_2=\mathsf H_{\rm res}(2),\qquad
 C_{11}=\ker(\partial:\mathbb Q^{30}\to\mathbb Q^{20}),
\]
其中每条有向边 e=(D_a,D_b) 满足 a<b，\(\partial e=-D_a+D_b\)。
HODGE 和 KUMMER_BOUNDARY 的实际 Nori 正合列扭转 (2) 后为
\[
 0\to A_2\to H^2_{\mathrm N}(U)(2)
       \xrightarrow{\rho}C_{11}\otimes\mathbb Q(0)\to0. \tag{1}
\]
记其扩张类为 \(\eta_{11}\)。对每个图循环 c，选择各 D_i 上的有理函数 \(f_i\)，使
\[
 \operatorname{div}_{D_i}f_i
  =-\sum_{e\ni i}\partial_{i,e}c_e\,p_e. \tag{2}
\]
每个 D_i 是有 k 有理点的光滑有理曲线，右端次数为零，故有理函数存在。有理系数 c 可先清分母再张量 Q。两支除子相消，因而
\[
 Z_c=\sum_i(D_i,f_i)\in\operatorname{CH}^2(X,1)_{\mathbb Q}.
\]

**比较命题。** 采用下述支持复形和局部留数约定，存在等式
\[
 \boxed{c^*\eta_{11}
 =r_{\mathrm N,res}(Z_c)
 \quad\text{in}\quad
 \operatorname{Ext}^1_{\mathrm N}(\mathbb Q(0),A_2).} \tag{3}
\]
这里右侧定义为 motivic 循环类经 Nori 实现、再经边界正交 Chow 投影子所得的类。应用 Hodge 实现，得到同样的 Hodge regulator 等式。式 (3) 是实际态射的比较，不是由 Ext 群的维数猜出一个同构。

所需实现框架是 Voevodsky 有理动机到 Nori 导出范畴、并到混合 Hodge 结构的实现，与六运算相容。使用 [Tubach，Theorem 0.1、4.5 和 Corollary 4.8](https://arxiv.org/html/2309.11999v4)。不使用其依赖 motivic t-structure 猜想的等价定理，也不假设 Hodge 实现对一般扩张全忠实。

## 二、支持复形给出循环的提升

令 \(i:D\hookrightarrow X\)、\(j:U\hookrightarrow X\)。支持上同调由
\(R\Gamma_D(X)=R\Gamma(X,i_*i^!\mathbf1_X)\) 定义。
因为无三重点，Mayer–Vietoris 与纯性给出
\[
 \bigoplus_e\mathbb Q(0)[-4]
 \xrightarrow{\partial_{\rm Gys}}
 \bigoplus_i R\Gamma(D_i)(1)[-2]
 \longrightarrow R\Gamma_D(X)(2)
 \longrightarrow\bigoplus_e\mathbb Q(0)[-3]. \tag{4}
\]
第一个箭头在有向边 e=(a,b) 上是“负的 a 支 Gysin，加上 b 支 Gysin”。所有扭转已经写在式内：点在曲面中的余维为二，曲线余维为一。

记 \(S=R\Gamma_D(X)(2)[3]\)。由 (4) 得到
\[
 H^{-1}(S)=\mathbb Q(1)^{20},\qquad
 H^0(S)=C_{11}\otimes\mathbb Q(0). \tag{5}
\]
第一项映到 \(H^2(X)(2)\) 的像是边界除子跨度；第二项由点至分量的关联核给出。不能把 S 整体误写为 C₁₁，因为它还有其他上同调。

审阅补充：关联余核还给出 H¹(S)=Q(0)，它对应整个连通边界的顶端支持类。第三节只使用 H^(−1) 和 H⁰ 的长正合列；并未截去 S 后假定它等于图循环商。Fib(g) 的目标 A₂[1] 由实际 Chow 投影子选出，故这一额外 H¹ 不改变式 (7)。

在 motivic cohomology 层面对 (4)[3] 取 \(\operatorname{Hom}(\mathbb Q(0),-)\)，相关部分为
\[
 \bigoplus_i k^\times\otimes\mathbb Q
 \longrightarrow H^3_{M,D}(X,\mathbb Q(2))
 \longrightarrow \mathbb Q^{30}
 \xrightarrow{\partial}\mathbb Q^{20}. \tag{6}
\]
这里最后的 Q²⁰ 是各 \(\operatorname{Pic}(D_i)_{\mathbb Q}\)，用次数识别。故 c∈ker∂ 可提升为
\[
 \alpha_c:\mathbb Q(0)\longrightarrow S_M.
\]
更具体地，在 Cone(∂) 的模型中，其循环条件是
\(d(f_i)+\partial c=0\)，恰为式 (2)。
将支持遗忘到 X，函数分量推出为 \(\sum(D_i,f_i)\)，故 \(\alpha_c\) 的像就是 Z_c 的 motivic 循环类。

这也说明为什么“除子相消”尚不证明类非零：支持遗忘映射可以有核。曲线函数描述和 tame 边界关系可对照 [Sato，§2.1](https://arxiv.org/html/2408.09102#S2.SS1)；此处额外保留了支持和图循环的提升信息。

## 三、同伦纤维给出实际扩张比较

现在全部对象取 Nori 实现。令
\[
 V=R\Gamma_{\mathrm N}(X)(2)[3],\qquad
 W=R\Gamma_{\mathrm N}(U)(2)[3].
\]
局部化给出三角 \(S\to V\to W\to S[1]\)。
由 [K3_KERNEL](K3_KERNEL.md) 的实际 Chow 投影子，V 的边界正交直和因子只有 H^(−1) 非零，故是 \(A_2[1]\)。记
\[
 g:S\longrightarrow V\longrightarrow A_2[1].
\]
它在 H^(−1) 上为零，因为边界类被投影子杀死。由同伦纤维长正合列得到
\[
 0\to A_2\to H^0\operatorname{Fib}(g)
       \to C_{11}\otimes\mathbb Q(0)\to0. \tag{7}
\]

自然映射
\(\operatorname{Fib}(S\to V)\to\operatorname{Fib}(g)\)
在 H⁰ 上给出
\[
 H^2_{\mathrm N}(U)(2)\longrightarrow H^0\operatorname{Fib}(g).
\]
在纯核上，它是由限制映射和 \(\pi_{\rm res}\) 固定的同构；在图循环商上为恒等。因此它是正合列 (1) 到 (7) 的同构。

将 \(\alpha_c\) 实现后复合 g，得到
\[
 g\alpha_c:\mathbb Q(0)\longrightarrow A_2[1].
\]
同伦纤维的自然性给出图
\[
\begin{array}{ccccccccc}
0&\to&A_2&\to&H^0\operatorname{Fib}(g\alpha_c)&\to&\mathbb Q(0)&\to&0\\
&&\Vert&&\downarrow&&\downarrow c\\
0&\to&A_2&\to&H^0\operatorname{Fib}(g)&\to&C_{11}\otimes\mathbb Q(0)&\to&0.
\end{array}
\]
该图是拉回图。因此上行的类是 \(c^*\eta_{11}\)。
另一方面，支持遗忘再投影所得 \(g\alpha_c\) 正是 \(r_{\mathrm N,res}(Z_c)\) 的定义，故得到式 (3)。整个论证先构造态射再使用实现；不从抽象 Hodge 同构倒推动机箭头。

**选择无关性。** 改变具有同一除子的 f_i 只会乘以 k× 常数，改变 Z_c 一个分解循环 \((D_i,a_i)\)。其投影 regulator 是单位类与被 \(\pi_{\rm res}\) 杀死的边界除子类之积，因此为零。即使 Z_c 在 Chow 群中依赖选择，式 (3) 的边界正交部分不依赖选择。

## 四、符号、Tate 扭转与物理下层的对应

在节点附近取坐标 x,y，使先后两个分量分别为 x=0、y=0。规定
\[
 \operatorname{Res}_{x=0}(d\log x\wedge d\log y)=d\log y,\qquad
 \operatorname{Res}_{y=0}(d\log x\wedge d\log y)=-d\log x.
\]
依次取点留数，第一支得 +1，第二支得 −1。这与式 (2) 的“第一支除子 +c_e，第二支 −c_e”一致，也与 KUMMER_BOUNDARY 的基 P 一致。
纯性使用 \(d\log/(2\pi i)\) 在正单位圆上积分为一的归一化；双留数不额外加入 2πi 或阶乘。Ext 类采用连接态射 f 与 \(\operatorname{Fib}(f)\) 对应的约定。

环面三列 Z 的实际杯积给出正合列 (1) 在 im Z 上的截面。取商并用 \(I-ZR\) 识别，正好得到
\[
 N(-1)=\operatorname{coker}(H^2_{\mathrm N}(T)(2)\to H^2_{\mathrm N}(U)(2)).
\]
对 \(c\in C=\ker\partial\cap\ker R\)，其拉回和 (1) 的拉回通过这个商映射在核、顶层均为恒等，故同构。若 η 是
\(0\to A=\mathsf H_{\rm res}(3)\to N\to C\otimes\mathbb Q(1)\to0\)
的类，则
\[
 \boxed{(c^*\eta)(-1)=r_{\mathrm N,res}(Z_c).} \tag{8}
\]
对物理同调对象使用已经固定的相对对偶同构，其 A 和 N 的识别与上述限制、连接一致。式 (8) 在 P 留数框架中书写；若改用一动机同调框架 −P，两侧对应方向同时取负，非零性判据不变。

对 \(m=(1,2,4,8)\)，先前具体函数构造的 c=3P₁+2P₂+P₃ 因而确实计算活动下层扩张，差别仅是选择活动单位生成元时明确的非零有理倍数。没有把上层 Kummer 的 log 2 本身等同于这个 regulator。

## 五、四个单位产生六个已分裂方向

U 上有四个全局单位
\[
 x_1=a/d,\quad x_2=b/d,\quad x_3=c/d,\quad h=M/d.
\]
h 的可逆性来自 \(MN=s\ne0\)。每个二次单位杯积给出实际 Nori 态射
\[
 \mathbb Q(0)\longrightarrow H^2_{\mathrm N}(U)(2).
\]
前三个坐标之间的杯积留数是 Z；另外三个 \(\{x_j,h\}\) 的留数矩阵记为 W。若 w_i=ord_(D_i)(M/d)，则
\[
 W_{(a,b),j}=v_a^j w_b-v_b^j w_a.
\]
有限计算证书证明
\[
 \operatorname{rank}(Z\mid W)=6.
\]
因此 (1) 沿这六维子空间全部分裂。这是六个实际杯积类给出的截面，不只是留数矩阵中出现六个线性无关向量。

去掉环境三列 Z 后，定义 \(\bar W=(I-ZR)W\)。在 P 框架中的坐标为
\[
 C_{\rm cup}=
 \begin{pmatrix}
 0&-1&-1\\
 1&0&1\\
 -1&-1&0\\
 -1&1&0\\
 0&1&1\\
 2&0&0\\
 1&0&1\\
 1&1&0
 \end{pmatrix}. \tag{9}
\]
对应的三个杯积先减去 \(RW\) 指定的环境杯积，再扭转 (1)，便给出实际子对象
\[
 J\simeq\mathbb Q(1)^3\subset N\subset E.
\]
其边界像为 \(B_{\rm cup}=\operatorname{im}C_{\rm cup}\)。
所以整个八维环境下层扩张映射满足
\[
 B_{\rm cup}\subset\ker L_{\eta,\rm env},\qquad
 \dim\operatorname{im}L_{\eta,\rm env}\le5. \tag{10}
\]
不声称这些是全部核；某些参数可能有额外分裂。

## 六、三个 Tate 方向可从环境对象直和剥离

质量旗标的支撑是 P 坐标前三维
\(B_{\rm mass}=\langle e_1,e_2,e_3\rangle\)。
精确核对给出
\[
 B_{\rm cup}\cap B_{\rm mass}=0.
\]
选择显式有理基
\[
 K=(e_1,e_2,e_3,e_4,e_5,C_{{\rm cup},1},
 C_{{\rm cup},2},C_{{\rm cup},3}),\qquad\det K=-2. \tag{11}
\]
令 \(B'=\langle e_1,\dots,e_5\rangle\)，则 \(B=B'\oplus B_{\rm cup}\)。
这里为有理分裂；行列式 −2 不给出整格直和分解。

取 N' 为 N 在 B' 上的拉回。J 的实际截面给出
\[
 N=N'\oplus J.
\]
将 \(0\to N\to E\to\mathbb Q(0)\to0\) 沿 N→J 推出。所得 Ext¹(Q(0),J) 类，就是实际边界 Kummer 类投影到 \(B_{\rm cup}\) 的结果；因全部质量列位于 B'，该类为零。其分裂唯一，因为 \(\operatorname{Hom}(\mathbb Q(0),J)=0\)。故 E→J 有一个限制在 J 上为恒等的实际退缩，得到
\[
 \boxed{E\simeq J\oplus E_{\rm red},\qquad
 J\simeq\mathbb Q(1)^3,\quad \operatorname{rank}E_{\rm red}=12.} \tag{12}
\]
E_red 的三个分级秩为 6、5、1。它是由所选补空间 B' 确定的约化对象，不声称该直和分裂不依赖补空间。

所有活动质量方向都在 B_mass⊂B'。退缩在 E_act 的下两层为零，故在 E_act 上因子化为 Q(0)→J；后者为零。因此
\[
 E_{\rm act}\subset E_{\rm red}.
\]
也可从实共轭验证物理链在 J 上投影为零；物理形式在 J 上的限制由 F³ 消失为零。此分解不再降低先前的秩八活动实现，但排除了环境三列作为独立 K3 扩张参数的可能。

## 七、coaction 约束与非零性的独立证据

选取与式 (12) 相容的 Betti、de Rham 框架，环境的比较矩阵和动机群表示均分块为三个纯 Tate 因子与 E_red。动机群矩阵的形状为
\[
 \begin{pmatrix}
 \chi_{\mathbb Q(1)}I_3&0\\
 0&G_{\rm red}
 \end{pmatrix}.
\]
因此不能把任意含八个自由 K3–Tate 列的块三角群直接识别为这个环境对象的实际群。余下五列整体独立性仍需证据；其中固定活动方向的非零性已由LOWER_REGULATOR补齐。

对目标 s=−64，式 (8) 和 MINIMALITY 的判据给出
\[
 \operatorname{rank}M(I^{\mathfrak m})=8
 \iff r_{\mathrm N,res}(Z_c)\ne0,\qquad c=3P_1+2P_2+P_3.
\]
Hodge regulator 非零是足够条件；Hodge regulator 为零不自动推出 Nori 扩张分裂。计算一个非零复代表也不够，必须排除其落入相应有理周期子群。

杯积证书明确证明 \(c\notin B_{\rm cup}\)，但这只排除了已知的三个分裂方向，不能由此推断 c 的 regulator 非零。DE_RHAM_MINIMALITY 已由显式质量提升和局部 Torelli 推导很一般参数非零；指定代数基点现由LOWER_REGULATOR的实际混合周期和实分裂障碍完成独立认证。

有限核对见 [regulator_cup_splitting.json](data/regulator_cup_splitting.json) 和 [verify_regulator_cups.py](code/verify_regulator_cups.py)。其责任是全部单位 valuation、六杯积秩、式 (9)–(11) 及质量支撑交集；式 (3)、(8)、(12) 的动机证明由本稿承担。
