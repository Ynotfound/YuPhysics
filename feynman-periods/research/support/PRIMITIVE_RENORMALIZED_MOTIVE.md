# 原始对数图减除振幅的实际相对动机

对于正质量欧氏数据、没有子发散的对数标量图，两运动学点减除可由一个明确的代数几何对及紧致相对链实现。构造保留减除点、符号和单位归一化，并在两个点重合时仍保持两个端点的标记。由此得到实际Kummer系数扩张、其推前的核商结构、运动学微分的秩界，以及共同尺度方向的形式coaction。

本稿处理的是规定范围内的减除振幅，不把它自动认作任意方案的有限部分。含子发散的全局森林、幂次反项及IR退化继续保留独立证明责任。参数对数减除的经典背景见Bloch–Kreimer[1]；Kummer周期与函数方程见Brown[2]。下面的几何标记和比较逐步给出，不宣称上述经典理论的原创性。

## 1. 精确输入与物理归一化

沿用CONTRACTION_FOREST_GEOMETRY的U、完整含质量F以及标准射影体积形式Ω。G连通，内边数N、圈数h≥1，n=N−1，D_0为正偶整数，传播子幂次a_e为正整数。给定
\[
 a_G=D_0h/2,\qquad
 d_S:=a_S-D_0h_S/2>0\quad(\varnothing\ne S\subsetneq E_G).       \tag{1}
\]
h_S按边子图的所有连通分支计算；不只检查连通子图。(1)是本稿采用的严格UV原始性条件。N>1时它排除桥和自环：删桥后的d为负；一个自环与其余边的d相加为零，不能同时正。

两组数据κ、κ_0均有每边质量平方μ_e>0和欧氏外动量，动量守恒；其系数属于固定k⊂R。在数域/代数参数纤维上得到k上的Nori对象；在参数族上先限制到上同调局部常值的光滑开集。κ_0可以与κ相同。

取D=D_0−2ε，欧氏动量测度每圈π^{−D/2}归一化，暂不含耦合、图对称因子及Minkowski的i因子。设
\[
 c_a=\frac1{\prod_e\Gamma(a_e)}\in\mathbb Q_{>0},\qquad
 \eta=\frac{\prod_e\alpha_e^{a_e-1}\Omega_G}{U_G^{D_0/2}},\qquad
 p_G=c_a\int_\sigma\eta.                                        \tag{2}
\]
σ取使Schwinger密度为正的取向；若用Σα=1且α_1,…,α_{N−1}递增的坐标，标准Ω拉回有(−1)^{N−1}号，须同时计入链取向。本稿不改动旧banana框架。

定义裸积分及减除
\[
 I_\epsilon(\kappa)=c_a\,\mu_R^{2h\epsilon}\Gamma(h\epsilon)
 \int_\sigma\eta\,U_G^{(h+1)\epsilon}F_\kappa^{-h\epsilon},\qquad
 R_G(\kappa,\kappa_0)=\lim_{\epsilon\to0}
                   (I_\epsilon(\kappa)-I_\epsilon(\kappa_0)).      \tag{3}
\]
减除顺序明确为当前点减参考点。某些文献振幅含相反整体符号；不能把其log比值直接替换(3)。

## 2. 收敛与取消极限有统一证明

取任意Hepp扇区，边序e_1,…,e_N满足α_{e_1}≤…≤α_{e_N}，设α_{e_N}=1及
\[
 \alpha_{e_i}=t_i\cdots t_{N-1},\quad 0\le t_i\le1,\quad
 S_j=\{e_1,\ldots,e_j\}.                                        \tag{4}
\]
Kruskal逐边扩展森林给一棵同时在每个S_j上取最大森林的生成树。故每个补树单项式的t_j次数至少h_{S_j}，且有一个单项式对所有j同时达到该下界。因此
\[
 U_G=\prod_{j=1}^{N-1}t_j^{h_{S_j}}\widehat U,
             1\le\widehat U\le\tau(G)\quad\hbox{在实单位立方体上}.           \tag{5}
\]
Jacobian及幂次分子给
\[
 \eta=\pm\prod_jt_j^{d_{S_j}-1}\widehat U^{-D_0/2}\,dt_1\cdots dt_{N-1}.
                                                                    \tag{6}
\]
由(1)，每个指数非负。于是p_G有限且正。

把任意二森林延伸为生成树，其补单项式是该树的补单项式乘一条边参数。故在α最大值为1的图上，Phi/U有有限上界。质量项则给F_κ≥μ_{e_N}U，因而
\[
 F_\kappa=\prod_jt_j^{h_{S_j}}\widehat F_\kappa,\qquad
       0<c\le\widehat F_\kappa,\widehat F_{\kappa_0}\le C<\infty.             \tag{7}
\]
常数可在任意给定的紧致正质量欧氏参数区域统一选取。

这里还有一个不依赖扇区的比较：正质量保证F的单项式支撑等于U·Σα的支撑；动量二森林项不增加新的单项式类型。因此F_κ和F_{κ0}系数逐项可比较，F_{κ0}/F_κ在整个正单纯形上有正的上下界。

**定理2.1。** (3)存在，且
\[
 R_G(\kappa,\kappa_0)=c_a\int_\sigma\eta
                            \log\frac{F_{\kappa_0}}{F_\kappa},
 \qquad \operatorname{Res}_{\epsilon=0}I_\epsilon=p_G/h.          \tag{8}
\]
**证明。** (5)–(7)把正规化后各端点幂改为d_{S_j}−1+h_{S_j}ε。对足够小的|Re ε|，所有积分及有限次ε导数由可积幂乘有限对数次幂支配。记F_0=F_{κ0}，则
\[
 \frac{F_\kappa^{-h\epsilon}-F_0^{-h\epsilon}}{h\epsilon}
                    \longrightarrow\log(F_0/F_\kappa).
\]
比值的对数有界，剩余公共因子也满足上述支配。又hεΓ(hε)→1，得到第一式。未作减除时Γ(hε)的简单极点给第二式。证毕。

这比在积分号内形式取消1/ε更强：给出了极限交换所需的端点控制。任意公共解析正规化因子A(ε)、A(0)=1，都不改变R_G；它可以改变单个有限部分。

## 3. 一个全局正链避让的模型

令X为P^{N−1}沿所有非空坐标中心作的标准扇区细分模型，坐标边界D_X正常交叉。等价地使用完整坐标旗标的光滑toric模型，(4)给其正扇区图。N=1时X为点。记Y_U、Y_F、Y_{F0}为相应非边界严格变换，置
\[
 V=X\setminus(Y_U\cup Y_F\cup Y_{F0}),\quad D=D_X\cap V,
                         r=F_0/F_\kappa.                         \tag{9}
\]
(5)–(7)证明提升后的正单纯形\(\widetilde\sigma\)是V内的紧致流形带角，其边界在D。F和F_0在每个坐标除子上具有相同阶数，故其法向单项式取消；r在V上是代数单位。η由(6)延伸为V上的正则n形式，且在D每个分量上的拉回为零。这些是全局构造所需的正链与形式，不只是扇区上的渐近记号。

V还删去F和F_0，可能承载超出UV周期的几何。因此不能仅根据p_G的最小对象判断H^n(V,D)的全部结构。UV周期本身来自不依赖κ、κ_0的较大开空间X\Y_U；限制到V是实际拉回，不改变给定η、σ的形式周期关系。

## 4. 分开端点的实际相对对象

若直接使用(Gm,{1,r})，r=1时两个端点合并，普通相对H¹由秩二降为秩一。这不是保持标记的Kummer家庭，不能用于无条件的对角特殊化。

引入额外坐标w，只用来保留端点标记：
\[
 Y=V\times\mathbb G_{m,z}\times\mathbb A^1_w,\qquad
 s_0(x)=(x,1,0),\quad s_1(x)=(x,r(x),1),
\]
\[
 B_R=(D\times\mathbb G_m\times\mathbb A^1)
                          \cup s_0(V)\cup s_1(V),\qquad
                M_R=H^{n+1}(Y,B_R).                              \tag{10}
\]
两个截面总不相交，包括r=1。B_R不必纯余维一；一般相对Nori上同调容许这种闭边界，不能为套用除子留数而假装它是光滑除子。

定义
\[
 \Omega_R=c_a\,\operatorname{pr}_V^*\eta\wedge\frac{dz}{z},
\quad\Gamma_R(x,t)=(x,1+t(r(x)-1),t),\quad
 (x,t)\in\widetilde\sigma\times[0,1].                            \tag{11}
\]
链取“基底在前，t在后”的取向。r>0保证z>0，所以Γ_R紧致且避开z=0；其端点落在s_0、s_1，其侧边落在D的逆像。

**定理4.1（实际减除周期）。** Ω_R定义M_R的de Rham类，Γ_R定义其有理Betti对偶标架，且
\[
 R_G^{\mathfrak m}(\kappa,\kappa_0)
           =[M_R,[\Gamma_R],[\Omega_R]]^{\mathfrak m},\qquad
          \operatorname{per}R_G^{\mathfrak m}=R_G(\kappa,\kappa_0).           \tag{12}
\]
**证明。** Ω_R闭；在s_0上拉回为零，在s_1上为η∧dlog r=0，因为V维数n；在D逆像上η已为零。在边界闭并的相对de Rham复形中取其余分量为零，得到实际相对类。链条件已由(11)验证。拉回dz/z时，基底微分项与η楔积为零，剩下
\[
 \Gamma_R^*\Omega_R=c_a\eta\wedge
                  \frac{r-1}{1+t(r-1)}dt.
\]
对t积分恰为log r，按(8)得到(12)。整个构造都在代数几何对上完成；Nori/Hodge比较使用实际六运算[3]，没有由数值周期相等倒推动机态射。证毕。

式(11)也把减除振幅写成收敛的有理形式积分。额外w坐标没有新增物理自由度，其作用是使端点为分离的代数截面。

## 5. Kummer系数及推前的精确核商

沿f:Y→V取相对端点上同调。纤维(Gm×A¹,{(1,0),(r,1)})的相对上同调只在一次非零。得到秩二Kummer系数对象\(\mathcal K_r\)，有
\[
 0\longrightarrow\mathbb Q(0)\longrightarrow\mathcal K_r
                      \longrightarrow\mathbb Q(-1)\longrightarrow0.       \tag{13}
\]
这是实际相对家庭的序列；并非凭扩张类猜测一个矩阵。设e_0=[dw]、e_1=[dz/z]为纤维相对de Rham基，Betti链取γ_0为(11)的纤维段，γ_1为正向z圆，则在正实支
\[
 \left(\int_{\gamma_j}e_i\right)_{i,j=0,1}
       =\begin{pmatrix}1&0\\\log r&2\pi i\end{pmatrix},\qquad
 \nabla e_0=0,\quad \nabla e_1=d\log r\,e_0.                    \tag{14}
\]
e_0虽然在绝对复形中是恰当形式，在端点相对复形中不是，因为w在两端取值0、1。Hodge层F¹由e_1生成，W_0由e_0生成。r为常数时dlog r=0仍可能有非平凡算术Kummer扩张；连接式不能取代比较条目log r。

记H^j=H^j(V,D)。相对推前与边界映射锥相容，给
\[
                  M_R\simeq H^n(V,D;\mathcal K_r).               \tag{15}
\]
可直接用Y的乘积投影及两个截面的Čech复形验证：相对纤维复形是K_r[−1]，限制到D给同一构造，再取相对纤维。此处不依赖未验证的一般非proper基变换。

令δ_j:H^j(-1)→H^{j+1}为(13)产生的连接，则
\[
 0\longrightarrow\operatorname{coker}\delta_{n-1}
       \longrightarrow M_R\longrightarrow\ker\delta_n\longrightarrow0.
                                                                    \tag{16}
\]
式(16)也可完全绕开系数对象的记号而从实际几何对导出。令A=D×Gm×A¹、C=s_0(V)∐s_1(V)。三重对和闭并的Mayer–Vietoris给
\[
 \cdots\to H^j(Y,A)\to H^j(C,C\cap A)
             \to H^{j+1}(Y,A\cup C)\to\cdots .
\]
以z=1分裂Gm的常值因子，Kunneth给H^j(Y,A)=H^j⊕H^{j-1}(-1)，而H^j(C,C∩A)=H^j⊕H^j。两个截面的拉回把(a,b)送至(a,a+δ_{j-1}b)，其中δ由r:V→Gm的类与b杯积得到，按“基底形式在前”固定符号。在j=n、n+1处分别取余核和核，便得到(16)。所用外积、截面拉回及正合列均为实际Nori态射。

这是需要保留的真实扩张结构，而不是无条件的M_R=H^n⊗某个固定秩二对象。其Betti连接由r的Kummer拓扑类杯积给出，de Rham中对应dlog r；连接之外还保留算术扩张数据。

**推论5.1。** 设b_n=dim H^n_B，则
\[
 \dim(M_R)_B=2b_n-\operatorname{rank}\delta_{n-1}
                         -\operatorname{rank}\delta_n\le2b_n.               \tag{17}
\]
M_R的简单组成来自H^n及H^n(-1)的子商。因此H^n若已证明mixed Tate，则M_R及减除周期的最小对象也mixed Tate；反向不成立。这里只给载体的精确核商和充分判据，没有声称由边表即可求得b_n或δ_j。

在r=1时，K_1实际分裂为Q(0)⊕Q(-1)：投影到(Gm,{1})分裂Q(-1)因子，投影到(A¹,{0,1})给Q(0)因子；两截面仍不同。更严格地，Γ_R由代数子对
\[
 (V\times\{1\}\times\mathbb A^1,
 D\times\{1\}\times\mathbb A^1\ \cup\ V\times\{1\}\times\{0,1\})
                       \longrightarrow(Y,B_R)
\]
中的相对链推前而来，Ω_R沿此态射的de Rham拉回恒为零。因此形式周期关系本身给R_G^m(κ_0,κ_0)=0，而非仅由数值积分为零倒推形式周期为零。对象没有因端点合并而降秩。

## 6. 减除点变换是形式周期恒等式

对第三组正质量欧氏数据κ_1，有
\[
 R_G^{\mathfrak m}(\kappa_2,\kappa_0)
 =R_G^{\mathfrak m}(\kappa_2,\kappa_1)
                         +R_G^{\mathfrak m}(\kappa_1,\kappa_0).              \tag{18}
\]
**证明。** 在同时删去三个F零集的共同基底上，r_{20}=r_{10}r_{21}。用标记为w=0、1、2的三个截面，分别取z=1、r_{10}、r_{20}。第一段是1→r_{10}；第二段由1→r_{21}经代数映射z↦r_{10}(x)z、w↦w+1得到。该映射将η∧dlog z拉回为自身，因为η∧dlog r_{10}=0。两段的和与1→r_{20}的正路径在z>0的纤维条带中相对同调，侧边在D。三个标记点的几何对及实际变换给(18)的形式周期关系，不需要假设period映射单射。证毕。

文献[2]给普通Kummer对数的同类乘法关系；本稿把端点标记、随基底变化的r以及η的高次形式一同保留。对固定参考点κ_0、κ_1，(18)的差与可变外数据κ_2无关，正是本原对数图的有限常数减除变化。

## 7. 运动学导数及微分系统上界

在共同光滑参数开集上，固定κ_0，η来自不变的U几何。(14)的相对de Rham模型给
\[
 \nabla_v[\Omega_R]=\iota\!\left[c_a\eta\,v\log r\right],\qquad
 vR_G=-c_a\int_\sigma\eta\,\frac{vF_\kappa}{F_\kappa},          \tag{19}
\]
其中\(\iota:H^n(V,D)\to M_R\)来自(13)的Q(0)项。第一式是类的等式：对ηe_1应用系数连接，η在参数方向水平；同一相对复形中得到η(vlog r)e_0。第二式再取实际链周期；第2节的正下界给微分与积分交换。

注意(19)的右端一般依赖运动学，不能当成常值UV局部投影。它消去一层对数，但留下F几何。

**定理7.1（单参数阶数界）。** 限制到一个光滑代数曲线参数s，设b_n为该族H^n(V,D)的秩。减除周期满足一个k(s)系数齐次微分方程，阶数至多b_n+1。

**证明。** \(\operatorname{im}\iota\)为水平子对象，所有∇_s^j[Ω_R]、j≥1均在其中。设这些导数生成的微分子模秩为d≤b_n，则连同[Ω_R]的循环子模秩至多d+1；d+2个连续导数必线性相关。对水平Betti链取周期即可。也可先找阶数≤d的P使P∂_sR=0，再写P∂_s。无须先算R的完整积分。证毕。

这是有明确定义的上界，不是声称最小阶数等于b_n+1，也没有给任意图的b_n求取算法。它与(17)的完整相对对象秩上界不同；一个标架生成的微分子模可以远小于整个动机载体。

## 8. 共同尺度方向、MS比较与coaction

令κ=λκ_0表示把所有质量平方及所有外动量平方不变量同时乘λ>0，保持其比例不变。此时F_κ=λF_0，r=λ^{-1}在α基底上为常数，故实际乘积对和Kunneth给
\[
 R_G^{\mathfrak m}(\lambda\kappa_0,\kappa_0)
                       =-p_G^{\mathfrak m}\log^{\mathfrak m}\lambda,
 \qquad \frac{dR_G}{d\log\lambda}=-p_G
                         =-h\operatorname{Res}_{\epsilon=0}I_\epsilon.      \tag{20}
\]
η、σ在共同基底上的UV周期与其较大常值模型由实际限制关系相同。式(20)不是保持所有质量固定而只改外动量的等式；后一种变化通常同时改变角变量。物理β函数还需耦合幂次、图权重、算符和RG时间约定，不直接等于p_G。

采用Brown[2]的右de Rham coaction约定，L^dr是Q(-1)的角色，
\[
 \Delta\log^{\mathfrak m}\lambda
   =\log^{\mathfrak m}\lambda\otimes L^{dr}
                         +1\otimes\log^{dr}\lambda,
\]
\[
 \Delta R_G^{\mathfrak m}
  =-\Delta p_G^{\mathfrak m}
       \big(\log^{\mathfrak m}\lambda\otimes L^{dr}
                              +1\otimes\log^{dr}\lambda\big).              \tag{21}
\]
这里p_G的完整coaction必须保留。若其算术对象不是Tate，不能因为λ依赖只有一个对数，就删去这些方向。对一般κ，(12)仍有实际规范coaction，但不由(21)提供完整显式条目。

作为函数，非零p_G使尺度振幅满足\((\lambda\partial_\lambda)^2R_G=0\)，且在k(λ)上最小齐次阶数恰二：围绕λ=0的解析延拓把log λ移为log λ+2πi，两个分支张成常数与log λ两个独立解。该参数微分阶数与p_G所需的固定算术对象秩不同。这里使用的是函数单值；没有把它升级成固定纤维Nori自同态。

最后，在同一正规化约定中写I_ε(κ)=p_G/(hε)+f(κ)+O(ε)。则
\[
 R_G(\kappa,\kappa_0)=f(\kappa)-f(\kappa_0),\qquad
                 C_{MS}(G)=-p_G/(h\epsilon).                    \tag{22}
\]
这是本原图的MS极点与点减除有限振幅的精确比较。未经额外归一化，f可带Euler常数或测度常数，不能从(12)推断每个未修正MS有限部分都是普通Nori周期。共同正规化常数在差中取消。多图耦合的有限方案变换也不能仅由(22)替代。

## 9. 一个实际物理反例和结论边界

四维单位幂次一圈massive bubble满足(1)，U=α_1+α_2，p_G=1。在Σα=1的正密度表示中
\[
 R_G(Q^2,Q_0^2)=\int_0^1
 \log\frac{m_1^2x+m_2^2(1-x)+Q_0^2x(1-x)}
          {m_1^2x+m_2^2(1-x)+Q^2x(1-x)}\,dx,\qquad
 \partial_{Q^2}R_G=-\int_0^1\frac{x(1-x)}{F_Q(x)}dx<0.            \tag{23}
\]
正质量欧氏域内分母正。所有这些数据具有同一UV周期和MS简单极点，但有限减除振幅随Q²严格变化。这是实际图的解析反例，说明UV数据不能重建完整有限振幅；没有为证明它求完整积分。

本批得到的是：(1)范围内全局实际相对周期、分离端点Kummer模型、推前核商、形式减除点恒等式、单参数微分上界，以及尺度射线的显式coaction和MS极点比较。它既不把有限振幅几何压缩成UV周期，也不把primitive模型冒充所有含重叠子发散图的方案。

下一步的真实缺口是让森林各项在同一几何对中抵消边界极点，并证明所用相对标架与指定全局减除方案相同；逐项对数若尚不收敛，不能分别套用本稿(12)后简单相加。一般图有效算法、物理基准余项及全文整合继续保留。

### 一手来源

[1] S. Bloch, D. Kreimer, *Mixed Hodge Structures and Renormalization in Physics*, arXiv:0804.4399v3，§1.1、§7.3–7.5，尤其Theorem7.3，https://arxiv.org/html/0804.4399 。参数对数减除和MOM/MS有限差别的先行框架；其整体符号不覆盖本文(3)的独立定义。

[2] F. Brown, *Notes on motivic periods*, §5.3，https://arxiv.org/html/1512.06410 。Kummer比较矩阵、函数方程和右de Rham作用。本文矩阵按“形式为行、链为列”，与文中显示矩阵转置；r=1改用分离端点模型。

[3] S. Tubach, *On the Nori and Hodge realisations of Voevodsky motives*, Theorem4.5、Corollary4.8，https://arxiv.org/html/2309.11999v4 。使用实际相对上同调与六运算相容，不依赖条件性动机t结构等价。

内部依赖：CONTRACTION_FOREST_GEOMETRY的UV法向计数；RG_CRITICAL_VISIBILITY的森林数据责任及归一化RG反演。所有新结论来自上述证明，没有新增数值积分、枚举或验证程序。
