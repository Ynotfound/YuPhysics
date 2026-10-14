# 实际割态射、留数过滤与二维 banana 最大割

本稿继续纯理论推导。目标是给一般图分层方法补上可证明的几何接口：先有带边界的几何对及实际留数态射，再讨论割周期的可见对象。文中不把图的字符、参数单值和固定纤维动机自同态混为一谈，没有新增程序、周期运输或验证计算。

## 1. 前批审阅与本稿的范围

DEGENERATE_EXTENSIONS中完整NS投影的下降使用有限扩域上的除子类、逆交叉矩阵和有理Chow群的限制单射。K3的Pic⁰=0使Galois对NS的作用确实作用于除子类；这不是从某个Hodge投影凭空得到代数对应。等质量顶层非零证明则先由Cartan式排除开集平坦提升，再对可数个有理提升取并。两个步骤均保留原范围，特别是后者不能决定任意指定代数参数。

下面采用特征零子域k⊂C及有理Nori动机，Hodge结论由实现得到。使用已核对的六运算和纯性；不假定Voevodsky动机存在未证明的动机t结构。[2] 为避免混淆，本稿中M是上同调对象，E=M∨是其同调对偶。

## 2. 一个光滑除子的实际留数

设X/k光滑，D⊂X为光滑除子，B⊂X为光滑闭边界或简单正常交叉边界；D与B的所有交叉层横截，且D不是B的分量。置
\[
 V=X\setminus D,\quad B_V=B\cap V,\quad B_D=B\cap D.
\]
边界为空亦允许。记相对复形为
\(R\Gamma(X,B)=\operatorname{Fib}(R\Gamma X\to R\Gamma B)\)。

**命题2.1。** 有实际Nori态射及正合段
\[
 H^n(X,B)\longrightarrow H^n(V,B_V)
 \xrightarrow{\ r_D\ }H^{n-1}(D,B_D)(-1)
 \longrightarrow H^{n+1}(X,B).                         \tag{2.1}
\]
特别地，留数像是原相对对象的商，也是目标割对象的子对象，但留数不必满射或单射。

**证明。** 对D⊂X使用局部化三角和余维一纯性
\[
 R\Gamma(D)(-1)[-2]\to R\Gamma(X)\to R\Gamma(V)\to.
\]
对B_D⊂B写相同三角，利用横截基变换连接二者。取相对纤维得到
\[
 R\Gamma(D,B_D)(-1)[-2]\to R\Gamma(X,B)
 \to R\Gamma(V,B_V)\to R\Gamma(D,B_D)(-1)[-1].         \tag{2.2}
\]
若B有多个正常交叉分量，在其有限交叉的Čech复形上逐项进行，再取总复形。横截性保证每一项纯性相容。取第n上同调得到(2.1)。所有箭头先在六运算范畴构造，再取Nori上同调，并非由Betti矩阵猜测。证毕。

在局部坐标z=0定义D时，闭对数形式可写为
\[
 \omega=\frac{dz}{z}\wedge\alpha+\beta,
 \qquad\operatorname{res}_D\omega=\alpha|_D.           \tag{2.3}
\]
相对形式及其边界分量使用同一映射锥约定。式(2.3)是未除以2πi的裸代数留数；(2.1)中的Tate扭转承担比较归一化。取正向小圆为法向取向，拓扑管映射满足
\[
 \int_{\operatorname{Tub}_D\gamma}\omega
       =2\pi i\int_\gamma\operatorname{res}_D\omega.   \tag{2.4}
\]
相对情形要求管的边界落在B_V中；横截管邻域可实现这一点。式(2.4)也由每个小圆上的Cauchy积分和相对Stokes公式直接给出。裸留数公式与含扭转的实际态射必须一同使用，不能省略Tate线后仍称比较矩阵不变。

## 3. 多重留数与真正的过滤

设D=⋃_{i=1}^m D_i为带标号的简单正常交叉除子，B与其所有交叉层横截。对有序集合C=(i₁,…,i_r)定义
\[
 Z_C=(D_{i_1}\cap\cdots\cap D_{i_r})
                         \setminus\bigcup_{j\notin C}D_j,
 \quad B_C=B\cap Z_C,
 \quad M=H^n(X\setminus D,B\setminus D).
\]
逐次应用(2.2)，得到
\[
 r_C:M\longrightarrow N_C:=H^{n-r}(Z_C,B_C)(-r).       \tag{3.1}
\]
每一步只暂时填回下一条待取留数的除子，未割分量仍从环境中删去。这个规则避免把交叉点误当作留数形式必须在其上正则的点。

局部按dz_{i₁}/z_{i₁}∧…∧dz_{i_r}/z_{i_r}的顺序取系数。交换两条除子使留数及有序法向管取向同时变号；固定一次顺序后
\[
 \int_{\operatorname{Tub}_C\gamma}\omega
       =(2\pi i)^r\int_\gamma\operatorname{res}_C\omega.\tag{3.2}
\]
这是经典迭代留数公式的相对版本，基本绝对公式见[1]的Theorems17–18。本稿不宣称留数定理的原创性。

令所有大小为r的割同时组成映射r^{(r)}=⊕r_C，并约定r^{(0)}=id_M。定义
\[
 K_r=\ker r^{(r)},\qquad E=M^\vee,\qquad
 T_r=\operatorname{im}\bigl((r^{(r)})^\vee\bigr).
                                                               \tag{3.3}
\]
空交叉或负上同调次数的目标记零。对r>min(n,dimX,m)，取K_r=M、T_r=0。

**定理3.1（固定模型的留数过滤）。**
\[
 0=K_0\subseteq K_1\subseteq K_2\subseteq\cdots\subseteq M,
 \qquad E=T_0\supseteq T_1\supseteq T_2\supseteq\cdots\supseteq0,
 \qquad T_r=K_r^\perp\simeq(M/K_r)^\vee.              \tag{3.4}
\]
这些都是实际子对象。

**证明。** 任意(r+1)重留数都分解为某个r重留数后再取一次留数。因此所有r重留数为零就迫使所有(r+1)重留数为零，得到K_r⊂K_{r+1}。中性Tannakian对偶正合，映射对偶的像为原核的湮灭子，给其余结论。顺序变号不改变核或像。证毕。

这给出一个确实具有包含关系的过滤，但其输入是固定的带标号正常交叉模型，而非笼统的“组合、几何、物理层”。目前没有证明它等于权重过滤、coradical过滤，或在加入任意解消例外分量后保持不变。它对保持标号及横截条件的拉回具有留数相容性；一般吹起的比较须另给几何图。

单除子时，(2.1)还明确
\[
 K_1=\operatorname{im}\{H^n(X,B)\to M\}.              \tag{3.5}
\]
所以所有管周期都看不见这部分延伸过D的类。仅凭割数据不能无条件恢复原周期。

## 4. 割周期的最小对象：等号需要哪些标记

固定C，令F_C=N_C∨，τ_C=r_C∨:F_C→E。取有理Betti链γ∈(F_C)_B和原读出ω∈E_dR∨。记ω_C=τ_C^*ω，其中Tate线已经包含在F_C内。则形式周期具有实际态射关系
\[
 [F_C,\gamma,\omega_C]^{\mathfrak m}
       =[E,\tau_C\gamma,\omega]^{\mathfrak m}.         \tag{4.1}
\]
使用裸积分书写时，(4.1)就是(3.2)，不能把右侧2πi因子漏掉。

**定理4.1。** 采用VISIBILITY_AND_CORADICAL的规范最小对象，有
\[
 \operatorname{Vis}(F_C;\gamma,\omega_C)
     \simeq\operatorname{Vis}(E;\tau_C\gamma,\omega).  \tag{4.2}
\]
若原物理链σ满足R_{τ_Cγ}⊂R_σ，则割的最小同调对象自然单射到原物理最小对象。尤其σ生成E时成立。

**证明。** τ_C把γ生成的最小子对象满射到τ_Cγ生成的最小子对象。其核被拉回读出消去；任何另一个被消去的子对象，像仍被ω消去。反向拉回目标最大不可见子对象，得到两者的最小商同构。若R_{τ_Cγ}⊂R_σ，两者使用同一个E内的最大不可见子对象I_ω，故
\[
 R_{τ_Cγ}/(R_{τ_Cγ}\cap I_\omega)
       \hookrightarrow R_\sigma/(R_\sigma\cap I_\omega).
\]
证毕。

若同时允许所有r重割的全部有理链，生成对象恰T_r：有限维Betti基张成各F_C，像之和正是(3.3)。因此在固定读出下
\[
 \mathcal V_r=T_r/(T_r\cap I_\omega),\qquad
 \mathcal V_{r+1}\hookrightarrow\mathcal V_r.          \tag{4.3}
\]
这是经过物理读出的实际割过滤。如果σ不是循环生成元，不得直接把V₀或全部V_r当作σ的最小对象；此时只有落在R_σ内的部分受其控制。

式(4.1)与规范coaction相容，因为矩阵系数关系本来就按实际态射取商。它不推出“所有图割恰好构成coaction第二因子的一个基”：后者还需完备性、线性独立性、对偶配对及收敛/重整化边界的额外证明。

## 5. 单值变分与留数态射必须分开

在一个已构造的解析族中，固定基点处T_ℓ是有理Betti单值。若局部拓扑已经给出
\[
 (T_\ell-1)\sigma=\sum_C a_C\operatorname{Tub}_C\gamma_C,
 \qquad a_C\in\mathbb Q,                              \tag{5.1}
\]
则积分线性性和(3.2)直接给
\[
 \operatorname{Disc}_\ell I
     =\sum_C a_C(2\pi i)^{|C|}
                      \int_{\gamma_C}\operatorname{res}_C\omega.
                                                               \tag{5.2}
\]
这里a_C包含实际交数与路径、方向信息，不默认为1。普通简单夹点时，这种表示由局部消失同调提供；但正能量δ_+表示还要求识别实际实支、正能量选择以及原链的夹挤。几何复杂性或Landau方程本身不提供这些标记。无穷远贡献、重整化边界及非孤立夹点也不能被(5.2)自动排除。

引用范围尤其需要准确：[1] Theorem41在该文的紧化约定下明确假定奇数D，并保留交数N；不能把该定理直接当作本文二维多圈情形的无条件Cutkosky证明。

**命题5.1。** 固定纤维的单值矩阵不一般是该纤维的Hodge或Nori自同态。

**证明及本项目实例。** 指定banana纤维的纯T满足End_Hdg(T)=Q，已由FIXED_PURE_GROUP证明。一个非平凡Picard–Lefschetz节点反射有一维−1空间和五维+1空间，因此不是标量，不可能属于End_Hdg(T)。若它来自Nori自同态，其Hodge实现也应是自同态，矛盾。证毕。

所以(5.1)是Betti向量等式；其两侧都可作为形式周期的链标记，但不能将T_ℓ−1未经证明升级为固定动机中的一个态射或幂等投影。近邻/消失循环的变分态射处在族或退化的范畴中，与固定纤维的断言有不同的源和靶。

## 6. 二维banana的最大割是一个明确代数簇

取h条内线的banana图，圈数h−1，二维Minkowski形式为q²=(q⁰)²−(q¹)²。在复数域取光锥坐标
\[
 x_i=q_i^0+q_i^1,\quad y_i=q_i^0-q_i^1,
 \quad x_iy_i=\mu_i\ne0,
 \quad\sum_i x_i=p^+,\quad\sum_i y_i=p^-.
\]
假设s=p^+p^-≠0。全部传播子在壳的复割簇为
\[
 Z_h=\{(x_i)\in\mathbb G_m^h:
             \sum_i x_i=p^+,\ \sum_i\mu_i/x_i=p^-\}.   \tag{6.1}
\]
其光滑点的复维数是h−2。

**命题6.1。** Z_h与下列射影环面超曲面同构：
\[
 \Big(\sum_i x_i\Big)\Big(\sum_i\mu_i/x_i\Big)=s,
 \quad[x_1:\cdots:x_h]\in\mathbb P^{h-1},\quad\prod x_i\ne0.
                                                               \tag{6.2}
\]
经x_i=μ_iα_i，这正是二维banana的第一、第二Symanzik组合
\[
 F_h=\Big(\sum_i\mu_i\alpha_i\Big)
               \Big(\sum_i\prod_{j\ne i}\alpha_j\Big)
                      -s\prod_i\alpha_i=0.             \tag{6.3}
\]

**证明。** (6.1)显然给(6.2)。反向，对(6.2)的任一射影点，s≠0保证Σx_i≠0；唯一缩放p^+/(Σx_i)把其和归一为p^+，第二个和随即为s/p^+=p^-。这给互逆的代数映射。代入x_i=μ_iα_i并乘Πx_i，公共因子Πμ_i消去即得(6.3)。证毕。

因此，h=2给一般两个点；h=3给一般光滑平面三次曲线，亏格由伴随公式为一；h=4给本项目的四次K3曲面模型，解消后为X_s。h=3的光滑性可直接检查：边界顶点的一个偏导为非零μ_i，边界非顶点处亦至少一个偏导非零；环面内排除临界值。s=0时三次方程退化为直线乘二次曲线，各不可约正规化分量亏格零。这已经说明，同一边表在不同运动学下并无单一“割曲线亏格”。h=4甚至根本不是曲线。

在正实物理支，x_i>0等价于未来向在壳q_i，且
\[
 s=\Big(\sum_i x_i\Big)\Big(\sum_i\mu_i/x_i\Big)
          \geq\Big(\sum_i\sqrt{\mu_i}\Big)^2.          \tag{6.4}
\]
这是Cauchy–Schwarz，等号当且仅当x_i∝√μ_i。因而s小于阈值时没有正能量实最大割，但复割簇仍存在。特别是s=−64的复K3周期不能被称为该欧氏点的正相空间积分；两者间还需选择解析延拓路径。

## 7. 四传播子的最大割形式：裸归一化的精确等式

现在h=4，用前三条内线作独立动量。置u_i=q_i^+、v_i=q_i^-（i=1,2,3），并定义
\[
 d_i=u_iv_i-\mu_i\quad(i\leq3),\qquad
 d_4=(p^+-u_1-u_2-u_3)(p^--v_1-v_2-v_3)-\mu_4.
\]
取有序裸六形式
\[
 \Theta=\frac{du_1\wedge dv_1\wedge du_2\wedge dv_2
                          \wedge du_3\wedge dv_3}
                         {d_1d_2d_3d_4}.              \tag{7.1}
\]
本节始终按d₁,d₂,d₃,d₄顺序取留数。

**定理7.1。** 在光滑最大割上，经命题6.1的同构，
\[
 \boxed{\operatorname{res}_{d_4}\operatorname{res}_{d_3}
 \operatorname{res}_{d_2}\operatorname{res}_{d_1}\Theta
       =\operatorname{Res}_{F_4=0}
                     \frac{\Omega_{\mathbb P^3,\alpha}}{F_4}
       =-\omega_s^{\rm abs}.}
                                                               \tag{7.2}
\]
这里ω_s^abs严格采用ABSOLUTE_COMPARISON式(1)的dlog(α₁/α₄)∧dlog(α₂/α₄)∧dlog(α₃/α₄)约定。下面定义的标准Ω在α₄=1图上等于−dα₁∧dα₂∧dα₃，所以本条有明确负号；不能把两种裸形式的取向差默默吸收进旧周期证书。

**证明。** du_i∧dv_i=du_i∧dd_i/u_i。把(dd₁,dd₂,dd₃)移至所有du之前需要1+2+3=6次交换，故前三次留数为
\[
 \frac{du_1\wedge du_2\wedge du_3}{u_1u_2u_3H},\quad
 H=x_4\Big(p^- -\sum_{i=1}^3\mu_i/u_i\Big)-\mu_4,
 \quad x_4=p^+-u_1-u_2-u_3.                           \tag{7.3}
\]
在Σx_i=p^+图上，令f=(Σx_i)(Σμ_i/x_i)。直接代数恒等式给
\[
 H=-\frac{x_4}{p^+}(f-s),\qquad
 \Omega_x=\sum_{i=1}^4(-1)^{i-1}x_i
             dx_1\wedge\cdots\widehat{dx_i}\cdots\wedge dx_4
            =-p^+du_1\wedge du_2\wedge du_3.           \tag{7.4}
\]
故(7.3)就是Ω_x/[Πx_i(f−s)]。最后一次留数可用H或f−s作局部方程，乘可逆函数不改变所得Poincaré留数。线性替换x_i=μ_iα_i使分子Ω_x和分母Πx_i(f−s)同时乘Πμ_i。最后在α₄=1图上，Ω_α/Πα_i=−dlogα₁∧dlogα₂∧dlogα₃，得到(7.2)及其相对于旧绝对框架的负号。证毕。

实际Lorentz测度还有
\[
 \bigwedge_{i=1}^3(dq_i^0\wedge dq_i^1)
       =-\tfrac18\bigwedge_{i=1}^3(du_i\wedge dv_i).  \tag{7.5}
\]
结合(7.2)，有序Lorentz测度的四重裸留数为+ω_s^abs/8。因此本稿没有擅自吸收传播子的i、每圈2π因子或交数。若采用(7.1)的光锥裸测度，四重管公式的因子就是(2πi)⁴；其他振幅约定必须另外代入(7.5)及其测度定义。旧Morse球取向、P6/P7/P8以及物理I的符号全部保持不变。

在上同调的Tate账目上，未归一化六维动量对象的留数靶为H²(U_s)(−4)。若每圈以一个Tate周期作归一化，整体扭转(3)后靶为H²(U_s)(−1)。其纯超越部分T(X_s)(−1)由杯积识别为
\[
                  T(X_s)(3)^\vee=A_T^\vee.            \tag{7.6}
\]
这解释纯核的权重为何相容。但本节只构造最大割簇与形式的同构，没有据此宣称已经构造完整动量相对对象与参数相对对象之间的全局同构；后者包括无穷远、相对链及推前的额外工作。

## 8. 最大割能证明什么，不能恢复什么

全纯K3二形式来自H²(X_s)，又被全部代数除子配对消去。因此(7.2)的上同调读出通过实际纯超越因子T(X_s)因子化，不需把整个混合H²(U_s)当作纯对象。

若所取割循环在T的对偶中非零，T的简单性使该形式周期的最小纯对象就是T，连同(7.6)所需的Tate扭转。在指定ρ16纤维，其秩为六；等质量很一般ρ19纤维，其秩为三。二者的幺幂深度均为零。已有实际消失球的非零性可以在指定例中提供所需链输入，无需重做运输。

这和原物理周期的秩八/深度二、等质量很一般秩四/深度一并不矛盾：最大割只看到纯核，不会自动恢复质量Kummer扩张、顶层相对扩张和它们的比较常数。相同纯K3模型可承载不同扩张；留数像的精确几何识别不等于完整混合coaction已由最大割唯一重建。

还有一个不需新微分计算的结果。设单参数纯Gauss–Manin微分模P的秩为r，几何单值绝对不可约，ω≠0是代数截面。则ω的所有协变导数生成P，否则其生成的真微分子模会给单值不变真子系统。相应纯周期满足r阶标量Picard–Fuchs方程；任何非零周期芽在k(s)上的最小齐次方程也恰为r阶：若存在更低阶方程，其所有单值分支仍为该方程的解，但不可约性使这些分支张成r维，与低阶解空间维数矛盾。

应用于已认证几何闭包O(Q₆)的指定质量射线，得到非零纯最大割周期的最小齐次Picard–Fuchs阶数为六。这使用既有单值不可约性，不需要再生成六阶算子；也不把完整混合物理积分的一参数微分阶数与其动机秩直接等同。

## 9. 一般图接口的理论验收条件

G1c的几何输入至少包括：基域与运动学、整数时空维数和环动量路由、分母及次数、分子形式、割集合、实际边界/无穷远模型，以及所讨论的链。只有边表不足以定义(3.1)。维数正规化的非整数D也不直接是这里的有限维代数簇。

在这些数据下，先判断割簇维数及横截性，再选合适不变量：零维割用有限代数与Artin对象；曲线先正规化并保留开点/相对边界，再区分几何亏格与完整相对H¹；高维割需要相应上同调及其实际因子。计算亏格并不计算扩张，也不证明PF不可约。

本稿完成的是实际留数/管态射、固定模型过滤、可见子商比较及二维banana的几何入口。未交付任意图解消、PF生成或不可约性判断的通用有效算法，更未给多项式复杂度结论。组合字符支撑仍由G1a处理，它与这里的动量壳上留数不是同一个操作。

## 一手来源及证明责任

[1] M. Mühlbauer, *Cutkosky’s Theorem for Massive One-Loop Feynman Integrals — Part 1*, arXiv:2206.08402v2，§2.3.1 Theorems17–18、§2.3.3及Theorem41：
https://arxiv.org/html/2206.08402v2#S2.SS3.SSS1
本批直接核对作者HTML中的裸留数、迭代管公式、交数和奇数D限制。一般相对Nori版本由本文局部化三角推导；没有用其一圈定理替代二维多圈物理切割证明。

[2] S. Tubach, *On the Nori and Hodge realisations of Voevodsky motives*, Theorem4.5、Corollary4.8：
https://arxiv.org/html/2309.11999v4
六运算/实现相容性沿用已核对的框架。本文不使用其后依赖动机t结构存在性的等价断言。

本稿的banana光锥消元和有序形式符号均在第6–7节逐式推导。K3解消及纯核引用本项目既有GEOMETRIC_COMPARISON、K3_KERNEL；指定不可约性与固定端同态引用GEOMETRIC_MONODROMY_GROUP、FIXED_PURE_GROUP。最小对象方法见VISIBILITY_AND_CORADICAL。此处的重新推导与整合不自动构成任何优先权声明。
