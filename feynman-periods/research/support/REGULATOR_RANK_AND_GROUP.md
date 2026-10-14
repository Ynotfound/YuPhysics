# 自同态域、regulator 秩与实际混合群：纯理论推导

本稿按用户最新要求，以范畴、表示论和 Hodge 变分推导为主，不新增数值计算或验证程序。旧计算证书只在已经明确的固定纤维结论中作为既有输入。

主要结果有二：
1. 对三个权重、纯核简单的实际扩张，幺幂根基的大小由下层扩张在纯核**自同态除环**上的秩决定；有理线性独立本身一般不够。
2. 在四参数 banana 族的很一般正质量欧氏纤维，得到秩十对象的最大27维幺幂根基，以及43维实际 Hodge/Nori 群的维数和其连通部分。这里“很一般”是排除可数个真解析条件，不是每个指定代数点。

## 1. 实际对象与假设

在有理中性 Tannakian 范畴中取实际三步过滤对象 E：
\[
 0\subset A\subset N\subset E,\qquad
 N/A=B\otimes V,\quad E/N=\mathbf1 ,
 \qquad V=\mathbb Q^\kappa.
\]
B是可逆的秩一对象；A简单，维数n>1；分级生成群R约化。令
\[
 S=A\otimes B^{-1},\quad D=\operatorname{End}(A)
   =\operatorname{End}(S).
\]
D是有限维有理除环。K3应用中D实际为数域。再假设S与B非同构；在本项目n=6与dim B=1中自动成立。

两个相邻扩张为
\[
 \eta\in\operatorname{Ext}^1(B\otimes V,A),\qquad
 \xi\in\operatorname{Ext}^1(\mathbf1,B\otimes V).
\]
它们已来自同一个E，不把 Yoneda 乘积消失当作凭空存在E的证明。选择V的有理基写
\(\eta=(\eta_1,\ldots,\eta_\kappa)\)、
\(\xi=(\xi_1,\ldots,\xi_\kappa)\)。
假定上层是活动的，即\(\xi_i\)在
\(\operatorname{Ext}^1(\mathbf1,B)\)中有理线性独立。

D通过推出来作用在下层Ext上，定义
\[
 r=\dim_D\langle\eta_1,\ldots,\eta_\kappa\rangle_D.        \tag{1}
\]
左D模的约定为 \(d\eta=d_*\eta\)。若D非交换，系数一律置于左侧。

**定理。** 若r>0，则E的幺幂根基U满足
\[
 \boxed{\dim U=nr+\kappa+n.}                             \tag{2}
\]
其相邻块是 \(X\oplus(B\otimes V)\)，其中
\(X\subset\operatorname{Hom}(B\otimes V,A)\simeq S^\kappa\)
作为R表示同构于 \(S^r\)；顶层到纯核的中心块为全部A。因此
\[
 U=U_{\max}
 \quad\Longleftrightarrow\quad r=\kappa,
 \qquad
 \dim U_{\max}=n\kappa+\kappa+n.                         \tag{3}
\]
此定理针对实际E及其生成群，不从一个允许的形式块矩阵推断群等号。

## 2. 下层为什么计D秩

记U_N为N在分级上作用恒等的群，嵌入向量群S^κ，其像记X。
在纤维函子下取权重分裂，幺幂作用的各列为 \(x_i\in A\)。

对 \(d=(d_1,\ldots,d_\kappa)\in D^\kappa\)，考虑Ext类
\(\eta_d=\sum_i(d_i)_*\eta_i\)。其对应的幺幂作用是
\[
                      (x_i)\longmapsto\sum_i d_i x_i .
\]
有
\[
 \eta_d=0
 \quad\Longleftrightarrow\quad
 \left.\sum_i d_i x_i\right|_X=0.                       \tag{4}
\]
正向来自分裂扩张的幺幂核为零。反向中，若该作用恒零，相关扩张的表示降到约化分级群；特征零表示半单，故扩张分裂。

所有对象与推出扩张均在N生成的张量子范畴内。这是充满子范畴，故“其中分裂”与原范畴中的分裂相同，不需要实现函子全忠实。

R约化，S简单，所以X是S^κ的半单子表示，必同构于S^q。限制映射
\[
 \operatorname{Hom}_R(S^\kappa,S)\longrightarrow
 \operatorname{Hom}_R(X,S)
\]
满射，因为X有R稳定补空间。左侧识别为D^κ，右侧的D维数为q。
式(4)说明其核恰是下层Ext关系的核，故q=r，dim X=nr。
这就是D秩而不是单纯Q秩出现的原因。

同理，上层\(\xi_i\)的有理独立性给出全部
\(Y=B\otimes V\)，dim Y=κ。

## 3. 联合像与全部中心

E到N及E/A的群态射是满射，幺幂根基亦满射到相应的核。
因此\(\operatorname{Lie}U\)在两个相邻块上的联合像投影满射到X与Y。

X只含简单因子S，Y只含简单因子B，且S不等于B。
R半单表示论排除两块之间的对角耦合，故联合像恰为X⊕Y。
这一步不能只用“两个投影分别非零”代替。

对严格上三角块记
\[
 n(x,y,z)=
 \begin{pmatrix}0&x&z\\0&0&y\\0&0&0\end{pmatrix}.
\]
交换子为
\[
 [n(x,y,z),n(x',y',z')]
       =n(0,0,xy'-x'y).                                  \tag{5}
\]
因r>0，X非零；Y是整个B⊗V，评价映射X⊗Y→A不为零。
其像为R子表示，而A简单，故像为全部A。
分别提升(x,0)和(0,y)便在Lie U中得到所有中心方向。

因此Lie U到相邻商的中心核维数n、相邻商维数nr+κ，得到式(2)。
特征零幺幂群由Lie代数确定。选Levi因子及与其相容的过滤分裂后，
\[
 U=\left\{\begin{pmatrix}
 I_A&x&z\\0&I_{B\otimes V}&y\\0&0&1
 \end{pmatrix}:x\in X,\ y\in B\otimes V,\ z\in A\right\}. \tag{6}
\]
式(6)按任意有理代数解释。

这里不把“顶层到纯核的中心块”与整个群中心混同。若
\(V_0=\ker[L_\eta:V\to\operatorname{Ext}^1(B,A)]\)，则
\[
 [U,U]=A,\qquad Z(U)=A\oplus(B\otimes V_0),
 \qquad\dim Z(U)=n+\dim_{\mathbb Q}V_0.
\]
确实，中心元素的x块必须与所有y对易，故x=0；其y块须被所有x消去。
由式(4)对有理系数的限制，这恰是B⊗V₀。最大情形r=κ时V₀=0，中心才恰为A。

r=0时本证明不强迫中心出现；还须分析独立的顶到底扩张。因此不把式(2)无条件延伸到r=0。

这是一种三步混合扩张的具体表示论论证。一般框架及相邻像/交换子的处理可比较
[Eskandari，Theorem A、Proposition 3.1.1](https://arxiv.org/html/2407.01379)。
本文不把其结果宣称为新的抽象Tannakian理论。

## 4. 最小周期与最大群不是同一个条件

若D=Q，下层映射
\[
 L_\eta:V\to\operatorname{Ext}^1(B,A)
\]
的秩就是r。既有MINIMALITY证明用其单射性判断物理框架的最小性。
所以当D=Q且上下层均活动时，
\[
 L_\eta\text{单射}\ \Longleftrightarrow\ U=U_{\max}.
\]
在这里两种结论由同一充分条件支持，但概念仍须区分。

若D严格大于Q，等价关系通常失效。取
\(\alpha\in D\setminus\mathbb Q\)、非零扩张η，并令
\[
                   \eta_1=\eta,\qquad \eta_2=\alpha\eta.
\]
两类在Q上线性独立：若aη+bαη=0而a+bα非零，该元素在D中可逆，迫使η=0。
它们却在D上只张成一维。因此对应两层扩张的幺幂核只含一个S，而非两个S。

这已给出“有理独立不保证最大下层幺幂核”的严格机制；不需要构造一个未经说明的三步动机来证明它。本项目原κ=1的十三维结论没有因此失效：一项非零扩张在任意除环上都具有秩一。κ>1的推广才必须追踪D。

## 5. 最大群时的框架最小性：一个无需数值的证明

设r=κ，取Betti向量σ在顶层的像非零，de Rham协向量ω在A上的限制非零。则该形式矩阵系数以E为最小对象。

先看含σ的子表示。中心平移σ产生整个A；全部Y平移再产生全部B⊗V；加上顶层非零向量，得到E。因此σ生成E。

再设J为被ω消去的子对象。若J∩A非零，A简单推出A⊂J，与ω|A非零矛盾。
若J在B⊗V上有非零像，全部Hom(B⊗V,A)作用把它送到非零A，仍矛盾。
若J在顶层有非零像，Y先产生非零中间层，回到前一种情形。
所以J=0，不存在承载同一形式矩阵系数的更小商。

这是纤维函子兼容的子对象论证，并不把Betti向量和de Rham协向量混为同一基中的向量；每个子对象在两个实现下对应。它不涉及周期求值映射是否单射，不声称复数数值的代数独立性，也不等同标量微分方程阶数。

## 6. 四参数banana族：很一般纯核没有额外自同态

采用MIXED_PARAMETER_CONNECTION的实际秩十对象，在其光滑带标记开集上讨论。纯核K秩六，极化签名(2,4)，Hodge数(1,4,1)。前稿DE_RHAM_MINIMALITY以十个固定节点、线丛截面和局部Torelli证明四维周期映射微分单射。

相应型IV周期域
\[
 \mathcal D=\{[\Omega]\in\mathbb P(K_\mathbb C):
 q(\Omega,\Omega)=0,\ q(\Omega,\overline\Omega)>0\}
\]
也是复四维。因此周期映射局部为开映射，而不只是非恒定。

在单连通标记图表固定任意非标量
\(\varphi\in\operatorname{End}_\mathbb Q(K_\mathbb Q)\)。
若φ为Hodge自同态，则必有
\(\varphi\Omega\in\mathbb C\Omega\)。
此条件定义真解析子集：否则一个周期域开集都由φ的特征线组成。
φ的特征线包含在有限个特征子空间的射影化之并中；非退化不可约二次曲面不包含于任何真线性子空间，也不能包含于它们的有限并。于是φ只能为标量，矛盾。

有理φ只有可数多个。排除这些真解析子集后，
\[
                       \operatorname{End}_{\rm Hdg}K=\mathbb Q.  \tag{7}
\]
另外，对每个非零有理v，额外(1,1)类条件q(v,Ω)=0也是一个真超平面截面。
排除其可数并后，边界正交核K等于完整超越部分T，Picard数为16。

这个证明依赖已经说明的四维局部Torelli与周期域维数；没有用指定s=−64的点数计算来冒充整个参数族的结论。它也不能直接套用任意一条特殊质量射线，因为一维周期像可能全落入某个额外自同态轨迹。

## 7. 很一般的三个下层regulator独立

对非零有理质量方向c=(c1,c2,c3)，置
\[
 \beta_c=\sum_{i=1}^3c_i\beta_i,\qquad
 \widetilde c=(c_1,c_2,c_3,0).
\]
后者不是常数四向量。沿固定质量只变化s，由已推导的Cartan恒等式
\[
 \nabla_s\beta_c=\sum_i\widetilde c_i\mu_i\nabla_{\mu_i}\omega .
                                                               \tag{8}
\]
局部Torelli给出：右侧模F²K的像，不在
\(\mathcal O\cdot\nabla_s\omega\)中。因为五参数空间中
\((\widetilde c_i\mu_i,0)\)不在共同缩放向量及s方向的张成空间内。

若此方向的下层扩张在某个参数开集有平坦有理Hodge截面q，
则留数与Hodge过滤迫使
\[
                           q=\beta_c+h\omega .
\]
沿s求导并模F²K，和式(8)矛盾。因此任意非零有理c都不能在非空参数开集给出分裂变分。

在单连通图表中，对固定c，有理平坦提升q仅有可数多个。
每个提升属于F²的条件为解析条件，且刚刚证明它不可能在开集恒成立。
对c、q及可数图表取并，得到可数个真解析子集。
在其外，三个下层扩张在Q上独立。结合式(7)，此时也就在D=Q上独立，r=3。

这不是把“无穷小障碍非零”直接当成某个指定纤维的非零。
无穷小论证先排除开集分裂，再通过可数提升论证得到很一般纤维结论。
任何给定代数点仍可能落入所排除的集合。

## 8. 上层独立与二十七维幺幂根基

上层一动机的实际单位为
\(r_i=\mu_i/\mu_4=x_i/x_4\)。
在正质量欧氏区域，任一有理关系可清分母写为
\[
                  \prod_{i=1}^3r_i^{a_i}=1,\qquad
                  (a_1,a_2,a_3)\in\mathbb Z^3\setminus\{0\}.
\]
每个关系是一个真解析（实际代数）超曲面；正实性排除了非平凡单位根歧义。
排除其可数并后，上层三个Kummer类独立，κ=3。

以上几个可数例外集不能包含欧氏实开集：它们是来自复参数开集的真解析条件，若在实开集恒为零，则全纯唯一性使其在复开集恒为零，与前面的证明矛盾。这里“很一般”可取在正质量欧氏实区域中。

于是n=6、κ=r=3，定理给出
\[
 \boxed{\dim R_uG_H(E_{\rm mass})=27,\qquad
        R_uG_H(E_{\rm mass})=U_{\max}.}                        \tag{9}
\]
其组成是18维下层、3维Kummer方向和6维交换子中心。

Hodge实现的忠实正合性给
\(G_H(E_{\rm mass})\subseteq G_N(E_{\rm mass})\)。
Nori群在分级上恒等的核又包含于同一U_max，因此它也恰为U_max。
分级对象A、B在Nori范畴中半单（A的Hodge实现简单，故A也简单），所以该核确为幺幂根基：
\[
                         R_uG_N(E_{\rm mass})=U_{\max}.        \tag{10}
\]
没有假设一般纯Nori范畴整体半单。

只要按前稿实际几何构造，在包含参数的实嵌入特征零基域上定义该纤维，以上Nori结论即成立。它不宣称所有有理参数点满足“很一般”条件。

## 9. 群维数、最小秩和实际coaction的范围

令S=T(1)、B=Q(1)，则A=S⊗B²。对式(7)的K3型结构，Zarhin定理给
\(G_H(S)=SO(q)\)。相关定理和假设见
[Charles §2.2，Theorem 13及第6页](https://www.imo.universite-paris-saclay.fr/~francois.charles/Picard.pdf)，原页已核对；这里取连通Hodge群的权重零表述。

SO(q)半单，无非平凡环面商，所以与独立的Tate角色合起来
\[
                    R_H=SO(q)\times\mathbb G_m,\qquad\dim R_H=16.
\]
实际Nori分级群满足
\[
 SO(q)\times\mathbb G_m\subseteq R_N\subseteq O(q)\times\mathbb G_m.
\]
因此
\[
 \boxed{\dim G_H(E_{\rm mass})=\dim G_N(E_{\rm mass})=43.}       \tag{11}
\]
两者幺幂根基均为式(9)；Hodge群连通，Nori群至多两个分支。
本稿未判定这个很一般纤维/基域上的行列式角色是否非平凡，故不自动写出两分支结论。指定Q纤维FIXED_PURE_GROUP的分支识别有单独证据，不能代用。

在选择权重分裂后，Hodge/Nori群分别是各自R_H/R_N在相应抛物群中的完整逆像。坐标元素为
\[
 g=\begin{pmatrix}H&X&z\\0&tI_3&y\\0&0&1\end{pmatrix},
 \quad X\in{\rm Mat}_{6\times3},\ y\in\mathbb Q^3,\ z\in\mathbb Q^6,
 \quad (H,t)\in R .
\]
所有X、y、z方向已经由实际扩张证明，故比较矩阵
\(\mathcal P=\left(\begin{smallmatrix}P_A&p&f\\0&\tau I_3&q\\0&0&1\end{smallmatrix}\right)\)
在既定右coaction约定下满足
\[
 \Delta p=P_A\otimes X+p\otimes t,\quad
 \Delta q=\tau\otimes y+q\otimes1,\quad
 \Delta f=P_A\otimes z+p\otimes y+f\otimes1 .                    \tag{12}
\]
式(12)给结构，不给未计算的P10绝对条目；p⊗y表示三个方向的收缩求和。

物理链在顶层非零，物理形式在A上非零。第五节遂给
\[
                 \boxed{\operatorname{rank}M(I^{\mathfrak m})=10}
\]
对同一很一般参数集合成立，且Hodge与实际Nori最小秩相同。
这不改变质量(1,2,4,8)、s=−64的已证最小秩八：其上层单位只有一个活动方向，不满足κ=3假设。

## 10. 证明责任与接下来的理论问题

本稿的依赖为实际相对对象及质量一动机比较、四维局部Torelli、Cartan导数恒等式，以及纯K3群的标准定理。下层D秩、联合像、中心、很一般纤维和框架最小性的逻辑分别写出；没有增加一份矩阵数值来替代这些证明。

明确的未完成内容仍包括：特殊参数的额外自同态/扩张退化分层；很一般Nori行列式角色的基域识别；不同物理操作间是否存在使“可见性”真正成为过滤的规范态射；一般图割几何/PF的有限适用域；连贯论文与文献优先权比较。代码化验证留在这些理论要求确实提出有限计算输入之后。

抽象幺幂群判据与现有混合扩张文献有重叠，不据独立推导宣称抽象新定理的优先权。该项目中的作用是补齐κ>1时的自同态域条件，并把实际banana变分与群结构的推理连成闭环。

