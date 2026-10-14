# 子图收缩的法向几何与森林留数

给定带质量、运动学及传播子幂次的图，收缩首先决定参数空间边界的法向首项，而不直接给出整个物理动机的直和分解。本稿从生成树展开推导这一首项，再区分相对边界连接与极点留数两个实际态射，推导对数森林的多尺度公式以及非对数发散所需的法向有限喷射。它为结构分类与局部RG之间提供可证明的接口；全文没有使用新的数值或程序计算。

## 1. 数据和经典结果的位置

设k⊂C为特征零域，G为连通图，内边数N、圈数h。每个顶点有外动量p_v，满足Σp_v=0；每条边有质量平方μ_e。先在欧氏正质量区域推导，再把有理恒等式延伸到其定义域。固定正整数传播子幂次a_e，记a_H=Σ_{e∈H}a_e。自环和平行边允许。外动量在收缩顶点处相加，保留边的质量不变。

使用补树约定
\[
 U_G=\sum_{T\text{生成树}}\prod_{e\notin T}\alpha_e,\qquad
 F_G=\Phi_G+U_G\sum_e\mu_e\alpha_e,\qquad
 \Phi_G=\sum_{T_1\sqcup T_2}(\sum_{v\in T_1}p_v)^2
                                  \prod_{e\notin T_1\sqcup T_2}\alpha_e.       \tag{1}
\]
二森林不计两个分支的交换重数。本文F是完整第二Symanzik多项式，对应Brown的Xi；不是该文仅含动量的Phi。

经典来源须明确：U的吹起边界分解见Bloch–Esnault–Kreimer[1, Proposition 3.5]；含质量的UV/IR分解及实际边界面关系见Brown[2, Proposition 2.2, Theorems 2.7, 6.5, 7.8]。下面给出适合本项目幂次、最小标架和RG接口的直接论证，不宣称这些经典分解的原创性。Bergbauer–Brunetti–Kreimer[3]提供位置空间嵌套集与重整化的先行框架；不把其位置空间模型直接等同本参数空间模型。

取真连通边子图γ，r=|E_γ|≥2，h_γ≥1，且G/γ至少有一条边。最初假定
\[
                         F_{G/\gamma}\not\equiv0.                  \tag{2}
\]
所有保留边质量正时自动满足(2)。无质量或退化运动学中必须独立检查。下文的“UV子图对象”会主动忘掉γ的质量及外动量，只保留U_γ；不能把它误认为γ在原物理数据下的完整动机。

## 2. 生成树给出的法向首项

在γ的边上置α_e=t x_e，其余边置α_e=y_e。记I_γ=(α_e:e∈γ)。有
\[
 U_G(tx,y)=t^{h_\gamma}\big(U_\gamma(x)U_{G/\gamma}(y)+tR_U(t,x,y)\big),
\]
\[
 F_G(tx,y)=t^{h_\gamma}\big(U_\gamma(x)F_{G/\gamma}(y)+tR_F(t,x,y)\big).       \tag{3}
\]
R_U、R_F是多项式；第二式并不包含F_γ作为UV主因子。

**证明。** 对任一生成树T，T∩γ是覆盖γ顶点的森林。设其分支数为c，则γ边变量在补树单项式中的次数为h_γ+c−1。因此最低次数为h_γ，且恰在T∩γ为γ的生成树时达到。这样的T与一对生成树(T_γ,T_{G/γ})一一对应：取逆像并并入T_γ即可，连通性及边数保证无圈。这证明第一式及余项界。

二森林中相同计数成立。达到h_γ时，γ完整位于同一分支；收缩不改变该分支总外动量，因此对应G/γ的二森林及其同一权重。故Phi的最低项为U_γ Phi_{G/γ}。质量和分成Σ_{e∉γ}μ_e y_e+tΣ_{e∈γ}μ_e x_e，后一项额外含t，得到第二式。U_γ U_{G/γ}非零；(2)保证F的指定首项也非零。证毕。

这也是I_γ-adic分级环中的等式。令L_γ={α_e=0:e∈γ}⊂P^{N−1}，吹起π:X=Bl_{L_γ}P^{N−1}→P^{N−1}。因N_{L/P}=O(1)^{⊕r}，例外除子
\[
                 E_\gamma\simeq P^{r-1}\times P^{N-r-1}.          \tag{4}
\]
设Y是V(U_G F_G)的严格变换之并。由(3)在E上去掉总变换的t幂后，得到
\[
 E_\gamma\setminus Y
  \simeq\big(P^{r-1}\setminus V(U_\gamma)\big)
       \times\big(P^{N-r-1}\setminus V(U_{G/\gamma}F_{G/\gamma})\big).
                                                                    \tag{5}
\]
这里用的是约化补空间，乘积中U_γ在总方程中出现几次不改变(5)。几何上，严格变换与E的交是相应超曲面沿L的射影法向锥；不是把t=0直接代入原U、F，后者会同时得到零。

若(2)失败，只能从(3)得ord_{Iγ}F_G≥h_γ+1，不能仍用(5)或下节的F阶数。一个明确例子是无质量一圈二点bubble，U=α_1+α_2、F=q²α_1α_2，缩放γ={1}。虽然h_γ=0，但G/γ为无尺度自环，F_{G/γ}=0；F实际一阶消失。这个例子说明退化分支必须由物理数据识别。它只用于首项失效，不满足本稿例外除子r≥2的假设。

## 3. 收缩面连接与UV留数是两个不同箭头

令\(\mathcal U=X\setminus Y\)，B为坐标超平面严格变换在\(\mathcal U\)中的并，不把E重复算入B。X光滑，E与坐标边界正常交叉；删去Y后仍成立。置E°=E∩\(\mathcal U\)，B_E=B∩E°，n=N−1。单次吹起时，(5)与边界一起是乘积对：
\[
 (E^\circ,B_E)=(V_\gamma\times V_c,
                    B_\gamma\times V_c\ \cup\ V_\gamma\times B_c),            \tag{6}
\]
其中V_γ=P^{r−1}\V(U_γ)，V_c=P^{N−r−1}\V(U_cF_c)，B_γ、B_c为各自坐标边界。定义
\[
 M_\gamma^{UV}=H^{r-1}(V_\gamma,B_\gamma),\quad
 M_c=H^{N-r-1}(V_c,B_c),\quad
 M_G^{\mathrm{face}}=H^n(\mathcal U,B\cup E^\circ).                \tag{7}
\]
若需要进一步解消，只能使用已经构造、且在E上与两个因子的解消相容的模型；不能假定任意新增吹起不改变所有标架。Brown的motic模型提供这样的系统性选择[2]。本节单次吹起的实际对已足够构造下列箭头。

**定理3.1（相对面连接）。** 外积与三重对的连接给实际Nori态射
\[
 M_\gamma^{UV}\otimes M_c
 \xrightarrow{\ \times\ }H^{n-1}(E^\circ,B_E)
 \xrightarrow{\ \partial_E\ }M_G^{\mathrm{face}}.                \tag{8}
\]
其中没有Tate扭转；不保证复合单射。本稿把连接号固定为Betti对偶提取E面的正乘积链：在正法向t≥0的领口中，通常向外边界取向与此相反，须同时调整连接的整体负号。这个选择不改变正合性。对提升后避开Y的正单纯形链σ_G，有
\[
 \langle\partial_E(\alpha\times\beta),\sigma_G\rangle
    =\langle\alpha,\sigma_\gamma\rangle
       \langle\beta,\sigma_c\rangle.                             \tag{9}
\]
如进一步边界仍碰到Y，(9)只在完成所需解消并得到实际相对链后使用。

**证明。** 三重对B⊂B∪E°⊂\(\mathcal U\)的纤维三角给
\[
 H^{n-1}(\mathcal U,B)\longrightarrow H^{n-1}(E^\circ,B_E)
 \xrightarrow{\partial_E}H^n(\mathcal U,B\cup E^\circ).            \tag{10}
\]
中间项由闭并的相对Mayer–Vietoris/切除得到。复合相对外积即(8)，不需要把全部H^{n−1}强行等同单个Kunneth分量。连接的Betti对偶按上述号约定取σ_G的正E面；正实吹起在该面记录γ参数之比及剩余参数，正好是σ_γ×σ_c。相对Stokes给(9)。这些步骤先在实际几何对与六运算中完成，再取Nori上同调；[4]保证与Hodge实现相容。证毕。

式(10)也给出精确可见性限制：外积类若从H^{n−1}(\(\mathcal U\),B)限制而来，则面连接消去它。另一方面，(9)右边非零时，连接类必非零。即使非零，它仍可能对应与原物理形式不同的de Rham类；“整个图对象含这个周期”不等于“给定振幅读出必须含它”。

**定理3.2（法向极点留数）。** 另有正合段
\[
 H^n(\mathcal U,B)\longrightarrow
 H^n(\mathcal U\setminus E^\circ,B\setminus E^\circ)
 \xrightarrow{\operatorname{Res}_E}
 H^{n-1}(E^\circ,B_E)(-1)
 \xrightarrow{\operatorname{Gys}_E}H^{n+1}(\mathcal U,B).          \tag{11}
\]
它来自余维一纯性及相对局部化，证明已在CUT_RESIDUE_VISIBILITY展开。因此留数像精确为Gysin核。未扭转的裸形式留数满足
\[
 \int_{\operatorname{Tub}_E\tau}\omega
                     =2\pi i\int_\tau\operatorname{res}_E\omega.              \tag{12}
\]
(8)从面对象进入相对对象；(11)从删去面的补空间出去并带(-1)。两者来源、方向、次数、Tate因子及物理链均不同。divergent正单纯形不能仅凭(8)当作删面补空间中的普通有限相对周期；需先给正规化和减除。

## 4. 幂次计数与对数边界形式

取偶整数D_0。去掉整体Gamma与测度常数后的标量射影形式为
\[
 \omega_G(D)=\Big(\prod_e\alpha_e^{a_e-1}\Big)
       \frac{U_G^{A_G(D)}}{F_G^{B_G(D)}}\Omega_G,
 \quad A_G=a_G-(h_G+1)D/2,\quad B_G=a_G-h_GD/2.                  \tag{13}
\]
取标准Ω_H=Σ_i(−1)^{i−1}α_i dα_1∧…∧省略dα_i∧…，子图和收缩图继承边顺序；留数按法向在前。D=D_0时是有理形式；D的非整数变化仅作解析正规化，并不称其为普通Nori对象。

置
\[
                         d_\gamma=a_\gamma-D_0h_\gamma/2.        \tag{14}
\]
在收缩图某一非零坐标设为1的图中，再用γ的一条非零方向坐标定义t，Jacobian是t^{r−1}；幂次分子额外给t^{aγ−r}，(3)的分母净贡献为t^{−D_0hγ/2}。故
\[
                     \operatorname{ord}_{E_\gamma}\pi^*\omega_G(D_0)
                                  =d_\gamma-1.                 \tag{15}
\]
等号指E的泛点，且(2)成立。a_e正整数保证单项式分子在方向泛点不消失；它不排除交叉层上的额外阶数。

**定理4.1（对数收缩留数）。** 若d_γ=0，则
\[
 \operatorname{res}_{E_\gamma}\pi^*\omega_G(D_0)
       =s_\gamma\eta_\gamma^{UV}\boxtimes\omega_{G/\gamma}(D_0),\qquad
 \eta_\gamma^{UV}=\frac{\prod_{e\in\gamma}x_e^{a_e-1}\Omega_\gamma}
                                      {U_\gamma^{D_0/2}}.        \tag{16}
\]
其中s_γ为体积形式的显式取向号。将γ的r条边移到最前、保持两块内部顺序，设该置换符号为sgn(π_γ)，则s_γ=(−1)^r sgn(π_γ)。采用相容正密度约定时可把这个号同时吸收入链取向，不能只改形式。

**证明。** 法向幂由(15)为dt/t。U_γ净指数为A_G−B_G=−D_0/2；且a_γ=D_0h_γ/2、h_G=h_γ+h_{G/γ}，给A_G=A_{G/γ}、B_G=B_{G/γ}。在γ已排最前、α_r=t、x_r=1、α_N=1的图上，直接楔积得到π*Ω_G=(−1)^r t^{r−1}dt∧Ω_γ∧Ω_c；恢复原边序再乘sgn(π_γ)。这给(16)。η_γ的总齐次次数为零，确实是γ方向射影形式。证毕。

若γ的UV方向积分收敛，定义其形式周期p_γ。它的几何对、形式与正链都不依赖原γ质量及外运动学。因此(16)给一个常值的算术对象乘上收缩图对象。这具体实现了局部性的一种来源：外运动学被忘掉，而p_γ的算术不必是Tate。在Tannakian层面，常值变分不等于Tate对象；不能把外动量导数为零当作混合Tate判据。

分子会改变判据。对齐次有理形式Ω_G P/(U_G^u F_G^v)，若P的γ法向次数为ν_γ(P)，則泛点阶数是
\[
              r-1+\nu_\gamma(P)-h_\gamma(u+v).                 \tag{17}
\]
其首项P_0一般是方向变量与收缩变量的和式，并不分解为一个乘积。因而它只给张量通道之和，且可使(16)的候选留数消失；单位标量的乘积结论不原样覆盖任意分子或Ward投影。

## 5. 极点由哪一阶法向数据决定

令D=D_0−2ε，先固定避开其余奇异面的紧致方向区域，或已经对其余面作出有效减除。则某个法向积分可写成
\[
 J(\epsilon)=\int_0^1 t^{d_\gamma-1+h_\gamma\epsilon}
                                A(t,\epsilon)\,dt,             \tag{18}
\]
A在t=0充分光滑，在ε=0解析，所需余项一致可积。这里Γ(B_G)、额外正规化因子和其他尺度的极点尚未乘入。

若d_γ=0，减加A(0,ε)立即得到
\[
 J(\epsilon)=\frac{A(0,\epsilon)}{h_\gamma\epsilon}
      +\int_0^1 t^{-1+h_\gamma\epsilon}(A(t,\epsilon)-A(0,\epsilon))dt,
 \qquad\operatorname{Res}_{\epsilon=0}J=A(0,0)/h_\gamma.          \tag{19}
\]
对数情况下，法向锥主项正好够用。但**这只是指定法向端点的裸积分极点**；反项采用C=−K时另有负号，整图的Gamma极点、森林和物理求和需按实际方案另行组合。

若d_γ=−k<0，展开A=Σ_{j=0}^k A_j(ε)t^j+t^{k+1}\widetilde A，逐项积分给
\[
 J(\epsilon)=\sum_{j=0}^k\frac{A_j(\epsilon)}{h_\gamma\epsilon+j-k}
                  +H(\epsilon),\qquad
 \operatorname{Res}_{\epsilon=0}J=A_k(0)/h_\gamma.                \tag{20}
\]
H在ε=0解析。首项A_0控制截断意义下最严重的幂次发散，但ε=0的留数由第k阶法向喷射给出。比如A=1+c t、d_γ=−1时，J=(h_γε−1)^{-1}+c/(h_γε)：同一首项可有不同的MS极点。这是直接反例，说明“只看法向锥即可给所有反项”不成立。

(20)是Schwinger尺度坐标中的公式，不直接等于外动量Taylor算符。两者在指定场论中的比较还需分子、局部算符插入和IR控制；不能把这里的k机械读成外动量微分阶数。一般正规化a_e(ε)=a_e+εr_e时，分母h_γ改为h_γ+Σ_{e∈γ}r_e；若它为零，这个ε方向根本不调节该端点。

## 6. 相容森林的多尺度因子化

采用明确而足够的相容条件：有限集合\(\mathscr F\)由真连通子图组成，任意两个或边集包含，或顶点集不交。对每个H∈\(\mathscr F\)∪{G}，令child(H)为其最大真子成员，定义
\[
                         H^\circ=H/\bigcup_{K\in child(H)}K.      \tag{21}
\]
要求每个H°至少一条边，根G°满足F_{G°}≠0。顶点不交条件保证子成员的圈数可加；可以放宽，但必须另证相应秩可加性。仅说“边不交”并不足以支持下面的证明。

在每个H°选一条标记边e_H。置
\[
 \alpha_e=u_e\prod_{\substack{\gamma\in\mathscr F\\e\in\gamma}}t_\gamma,
                  \qquad u_{e_H}=1.                             \tag{22}
\]
根标记固定射影尺度。每层标记先定义该层相对父层的尺度，故(22)在非零标记图上可逆；t_γ=0是相应迭代吹起图的边界。体积形式的Jacobian因子为∏_γ t_γ^{|E_γ|−1}，可按最外层到最内层逐次计算；嵌套子成员的所有边都计入外层指数。

**定理6.1。** 去掉法向单项式后在t_γ=0对所有γ取值，有
\[
 U_G=\Big(\prod_\gamma t_\gamma^{h_\gamma}\Big)
       \left(\prod_{H\in\mathscr F\cup\{G\}}U_{H^\circ}+R\right),
\]
\[
 F_G=\Big(\prod_\gamma t_\gamma^{h_\gamma}\Big)
       \left(F_{G^\circ}\prod_{\gamma\in\mathscr F}U_{\gamma^\circ}
                                                       +R'\right),          \tag{23}
\]
R、R'属于由全部t_γ生成的理想。若所有d_γ=0，则
\[
 \operatorname{res}_{\mathscr F}\pi^*\omega_G(D_0)
        =s_{\mathscr F}\mathop{\boxtimes}_{\gamma\in\mathscr F}\eta_{\gamma^\circ}^{UV}
                                                \boxtimes\omega_{G^\circ}.
                                                                    \tag{24}
\]
这里s_F由选定消去次序中(16)的各取向号及把法向移到统一顺序的Koszul号相乘得到。法向顺序及因子顺序同时固定；交换法向次序的反号不被省略。

**证明。** 最内层生成树可先在顶点不交的子成员中独立选取。其并是森林，收缩后选父商图生成树，再提升，逐层得到同时适应全部成员的树。反过来，每棵适应树唯一恢复这些商图树。对每个γ，非适应树的补树单项式中γ次数至少h_γ+1；所有树的该次数至少h_γ。故可同时提出∏t_γ^{h_γ}，常数项精确为各商图U之积。根用二森林并保留外运动量及剩余质量，给F公式；子成员内部质量每次进入余项。这个论证也证明首项不依赖选取哪一个包含关系的线性延伸次序。

圈数和a在收缩下可加；全部d_γ=0蕴含每个非根商图γ°仍对数。将(23)及Jacobian代入(13)，每个法向给dt_γ/t_γ，各非根U指数为−D_0/2，根指数与其自身(13)相同，证明(24)。证毕。

在多尺度A(t,ε)解析且余下积分一致可积的单一乘积图上，(24)还给
\[
 [\epsilon^{-|\mathscr F|}]\int\prod_{\gamma\in\mathscr F}
                t_\gamma^{h_\gamma\epsilon-1}A(t,\epsilon)\,dt
                       =\frac{A(0,0)}{\prod_{\gamma\in\mathscr F}h_\gamma}.
                                                                    \tag{25}
\]
若还有方向变量，应对右侧继续积分。证明是对各变量反复使用(19)；含至少一个t_γ的余项最多有少一个法向极点。式(25)不含整体Gamma、其他图、其他扇区及其重叠计数，因此不是未经拼接的整图最高极点公式。

**最高留数不足以重建RG简单极点。** 两个相容尺度就有反例：
\[
 \int_0^1\!\int_0^1
 t_1^{h_1\epsilon-1}t_2^{h_2\epsilon-1}(1+c t_1)\,dt_1dt_2
 =\frac1{h_1h_2\epsilon^2}
                  +\frac{c}{h_2\epsilon(1+h_1\epsilon)}.         \tag{25a}
\]
所有c的最深角点首项都等于1，但简单极点相差c/h_2。即使所有端点均对数，也必须保留较浅面上的数据及其余变量积分；一般还须保留其ε展开。这个有理积分反例与RG_CRITICAL_VISIBILITY中“有限项先插入后可进入极点”的代数反例相容。它不声称任意c对应实际Feynman图；它严格排除了仅凭最深森林乘积首项推出全部局部RG系数的形式论证。

## 7. 重叠子图为什么不能同时作为独立残数

令A、C为两个不可比边集，L_A、L_C为坐标中心。若A∪C≠E_G，则交为L_{A∪C}。先吹起该交后，L_A、L_C的严格变换不交：它们在交中心法向空间中的方向子空间交为零，射影化后没有公共点。后续吹起二者产生的例外除子也无共同层。若A∪C=E_G，则原射影中心已经不交。

这说明在包含该交中心的模型里，不能给这两个原标签捏造共同的二重留数面。对于相交边集，这是处理重叠的一条具体规则；对边不交的横截中心，也可以选择再吹起其交，此时模型把同时面细分为旗标面。因此**边界层由选定building set决定，不能无条件说“所有Zimmermann森林就是某一个模型的全部层”**。不同细分必须通过实际拉回、链分解和取向比较。

例：四维三线sunset的三个二边bubble均对数发散，但对应P²中的三个不同坐标点。吹起后所得三个例外曲线两两不交。因此一个局部图最多含其中一个proper-bubble留数；三项须在不同面/区域求和，不能把三bubble相乘当作共同三重极点。整体图另有幂次发散及整体尺度，故这里也没有断言整图只有简单ε极点。

在已经得到的正常交叉对数图上，局部Taylor零阶算符\(T_i A=A|_{t_i=0}\)彼此交换，
\[
 \prod_i(1-T_i)A=\sum_{S}(-1)^{|S|}\prod_{i\in S}T_i A          \tag{26}
\]
消去所有对应端点的对数不可积性。更高幂次用第5节所需有限阶Taylor算符。这个局部几何包含–排除式解释森林项的来源；要等同某个全局物理BPHZ/MS方案，仍需给全局拼接、重叠计数、外动量Taylor/局部插入比较及IR处理。[3]在其位置空间、至多对数的主范围内完成了相应方案理论，不能未经比较就移植为本项目所有参数积分的定理。

## 8. 本项目banana的检验性推论

二维h线banana，单位传播子、正质量，任取k条边组成的真连通子图γ，有h_γ=k−1、a_γ=k，因此
\[
                             d_\gamma=k-(k-1)=1.                \tag{27}
\]
每个这样的UV面上参数形式泛点阶数为零，**没有对数UV留数**。整体也有a_G−D_0h_G/2=1，无整体UV Gamma极点。四线情形(13)恰为Ω_G/F_G；其K3及已证混合周期来自完整几何与边界扩张，而不是UV极点。

所以“UV筛选为空”绝不推出原物理最小对象是Tate或零。本项目指定四线对象的纯K3核与秩八最小性是既有独立证据；(27)只用圈数和法向计数，未重算这些证据。

另一方面，若一个对数子图的p_γ已经在实际形式周期中证明含某算术方向，(16)只证明该方向可沿对应收缩通道出现。完整森林与物理求和仍可能消去它。只有它在规定MS型的完整局部RG数据中非零后，RG_CRITICAL_VISIBILITY的反演才保证它不会在所有临界系数中同时消失。这样连接三层时，每一次“可见”都有自己的实际映射与非零条件。

## 9. 已解决的接口与剩余证明责任

本稿给出：带幂次的法向阶数；正质量UV面乘积对；实际Nori面连接与留数/Gysin两套箭头；对数收缩和相容森林的形式因子化；幂次发散的有限喷射规则；重叠层不能任意相乘的几何反例；以及二维banana无UV留数但有非Tate物理周期的区分。

仍需独立完成的接口包括：无质量或特殊运动学下退化F首项与IR层；任意分子及局部算符插入；具体物理方案的全局森林拼接；减除后形式和物理链的实际Nori标架；它们与RG提取、量子场论恒等式及coaction的比较。没有从(23)宣称完整M_G等于所有商图动机的张量积，没有由局部Laurent计算构造全部重整化振幅的动机，也没有把普通生成树分解称为新分类定理。

### 一手来源与阅读范围

[1] S. Bloch, H. Esnault, D. Kreimer, *On Motives Associated to Graph Polynomials*, Proposition 3.5，https://arxiv.org/html/math/0510011 。原始HTML的法向锥及例外除子公式已核对；未把抽取的树边数不等式方向当作证明，本文直接按分支数计数。

[2] F. Brown, *Feynman amplitudes, coaction principle, and cosmic Galois group*, arXiv:1512.06409v3，§2.1–2.3、§6.4、§7.5，https://arxiv.org/html/1512.06409 。本文(3)、(5)、(8)分别对应经典分解与面关系。该文主要在实现范畴表述；本稿Nori箭头由实际对及[4]构造，未把实现范畴和Nori范畴等同。

[3] C. Bergbauer, R. Brunetti, D. Kreimer, *Renormalization and resolution of singularities*, arXiv:0908.0633v2，§4.4、§5、§6.4，https://arxiv.org/html/0908.0633 。嵌套集依building set、适应树和正常交叉Laurent结构的背景；主范围为位置空间、至多对数奇性。

[4] S. Tubach, *On the Nori and Hodge realisations of Voevodsky motives*, Theorem4.5、Corollary4.8，https://arxiv.org/html/2309.11999v4 。使用六运算相容，不使用条件性动机t结构等价。

内部比较：CUT_RESIDUE_VISIBILITY的相对纯性与管取向；RG_CRITICAL_VISIBILITY的森林自然性、有限项插入与临界反演；THEORY_MANUSCRIPT的实际最小标架。原始资料及所有旧K3/周期数值证书保留不变。
