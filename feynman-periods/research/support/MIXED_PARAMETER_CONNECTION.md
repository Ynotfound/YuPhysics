# 四参数实际混合连接：三条质量方向与秩十相对对象

本稿在既定二十曲线解消和相对复形约定下，构造保留全部三个质量 Kummer 方向的实际秩十对象，给出四个完整有理连接矩阵。新增下层耦合具有非常简短的仿射多项式表达式。几何识别由实际拉回、Cartan 公式、留数和唯一性证明承担；六百个曲率恒等式由两套算法精确核对。

这是四参数混合连接的结果。指定 \(m=(1,2,4,8),s=-64\) 的最小活动对象仍为已经认证的秩八对象，其 P8 和物理读出不变。本稿没有计算秩十绝对比较矩阵，也没有识别固定纯动机群或饱和整格标记。

## 1. 参数、形式与实际对象

令 \(x_i=\mu_i/s\)，\(D_i=x_i\partial_{x_i}\)，\(\Theta=\sum_iD_i\)。在环面图表中使用另一组坐标
\[
 \alpha=(a,b,c,1),\qquad N=\sum_i\alpha_i^{-1},\quad
 M=\sum_i x_i\alpha_i,\quad f=MN,\quad F=f-1.
\]
设 \(U_x=\{F=0\}\subset T=\mathbb G_m^3\)，
\[
 \eta=d\log a\wedge d\log b\wedge d\log c,\quad
 \Omega=\operatorname{Res}_{F=0}\frac{\eta}{F},\quad
 V_i=\nabla_{D_i}\Omega .
\]
在恢复的五参数记号中 \(\Omega=s\omega\)。下文全部形式在未扭转上同调中；最后扭转 (3) 才得到此前的物理对象。

取保证二十曲线标记及光滑解消成立的参数开集，再限制到下述导数框架非退化的开集。纯核 \(K\) 是边界正交因子，秩六；其 Hodge 维数为 \(1,5,6\)。不要求每个纤维的 Picard 数都等于十六：只有满足该条件的纤维才把 K 称为完整超越部分。

[KUMMER_BOUNDARY](KUMMER_BOUNDARY.md) 的局部化、杯积与同伦纤维构造在本参数族中同样成立：曲线、交点和旗标有标记，质量比为开集上的可逆函数。其实际一动机商在物理扭转约定下为
\[
 \mathsf K_{\rm bd}
 =[\mathbb Z\longrightarrow\mathbb G_m^8]_{\mathbb Q},
 \qquad1\longmapsto(r_1,r_2,r_3,1,1,1,1,1),
 \quad r_i=x_i/x_4 .
\]
这里一动机余字符基是此前的 \(-P\)；de Rham 留数基仍为 P，不能省略这一区别。

前三个环面坐标给出实际子一动机
\[
 \mathsf K_{\rm mass}
 =[\mathbb Z\longrightarrow\mathbb G_m^3]_{\mathbb Q},
 \qquad1\longmapsto(r_1,r_2,r_3).
\]
由 \(\mathsf E=H^3_{\rm N}(T,U_x)(3)\to\mathsf K_{\rm bd}\) 拉回它，定义
\[
 \mathsf E_{\rm mass}
   =\mathsf E\times_{\mathsf K_{\rm bd}}\mathsf K_{\rm mass},
 \qquad
 \mathsf E/\mathsf E_{\rm mass}\simeq\mathbb Q(1)^5.                 \tag{1}
\]
这是先在实际对象中取拉回，再取实现；没有从一个平坦矩阵猜测动机。可在函数域 \(\mathbb Q(x_1,\ldots,x_4)\) 采用原有有理 Nori 构造，并由同一几何态射在适当开集上给出相对实现与变分。此处不额外假设一般三角动机范畴的 motivic t-structure。

记未扭转对象 \(E_{\rm mass}=\mathsf E_{\rm mass}(-3)\)。它的权重分级秩为
\[
 \operatorname{Gr}^W_2=K,\qquad
 \operatorname{Gr}^W_4=\mathbb Q(-2)^3,\qquad
 \operatorname{Gr}^W_6=\mathbb Q(-3),
 \qquad\operatorname{rk}E_{\rm mass}=10.                           \tag{2}
\]
正质量欧氏纤维的物理链在扭转后属于 \(\mathsf E_{\rm mass}\)：沿用已有相对对偶识别，物理链共轭不变，而式 (1) 的分裂实商上共轭为负。因此商中的物理像为零。这不声称所有特殊纤维都以秩十为最小实现。

## 2. 实际 Hodge 提升与连接约束

沿用 [DE_RHAM_MINIMALITY](DE_RHAM_MINIMALITY.md) 中已经核对的三十处双留数，定义
\[
 \beta_i=-x_i\alpha_iN\,\Omega,\qquad
 \sum_{i=1}^4\beta_i=-\Omega.                                      \tag{3}
\]
前三个 \(\beta_i\) 的留数分别为 \(P_1,P_2,P_3\)，第四个为它们的负和。所有等式在 \(H^2(U_x)/\operatorname{im}H^2(T)\) 中理解，继而用连接映射识别到相对对象；环境杯积方向已按既定 \(I-ZR\) 投影处理。

相对复形及符号采用 [RELATIVE_TOP_COMPARISON](RELATIVE_TOP_COMPARISON.md)：
\[
 d(\alpha,\beta)=(d\alpha,-i^*\alpha-d\beta),\qquad e=(\eta,0),
 \qquad\nabla_{D_i}e=\delta\beta_i.                                \tag{4}
\]
\(e\) 属于 F³；式 (1) 的未扭转商只有型 (2,2)，故 e 在 \(E_{\rm mass}\) 中。再由连接稳定性，四个 \(\delta\beta_i\) 也在其中。此后略去注入符号 \(\delta\)。

令
\[
 C_{ij}=\nabla_{D_j}\beta_i.
\]
边界留数为常数，所以 \(C_{ij}\in K\)；Griffiths 横截性又给出 \(C_{ij}\in F^1K\)。由式 (4) 的平坦性以及式 (3)，
\[
 C_{ij}=C_{ji},\qquad \sum_j C_{ij}=-V_i.                         \tag{5}
\]
因此唯一地写成
\[
 C_{ij}=a_{ij}\Omega+\sum_{r=1}^4b_{ij}^{\,r}V_r.                 \tag{6}
\]
这里使用已经证明的局部 Torelli 性，\(\Omega,V_1,\ldots,V_4\) 是 F¹K 的一组基。

## 3. 环面分部积分确定对角项

令 \(A_i=\alpha_iN\)，以免与连接矩阵混淆，暂记环面向量场为
\[
 T_i=\alpha_i\partial_{\alpha_i}\quad (i\le3),\qquad
 T_4=-T_1-T_2-T_3.
\]
在 \(\alpha_4=1\) 图表，这四个向量场都保持 \(\eta\)。直接计算
\[
 T_iA_j=\delta_{ij}A_j-\alpha_j/\alpha_i,\qquad
 T_if=x_iA_i-M/\alpha_i,
\]
以及恒等式
\[
 \partial_{x_j}\!\left(\frac{-x_iA_i}{F}\right)
       +T_i\!\left(\frac{A_j}{F}\right)
       =\frac{\alpha_j/\alpha_i}{F^2}.                            \tag{7}
\]
乘以 \(\eta\) 后，第二项为
\(d\,\iota_{T_i}((A_j/F)\eta)\)，故其 de Rham 类为零。取留数与 Gauss–Manin 导数相容，得到
\[
 \nabla_{\partial_{x_j}}\beta_i
   =\operatorname{Res}\frac{(\alpha_j/\alpha_i)\eta}{F^2}.
\]
又因
\[
 \Theta\Omega=-\operatorname{Res}\frac{f\eta}{F^2}
             =-\Omega-\operatorname{Res}\frac{\eta}{F^2},
\]
故
\[
 \boxed{C_{ii}=-x_i(\Omega+V_1+V_2+V_3+V_4).}                    \tag{8}
\]
独立验证器在取留数之前逐项检查了全部十六个式 (7)，包括第四个射影方向。此证明不需要对发散的实积分作形式分部积分；它是在开代数环面的有理 de Rham 复形中进行。

## 4. 全部混合耦合的显式公式

对 \(i\ne j\)，令 \(\{k,l\}=\{1,2,3,4\}\setminus\{i,j\}\)。则
\[
 \boxed{\begin{aligned}
 a_{ij}&=\frac{2x_i+2x_j-x_k-x_l}{6},\\
 b_{ij}^{\,i}&=\frac{x_i-1}{3},&
 b_{ij}^{\,j}&=\frac{x_j-1}{3},\\
 b_{ij}^{\,k}&=\frac{3x_i+3x_j-x_k-3x_l+1}{6},&
 b_{ij}^{\,l}&=\frac{3x_i+3x_j-3x_k-x_l+1}{6}.
 \end{aligned}}                                                  \tag{9}
\]
交换 k 与 l 只交换相应的系数。结合式 (8)，这给出所有 \(C_{ij}\)。系数均为仿射多项式，最后一个纯框架方向的系数为零。

下面证明：这些公式确实是实际质量提升的连接，而不只是某个平坦扩张的系数。

### 4.1 纯配对与对称三张量

令 \(Y_{rl}=\langle V_r,V_l\rangle_{\rm pub}\)。[FOUR_PARAMETER_CONNECTION](FOUR_PARAMETER_CONNECTION.md) 的四阶紧凑矩阵记为 \(C^{\rm met}\)，则
\[
 Y^{-1}=C^{\rm met}/3,\qquad Y=3(C^{\rm met})^{-1}.                \tag{10}
\]
这已经包含从普通导数到对数导数的换基；不能直接把普通导数的逆配对代入。实际几何杯积比公开配对多一个整体 \(-2\)，它在以下配对与逆配对的复合中抵消。

设
\[
 T_{ijl}=\langle C_{ij},V_l\rangle=\sum_r b_{ij}^{\,r}Y_{rl},
 \qquad y_l=\sum_rY_{rl}.
\]
由式 (5)，T 在前两个指标对称。对
\(\nabla_{D_l}C_{ij}=\nabla_{D_j}C_{il}\) 配对 \(\Omega\)，利用
\(\langle F^1K,F^2K\rangle=0\)，得到
\[
 -T_{ijl}+T_{ilj}=0.
\]
故 T 为完全对称的三张量。式 (8) 和式 (5) 分别给
\[
 T_{iil}=-x_i y_l,\qquad \sum_jT_{ijl}=-Y_{il}.                   \tag{11}
\]

只有 \(t_{123},t_{124},t_{134},t_{234}\) 四个全异指标项尚未确定。定义
\[
 R_{il}=-Y_{il}+x_i y_l+x_l y_i.
\]
对四个指标对 (12),(13),(14),(23)，式 (11) 给出的系数矩阵为
\[
 \begin{pmatrix}
 1&1&0&0\\
 1&0&1&0\\
 0&1&1&0\\
 1&0&0&1
 \end{pmatrix},\qquad \det=-2.                                  \tag{12}
\]
从而
\[
 \begin{aligned}
 t_{123}&=(R_{12}+R_{13}-R_{14})/2,\\
 t_{124}&=(R_{12}-R_{13}+R_{14})/2,\\
 t_{134}&=(-R_{12}+R_{13}+R_{14})/2,\\
 t_{234}&=R_{23}-t_{123}.
 \end{aligned}
\]
再右乘 \(Y^{-1}\) 唯一恢复全部 \(b_{ij}^{\,r}\)。代入式 (10) 精确化简，所得恰为式 (8)–(9) 的 b 系数；剩下两组指标对方程也成立。

### 4.2 全纯分量的唯一性

纯框架取
\[
 L=(\Omega,V_1,V_2,V_3,V_4,\partial_{x_1}^2\Omega)^{\mathsf t},
 \qquad \nabla_{D_j}L=\mathcal A_jL .
\]
第六项是普通二阶导数，不能混同旧极化框架的各向同性补足。用从 0 开始的矩阵下标，对任意 \(l\ne j\)，在
\(\nabla_{D_l}C_{ij}-\nabla_{D_j}C_{il}=0\) 中取 \(V_l\) 系数，得到
\[
 a_{ij}=-D_l b_{ij}^{\,l}+D_j b_{il}^{\,l}
       -\sum_r b_{ij}^{\,r}(\mathcal A_l)_{r,l}
       +\sum_r b_{il}^{\,r}(\mathcal A_j)_{r,l},\quad r=1,\ldots,4.
                                                                    \tag{13}
\]
程序中的纯 V 行列为 1,…,4，与此公式一致。右侧没有未知 a；精确化简给出式 (8)–(9) 的 a。

等价地，若两个实际候选已具有相同 b，其差只能为 \(h_{ij}\Omega\)。同一平坦性恒等式的 \(V_l\) 系数在 \(l\ne j\) 时即为 \(h_{ij}\)，所以所有差都为零。至此实际质量提升的系数已被几何约束唯一确定。程序的六百个平坦性检验又确认写出的完整公式满足全部方程，避免只解一个欠定子系统。

## 5. 四个秩十连接矩阵及 Hodge 过滤

前批普通导数纯框架为
\[
 W=(\Omega,\partial_1\Omega,\ldots,\partial_4\Omega,\partial_1^2\Omega)^{\mathsf t},
 \quad \nabla_{\partial_j}W=A_jW.
\]
取 \(S=\operatorname{diag}(1,x_1,x_2,x_3,x_4,1)\)，则
\[
 \mathcal A_j=(D_jS)S^{-1}+x_jSA_jS^{-1}.                         \tag{14}
\]
完整混合框架为
\[
 \mathcal B=(L,\beta_1,\beta_2,\beta_3,e)^{\mathsf t}.
\]
其连接 \(\nabla_{D_j}\mathcal B=\mathbb A_j\mathcal B\) 由下列三行规则完全指定：
\[
 \nabla_{D_j}L=\mathcal A_jL,\qquad
 \nabla_{D_j}\beta_i=C_{ij}\ (i=1,2,3),\qquad
 \nabla_{D_j}e=
 \begin{cases}
 \beta_j,&j=1,2,3,\\
 -\Omega-\beta_1-\beta_2-\beta_3,&j=4.
 \end{cases}                                                     \tag{15}
\]
所有未写出的块为零。因此只需已有纯6×6块和式 (9) 即可重建四个10×10矩阵。

在当前框架开集上，
\[
 \begin{aligned}
 F^3E_{\rm mass}&=\langle e\rangle,\\
 F^2E_{\rm mass}&=\langle\Omega,\beta_1,\beta_2,\beta_3,e\rangle,\\
 F^1E_{\rm mass}&=\langle\Omega,V_1,V_2,V_3,V_4,\beta_1,\beta_2,\beta_3,e\rangle,\\
 F^0E_{\rm mass}&=E_{\rm mass}.
 \end{aligned}                                                     \tag{16}
\]
其维数为 \(1,5,9,10\)。扭转 (3) 后过滤整体平移，不能把未扭转的 F³写入物理扭转对象。

顶层连接还有直观写法
\[
 \nabla e=-\Omega\,d\log x_4
       +\sum_{i=1}^3\beta_i\,d\log(x_i/x_4).                      \tag{17}
\]
模去纯核后，这正是三个实际 Kummer 单位的连接。沿用
\(\nu=e_F-e_B\) 和 P 留数基，比较差为 \(+\log r_i\)；式 (17) 没有改变此前 P 与 −P 的校准。对三个坐标除子 \(x_1=0,x_2=0,x_3=0\) 取赋值，三个 \(r_i\) 的赋值矩阵为单位阵，证明它们在函数域模常数乘法独立。

## 6. 精确证书及独立复核

所有新数据在 [data/banana_mixed_parameter](data/banana_mixed_parameter)。无数值抽样承担有理恒等式证明。

- reconstruction：从纯配对的对称三张量恢复 b，再由纯连接恢复 a，写出四个混合矩阵。
- flatness：清分母后检查六对参数的全部 600 个曲率多项式系数为零。
- replay：不调用重建器，从式 (8)–(9) 的紧凑仿射公式独立生成耦合；检查 16 个取留数前的环面恒等式、96 个保存耦合条目、144 个纯换基条目、全部混合块和 600 个曲率条目。
- 独立微分按每个多项式单项式的指数手写 Euler 导数，矩阵乘法用外积累加；人工改动一个下层耦合符号/常数的负对照被拒绝。两套多项式运算共同依赖 FLINT 的精确有理算术，不能称为两个完全不同的算术库。

对列向量形式框架，程序使用
\[
 D_j\mathbb A_i-D_i\mathbb A_j
       +\mathbb A_i\mathbb A_j-\mathbb A_j\mathbb A_i=0.
\]
这是先对框架应用两次连接所得的约定；与把整个等式取负的写法相同。

| 文件 | SHA256 |
|---|---|
| Alog1.json | 93e9fe3b9159f1fe01b4bd63aa6aaa4e84c68b62dd834bd5315395b670922d21 |
| Alog2.json | 48c92a336a9b02697f1f43fb5e19466e3d48c8a31886639477d6ae7cdece436a |
| Alog3.json | 9a14e9a5dcf1f63adbe768bbab47e0bb78db5846e20e78864c4beb9b2e2ae419 |
| Alog4.json | cdce2444925b8bb5d1035861f3507bbafdb119e5cb965c2910860b8b57f1e338 |
| flatness.json | 703a7ecd14742715ceb49a6dbfc41ba3e58a0d7fca54d779a20646f7faf79327 |

完整终态和独立验证器哈希见 replay.json；final_audit.json 另固定本证明与全部旧证书哈希。reconstruction/flatness 中未完成独立复核的布尔字段记录的是生成瞬间，终态由 replay 和 final_audit 给出，不回写旧文件。

## 7. 指定质量射线与既有秩八对象

沿 \(x_i=\mu_i/s\)，有 \(\nabla_s=-\Theta/s\)。由式 (3)、(5)，
\[
 \nabla_s e=\Omega/s=\omega,\qquad
 \nabla_s\beta_i=V_i/s=\mu_i\nabla_{\mu_i}\omega.                  \tag{18}
\]
因而 \(K,\beta_c,e\)，其中 \(\beta_c=3\beta_1+2\beta_2+\beta_3\)，沿质量 \((1,2,4,8)\) 的射线构成既有秩八稳定子对象的框架。这里边界单位关系已降为单一 \(\log2\)，不同于四参数函数域的三个独立方向。

式 (18) 逐项恢复此前的纯六维、活动下层七维及顶层八维连接；既有 P6、P7、P8 的数值证书和基点几何证明保持原哈希。本稿没有重复纯运输，也没有把新的十维环境子对象误称为这个特殊纤维的最小对象。

## 8. 后续范围与文献归属

新完成的是实际质量子对象、完整混合有理连接、Hodge 框架及其精确平坦性。后续可以复用已保存纯多项式，分别积分三个 \(\beta_i\) 的行，构造秩九下层和秩十相对比较；这仍须给出三个边界 Betti 列及顶层填充的几何标记，不能仅把已有 P8 补零。

质量方向的整闭路或边界粘合还须区分三个饱和候选。固定纯 \(R_H,R_N\)、有界一般图算法、另设物理基准及连贯论文仍为完整总目标；这里未据连接平坦性识别这些对象。

纯微分系统和公开配对来自 [Duhr–Maggio–Porkert–Semper–Stawinski](https://arxiv.org/abs/2507.23061) 及 [Maggio 的 Bonn 原始补充数据](https://doi.org/10.60507/FK2/J4OPZY)，归属原作者。相对及支持态射使用前稿所列 [Tubach 的实现与六运算](https://arxiv.org/html/2309.11999v4) 和边界一动机构造。本文新增的混合仿射公式通过式 (7)–(13) 自行推导；新颖性仍待完整文献比较，不仅凭这里的独立推导作优先权声明。

