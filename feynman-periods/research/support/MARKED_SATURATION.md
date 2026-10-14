# 由二阶判别式型唯一确定实际饱和半循环

固定 \(m=(1,2,4,8),s=-64\)。这里的 \(\delta_1,\ldots,\delta_6\) 按既有次序指节点 \(1,9,25,49,81\) 的球及 \(T_0\delta_1\)，并保持 Q6 证书的取向。

本稿证明
\[
 \boxed{w=\frac{\delta_1+\delta_2+\delta_3+\delta_4}{2}\in T(X,\mathbb Z),
 \qquad T=L+\mathbb Zw,\quad L=\sum_i\mathbb Z\delta_i.}       \tag{1}
\]
所以 \(\mathcal B=(w,\delta_2,\delta_3,\delta_4,\delta_5,\delta_6)\) 是完整整 Betti 基。结论由已经证明的几何输入及有限判别式型计算推出，不需要将四个旧质量矩阵认作实际闭路。

## 1. 输入与之前限制的区别

使用下列已保存证明：

1. Q6_GEOMETRIC_CERTIFICATE：六个规定整循环确在完整超越格 T 中，实际 Gram 为 Q6，行列式192。
2. LATTICE_SATURATION：目标纤维的完整整格为
   \[
   T_{\rm std}=U\oplus U(2)\oplus A_2(-2),\qquad\det T=48.
   \]
   该结论使用 Hessian 格定理及本原性在带标记光滑族中的保持，已结合目标点 Picard 数16证明。
3. 因此 \([T:L]=2\)，T 为偶格。

本轮再次查看 Dardanelli–van Geemen §1.6 原 PDF 第3页，公式确为上述整格，已在 LATTICE_SATURATION 中说明其在本纤维的适用性。仅凭秩、行列式、Smith 形不足以完成本稿；需要实际已证明的判别式二次型。

MARKED_LATTICE 关于“全射线单值保持三种补全”的结论仍然正确。它只检验保持性。本稿加入完整 T 的判别式型，排除其中两种，并同时排除十一种初始候选中的其余所有类。之前“必须新增质量闭路或边界粘合才能区分”的下一步建议因而被更强的现有输入替代。

## 2. 二阶判别式型的不变量

对非退化偶整格 K，记
\[
 A_K=K^\vee/K,\qquad
 q_K(\bar x)=(x,x)\pmod{2\mathbb Z}.
\]
考虑性质
\[
  \mathcal P(K):\quad q_K(A_K[2])\subset\mathbb Z/2\mathbb Z.       \tag{2}
\]
它只要求二阶挠元上的二次型值为整数模2；不要求整个判别群的值都为整数。特别是本例的三阶部分不在断言中。

式 (2) 在整格同构下不变。选 Gram 矩阵 G 后，每个二阶挠元唯一写为 \(\varepsilon/2\) 模 K，其中 \(\varepsilon\in\{0,1\}^{\operatorname{rk}K}\)，并满足
\[
 G\varepsilon\equiv0\pmod2,\qquad
 q_K(\varepsilon/2)=\frac{\varepsilon^{\mathsf t}G\varepsilon}{4}
                      \pmod{2\mathbb Z}.                         \tag{3}
\]
所以这是至多 \(2^{\operatorname{rk}K}\) 次整系数判定，不需要数值运算或格同构搜索。

**引理。** \(\mathcal P(T_{\rm std})\) 成立。

**证明。** U 是幺模格，没有非零判别元。U(2) 的 Gram 为
\(\left(\begin{smallmatrix}0&2\\2&0\end{smallmatrix}\right)\)，故
\(q((a,b)/2)=ab\)，为整数。A2(−2) 的 Gram 为
\(\left(\begin{smallmatrix}-4&2\\2&-4\end{smallmatrix}\right)\)，其所有二阶类由 \((c,d)/2\) 表示，二次型值
\[
 -c^2+cd-d^2
\]
也为整数。直和给出 (2)。其二阶子群共16个元素。证毕。

## 3. 穷尽所有指数二补全

由 \([T:L]=2\)，有 \(2T\subset L\)。所以任一可能的 T 都形如
\[
 K_v=L+\mathbb Z(v/2),\qquad
 0\ne v\in\{0,1\}^6.                                           \tag{4}
\]
整配对与偶性分别等价于
\[
 Q_6v\equiv0\pmod2,\qquad v^{\mathsf t}Q_6v\equiv0\pmod8.         \tag{5}
\]
第一条件因 Q6 所有条目为偶数而自动满足。六十三个非零 v 中，恰十一种满足第二条件。

对每个 v，令 \(p(v)\) 为首个非零坐标，取基矩阵 B_v 为单位阵把第 p(v) 列替换为 v/2，令 \(G_v=B_v^{\mathsf t}Q_6B_v\)。以下 \(\varepsilon=e_j\) 均是这个特定 B_v 基中的坐标；并不是原六循环基的第 j 个半球。

| v（按六个坐标连写） | 判别元见证 \(\varepsilon\) | \(G_v\varepsilon\) 模2 | \(q_{K_v}(\varepsilon/2)\) 模2 |
|---|---|---|---|
| 001100 | e6 | 0 | 3/2 |
| 010100 | e6 | 0 | 3/2 |
| 011011 | e4 | 0 | 3/2 |
| 011101 | e5 | 0 | 3/2 |
| 011110 | e6 | 0 | 3/2 |
| 100001 | e6 | 0 | 3/2 |
| 100010 | e5 | 0 | 3/2 |
| 100111 | e6 | 0 | 3/2 |
| 101000 | e3 | 0 | 3/2 |
| 110000 | e2 | 0 | 3/2 |
| 111100 | 无障碍；16个二阶值全部整数 | — | 0 或 1 |

前十行违反 \(\mathcal P(T)\)，不可能是实际 T。实际补全必然存在，故唯一剩余 \(v=(1,1,1,1,0,0)\) 就是其类。这证明式 (1)。

计算还显示，被排除的每个候选都有8个非整数二阶值；证明只需表中一个见证。所有十一候选的行列式均为48，说明单用行列式确实看不出这个区别。

## 4. 标记整格和标准基

以 \(\mathcal B\) 为整基，实际 Gram 为
\[
G=
\begin{pmatrix}
0&-1&-1&3&-1&-1\\
-1&-2&0&2&0&-4\\
-1&0&-2&2&0&-4\\
3&2&2&-2&0&8\\
-1&0&0&0&-2&-4\\
-1&-4&-4&8&-4&-2
\end{pmatrix}.                                                   \tag{6}
\]
取
\[
 R=\begin{pmatrix}
1&1&4&4&6&4\\
0&-1&1&-1&-1&-1\\
0&0&-2&-1&-2&-2\\
0&0&-1&-1&-2&-1\\
0&0&0&0&-1&1\\
0&0&-2&-1&-2&-1
\end{pmatrix}.
\]
精确计算给
\[
 \det R=1,\qquad
 R^{\mathsf t}GR=
 \begin{pmatrix}0&1\\1&0\end{pmatrix}
 \oplus\begin{pmatrix}0&2\\2&0\end{pmatrix}
 \oplus\begin{pmatrix}-4&2\\2&-4\end{pmatrix}.                       \tag{7}
\]
R 原来只是旧候选稿中的代数同构矩阵。现在 (1) 独立确定了其实际输入整基，再加式 (7) 的精确核验，故 \(\mathcal B R\) 真正是规定几何循环标记下的标准整基。

这一步证明的是整同调类的归属和基，不另外画出 w 的嵌入曲面代表，也不需要这样的代表来验证整格等式。

## 5. 完整整 Betti 基的绝对比较

令 P 为 ABSOLUTE_COMPARISON 已认证的裸留数 B 框架、六循环有理基的周期矩阵。B=B_(111100)。则
\[
 P_{\rm sat}=PB,\qquad P_{\rm std}=PBR.                            \tag{8}
\]
这两个矩阵现在都采用完整整 Betti 基。第一列分别是已明确整循环的周期，不再只是一个未经证明可整除的半和。

原 P 的条目中心为 Gaussian 整数除 \(2^{320}\)，统一矩阵无穷范数误差为已保存 \(e_P\)。B、BR 的条目在半整数中，所以式 (8) 的中心可精确储存在分母 \(2^{321}\) 上，没有新增舍入。又
\[
 \|B\|_\infty=3/2,\qquad \|BR\|_\infty=10,
\]
于是
\[
 \|P_{\rm sat}-\widehat P_{\rm sat}\|_\infty
       \le\tfrac32 e_P<2^{-84},\qquad
 \|P_{\rm std}-\widehat P_{\rm std}\|_\infty
       \le10e_P<2^{-81}.                                        \tag{9}
\]
对应的36+36个条目及精确误差存于 data/marked_saturation/periods.json。两者仍使用旧的明确有理 de Rham B 框架；这里没有声称恢复旧标量 theta 框架。

实际配对关系随换基变为
\[
 P_{\rm sat}^{\mathsf t}K_BP_{\rm sat}=8\pi^2G,\qquad
 P_{\rm std}^{\mathsf t}K_BP_{\rm std}=8\pi^2Q_{\rm std}.
\]
无需重跑 ODE、Morse 初值、闭路或整数迹。

## 6. 独立复核、前批审阅和剩余范围

生成入口 code/certify_marked_saturation.py；独立入口 code/verify_marked_saturation.py。独立验证不导入生成器或 SymPy，用整数、Fraction、直接双线性和及 Leibniz 行列式完成：

- 穷尽63个非零模二类，核对全部11个偶补全；
- 为每个补全枚举全部二阶判别元，验证10份排除见证及唯一剩余类；
- 核对 \(\det R=1\) 和完整标准 Gram；
- 直接重算两份整基周期的72个条目及误差传播；
- 拒绝错误半循环 (100001)/2，检查所有输入哈希。

certificate.json SHA c486af51793ea8e675c3b471c3e8bb278e742e4d12bd77e2abea4209b16f3f7b；
periods.json SHA 39f550204a9ca54d90e63ffde6d8bd3f8c5f775947f80be7bc443c0d2443e3d2。
几何输入的证明由 Q6_GEOMETRIC_CERTIFICATE 和 LATTICE_SATURATION 承担，代码不把无标记候选 Q6 自行升级成几何循环。

本批同时审阅 FIXED_PURE_GROUP 的用途：实乘/CM排除使用固定 p=11、辅助系数素数29和基域Q；Nori的额外分支来自有理比较行列式而非射线群；de Rham纯关系使用实际杯积。没有改变其证明或证书。初步考察复共轭及候选格的判别群作用并没有选出半循环；该探索不作为本稿认证输入。真正完成选择的是二次型，不只是判别群或群作用。

因此，指定目标的实际饱和标记与完整整基绝对纯比较现已完成。四个旧质量方向矩阵的几何真实性、旧theta框架、秩十环境比较仍未由本稿认证。P7/P8仍是原明确有理基比较，不能把纯格整基结果直接升级为整个混合整格基。一般图算法、独立物理基准和完整论文继续，不缩小总目标。

