# 固定 banana 纤维的纯 Hodge/Nori 群及完整活动群

本稿只固定
\[
 X=X_{(1,2,4,8),-64}/\mathbb Q,\qquad T=T(X)_{\mathbb Q},
 \qquad S=T(1),\quad A=T(3),\quad B=\mathbb Q(1).
\]
系数范畴为此前使用的有理 Nori 动机，基域明确为 **\(\mathbb Q\)**；群默认采用 Betti 纤维函子。已有 \(\rho(X)=16\)、实际六循环 Gram \(Q_6\)、\(\det Q_6=192\)、指数二、绝对 P6 和活动秩八对象 E 均作为已保存输入。

**结论。**
\[
 \boxed{\operatorname{End}_{\rm Hdg}T=\mathbb Q,\qquad
        G_H(S)=\operatorname{SO}(T,q),\qquad
        G_{\rm N,\mathbb Q}(S)=\operatorname O(T,q).}             \tag{1}
\]
取 \(R_H=G_H(A\oplus B)\)、\(R_N=G_{\rm N,\mathbb Q}(A\oplus B)\)，有
\[
 \boxed{R_H\simeq SO(q)\times\mathbb G_m,\qquad
        R_N\simeq O(q)\times\mathbb G_m.}                        \tag{2}
\]
在 A、B 上的作用分别为 \((g,t)\mapsto(t^2g,t)\)。

结合既有完整十三维幺幂根基，实际活动对象的 Hodge/Nori 群均成为明确矩阵群；它们维数都为29，前者连通，后者有两个几何连通分支。这是固定纤维的证明，没有使用参数射线的几何单值 O(Q6) 来替代它。

## 1. 旧证书中本次需要的信息

[PICARD_CERTIFICATE](PICARD_CERTIFICATE.md) 已证明 p=11 良约化且几何 Picard 数仍为16。其残余六维 Frobenius 多项式为
\[
 P(Z)=Z^6+10Z^5+55Z^4-484Z^3+6655Z^2+146410Z+1771561.             \tag{3}
\]
旧点数为 \(288,16568,1794960\)，残余幂迹 \((-10,-10,2102)\)。本次只从保存的点数重算 Newton 恒等式，不再枚举有限域点。

P 模5满足完整六次 Rabin 不可约性判据，故在 \(\mathbb Q\) 上不可约。归一化 \(11^{-6}P(11Z)\) 的系数不全为整数，故 P 没有11倍单位根根。另有
\[
 P(Z)\equiv Z^5(Z-1)\pmod {11},\qquad P'(1)\equiv1\pmod {11}.     \tag{4}
\]
因此 P 在 \(\overline{\mathbb Q}_{11}\) 中恰有一个单位根，此处“单位根”指 **11-adic unit root**，不是有限阶的 root of unity。下文称它为 p-adic 单位特征根，以避免混淆。

实际六循环形成有理基。它们的格虽未饱和，但
\[
 \det q\equiv\det Q_6=192\equiv3
 \quad\text{于 }\mathbb Q^\times/\mathbb Q^{\times2}.             \tag{5}
\]
指数二只改变判别式的平方因子，不影响式 (5)。

## 2. 排除实乘

令 \(F=\operatorname{End}_{\rm Hdg}T\)。K3 型的标准定理给出：F 是全实域或 CM 域，其对极化的伴随分别为恒等或复共轭，且 \([F:\mathbb Q]\mid6\)。采用 [Charles §2.2 所列 Zarhin 定理](https://www.imo.universite-paris-saclay.fr/~francois.charles/Picard.pdf)，其原页第6页已核对。

先说明 F 全实时 \(\dim_FT\ge3\)，以免遗漏三次域。

- 若 \(\dim_FT=1\)，区别嵌入对应的复特征空间只有一维，却必须同时容纳 \(H^{2,0}\) 与 \(H^{0,2}\)，矛盾。
- 若 \(\dim_FT=2\)，写极化为一个 F 上二元对称型的迹，F 上对角化为 \(\langle a,b\rangle\)。每个实嵌入处该二元型均定号：区别嵌入为正定，其余为负定。矩阵
  \[
  J=\begin{pmatrix}0&-b\\a&0\end{pmatrix},\qquad J^2=-ab
  \]
  在区别嵌入处保持两条复各向同性线，即 \(H^{2,0}\)、\(H^{0,2}\)；其余嵌入全部为 (1,1)。故 J 是有理 Hodge 自同态。又因 \(-ab\) 全负，J 不在全实标量域 F 中，矛盾。

所以非平凡全实 F 唯一可能为二次域，\(\dim_FT=3\) 为奇数。Charles Theorem 1(2) 在 p=11 给出
\[
 \rho(X_{\overline{\mathbb F}_{11}})\ge16+2=18,
\]
与已认证的16矛盾。原定理的特征至少5条件满足。此处采用的是固定良约化，而非很一般参数的特化。

亦可用其 Lemma 16 的群表述理解这个矛盾：某个归一化 Frobenius 幂属于三个 F 维度的特殊正交群，必须有特征值1；这会给 P 的11倍单位根特征值，与第一节矛盾。

因此 F 全实时只能为 \(\mathbb Q\)。

## 3. CM 候选必含 \(\mathbb Q(\sqrt{-3})\)

F 若为 CM 域，其次数只能是2或6。必须谨慎处理次数6的非 Galois 可能性，不能未经证明假定所有 CM 域都含虚二次子域。

[Elsenhans–Jahnel Theorem 7.5](https://wwwuser.gwdguser.de/~jjahnel/Arbeiten/Komponentengruppe_v4f.pdf) 证明：若 \(X\times X\) 的 Hodge 猜想成立，基域 K 与端同态域的合成 KF 在 K 上正规。同文 Theorem 7.7 的证明明确援引 CM 情形下该 Hodge 猜想已知成立；这里 CM 按其 §5.1 指端同态域为 CM 域，不要求 \(\dim_FT=1\)。原文第21、22页均已查看。

仅在当前反证的 CM 情形应用此已知结论。因为 \(X/\mathbb Q\)，F 在 \(\mathbb Q\) 上 Galois。次数6的 Galois CM 域具有中央复共轭，故其六阶 Galois 群不能是中心平凡的 \(S_3\)，只能为 \(C_6\)。其三阶子群固定一个虚二次域 \(K_0\)。次数2时直接取 \(K_0=F\)。

写 \(K_0=\mathbb Q(\sqrt{-d})\)，d 为正无平方因子整数。T 是三维 \(K_0\) 向量空间，伴随在 \(K_0\) 上为复共轭，故极化是一个 Hermitian 型的迹。Hermitian 型可在 \(K_0\) 上对角化，三个对角系数 \(h_j\) 均在 \(\mathbb Q^\times\) 中。相对于每个 \(K_0\) 基向量及其 \(\sqrt{-d}\) 倍，迹型的二阶块为
\[
 \begin{pmatrix}2h_j&0\\0&2dh_j\end{pmatrix}.
\]
于是
\[
 \det q=2^6d^3(h_1h_2h_3)^2\equiv d
          \pmod{\mathbb Q^{\times2}}.
\]
结合式 (5)，d=3。故每个尚存 CM 候选都包含
\[
                         K_0=\mathbb Q(\sqrt{-3}).              \tag{6}
\]
这里用的二次域来自极化判别式，不是根据某个数值周期猜测。

## 4. 一个辅助29-adic特征线排除全部 CM 候选

### 引理：P 的任意正整数幂谱都无重复

设 \(\alpha\ne\beta\) 是 P 的根。若 \(\alpha/\beta\) 为有限阶单位根，因为 P 不可约，可以在其分裂域中取一个 \(\mathbb Q\) 自同构，将 \(\alpha\) 送到第一节唯一的11-adic单位特征根。它将 \(\beta\) 送到另一个与它同11-adic赋值的根，也必须为单位根，因唯一性两根相等，矛盾。

所以任意不同根的比均不是有限阶单位根。特别地，对每个 \(n\ge1\)，所有 \(\alpha^n\) 互不相同。这个论证实际处理了“需要先取 Frobenius 的某个有限幂”的问题；只检查 P 本身无重根并不够。

### 29-adic 一维不变空间

直接模29计算
\[
 P(Z)\equiv(Z-1)(Z-5)(Z^4-13Z^3+Z^2-7Z-4)\pmod {29},            \tag{7}
\]
且 \(P'(1)\equiv1\)。Hensel 引理给出一个根 \(\lambda\in\mathbb Q_{29}\)。式 (3) 是未扭转 \(T_{29}\) 上 p=11 几何 Frobenius \(\Phi\) 的特征多项式。因此 \(\Phi\) 有一个一维 \(\mathbb Q_{29}\) 特征空间。前述引理保证，这个空间对每个 \(\Phi^n\) 仍然恰为一维。

K3 的 Hodge 自同态具有相容的29-adic实现，并与 Galois 群的某个有限指数子群交换；采用 Charles §2.2 第6页的绝对 Hodge 实现说明。因此存在 n 使 \(\Phi^n\) 与 \(F\otimes\mathbb Q_{29}\) 交换。若式 (6) 成立，则其一维特征线成为
\[
 K_0\otimes\mathbb Q_{29}=\mathbb Q_{29}(\sqrt{-3})
\]
的含幺模。但
\[
 (-3)^{14}\equiv-1\pmod {29},
\]
故这是一个二次域扩张，不能在一维 \(\mathbb Q_{29}\) 空间上作含幺作用。矛盾。

29只作为系数素数 \(\ell\)；此步骤不要求 X 在29处良约化，也没有新算29处点数。良约化 Frobenius 始终是已认证的 p=11。

由第二至四节得到
\[
                         \operatorname{End}_{\rm Hdg}T=\mathbb Q.
                                                                    \tag{8}
\]

## 5. 固定 Hodge 群

Zarhin 的 K3 群描述在权重零形式上给出
\[
                         G_H(S)=SO(q).                           \tag{9}
\]
可采用 Charles 的 K3 Mumford–Tate 定理及 Lemma 16 的连通群表述：E=\(\mathbb Q\) 时中心化群的恒等分支是 SO。Mumford–Tate 群连通，不能把未标恒等分支的 GO/O 符号机械照抄成不连通的 Hodge 群。其第8页的 \(T(1)\) 约定已经核对。

A=S⊗B²。S 的群 SO 是连通半单群，B 的群是 \(\mathbb G_m\)。两者没有非平凡共同代数群商，故
\[
 G_H(S\oplus B)=SO(q)\times\mathbb G_m.
\]
变回 A、B 的框架即得到式 (2) 的 \(R_H\)。这里证明的是固定点，不使用射线的几何群。

## 6. 比较行列式检测 Nori 群的另一分支

Hodge 实现以及实际纯核的杯积先给出夹逼
\[
 SO(q)\times\mathbb G_m
   \subseteq R_N
   \subseteq O(q)\times\mathbb G_m.                               \tag{10}
\]
两端仅差一个行列式符号，故只须判断一个实际一维对象
\[
                              D=\det S
\]
是否为单位对象。杯积在 Nori 范畴中给出 \(D^{\otimes2}\simeq\mathbb Q(0)\)。

用已有六循环有理 Betti 基和裸留数 B de Rham 框架记纯周期矩阵为 P。令 \(z=2\pi i\)，则 S 的比较矩阵为 \(C_S=z^{-1}P^{\mathsf t}\)，故 D 的比较系数为
\[
                           c_D=\det(P)/z^6.
\]
[ABSOLUTE_COMPARISON](ABSOLUTE_COMPARISON.md) 的实际定标给出
\[
                          P^{\mathsf t}K_BP=8\pi^2Q_6.
\]
独立精确代数在 s=−64 求得
\[
 \det K_B=\frac{353067798239177216000}{3},
\]
从而
\[
 \boxed{c_D^2=\frac{64\det Q_6}{\det K_B}
        =\frac9{86198192929486625}
        =\frac9{85^2\cdot11930545734185}.}                        \tag{11}
\]
其中
\[
 11930545734185=5\cdot13\cdot29\cdot37\cdot73\cdot89\cdot113\cdot233
\]
为无平方因子整数。因此 \(c_D\notin\mathbb Q\)。

若 D 在基域 \(\mathbb Q\) 的 Nori 范畴中为单位对象，任何有理 Betti/dR 基之间的比较系数都在 \(\mathbb Q^\times\)，与式 (11) 矛盾。这只用“实际单位对象必有有理比较”这一方向，没有使用周期映射的单射猜想，也没有从一般数值关系反推出动机同构。

故行列式角色在 R_N 上非平凡。因其值域为 \(\mu_2\)，式 (10) 的右端被全部取到，即
\[
                         R_N=O(q)\times\mathbb G_m.              \tag{12}
\]
同样得到 \(G_{\rm N,\mathbb Q}(S)=O(q)\)。

基域不能省略：这是 \(\mathbb Q\) 上的 Nori 群；本稿没有声称任意扩大基域后仍保留同一个非平凡分支。式 (11) 给出一个非平凡二次比较角色，但本稿也未仅凭周期数值就宣称 D 与某个指定 Artin 对象同构。

## 7. 完整固定活动群与 coaction 关系

[MIXED_UNIPOTENT_GROUP](MIXED_UNIPOTENT_GROUP.md) 已证明固定活动 E 的幺幂核为
\[
 U=\left\{\begin{pmatrix}I_6&u&w\\0&1&v\\0&0&1\end{pmatrix}\right\},
 \qquad \dim U=13.
\]
它同时证明实际群是其分级群在权重抛物群中的完整逆像。现代入式 (2)，在选定有理权重分裂后：
\[
 \boxed{
 G_{\rm N,\mathbb Q}(E)=
 \left\{\begin{pmatrix}t^2g&u&w\\0&t&v\\0&0&1\end{pmatrix}:
     g\in O(q),\ t\ne0,\ u,w\in\mathbb Q^6,\ v\in\mathbb Q\right\}.
 }                                                               \tag{13}
\]
式 (13) 按任意 \(\mathbb Q\) 代数解释为群概形；Hodge 群由再加 \(\det g=1\) 得到。这些群在适当分裂下分别为 U 与 \(O(q)\times\mathbb G_m\)、\(SO(q)\times\mathbb G_m\) 的半直积。两者维数 \(13+15+1=29\)，Nori 群两个分支，Hodge 群一个。

**基的提醒。** 原六循环的 q 为 Q6；若使用它们的对偶上同调 Betti 基，则 q 为 \(Q_6^{-1}\)，两者由有理杯积同构联系。de Rham 群须用实际 de Rham 杯积，而不能直接用 Q6 替代。旧 B 框架下为
\[
                         q_{\rm dR}=-2K_B^{-1}.
\]
物理 A=T(3) 上的 Tate 值域为 B⁴，因此 coaction 纯块 H、Tate 块 t 的全部关系可写为
\[
                  H^{\mathsf t}q_{\rm dR}H=t^4q_{\rm dR},
                  \qquad t\ne0.                                \tag{14}
\]
原来尚未识别的固定纯关系在此明确为正交相似关系；不额外强加 \(\det H=t^{12}\)，那会错误地删去 Nori 的另一分支。Hodge 群才增加该行列式等式。

式 (14) 连同自由的 u、v、w 和既有矩阵乘法余乘法，给出实际活动对象的完整 Hopf 描述。物理对偶对象仍取 \(C_M=C_E^{-\mathsf t}\)，既有
\[
 I^{\mathfrak m}=2(C_M^{\mathfrak m})_{8,6},\qquad
 \Delta I^{\mathfrak m}
 =2\sum_{j=1}^8(C_M^{\mathfrak m})_{8,j}
       \otimes(g_{E,\rm dR}^{-\mathsf t})_{j,6}
\]
不变，且其作用群现在已在指定纤维和基域下识别。这不意味着复周期值已经证明代数独立，也不把数值包围改称闭合特殊函数表达式。

## 8. 验证责任、前批审阅与范围

算术证书为 [certificate.json](data/fixed_pure_group/certificate.json)，独立复核为 [replay.json](data/fixed_pure_group/replay.json)。

生成器从旧点数与旧极化配对恢复本次见证。独立验证器不使用 SymPy 或生成器，采用整数多项式的模幂/gcd、独立有理表达式解析器、Fraction 和 Leibniz 行列式：重新验证模5不可约、模11唯一单位根、模29简单根及−3不为平方、det Q6、det K_B 和式 (11) 的全部素因子。改坏常数项的负对照被拒绝。没有新点数和纯运输。

本稿另外审阅 MIXED_PARAMETER_CONNECTION 的实际拉回和唯一性：子一动机沿前三坐标嵌入、F³顶层落在拉回、常数留数使导数进入纯核、配对三张量的对称性及 det−2 唯一性均相容；未改变其证明与600条目证书哈希。这项审阅不等于对全部早期几何稿的外部同行评审。

Charles 原 PDF 第3、6、8页及 Elsenhans–Jahnel 第21、22页已核对，下载哈希及原页图保存在 external/fixed-pure-group。实乘排除、CM 正规性、相容实现和群描述分别有明确条件，计算代码不承担这些理论定理。

本批只关闭指定 \(\mathbb Q\) 纤维的纯群及其活动群关系缺口。三种饱和整格类的实际标记、秩十环境绝对比较、有界一般图算法、独立物理基准和连贯论文仍继续。其他质量、其他 s 及扩大基域后的群均不作未经证明的推广。

