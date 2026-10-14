"""Save the complete relative comparison checkpoint without modifying old evidence."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[1]
def read(n):return (ROOT/n).read_text(encoding='utf-8')
def write(n,s):(ROOT/n).write_text(s,encoding='utf-8')
def append(n,s):
    old=read(n)
    if s.strip() not in old:write(n,old.rstrip()+'\n\n'+s.strip()+'\n')
def main():
    checkpoint="""### 最新批次：完整秩八相对比较、物理链标记与物理积分

2026-10-13 07:38（Asia/Shanghai）。新顶层生成15774、独立复核8207均exit0；物理读出生成也已命令内exit0。无本项目待观察计算。原纯运输、P7、下层regulator及群论证明哈希均保留。证明为RELATIVE_TOP_COMPARISON.md；新证书在data/banana_pure_ray/relative_top。

1. 相对复形d(alpha,beta)=(dalpha,−i*alpha−dbeta)和配对int_Gamma alpha+int_boundary beta固定后，顶层e=(eta,0)、eta=dlogx∧dlogy∧dlogz满足nabla_s e=delta omega，mu_i nabla_mu_i e=delta beta_i。这是实际相对Cartan公式，不是自由添原函数。共同缩放由s omega+sum beta_i=0消去。完整四参数纯平坦性仍需独立核对。

2. 五个实际消失膜Gamma_a以收缩零初值唯一固定q_a=int eta，其导数为omega球周期。精确积分旧320个纯节点步。零点顶层耦合g_eta to_D也为常数，复用旧128次规范和64步，q6=q1+j_eta,* W^-1(v6−v1)，明确含相对顶层单值。没有重跑任何纯ODE。

3. 边界环面先用e_(E12,F3)^t(I−ZR)修正，使它在H2(T)中为零且不改变P7的边界列，然后取有理填充与其复共轭的平均Gamma7+。质量旗标给q7/(2pi i)^2=−2log2模2pi i Q；实提升消去纯虚歧义，故q7=8pi²log2。Gamma8=T³的eta周期为(2pi i)^3。

4. 完整未扭转P8=[[P6,0,0],[b,(2pi i)^2,0],[q,−2(2pi i)^2log2,(2pi i)^3]]的全部64项已独立认证，矩阵误差<2^(-85)。这是E_act(-3)的明确有理基比较；实际dR→Betti算子C_E=(2pi i)^(-3)P8^t。特定饱和整基仍未给出。矩阵SHA=0eeca4619c39f5a3c3156dd1d31beeb43585e838bb824d8565037fb14c46c5b2。

5. 正实坐标上f≥(sum m)^2=225，而前六消失膜的f路径避开[225,∞)，所以它们与正物理链不相交。Gamma7+与正链均共轭不变，ambient复三维共轭反定向，交数亦零。正链在有序T³上的交数按显式坐标为−1；沿用顶层+1归一化取负的有序相对交数，故sigma_B=(0^7,1)。这个拓扑标记由证明承担，不由数值猜整数。

6. nu=e_F−sigma满足P7^t nu=(q1..q7)^t，nu_beta=−2log2。局部环带和上述整体对偶符号给物理读出−C_geom(omega,alpha)，结合C_geom(omega,w_published)=−2及F过滤，I=−2nu_w。平衡D框架的精确Neumann逆/输入传播认证I(-64)=0.735547022661485632257667775451附近，误差<2^(-82)，粗严格区间(0.73554702266148,0.73554702266149)。readout SHA=8cb313ee98adccc4811b116f5c90d8068473ce862f52eebd647883204e2d33e3；独立replay SHA=ca66d4beb668ac4c9f6a9febe15385466ff843fb33109350ed2f4bae54c10307。原参数24/40/64阶求积一致，仅诊断，不承担认证。

7. 全部实际比较块已定，物理M=E∨有C_M=C_E^-t，I^m=2(C_M^m)_(8,6)，右coaction=2 sum_j(C_M^m)_(8,j) tensor(g_E^-t)_(j,6)。既有完整13维幺幂核仍适用；固定纯R_N未知，不能以射线O(Q6)替换。不能再机械列顶层比较或指定物理读出为缺失，也不能把数值包围称为全部闭合特殊函数表达式。

8. 下一项先审阅本稿实际膜、共轭填充、整体对偶符号与物理标记的几何证明，随后继续完整四参数平坦性/质量方向标记、固定纯群和主稿整合。有界一般图算法、另行定义的物理基准及全文仍是完整总目标。P8和I不等于所有物理基准完成。单一heartbeat继续；没有缩小目标或标记整体完成。

"""
    s=read('CONTINUATION.md');marker='### 最新批次：完整秩八相对比较、物理链标记与物理积分'
    if marker not in s:
        s=s.replace('更新于 2026-10-13 07:20','更新于 2026-10-13 07:38',1)
        old='### 最新批次：完整下层七维绝对比较';assert old in s
        s=s.replace(old,checkpoint+'### 前一批次：完整下层七维绝对比较',1)
        start=s.index('下一批优先（以开头最新批次为准）：');end=s.index('已读取参考：',start)
        s=s[:start]+"""下一批优先（以开头最新批次为准）：

1. 审阅RELATIVE_TOP_COMPARISON的膜扫迹避开物理阈值、边界投影/共轭填充、相对复形与整体对偶符号。P8数值和物理读出已独立重放；几何标记靠文中证明，不应拿浮点吻合作替代。不重复纯运输、P7或新顶层生成。
2. 完整四参数纯连接平坦性以及质量方向的几何标记继续；顶层公式nabla e=delta omega ds+sum delta beta_i dlog mu_i现已有，可作多参数相容性检查。
3. 固定纯R_H/R_N与三种饱和类的几何区分仍待独立推进。已知射线O(Q6)和13维混合幺幂核不能代替这两项。
4. 推进有界一般图算法、范围明确的物理基准、文献差异与全文整合。指定banana的认证I是一个物理检查点，不自动完成massive phi4等另设基准。

"""+s[end:]
        write('CONTINUATION.md',s)
    append('TASKS.md',"""## 完整八维比较之后的当前任务

H3/M1的指定纤维比较部分已由RELATIVE_TOP_COMPARISON完成：实际P8、正物理链标记、nu_E全部坐标和I(-64)的严格读出均有证明及独立证书。P8误差<2^(-85)，I误差<2^(-82)。M1剩余固定纯分级群R_N的识别及其关系；H5d的全四参数平坦性、质量方向的饱和标记、G1一般图、F1另设物理基准、P2全文和文献审阅仍继续。首先审阅新顶层几何/对偶符号，然后推进这些未完成项。旧批次的“完整八维比较未算”已被本批更新。""")
    s=read('TASKS.md').replace('环境和活动 Yoneda 消失已推导，nu_E 的纯核提升待计算','环境和活动Yoneda消失已推导；指定s=−64的nu_E纯核提升与物理标记现已由RELATIVE_TOP_COMPARISON认证')
    s=s.replace('固定纯群和完整混合比较仍未识别','完整混合P8及实际物理coaction系数已认证；固定纯分级群仍未识别')
    write('TASKS.md',s)
    s=read('CLAIMS.md')
    s=s.replace('| 全部秩八比较条目已求值 | 未完成 | LOWER_COMPARISON已完成P7，包括第六混合条目和边界Tate列；顶层相对条目仍待算 |',
      '| 全部秩八比较条目已求值 | RELATIVE_TOP_COMPARISON已认证P8全部64项 | 指定s=−64、明确有理几何基；不是饱和整基或闭合特殊函数表达式 |')
    s=s.replace('| 完整秩八比较 | 仍未完成 | 剩余顶层相对框架和比较条目；下层不再缺项 |',
      '| 完整秩八比较 | RELATIVE_TOP_COMPARISON已完成 | 消失膜、实边界填充及T³给明确Betti基；需按文中几何证明理解其标记 |')
    write('CLAIMS.md',s)
    append('CLAIMS.md',"""## 顶层相对比较与物理标记

| 主张 | 状态 | 范围 |
|---|---|---|
| 顶层相对连接nabla_s e=delta omega | 实际Cartan公式与相对复形符号证明 | 同时得到mu_i nabla_mu_i e=delta beta_i；不代替全多参数纯平坦性 |
| 实际消失膜的六个eta周期 | 收缩初值、旧纯多项式积分及完整T0增量独立认证 | 不重跑纯ODE |
| 边界顶层项q7=−2(2pi i)^2log2 | 环境投影后取实共轭平均填充，Kummer类固定余下常数 | 原环面不能未经环境零同调检验就填充 |
| 完整P8 | 全64项矩阵误差<2^(-85)，独立重放 | 实际E_act(-3)的有理基；非饱和整格声明 |
| 正物理链在对偶相对基为(0^7,1) | 阈值225的支撑避让、边界膜共轭、T³交数归一化 | 拓扑证明，不是数值取整 |
| 指定物理I(-64) | I=−2nu_w，误差<2^(-82)，严格区间(0.73554702266148,0.73554702266149) | 无额外动量测度因子的正欧氏参数积分 |
| 实际物理coaction的全部比较系数 | C_M=C_E^-t，I^m=2(C_M^m)_(8,6)，八项基公式已指定 | 固定纯R_N关系仍未知，不能称整个绝对动机群已识别 |""")
    s=read('PAPER.md')
    s=s.replace('The top relative comparison entries and the fixed pure motivic quotient remain to be determined.',
      'Collapsing relative thimbles and a real boundary filling now determine the full rank-eight comparison. A geometric intersection argument identifies the positive physical chain in this frame and yields a certified Euclidean integral at s=-64, with error below 2^(-82). The fixed pure motivic quotient, marked integral saturation, and full multiparameter verification remain separate problems.')
    write('PAPER.md',s)
    append('PAPER.md',r"""## 可纳入正文的完整相对比较与物理积分定理

**Theorem.** At m=(1,2,4,8) and s=-64, let e=(dlog x∧dlog y∧dlog z,0) be the top relative Hodge lift. In the form basis (B,β_c,e), use the six collapsing relative thimbles Γ₁,…,Γ₅,Γ₆=T₀Γ₁, the real conjugation-invariant filling of the ambient-corrected boundary torus, and T³ as the rational relative Betti basis. Then the full comparison of E_act(-3) is

\[
P_8=\begin{pmatrix}
P_6&0&0\\
\boldsymbol b&(2\pi i)^2&0\\
\boldsymbol q&-2(2\pi i)^2\log2&(2\pi i)^3
\end{pmatrix}.
\]

All entries admit an independently replayed rational enclosure with matrix error below 2⁻⁸⁵. In the normalized relative duality, the positive physical chain has coordinates (0⁷,1). Writing ν=e_F−σ in the B,β frame gives ν_β=−2log2 and I(-64)=−2ν_w. The resulting physical parameter integral has error below 2⁻⁸² and lies strictly between 0.73554702266148 and 0.73554702266149.

**Proof.** Cartan's formula in the fixed relative cone gives ∇_s e=δω and μ_i∇_μ_i e=δβ_i. Collapsing thimbles fix the six scalar integration constants; a constant coupling in the zero-point regular gauge permits the actual relative T₀ increment to be recovered from archived pure certificates. The boundary torus is first projected by I−ZR to kill its ambient homology class, then filled and averaged under real conjugation. The existing Kummer class and reality fix its η-period exactly. The positive simplex satisfies f≥225, whereas the first six thimble paths avoid [225,∞); their intersection numbers vanish. The boundary filling has vanishing intersection by conjugation, and the ordered T³ intersection fixes the remaining sign. A local annulus computation combines this normalization with the independently calibrated pure cup row to give I=−2ν_w. Exact residual bounds and a Neumann inverse certify the values. See RELATIVE_TOP_COMPARISON for the complete geometry, sign conventions, bounds and hashes.

This supplies the actual eight-term coaction in the dual comparison frame: C_M=C_E^-t, I^m=2(C_M^m)_(8,6), and ΔI^m=2Σ_j(C_M^m)_(8,j)⊗(g_E^-t)_(j,6). The fixed pure quotient of the motivic group is retained as R_N; its relations have not been replaced by those of the geometric ray monodromy group. The complete comparison theorem is specific to this fibre and these rational geometric bases.""")
    append('README.md',"""## 当前结果：完整相对比较与物理积分

[RELATIVE_TOP_COMPARISON](RELATIVE_TOP_COMPARISON.md)完成指定纤维的完整八维比较P8、正物理链标记与nu_E坐标。P8误差<2^(-85)，物理参数积分I(-64)约0.7355470226614856，认证误差<2^(-82)，独立复核通过。旧证书和P7保持不变；下一项是新几何证明审阅、固定纯群、多参数平坦性、饱和标记、一般图算法、另设物理基准与全文整合。""")
    append('SOURCE_AUDIT.md',"""## 顶层比较及物理符号审计

RELATIVE_TOP_COMPARISON直接从相对复形和Cartan公式推出顶层连接。BKV §6.2原始公式6.2.3–6.2.5已在原HTML核对：其物理差为sigma−s_F，与本项目nu=e_F−sigma相反。膜的零初值、边界环境投影及实填充、阈值225的交数避让由本文单独证明。按正实log坐标后接单位圆坐标的排列有三次交换，明确T³交数为−1；归一化顶层+1后，物理读出取负的几何留数杯积。结合旧C_geom(omega,w)=−2得到I=−2nu_w。粗非严格求积只作诊断，没有取代这些符号证明和严格区间。""")
    p=ROOT/'data/input_manifest.json';d=json.loads(p.read_text())
    d.update(full_active_eight_dimensional_comparison_certified=True,
      full_relative_comparison_proof='../RELATIVE_TOP_COMPARISON.md',
      full_relative_comparison_certificate='banana_pure_ray/relative_top/period_matrix.json',
      full_relative_comparison_replay='banana_pure_ray/relative_top/replay.json',
      full_relative_comparison_uniform_error='2^(-85)',
      physical_positive_chain_Betti_coordinates=[0]*7+[1],
      physical_readout_certificate='banana_pure_ray/relative_top/physical_readout.json',
      physical_readout_uniform_error='2^(-82)',
      physical_integral_strict_interval=['0.73554702266148','0.73554702266149'],
      physical_integral_scope='2D banana m=(1,2,4,8), s=-64, no extra momentum-measure normalization',
      explicit_physical_coaction_framing_certified=True,
      full_mixed_comparison_special_function_closed_forms_complete=False,
      absolute_motivic_group_identified=False)
    p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('Relative-top evidence reconciled; full goal retains all remaining work.')
if __name__=='__main__':main()

