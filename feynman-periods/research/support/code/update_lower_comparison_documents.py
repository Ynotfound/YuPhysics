"""Reconcile the independently verified seven-dimensional comparison."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[1]
def read(n):return (ROOT/n).read_text(encoding='utf-8')
def write(n,s):(ROOT/n).write_text(s,encoding='utf-8')
def append(n,s):
    old=read(n)
    if s.strip() not in old:write(n,old.rstrip()+'\n\n'+s.strip()+'\n')
def main():
    marker='### 最新批次：完整下层七维绝对比较'
    checkpoint="""### 最新批次：完整下层七维绝对比较

2026-10-13 07:20（Asia/Shanghai）。生成与独立验证均已在命令内exit0，没有新的长运行会话。旧纯运输、局部种子、下层五球和所有旧证书未重算；全部旧哈希保留。证明为LOWER_COMPARISON.md，计算在data/banana_pure_ray/lower_comparison。

1. 已审阅LOWER_REGULATOR的实分裂必要条件及MIXED_UNIPOTENT_GROUP的相邻简单块/交换子/实现夹逼论证。指定纤维非零、最小秩八、完整13维幺幂核的范围保持。固定纯群仍未知。

2. 精确计算g_D to_D(s)，得到常数行g_0，无1/s项。用旧128次零点解析规范G及N_0²=0，构造原函数j(t)=h sum t^(n+1) g_0 G_n[I/(n+1)−N_0/(n+1)^2]。尾界≤|h| ||g_0||(1+||N_0||)epsilon_G；独立验证129个系数恒等式和尾界。

3. 只对既有64个零点纯多项式作精确积分，得到共同高基点j_*。实际第六球的beta值用b6=b1+j_* W^-1(v6−v1)，明确含混合单值增量。反向干路及Morse幅度误差均计入。s=−64的beta6约82.17960185016268−68.32617263832317i，严格误差<2^(-88)。

4. 边界点E12∩F3按该次序的小整环面tau双留数为+1，纯K类积分为0，beta积分精确(2pi i)^2。故P7=[[P6,0],[b1..b6,(2pi i)^2]]是完整V_c的有理Betti比较；49项矩阵无穷范数误差<2^(-85)，独立标量重放通过。F²在Betti坐标由omega、beta两行张成，F¹由omega、v1..v4、beta六行张成。物理下层N=V_c(3)。不是完整整格基，也不等于完整秩八E的比较。

5. sixth_period.json SHA ca3a03d8224639dc0c995129b9c12181f1842a292b275a8000d28cd68c061bbf；period_matrix.json SHA 264ddd733392c8b4a08afefebc88a29ac3a29d956d88698aff1b20467149f939。独立verify_lower_comparison.py重算所有新积分、框架作用、取向、Tate列、误差和49条目，拒绝零尾界并排除漏掉混合增量；replay.json已保存。

6. 下一理论项是完整秩八对象的顶层相对de Rham/Betti框架、实际链与其Hodge提升的比较；已有P7可直接作为下层输入。不要重算第六beta或边界Tate列。还须固定纯R_N识别、全四参数平坦性、三种饱和类的质量方向/边界几何区分、有界一般图算法、物理基准和主稿。整体目标未完成，单一heartbeat继续；不新增任务或聊天。

"""
    s=read('CONTINUATION.md')
    if marker not in s:
        s=s.replace('更新于 2026-10-13 07:05','更新于 2026-10-13 07:20',1)
        old='### 最新批次：固定基点 regulator 非零、秩八最小性及十三维幺幂群'
        assert old in s
        s=s.replace(old,checkpoint+'### 前一批次：固定基点 regulator 非零、秩八最小性及十三维幺幂群',1)
        start=s.index('下一批优先（以开头最新批次为准）：');end=s.index('已读取参考：',start)
        s=s[:start]+"""下一批优先（以开头最新批次为准）：

1. 完整下层P7已完成，审阅LOWER_COMPARISON后继续顶层相对比较：选择与当前下层七维基相容的相对de Rham提升及有理Betti链，固定Kummer单位与nu_E的符号，再计算缺少的顶层条目。纯P6、五球beta、第六beta、边界Tate列和旧运输均不重跑。
2. 固定纯群R_H/R_N与四参数平坦性保持独立任务；几何射线O(Q6)不能替代固定纯群。已证完整13维幺幂核不需反复重证。
3. 三种饱和类需质量方向或边界粘合；已证明全射线单值都保持三种，不能靠提高同一射线精度选出唯一类。
4. 将已完成的相对几何、非零性、秩八最小性、P7与群结构整合到论文；一般图算法及物理基准仍属完整总目标。

"""+s[end:]
        write('CONTINUATION.md',s)
    for name in ['TASKS.md','CLAIMS.md']:
        s=read(name)
        s=s.replace('当前有纯P6及五个球上的beta周期；第六混合条目和顶层仍待算','LOWER_COMPARISON已完成P7，包括第六混合条目和边界Tate列；顶层相对条目仍待算')
        write(name,s)
    append('TASKS.md',"""## 完整下层比较之后的当前任务

LOWER_COMPARISON已审阅前批非零/群论证明，并认证完整七维P7。第六球T0delta1的beta周期、边界Tate列及全部49项的有理包围已完成，矩阵误差<2^(-85)。H3/M1的下一项是秩八顶层相对比较、相容框架和实际链；固定纯群R_N、全四参数核对、饱和类几何、一般图、物理基准与论文整合仍继续。旧批次中的“第六beta待算”现仅为历史状态。""")
    append('CLAIMS.md',"""## 完整下层七维比较

| 主张 | 状态 | 范围 |
|---|---|---|
| 第六球T0delta1的beta周期 | 已由零点解析原函数和旧64步积分认证，误差<2^(-88) | 必须包含j_*W^-1(v6−v1)混合单值增量 |
| 边界Tate列 | E12∩F3双留数+1的小环面给(0,…,0,(2pi i)^2) | 双圆周次序和符号已固定 |
| 完整下层V_c七维比较P7 | 全49项及矩阵误差<2^(-85)，独立重放通过 | 明确有理de Rham/Betti基；不称完整整格基；物理下层须扭转(3) |
| 完整秩八比较 | 仍未完成 | 剩余顶层相对框架和比较条目；下层不再缺项 |""")
    s=read('PAPER.md')
    s=s.replace('Complete mixed comparison entries and the fixed pure motivic quotient remain to be determined.',
        'An explicit nilpotent primitive computes the sixth mixed period, and a boundary double-residue torus fixes the Tate column. Together these yield the complete rank-seven lower comparison at s=-64, with matrix error below 2^(-85). The top relative comparison entries and the fixed pure motivic quotient remain to be determined.')
    s=s.replace('This does not determine their absolute period amplitudes or the marked saturation coset.',
        'Morse normalization separately determines the absolute pure period amplitudes; the marked saturation coset remains undetermined.')
    write('PAPER.md',s)
    append('PAPER.md',r"""## 可纳入正文的完整下层比较定理

**Theorem (the lower comparison at a prescribed algebraic fibre).** Let m=(1,2,4,8), s=-64, and let V_c be the extension of Q(-2) by T(X) selected by c=3P₁+2P₂+P₃. Fix the bare Poincare residue, the rational pure frame B, and the oriented vanishing cycles δ₁,…,δ₅, δ₆=T₀δ₁. Let τ be the ordered boundary torus at E12∩F3. In the form basis (B,β_c) and the cycle basis (δ₁,…,δ₆,τ), the full comparison matrix is

\[
P_7=\begin{pmatrix}P_6&0\\ \boldsymbol b&(2\pi i)^2\end{pmatrix},
\qquad \boldsymbol b=(b_1,\ldots,b_6).
\]

Every entry has a certified rational enclosure; the matrix infinity-norm error is less than 2⁻⁸⁵. The sixth mixed entry has error less than 2⁻⁸⁸. These cycles form a rational Betti basis for the dual of V_c. No assertion of integral saturation is needed.

**Proof.** The first five mixed entries are fixed by collapse of the actual vanishing spheres and the exact integral of the certified pure transport. In the zero-point regular frame the lower coupling is a constant row g₀. For the saved analytic gauge G(t) and N₀²=0, the row j(t)=hΣtⁿ⁺¹g₀Gₙ[I/(n+1)−N₀/(n+1)²] satisfies j′+jN₀/t=hg₀G. Its tail is bounded by integrating the certified gauge error against 1+|log t|‖N₀‖. Transport gives j_* at the common basepoint, and the complete mixed monodromy gives b₆=b₁+j_*W⁻¹(v₆−v₁). Both the previously certified pure frame and the new row integral carry their input, tail and rounding errors. At the boundary flag the double residue of β is +1, so τ has β-period (2πi)² and all pure periods zero. Nondegeneracy follows from det P₇=(2πi)² det P₆. Independent scalar replay verifies all new coefficients, quadratures, orientations and error estimates. Details and reproducible certificates are in LOWER_COMPARISON.

This theorem completes the lower comparison only. For N=V_c(3), apply the explicit Tate factor (2πi)⁻³ in the corresponding tensor bases. The remaining top relative column of the rank-eight object, its compatible framings, and the fixed pure motivic quotient require separate work.

当前正文整合以本段和开头更新后的摘要为准；此前各计算批次段落记录了递进证据，不能把其中旧的“尚未完成”再次写入当前定理范围。""")
    append('README.md',"""## 当前结果：完整下层比较

[LOWER_COMPARISON](LOWER_COMPARISON.md)已补齐第六混合周期及边界Tate列，完成s=−64下层七维绝对比较P7，全部49项误差在矩阵无穷范数下<2^(-85)。生成和独立复核均完成，旧纯证书未修改。第六周期使用完整混合单值增量；边界列来自双留数+1的小整环面。下一项转入秩八顶层相对比较，保持固定纯群、四参数、饱和类、一般图、物理基准和主稿的完整范围。""")
    append('SOURCE_AUDIT.md',"""## 完整下层七维比较的证据审计

新增LOWER_COMPARISON及两份计算/验证入口。g_D to_D为常数的恒等式直接来自保存的几何耦合及精确规范；原函数由N₀²=0逐项推导。独立检查129个系数、64个已有纯步的新积分行、实际混合单值增量、前五球相对于Q6列的+取向以及边界Tate负实值。原有纯证书和五球证书保持哈希。Sato §2.2只用于核对一般regulator的Hodge商约定，不将其非消失定理移植到banana。完整P7误差<2^(-85)，顶层相对比较及固定纯群仍未完成。""")
    p=ROOT/'data/input_manifest.json';d=json.loads(p.read_text())
    d.update(full_lower_seven_dimensional_comparison_certified=True,
      lower_comparison_proof='../LOWER_COMPARISON.md',
      lower_seven_comparison_certificate='banana_pure_ray/lower_comparison/period_matrix.json',
      lower_seven_comparison_replay='banana_pure_ray/lower_comparison/replay.json',
      lower_seven_comparison_uniform_error='2^(-85)',
      sixth_lower_beta_period_certificate='banana_pure_ray/lower_comparison/sixth_period.json',
      sixth_lower_beta_period_uniform_error='2^(-88)',
      lower_boundary_Tate_flag=['E12','F3'],
      lower_boundary_Tate_period='(2*pi*i)^2',
      full_active_eight_dimensional_comparison_certified=False)
    p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    archive=json.loads(read('data/banana_pure_ray/lower_regulator/final_audit.json'))
    for f,h in archive['old_hashes_preserved'].items():
        assert hashlib.sha256((ROOT/'data/banana_pure_ray'/f).read_bytes()).hexdigest()==h
    for f,h in archive['new_proof_sha256'].items():
        assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h
    print('Documents reconciled; old pure and preceding proof hashes preserved.')
if __name__=='__main__':main()

