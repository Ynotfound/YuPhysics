"""Save marked saturation, integral comparisons, and a nonduplicating checkpoint."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import json,hashlib,re
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/marked_saturation'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rd(p):return (ROOT/p).read_text(encoding='utf-8')
def wr(p,s):(ROOT/p).write_text(s,encoding='utf-8')
def append(p,s):
    old=rd(p)
    if s.strip() not in old:wr(p,old.rstrip()+'\n\n'+s.strip()+'\n')
def main():
    c=json.loads(rd('data/marked_saturation/certificate.json'));v=json.loads(rd('data/marked_saturation/replay.json'))
    assert v['status']=='marked_saturation_independently_replayed'
    assert v['certificate_sha256']==sha(OUT/'certificate.json') and v['periods_sha256']==sha(OUT/'periods.json')
    assert v['verifier_sha256']==sha(ROOT/'code/verify_marked_saturation.py')
    prior=json.loads(rd('data/fixed_pure_group/final_audit.json'))
    old=dict(prior['preserved_old_hashes'])
    for k,h in c['inputs_sha256'].items():old[k]=h
    old['FIXED_PURE_GROUP.md']=prior['proof_sha256']
    for k,h in old.items():assert sha(ROOT/k)==h
    stamp=datetime.now(timezone.utc).astimezone(timezone(timedelta(hours=8))).strftime('%Y-%m-%d %H:%M')
    latest=f"""### 最新批次：判别式二次型唯一标记饱和格与整基绝对比较

{stamp}（Asia/Shanghai）。MARKED_SATURATION.md为新证明，data/marked_saturation保存精确证书。生成、独立整数/Fraction验证和文档更新均命令内exit0，无待观察进程。旧P6/P7/P8、固定群和全部既有证书哈希保持，没有新点数、闭路或ODE运输。

1. 审阅FIXED_PURE_GROUP的定理条件、基域Q和deRham杯积约定，旧证明未改。先探索复共轭及判别群作用，未因此选出半循环；probe_lattice_discriminant_action.py仅为探索，不是新证书输入。
2. 实际L⊂T、Q6及[T:L]=2已有证据。任何补全形如L+Z(v/2)，v非零模2，整性/偶性为Q6v=0模2及v^tQ6v=0模8，共11类。完整已证T=U⊕U(2)⊕A2(-2)的判别式二次型在A_T[2]上全取整数模2。这比行列式、Smith形和单值保持性更强。
3. 全部11类的二阶判别元逐项枚举，其中10类各有8个非整数值，并给出一个值3/2的明确反例；唯一余下v=(1,1,1,1,0,0)。故w=(delta1+delta2+delta3+delta4)/2是真实整类，(w,delta2,...,delta6)为完整整基。无需再认证四个旧质量矩阵来证明这一点；之前三种歧义仅在“射线单值保持性”准则下仍成立。
4. 旧标准同构R现通过detR=1及R^t G R=U⊕U(2)⊕A2(-2)的独立精确复核，成为实际循环标记下的标准整基。完整T与Q6的几何识别、半循环及标准变换均完成。
5. 原绝对P6直接右乘B和BR，Gaussian整数中心改用分母2^321，无新舍入。两个完整整Betti基比较矩阵误差分别<2^-84、<2^-81；72项独立重放。仍采用旧明确有理dR B框架，不是旧theta基；不把P7/P8自动升级为混合整格基。
6. certificate SHA {sha(OUT/'certificate.json')}；periods SHA {sha(OUT/'periods.json')}；replay SHA {sha(OUT/'replay.json')}。原Hessian论文第3页整格公式已再次查看。所有旧证书保留，完整目标未完成。
7. 下一项优先G1有界一般图算法与F1独立物理基准，再把积累证明整理成连贯全文。不要继续把实际饱和半循环列为缺失，也不必为此启动质量闭路。四个旧质量矩阵和旧theta框架仍未比较；如继续须独立范围。秩十P10可另行拓展但不替代一般图、物理基准和全文。

"""
    s=rd('CONTINUATION.md');marker='### 最新批次：指定Q纤维纯群与完整活动群'
    if '### 最新批次：判别式二次型唯一标记饱和格与整基绝对比较' not in s:
        assert marker in s
        s=re.sub(r'更新于 .*?（Asia/Shanghai）',f'更新于 {stamp}（Asia/Shanghai）',s,count=1)
        s=s.replace(marker,latest+'### 前一批次：指定Q纤维纯群与完整活动群',1)
        a=s.index('下一批优先（以开头最新批次为准）：');b=s.index('已读取参考：',a)
        s=s[:a]+"""下一批优先（以开头最新批次为准）：

1. 审阅MARKED_SATURATION的二阶判别式型与十一类穷尽证明。指定实际饱和半循环、标准整基及两份整基绝对P6已完成；不重跑运输、不继续为已消除的三种格歧义启动质量闭路。
2. G1下一批主线：从原稿恢复组合定义与可判定任务，规定有限输入/输出/边界，执行可复核的有界一般图算法；不把字符支撑等同动机分类，不把一个coNP命题外推至全问题。先核对原页和一手文献。
3. F1并行后续：定义并完成独立物理基准的有限范围与归一化，banana单点积分不替代。实标量phi4无需虚构规范Ward恒等式。
4. P1/P2：将现有证明集合整合连贯正文，清理过时待证语句与来源归属；一般可见性过滤若无嵌套对象证明继续标条件。
5. 旧theta框架、四个旧质量矩阵仍未几何比较，但不再是实际整格认证的依赖。P7/P8仍为明确有理混合基，纯格整基不自动认证混合整基。秩十P10/群若拓展须另证边界列与填充，不能取代核心未执行任务。
6. 本批全部命令已结束，无待观察项目计算；旧哈希保留。单一heartbeat继续，get_goal旧blocked状态不代表当前无进展，不创建重复目标或自动任务。

"""+s[b:];wr('CONTINUATION.md',s)
    s=rd('TASKS.md')
    replacements={
      'H4e':'三个样本完整T已证；MARKED_SATURATION由二阶判别式型唯一认证实际w及标准整基，无需新增质量闭路',
      'Q4':'21迹、实际Q6、指数二及实际T0已证；十一补全由判别式二次型唯一选出w=(delta1+delta2+delta3+delta4)/2',
      'Q5':'原六循环绝对P6已认证；完整整基及标准整基P6现误差<2^-84、<2^-81；旧theta框架未比较'}
    lines=[]
    for line in s.splitlines():
        for key,status in replacements.items():
            if line.startswith('| '+key+' |'):
                parts=line.split('|');parts[-2]=' '+status+' ';line='|'.join(parts)
        lines.append(line)
    s='\n'.join(lines)+'\n'
    a=s.index('## 当前理论批次');b=s.index('## 决策门槛',a)
    s=s[:a]+"""## 当前理论批次

1. H4e/Q4/Q5：MARKED_SATURATION完成实际半循环、完整整基、标准格变换及整基绝对P6；十一类判别式型与72个周期条目独立复核。审阅新证明后转入尚未执行任务，不为已解的格歧义重复运输。
2. G1：下一主线为有界一般图算法。恢复定义、确定有限输入输出、证明终止和复杂度、执行有独立证书的小型覆盖。不能自动称任意图动机分类。
3. F1：实施范围明确的独立物理基准，写清方案、可观测量和检查内容；banana单点读出不替代此项。
4. P1/P2：从证明包改写连贯全文，完成引用、新颖性和状态审计。M1的固定Q纤维纯/混合群及完整有理比较已经明确，不重复列缺口。
5. 旧theta框架、四个质量矩阵尚未比较；秩十环境绝对比较及群可另作拓展。新纯整基不自动构成完整混合整基。下方历史段落由本节和最新检查点覆盖。

"""+s[b:];wr('TASKS.md',s)
    s=rd('CLAIMS.md')
    s=s.replace('detQ6=192及detT=48已证；指定补全半循环及四质量矩阵的几何标记仍待认证',
      'detQ6=192及detT=48已证；MARKED_SATURATION现认证指定半循环；四质量矩阵另未认证')
    wr('CLAIMS.md',s)
    append('CLAIMS.md',"""## 实际饱和半循环与完整整基比较（更新先前补全未证状态）

| 主张 | 证据 | 范围 |
|---|---|---|
| 实际w=(delta1+delta2+delta3+delta4)/2为整类 | 已证L⊂T与指数二；完整T的二阶判别式型为整数值；10份非整数反例排除全部其他偶补全 | 唯一格标记，不依赖四个候选质量矩阵 |
| 完整整基与标准同构 | B替换首列，detR=1及R^t(B^tQ6B)R=标准T | 原R虽来自候选稿，现由实际输入基和精确等式赋予几何标记 |
| 完整整基绝对P6 | P_sat=P6 B、P_std=P6 BR，中心精确分母2^321 | 误差<2^-84、<2^-81，旧有理dR B框架；不是旧theta框架 |
| 三种歧义是否仍阻塞格认证 | 不再阻塞；保持性有三种，但完整判别式二次型只容许一种 | 旧射线保持性限制仍正确，不能再泛称所有现有证据有三解 |
| 混合整基与旧质量矩阵 | 本批未认证 | 纯整格结果不自动认证P7/P8整性或未运输质量闭路 |
""")
    s=rd('PAPER.md').replace('Q6 认证与条件性指数二补全','Q6 认证与实际指数二补全')
    s=s.replace('Marked integral saturation, the broader graph algorithm, an independent physical benchmark and the optional full rank-ten comparison remain separate tasks.',
      'The two-torsion discriminant quadratic form uniquely identifies the marked index-two saturation, and exact basis changes give the full integral pure comparison. The broader graph algorithm, an independent physical benchmark and the optional full rank-ten comparison remain separate tasks.')
    s=s.replace('the marked saturation coset remains undetermined.',
      'the two-torsion discriminant form selects the unique actual saturation coset (1,1,1,1,0,0)/2.')
    s=s.replace('固定纯群仍待识别；','固定Q纤维纯群现由FIXED_PURE_GROUP识别；')
    s=s.replace('特定补全半循环仍待标记；','特定补全半循环现由MARKED_SATURATION认证；')
    s=s.replace('具体指数二补全类、绝对周期矩阵、质量方向标记和固定s=−64的下层regulator非零仍需各自的证据。',
      '具体补全、绝对比较与固定regulator非零分别见后续证明；四个旧质量矩阵的几何标记另未认证。')
    s=s.replace('完整整基的具体补全与旧theta框架比较另列。',
      '完整整基补全由MARKED_SATURATION证明，旧theta框架比较另列。')
    wr('PAPER.md',s)
    append('PAPER.md',r"""## 实际饱和标记及整基绝对比较定理

**Theorem.** At the prescribed fibre, retain the six geometrically certified oriented cycles and their Gram matrix \(Q_6\). Then
\[
 T(X,\mathbb Z)=L+\mathbb Z\frac{\delta_1+\delta_2+\delta_3+\delta_4}{2},
 \qquad L=\sum_i\mathbb Z\delta_i.
\]
Replacing the first cycle by this half-sum gives a complete integral basis. The matrix \(R\) in MARKED_SATURATION has determinant one and sends its Gram matrix to \(U\oplus U(2)\oplus A_2(-2)\). The absolute comparison matrices in these two integral Betti bases have infinity-norm errors below \(2^{-84}\) and \(2^{-81}\), respectively.

**Proof.** The actual inclusion and index two restrict all possible lattices to the eleven even overlattices \(K_v=L+\mathbb Z(v/2)\). In the known full transcendental lattice, the discriminant quadratic form takes integral values modulo two on all two-torsion elements: on \(U(2)\) it is \(ab\), and on \(A_2(-2)\) it is \(-c^2+cd-d^2\). For ten of the eleven candidates, the explicit witness table in MARKED_SATURATION supplies a dual half-vector of square \(3/2\) modulo two. The sole remaining class is \(v=(1,1,1,1,0,0)\); existence of the actual lattice proves that it is the required class. Direct integer multiplication then verifies the stated standard isometry.

Let B replace the first identity column by v/2. Right multiplication of the existing period matrix by B and BR gives the two integral comparisons. Since their infinity norms are \(3/2\) and 10, the previously certified error propagates to the stated bounds. All centres are transformed exactly on the dyadic grid with denominator \(2^{321}\). No new analytic continuation is needed.

This proof uses the full geometric lattice identification, not just its determinant or discriminant group. It does not authenticate the old mass-direction matrices or identify an integral basis for the mixed extensions.""")
    append('README.md',"""## 当前结果：实际饱和格与整基纯比较

[MARKED_SATURATION](MARKED_SATURATION.md)用二阶判别式二次型唯一确定w=(delta1+delta2+delta3+delta4)/2为实际整类，补齐完整整基和标准格同构。十种其他候选各有明确排除见证。完整整Betti基及标准基的绝对P6误差分别<2^-84、<2^-81，72项独立精确重放，旧证书未改。三种射线保持性候选不再阻塞格标记。下一重点转入有界一般图算法、独立物理基准和论文全文。""")
    append('SOURCE_AUDIT.md',"""## 饱和标记的二阶判别式输入

MARKED_SATURATION再次核对Dardanelli–van Geemen §1.6原PDF第3页的完整整格公式；适用指定纤维所需Picard数和本原性仍引用LATTICE_SATURATION已有证明。本批新增的是用q_K在A_K[2]上的整数值性质排除十种其他补全，并非新格分类定理。全部候选与半向量反例由两个独立精确程序重放。旧标准同构R只作为待核验矩阵输入，detR与Gram等式再次认证；没有使用四个旧质量矩阵的真实性。探索性复共轭计算未构成认证输入。固定纯群前批证明的条件和基域已复查，原哈希不变。""")
    m=json.loads(rd('data/input_manifest.json'))
    m.update({
      'comparison_status':'Specified-fibre pure and mixed rational comparisons, fixed-Q active groups, actual marked saturation and complete integral pure comparisons proved. Bounded graph algorithm, independent physical benchmark and full manuscript remain.',
      'specific_half_cycle_completion_certified':True,'specific_half_cycle_numerator':[1,1,1,1,0,0],
      'full_integral_Betti_basis_certified':True,'full_integral_Betti_basis_scope':'Pure T only; (w,delta2,...,delta6) and its explicit standard basis.',
      'geometric_marking_scope':'Actual six cycles, Gram, index, unique saturation half-class and full integral standard lattice basis.',
      'absolute_comparison_scope':'Bare residue rational B de Rham frame; original six-cycle rational matrix retained, plus certified full integral saturated and standard Betti bases.',
      'saturation_method':'Two-torsion discriminant quadratic form excludes ten of eleven even index-two completions.',
      'saturation_proof':'../MARKED_SATURATION.md','saturation_certificate':'marked_saturation/certificate.json',
      'saturation_replay':'marked_saturation/replay.json','integral_pure_periods':'marked_saturation/periods.json',
      'integral_pure_period_errors':['<2^-84','<2^-81'],
      'mixed_integral_betti_basis_certified':False,'old_mass_matrices_geometrically_certified':False})
    wr('data/input_manifest.json',json.dumps(m,ensure_ascii=False,indent=2)+'\n')
    audit={'status':'actual_marked_saturation_and_integral_pure_comparisons_complete','timestamp':stamp,
      'proof_sha256':sha(ROOT/'MARKED_SATURATION.md'),'certificate_sha256':sha(OUT/'certificate.json'),
      'periods_sha256':sha(OUT/'periods.json'),'replay_sha256':sha(OUT/'replay.json'),
      'old_hashes_preserved':old,'new_ODE_transports':0,'all_started_calculations_ended':True,
      'fixed_pure_proof_reviewed_without_change':True,'mass_loop_matrices_certified':False,
      'mixed_integral_basis_certified':False,'full_goal_complete':False}
    wr('data/marked_saturation/final_audit.json',json.dumps(audit,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in audit.items() if k!='old_hashes_preserved'},ensure_ascii=False))
if __name__=='__main__':main()

