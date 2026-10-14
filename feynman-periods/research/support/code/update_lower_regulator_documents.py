"""Reconcile the fixed-fibre proof into project checkpoints; preserve certificates."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]

def edit(name,changes):
    p=ROOT/name;s=p.read_text(encoding='utf-8')
    for old,new in changes:
        assert old in s,(name,old[:80])
        s=s.replace(old,new)
    p.write_text(s,encoding='utf-8')

def append(name,text):
    p=ROOT/name;s=p.read_text(encoding='utf-8');assert text.strip() not in s
    p.write_text(s.rstrip()+'\n\n'+text.strip()+'\n',encoding='utf-8')

def main():
    checkpoint="""### 最新批次：固定基点 regulator 非零、秩八最小性及十三维幺幂群

2026-10-13 07:05（Asia/Shanghai）。新增生成会话55106、独立复核34451均exit0；无待观察计算。原有纯运输没有重启，旧证书哈希保持不变。证明见LOWER_REGULATOR.md与MIXED_UNIPOTENT_GROUP.md。

1. 已审阅Morse/Tate约定并构造实际混合周期入口。beta_c=−s(3alpha1+8alpha2+16alpha3)/M·omega，在平衡D框架满足beta'=g_D z，g_D=(0,−3/64,−2/64,−1/64,0,0)。小消失球的beta周期收缩为零，故局部初值唯一为h g_D Q sum P_n/(n+1/2)，尾界≤2|h| ||g_DQ|| epsilon_Q。五节点首项beta/omega均与−S(3r1+2r2+r3)完全一致。

2. 新增certify_lower_regulator.py与独立verify_lower_regulator.py，只对已认证纯多项式作精确积分。行R=h g sum P_n/(n+1)，行误差≤|h| ||g|| epsilon，累计包括输入、行误差、舍入及反向干路减法。首个旧步的变量是s+64，明确乘i^n后积分。五个节点1/9/25/49/81在s=−64的omega和beta周期均认证；omega误差<1.2e−30，beta误差<7.8e−28。所有相位按实际Morse球确定，可与旧Q6球相差取向符号，不影响非零性。

3. 若下层扩张实Hodge分裂，则某个g∈C使beta+g omega在每个整球上的周期为实数，因为Q(−2)的比较常数(2pi i)^2为实。故三列(Re p,Im p,Im b)应秩≤2。节点1、9、25的三阶行列式中心8.774987805168132，严格误差<1.2e−27，精确区间(877/100,439/50)排除零；其余九子式也排除零。独立标量积分与七误差子集算法重放通过，零尾界负对照拒绝。

4. 因此指定m=(1,2,4,8)、s=−64的下层实Hodge regulator非零，进而有理Hodge/Nori扩张非零。MINIMALITY的链生成性和κ=1判据现在给出固定纤维最小秩八。显式Z_(3P1+2P2+P3)在CH²(X,1)_Q不可分解商中亦非零。不要再把该指定基点非零性或秩八最小性列为待证，也不要外推到所有代数参数。

5. data/banana_pure_ray/lower_regulator/obstruction.json SHA=f6cfda8716d7d963e5ae4563d5c164b611359b30e9b9dbe228bfb5035929c1c8；trunk_integrals.json SHA=f273bb2f0e2b01c0f1477ebe46a323c1d19ce4da70a3b21b7cea21186760f6ec；五份node_a.json及replay.json已保存。既有21迹、纯群、T0、绝对P6哈希保持不变。

6. 由新下层非零、上层Kummer log2非零及纯A=T(3)简单性，Hodge群相邻幺幂块分别为Hom(B,A)六维和Hom(C,B)一维；两表示非同构简单，所以联合像为全部七维。交换子xy'−x'y迫使中心Hom(C,A)全部六维，故幺幂根基为完整13维二步群。忠实Hodge→Nori实现夹逼给Nori群同样的全部幺幂根基。两个实际群分别是其各自实际分级群R_H、R_N在权重抛物群中的完整逆像。固定纯群R_N仍未识别，不能用射线O(Q6)替换。

7. MIXED_UNIPOTENT_GROUP给出了实际群的坐标余乘法及右motivic coaction结构公式；未声称所有比较条目已求值。当前下层只有五个球的beta条目，还需T0delta1的第六条目及边界Tate列的具体基；顶层相对比较仍未完成。

8. 下一项：审阅新实分裂障碍和13维群论证明并整合主稿；继续第六下层混合周期、完整秩八比较/顶层条目、固定纯群与四参数核对。三种饱和类仍需质量方向或边界粘合。一般图算法、物理基准、完整显式coaction/绝对群及论文目标全部保留。get_goal旧blocked状态未被重复创建或改为完成；单一续跑继续。

"""
    p=ROOT/'CONTINUATION.md';s=p.read_text(encoding='utf-8')
    s=s.replace('更新于 2026-10-13 06:45','更新于 2026-10-13 07:05',1)
    marker='### 最新批次：实际T0标记、全部节点、补全歧义与绝对P6'
    assert marker in s;s=s.replace(marker,checkpoint+'### 前一批次：实际T0标记、全部节点、补全歧义与绝对P6',1)
    start=s.index('下一批优先（以开头最新批次为准）：')
    end=s.index('已读取参考：',start)
    s=s[:start]+"""下一批优先（以开头最新批次为准）：

1. 审阅LOWER_REGULATOR的实分裂必要条件、实际球初值及误差链；固定s=−64非零和最小秩八已有完整证据，不再机械标为待证。审阅MIXED_UNIPOTENT_GROUP的13维核和实际群相对于分级群的逆像描述。
2. 继续第六个球T0delta1的beta周期及完整下层比较。可在零点局部幂零框架上积分g_D z，或复用旧纯证书构造完整七维局部规范；不能只作用纯T0并遗漏beta单值增量。顶层相对初值和秩八比较仍需另做。
3. 全四参数平坦性、固定纤维纯群R_N、三种饱和类的几何区分继续；不要把射线O(Q6)替换固定纯群。全部有限纯闭路和P6证书不重跑。
4. 主稿整合固定点非零、最小性和混合群；一般图算法及物理基准保持完整总目标的独立验收。

"""+s[end:]
    s=s.replace('指定饱和半循环、混合比较和固定regulator非零继续。','指定饱和半循环、完整混合比较和固定纯群识别继续；指定regulator非零已在新批次完成。')
    p.write_text(s,encoding='utf-8')

    edit('TASKS.md',[
      ('因子化、链生成性、扩张判据及两条射线很一般参数最小性已有证明稿；固定基点待证','因子化与链生成性已有证明；LOWER_REGULATOR认证固定s=−64实Hodge非分裂，最小秩八成立；其他参数保留原范围'),
      ('DE_RHAM_MINIMALITY 已推导很一般参数非零；s=−64 仍待证明','很一般参数和固定s=−64现分别有证明；固定点三阶行列式严格处于(8.77,8.78)，独立重放通过'),
      ('实际对象及结构约束已有；指定质量射线纯几何单值Zariski闭包=O(Q6)已证，15维见证独立核对；显式比较条目、相对混合及绝对动机群仍未完成','纯射线O(Q6)已证；固定活动对象Hodge/Nori幺幂核均为完整13维，实际群为各自分级群的逆像；固定纯群和完整混合比较仍未识别'),
      ('继续固定 s=−64 的比较值/非零性；不把很一般参数定理直接套到代数基点。','LOWER_REGULATOR已独立证明固定s=−64非零；继续第六下层和顶层比较，不外推其他代数点。')])
    append('TASKS.md',"""## 固定纤维非零与实际幺幂群之后

H5/H5b在指定m=(1,2,4,8)、s=−64的目标已通过；五个混合球周期及实分裂障碍有独立证书。M1的13维幺幂部分与结构coaction已推导，剩余是固定纯分级群R_N及完整比较条目。下一项优先第六球T0delta1的beta周期、边界Tate列及顶层相对比较；全四参数核对、饱和类标记、一般图、物理基准、全文整合仍独立继续。上方旧批次顺序为历史记录，以本节及CONTINUATION开头为准。""")

    edit('CLAIMS.md',[
      ('| s=−64 的秩八形式周期最小性 | 已等价归约到 eta_act 非零；非零性未证 | 需下层 Ext¹(Q(1),T(X)(3)) 非零，不能用上层 log 2 代替 |','| s=−64 的秩八形式周期最小性 | LOWER_REGULATOR已认证下层实Hodge非分裂，结合MINIMALITY判据得最小秩八 | 实际混合周期的三阶障碍；不由上层log2或很一般定理代替 |'),
      ('下层 regulator 比较已在 REGULATOR_COMPARISON 推导；固定基点非零性待证','比较见REGULATOR_COMPARISON；固定纤维LOWER_REGULATOR另证非零及不可分解性'),
      ('| 很一般参数的最小性证明覆盖 s=−64 | 不成立 | 可数例外集可能包含给定代数点；固定纤维须另证 |','| 很一般参数的最小性证明覆盖 s=−64 | 不成立 | 固定点现由LOWER_REGULATOR独立认证；不是很一般结论的直接特化 |'),
      ('| 块三角矩阵给出完整绝对动机群及显式 coaction | 未证明 | 指定基、完整关系、纯 K3 与质量单位比较，以及群等号的独立证据 |','| 活动对象的实际块三角群与coaction | MIXED_UNIPOTENT_GROUP证明全部13维幺幂核及相对于实际分级群的完整逆像 | 仍须识别固定纯群R_N和求值所有比较条目；不能预设纯正交群 |'),
      ('| 纯绝对P6认证推出固定regulator非零 | 不成立 | 仍需下层混合比较和有理周期商中的非零论证 |','| 纯绝对P6认证推出固定regulator非零 | 不成立 | 固定非零现由另算的beta混合周期和实分裂障碍证明 |')])
    append('CLAIMS.md',"""## 固定 regulator 与混合群新增结论

| 主张 | 状态 | 范围 |
|---|---|---|
| 固定s=−64活动下层扩张实Hodge非零 | 五球混合周期及独立精确复核；节点1/9/25的(Re p,Im p,Im b)行列式在(8.77,8.78) | LOWER_REGULATOR；不需要排除无穷多有理候选，已排除实分裂 |
| 显式Z_c在固定曲面不可分解商中非零 | 完整超越regulator非零，分解循环投影为零 | 使用已证rho16；CH²(X,1)_Q |
| 固定活动对象Hodge/Nori幺幂根基维数13 | 两个相邻非零扩张、非同构简单块及交换子给完整核 | MIXED_UNIPOTENT_GROUP；不是参数单值群 |
| 实际混合群相对于分级群的关系 | G_H(E)、G_N(E)分别是R_H、R_N在权重抛物群中的完整逆像 | 固定纯R_N本身未识别 |
| 全部秩八比较条目已求值 | 未完成 | 当前有纯P6及五个球上的beta周期；第六混合条目和顶层仍待算 |""")

    edit('MINIMALITY.md',[
      ('固定 s=−64 的非零性仍未证明。','固定 s=−64 的非零性现由 LOWER_REGULATOR.md 的实分裂障碍另证，故该点最小秩八。'),
      ('当前尚未确定这个非零性，不能写“秩八最小性已经完成”。','LOWER_REGULATOR现已用实际混合周期认证该类非零，所以本目标纤维的秩八最小性成立。'),
      ('固定代数基点的非零性仍是独立步骤。','固定代数基点的非零性由LOWER_REGULATOR独立完成。'),
      ('的一个循环代表；并未证明其 Chow 类非零。','的一个循环代表；本节的除子检验本身不证明非零，LOWER_REGULATOR现另证固定纤维的不可分解类非零。')])
    edit('REGULATOR_COMPARISON.md',[
      ('固定 s=−64 的非零性仍未证明，不能据本稿的比较宣布该基点最小性完成。','固定 s=−64 的非零性现由LOWER_REGULATOR独立认证，结合MINIMALITY得到最小秩八；该新证据不能省去。'),
      ('## 七、coaction 约束与尚未完成的非零性','## 七、coaction 约束与非零性的独立证据'),
      ('余下五列以及活动子空间中的列是否非零、是否独立，仍需证据。','余下五列整体独立性仍需证据；其中固定活动方向的非零性已由LOWER_REGULATOR补齐。'),
      ('指定代数基点仍需独立计算或非消失论证。','指定代数基点现由LOWER_REGULATOR的实际混合周期和实分裂障碍完成独立认证。')])
    for name in ['DE_RHAM_MINIMALITY.md','ABSOLUTE_COMPARISON.md','PURE_RAY_SYSTEM.md','HODGE.md']:
        p=ROOT/name;s=p.read_text(encoding='utf-8');first=s.index('\n')+1
        note='\n当前状态更新（2026-10-13 07:05）：LOWER_REGULATOR已独立认证指定m=(1,2,4,8)、s=−64的下层实Hodge非零和最小秩八；MIXED_UNIPOTENT_GROUP给出活动对象的完整13维幺幂核。下文较早“固定非零未证”的文字仅记录本稿原范围；完整混合比较和固定纯群仍待完成。\n'
        s=s[:first]+note+s[first:];p.write_text(s,encoding='utf-8')

    edit('PAPER.md',[
      ('这个判据已有证明，非零性仍保留为明确问题。证明分别见 PICARD_CERTIFICATE.md 和 MINIMALITY.md，不能将“得到判据”改写成“已证最小”。','判据见MINIMALITY，指定点的非零性现由LOWER_REGULATOR的严格实分裂障碍补齐，故最小秩八成立。'),
      ('指定 s=−64 的非消失仍不能据此宣布完成。','指定s=−64的非消失另由LOWER_REGULATOR认证，未使用很一般结论的直接特化。'),
      ('Nonvanishing at prescribed algebraic fibres, mixed comparison entries, and identification of the full motivic group remain separate questions.','At masses (1,2,4,8) and s=-64, certified mixed periods give a nonzero real Hodge splitting obstruction, proving minimal rank eight at this prescribed algebraic fibre. The two nonzero adjacent extensions force the full thirteen-dimensional unipotent radical, both in the Hodge and Nori settings. Complete mixed comparison entries and the fixed pure motivic quotient remain to be determined.'),
      ('指定 s=−64 的最小性及完整群等号仍待证明；','指定s=−64最小性及混合群相对于实际分级群的完整逆像现已证明；固定纯群仍待识别；'),
      ('特定补全半循环仍待标记，绝对比较矩阵另需Q5。','特定补全半循环仍待标记；Q5已完成明确有理框架的六循环绝对比较。')])
    append('PAPER.md',r"""## 固定纤维非消失和活动混合群定理

**定理（指定代数纤维）。** 对二维四传播子banana、m=(1,2,4,8)、s=−64，活动方向c=3P₁+2P₂+P₃的下层扩张在实混合Hodge范畴中非零。显式高阶Chow循环Z_c的超越regulator非零，物理形式周期的最小Hodge及Nori对象秩均为八。

证明从实际提升beta_c入手：消失球收缩令其周期零初值唯一；几何质量导数恒等式给beta'=g_D z。对既有纯运输多项式作精确积分并含输入/尾界/尺度误差，得到五个实际球的混合周期。实Hodge分裂会迫使某个复g使所有beta+g omega周期为实数。节点1、9、25的矩阵(Re p,Im p,Im b)行列式严格位于(8.77,8.78)，矛盾。此判据与有理周期商相容，未把一个非零复代表直接等同非零扩张。详细证明和独立证书见LOWER_REGULATOR。

**定理（实际幺幂群）。** 该活动对象的Hodge群和Nori群在权重分级上的核均为完整十三维二步幺幂群，维数来自相邻六维K3–Tate块、一维Kummer块和六维中心。两个相邻扩张非零，相关表示非同构简单；交换子迫使中心全部出现。两群分别是其实际分级群R_H、R_N在权重抛物群中的完整逆像。MIXED_UNIPOTENT_GROUP给出证明和明确Hopf/coaction坐标式。

本文仍不识别固定纯群R_N，也不声称全部秩八比较条目已经求值。当前五个beta球周期还需第六混合条目和顶层相对比较才能形成完整比较矩阵；指定饱和类、四参数核对及原有总目标其余部分继续。""")

    append('README.md',"""## 当前结果：固定基点最小性与实际混合群

[LOWER_REGULATOR](LOWER_REGULATOR.md)已完成指定m=(1,2,4,8)、s=−64的下层实Hodge非零和最小秩八证明：三个实际球的混合周期给出严格处于(8.77,8.78)的分裂障碍行列式，独立复核通过。[MIXED_UNIPOTENT_GROUP](MIXED_UNIPOTENT_GROUP.md)进一步确定Hodge/Nori群的完整13维幺幂根基，并给出相对于实际纯分级群的群等式及coaction结构。上方旧进度中的固定非零缺口已由本批补齐。

下一步是第六下层比较条目、顶层相对比较、固定纯群R_N和四参数核对；饱和类标记、一般图算法、物理基准及主稿仍保留。""")
    append('SOURCE_AUDIT.md',"""## 固定点 regulator 证据与三步扩张群

2026-10-13 07:05：新增LOWER_REGULATOR、MIXED_UNIPOTENT_GROUP。固定s=−64非零由本项目实际混合球周期严格包围证明，不使用Sato特殊K3的非消失定理，也不使用很一般参数特化。Sato arXiv:2408.09102 §2.2–2.4本轮核对了原始HTML的regulator商描述。

Eskandari arXiv:2407.01379v3（2026-08-05）原始HTML Theorem A给出三步扩张幺幂核的一般框架；本稿以非同构简单相邻表示和直接矩阵交换子给出指定对象的13维结论。固定纯群仍未知，未把纯几何单值O(Q6)替换为它。""")

    p=ROOT/'data/input_manifest.json';d=json.loads(p.read_text(encoding='utf-8'))
    d.update(status='fixed_fibre_rank_eight_minimality_and_full_mixed_unipotent_kernel_proved',
      comparison_status='Fixed s=-64 real Hodge lower extension nonzero by certified mixed periods; complete mixed comparison and fixed pure group remain open.',
      rank_eight_minimality_certified=True,
      fixed_lower_real_Hodge_regulator_nonzero_certified=True,
      fixed_regulator_proof='../LOWER_REGULATOR.md',
      fixed_regulator_certificate='banana_pure_ray/lower_regulator/obstruction.json',
      fixed_regulator_replay='banana_pure_ray/lower_regulator/replay.json',
      fixed_regulator_determinant_interval=['877/100','439/50'],
      fixed_regulator_primary_nodes=[1,9,25],
      fixed_explicit_Chow_cycle_indecomposable_over_Q=True,
      certified_lower_mixed_period_nodes=[1,9,25,49,81],
      full_lower_seven_dimensional_comparison_certified=False,
      full_active_eight_dimensional_comparison_certified=False,
      fixed_Hodge_unipotent_radical_dimension=13,
      fixed_Nori_unipotent_radical_dimension=13,
      full_mixed_unipotent_kernel_proof='../MIXED_UNIPOTENT_GROUP.md',
      actual_mixed_group_relative_to_graded_group='Full inverse image of the actual graded group; fixed graded group itself remains unidentified.',
      absolute_motivic_group_identified=False)
    p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('Proofs, tasks, claims, main manuscript and continuation reconciled.')

if __name__=='__main__':main()
