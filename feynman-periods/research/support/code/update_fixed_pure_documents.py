"""Integrate fixed-Q pure group proof after independent arithmetic replay."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import json,hashlib,re
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/fixed_pure_group'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rd(p):return (ROOT/p).read_text(encoding='utf-8')
def wr(p,s):(ROOT/p).write_text(s,encoding='utf-8')
def append(p,s):
    old=rd(p)
    if s.strip() not in old:wr(p,old.rstrip()+'\n\n'+s.strip()+'\n')
def main():
    c=json.loads(rd('data/fixed_pure_group/certificate.json'))
    v=json.loads(rd('data/fixed_pure_group/replay.json'))
    assert v['status']=='fixed_pure_arithmetic_independently_replayed'
    assert v['certificate_sha256']==sha(OUT/'certificate.json')
    assert v['verifier_sha256']==sha(ROOT/'code/verify_fixed_pure_group.py')
    assert c['code_sha256']==sha(ROOT/'code/certify_fixed_pure_group.py')
    prior=json.loads(rd('data/banana_mixed_parameter/final_audit.json'))
    preserved={}
    for field,base in [('old_ray_hashes_preserved','data/banana_pure_ray'),
                       ('old_four_parameter_hashes_preserved','data/banana_four_parameter'),
                       ('old_proof_hashes_preserved','')]:
        for name,h in prior[field].items():
            p=Path(base)/name
            assert sha(ROOT/p)==h
            preserved[p.as_posix()]=h
    for name,h in prior['matrix_sha256'].items():
        p='data/banana_mixed_parameter/'+name
        assert sha(ROOT/p)==h;preserved[p]=h
    for name,h in [('MIXED_PARAMETER_CONNECTION.md',prior['proof_sha256']),
                   ('data/banana_mixed_parameter/flatness.json',prior['flatness_sha256']),
                   ('data/banana_mixed_parameter/replay.json',prior['replay_sha256'])]:
        assert sha(ROOT/name)==h;preserved[name]=h
    sources=json.loads(rd('external/fixed-pure-group/manifest.json'))
    for source in sources:
        assert sha(ROOT/'external/fixed-pure-group'/source['name'])==source['sha256']
    proof=rd('FIXED_PURE_GROUP.md')
    assert all(ch not in proof for ch in ['\x08','\x0c'])
    stamp=datetime.now(timezone.utc).astimezone(timezone(timedelta(hours=8))).strftime('%Y-%m-%d %H:%M')
    checkpoint=f"""### 最新批次：指定Q纤维纯群与完整活动群

{stamp}（Asia/Shanghai）。FIXED_PURE_GROUP.md给出固定m=(1,2,4,8)、s=−64、基域Q的证明。生成器及独立整数/Fraction验证器均已命令内exit0，没有待观察计算；初次生成仅在JSON整型序列化处失败，修正后完成，未重复点数或运输。全部旧P6/P7/P8、纯/混合四参数证书及前批证明哈希保持。

1. 已审阅MIXED_PARAMETER_CONNECTION的实际一动机拉回、常数留数、IBP和三张量唯一性，原600条目证书未改。本批转向固定纯群，不扩张比较矩阵规模。
2. 旧p=11残余P=Z6+10Z5+55Z4−484Z3+6655Z2+146410Z+1771561模5不可约、模11有唯一p-adic单位根。不可约性加唯一单位根证明任意Frobenius幂仍有简单谱。rho_p=16与Charles排除全实二次乘法；全实三次情形的二维F空间会强迫额外CM自同态，也排除。
3. CM候选由Elsenhans–Jahnel的正规性定理（仅在CM已知Hodge猜想的条件下应用）及q判别式平方类3，必含Q(sqrt(−3))。模29简单特征根给一维Q29特征线，而−3在Q29非平方；绝对Hodge实现使某个Frobenius幂与CM域交换，矛盾。故End_Hdg T=Q。原PDF的Charles第3/6/8页和EJ第21/22页已核对。
4. S=T(1)的Hodge群为SO(q)。已有绝对比较给det(C_S)^2=9/(85²·11930545734185)，后者无平方因子核非平凡；若detS是Q上的Nori单位对象，比较系数必为有理数，故排除。结合Hodge下界和杯积上界，G_N,Q(S)=O(q)。这一步不使用周期映射单射猜想，也不以射线几何群代替固定群。
5. A=T(3)、B=Q(1)的实际分级群分别为SO(q)×Gm、O(q)×Gm，作用(t²g,t)。结合旧完整U13，指定秩八活动对象的两个实际群维数均29，Hodge连通、Nori两个分支。deRham纯块关系用q_dR=−2K_B^-1，H^t q_dR H=t⁴q_dR；实际八项物理coaction不变，其纯关系现已明确。不要将此结果无条件推广其他纤维、扩大基域或秩十环境。
6. 算术certificate SHA {sha(OUT/'certificate.json')}；replay SHA {sha(OUT/'replay.json')}。独立算法不用SymPy或生成器，重新验证模多项式、Newton迹、判别式和比较行列式，无新点数。生成文件中的independent_replay_completed=false仅为生成时记录，终态读取replay及final_audit。
7. 下一项优先实际饱和半循环（三种候选仍未区分）、G1有界一般图算法、F1独立物理基准及P1/P2全文整合。秩十绝对比较可择时推进，需三个真实边界列；现有P8完整。指定Q纤维固定纯群不再列为缺失；全目标仍未完成，单一heartbeat保留。

"""
    s=rd('CONTINUATION.md')
    marker='### 最新批次：实际秩十多参数混合对象与完整平坦连接'
    if '### 最新批次：指定Q纤维纯群与完整活动群' not in s:
        assert marker in s
        s=re.sub(r'更新于 .*?（Asia/Shanghai）',f'更新于 {stamp}（Asia/Shanghai）',s,count=1)
        s=s.replace(marker,checkpoint+'### 前一批次：实际秩十多参数混合对象与完整平坦连接',1)
        a=s.index('下一批优先（以开头最新批次为准）：');b=s.index('已读取参考：',a)
        s=s[:a]+"""下一批优先（以开头最新批次为准）：

1. 审阅FIXED_PURE_GROUP的实乘/CM排除与比较行列式论证，特别区分Q上的Nori群与扩大基域后的群。指定纤维End_Hdg T=Q、SO/O及完整29维活动群已有证明与独立算术见证；不要重新计算点数、纯运输或旧P8。
2. H4e：三种偶指数二饱和类须实际质量方向闭路或边界粘合区分。群的Zariski闭包和更高精度不选出实际半循环。四个旧质量矩阵仍只是候选。
3. G1：设计明确有限输入、输出和适用范围的有界一般图算法；先审计原稿组合判据及复杂度含义，不能宣称自动分类任意图动机。F1：执行另外定义的物理基准，banana单点读出不代替它。无规范对称性的标量phi4无需虚构Ward恒等式。
4. P1/P2：把已有定理、证明和证书改写成连贯全文，清理历史段落的过时“待证”及来源归属；一般可见性过滤若无嵌套对象证明仍须条件表述。
5. 秩十环境的绝对比较仍可复用旧纯多项式积分三个beta_i；先构造实际三个边界环面、I−ZR校正及实填充，第六球计入混合T0增量。不能P8补零，亦不把该目标取代一般图和独立基准。
6. 所有本批计算均已结束，旧哈希保留。无本项目待观察进程；不创建重复聊天或自动任务。get_goal仍显示旧blocked且工具不能恢复，用户续跑授权有效，已有实质推进，不新建目标。

"""+s[b:]
        wr('CONTINUATION.md',s)
    s=rd('TASKS.md')
    s=re.sub(r'^\| M1 \|.*$', '| M1 | 显式 coaction 与实际群比较 | H3、H4a 至 H4e、H5、Q5 | 基、扩张条目、关系和群等号分别给出证明 | 指定Q纤维P8与实际coaction已认证；FIXED_PURE_GROUP识别End_Hdg=Q，分级群SO×Gm/O×Gm，结合U13得到完整29维活动群；其他纤维与秩十环境群不外推 |',s,flags=re.M)
    a=s.index('## 当前理论批次');b=s.index('## 决策门槛',a)
    s=s[:a]+"""## 当前理论批次

1. M1：FIXED_PURE_GROUP已关闭指定m=(1,2,4,8)、s=−64、基域Q的固定纯群及完整秩八活动群关系缺口。先审阅新的算术与文献定理条件，旧P8和实际coaction保持。所有已求条目为严格包围，不声称闭合特殊函数表达式或周期代数独立。
2. H4e：实际六循环Q6与指数二已证；三个饱和类仍须质量方向或边界粘合几何区分。
3. G1/F1：执行有界一般图算法及独立物理基准。必须明确输入、输出、范围和可重放验收，不能以banana单点积分替代。
4. P1/P2：从现有证明稿整合连贯全文并作来源/新颖性审计。下方小节保留历史，以本节及CONTINUATION最新批次为准。
5. H5e实际秩十连接已完成；其P10若推进须另证三个边界Betti列与实填充。指定纤维最小活动对象仍是秩八。

"""+s[b:];wr('TASKS.md',s)
    append('CLAIMS.md',"""## 指定Q纤维的固定纯群与完整活动群（覆盖前述旧“群未知”状态）

| 主张 | 证据 | 适用范围 |
|---|---|---|
| End_Hdg T=Q | Charles排除实乘；CM正规性与detq平方类3；p11幂谱分离及Q29一维特征线排除sqrt(−3) | m=(1,2,4,8)、s=−64，不用一般参数特化 |
| G_H(T(1))=SO(q) | Zarhin/Charles连通K3群描述 | 固定有理Hodge结构 |
| G_N,Q(T(1))=O(q) | Hodge下界、杯积上界、实际比较行列式平方9/86198192929486625非平方 | 基域Q；不宣称扩域后分支不变 |
| R_H=SO(q)×Gm，R_N=O(q)×Gm | S与Tate的群夹逼，SO无非平凡交换商 | 在A=T(3)、B=Q(1)上的作用为(t²g,t) |
| 固定活动群维数29 | 已证明U13及完整逆像，加本批纯群 | Hodge连通；Nori两个分支；非秩十环境群 |
| 物理coaction的纯关系已明确 | q_dR=−2K_B^-1，H^t q_dR H=t⁴q_dR | 采用deRham杯积，不能替换为未换基Q6；不增加detH=t¹²约束 |
| 计算与理论分工 | 独立整数/Fraction算术重放；外部定理逐条核对条件 | 完整证明见FIXED_PURE_GROUP；不是仅凭浮点周期猜群 |
""")
    # Replace live summaries, while keeping earlier proof artifacts immutable.
    s=rd('PAPER.md')
    s=s.replace('最小性与绝对群识别仍是独立问题。','指定纤维最小性和基域Q上的绝对活动群现分别由LOWER_REGULATOR与FIXED_PURE_GROUP证明；其他纤维仍须独立论证。')
    s=s.replace('The fixed pure motivic quotient, marked integral saturation and full rank-ten absolute comparison remain separate problems.',
      'At the prescribed fibre over Q, an arithmetic exclusion of real and complex multiplication proves End_Hdg(T)=Q. The normalized pure Hodge group is SO(q); an irrational determinant comparison period proves that the corresponding Nori group over Q is O(q). Together with the full unipotent radical this identifies both rank-eight active groups of dimension twenty-nine and all their coaction relations. Marked integral saturation, the broader graph algorithm, an independent physical benchmark and the optional full rank-ten comparison remain separate tasks.')
    s=s.replace('The remaining top relative column of the rank-eight object, its compatible framings, and the fixed pure motivic quotient require separate work.',
      'The top relative comparison is supplied in the following theorem; FIXED_PURE_GROUP subsequently identifies the fixed pure groups over Q.')
    s=s.replace('The fixed pure quotient of the motivic group is retained as R_N; its relations have not been replaced by those of the geometric ray monodromy group.',
      'The fixed graded quotient R_N is identified over Q in FIXED_PURE_GROUP by arithmetic and determinant comparison, independently of the geometric ray group.')
    wr('PAPER.md',s)
    append('PAPER.md',r"""## 固定纤维的纯群与完整活动群定理

**Theorem.** Let \(X=X_{(1,2,4,8),-64}/\mathbb Q\), \(S=T(X)(1)\), and \(E\) be the rank-eight active object constructed above. Then
\[
 \operatorname{End}_{\mathrm{Hdg}}T=\mathbb Q,\qquad
 G_H(S)=SO(q),\qquad G_{\mathrm{N},\mathbb Q}(S)=O(q).
\]
With \(A=T(3)\) and \(B=\mathbb Q(1)\), their graded groups are \(SO(q)\times\mathbb G_m\) and \(O(q)\times\mathbb G_m\), acting by \((t^2g,t)\). The full active groups are their complete inverse images in the weight parabolic, with unipotent radical of dimension thirteen. Both groups have dimension twenty-nine; the Hodge group is connected and the Nori group over \(\mathbb Q\) has two geometric components.

**Proof.** The already certified reduction at 11 has residual polynomial
\[
 P(Z)=Z^6+10Z^5+55Z^4-484Z^3+6655Z^2+146410Z+1771561.
\]
It is irreducible modulo 5 and has exactly one 11-adic unit root. Transitivity therefore rules out root-of-unity ratios between distinct roots, so every Frobenius power has simple spectrum. Charles's Picard-jump bound excludes real quadratic multiplication. A two-dimensional space over a totally real endomorphism field would admit an additional imaginary Hodge endomorphism from the binary polarization, excluding real cubic multiplication as well.

In the CM case, the known Hodge conjecture for the self-product and Elsenhans–Jahnel's normality theorem force the endomorphism field over \(\mathbb Q\) to be Galois. Its degree is two or six. It consequently contains an imaginary quadratic field, whose trace-Hermitian determinant and the square class \(\det q=3\) force \(\mathbb Q(\sqrt{-3})\). But \(P\) has a simple root over \(\mathbb Q_{29}\), while \(-3\) is not a square there. Absolute Hodge realization makes a finite Frobenius power commute with this quadratic field, contradicting its action on the one-dimensional eigenline. This proves the endomorphism assertion; the K3 Hodge group is then \(SO(q)\).

The actual cup product bounds the Nori group by \(O(q)\). Its determinant object has comparison coefficient
\[
 c_D=\det(P_6)/(2\pi i)^6,\qquad
 c_D^2=\frac9{85^2\cdot11930545734185},
\]
where \(11930545734185=5\cdot13\cdot29\cdot37\cdot73\cdot89\cdot113\cdot233\) is squarefree. A unit object over \(\mathbb Q\) has rational comparison in rational bases, so this determinant object is nontrivial. Together with the Hodge lower bound this proves the full orthogonal group. The Tate factor is independent because \(SO(q)\) has no nontrivial commutative quotient. Finally apply the already established full-unipotent-kernel theorem.

In the explicit de Rham frame, the pure block relation is \(H^{\mathsf t}q_{\rm dR}H=t^4q_{\rm dR}\), with \(q_{\rm dR}=-2K_B^{-1}\). Only the Hodge subgroup imposes \(\det H=t^{12}\). These relations complete the actual Hopf description underlying the eight-term physical coaction. This is a statement over the specified base field and at the specified fibre, not an identification after arbitrary field extension or an algebraic-independence statement about numerical periods. Full details, source conditions and independent arithmetic certificates are in FIXED_PURE_GROUP.""")
    append('README.md',"""## 当前结果：指定Q纤维的实际群已识别

[FIXED_PURE_GROUP](FIXED_PURE_GROUP.md)以旧p=11证书、极化判别式及29-adic特征线证明End_Hdg T=Q；绝对比较行列式区分SO与O。指定m=(1,2,4,8)、s=−64的Hodge与Q上Nori活动群均为明确29维群，后者有两个分支，完整coaction关系已给出。独立算术证书在data/fixed_pure_group；没有重算点数或运输。此前“固定纯群未知”的条目在该点及基域Q的范围内由本稿更新。下一重点为饱和格标记、有界一般图算法、独立物理基准及连贯论文。""")
    append('SOURCE_AUDIT.md',"""## 固定纯群的外部定理与新算术应用

FIXED_PURE_GROUP使用Charles《On the Picard number of K3 surfaces over number fields》第3、6、8页：全实奇维Picard跳跃、K3端同态域、绝对Hodge的有限指数交换性、归一化SO恒等分支。特征p=11满足定理条件。EJ《On the component group of the algebraic monodromy group of a K3 surface》Theorem7.5及7.7证明（第21、22页）给端同态域正规性，并在CM情况下援引已知自积Hodge猜想。此处只在CM反证分支调用，未无条件假设一般Hodge猜想。原PDF及已核对页图保存于external/fixed-pure-group，manifest含源URL与哈希。全实二维端同态空间排除、模29特征线矛盾、比较行列式分支判定在本稿展开。计算见证独立重放，旧点数和绝对杯积为明确输入。群同构限定于指定Q纤维，不从参数几何单值偷换，也不声称未经检索的文献首创性。""")
    m=json.loads(rd('data/input_manifest.json'))
    m.update({
      'comparison_status':'Specified-fibre P8, physical coaction and fixed-Q pure/active groups identified; actual rank-ten family mixed connection certified. Marked saturation, bounded general-graph algorithm, independent physical benchmark and manuscript remain unfinished.',
      'fixed_pure_group_identified':True,'fixed_pure_group_scope':{'base_field':'Q','masses':[1,2,4,8],'s':-64,'normalization':'S=T(1)'},
      'fixed_hodge_endomorphism_field':'Q','fixed_pure_hodge_group':'SO(q)','fixed_pure_nori_group_over_Q':'O(q)',
      'fixed_graded_hodge_group':'SO(q) x Gm','fixed_graded_nori_group_over_Q':'O(q) x Gm',
      'actual_mixed_group_relative_to_graded_group':'Full inverse image; graded groups identified for specified Q-fibre in FIXED_PURE_GROUP.',
      'fixed_active_group_dimension':29,'fixed_active_nori_geometric_components':2,
      'fixed_pure_group_proof':'../FIXED_PURE_GROUP.md','fixed_pure_group_certificate':'fixed_pure_group/certificate.json',
      'fixed_pure_group_replay':'fixed_pure_group/replay.json','rank_ten_group_identified':False})
    wr('data/input_manifest.json',json.dumps(m,ensure_ascii=False,indent=2)+'\n')
    audit={'status':'fixed_Q_fibre_pure_and_active_group_batch_complete','timestamp':stamp,
      'proof_sha256':sha(ROOT/'FIXED_PURE_GROUP.md'),'certificate_sha256':sha(OUT/'certificate.json'),
      'replay_sha256':sha(OUT/'replay.json'),'generator_sha256':c['code_sha256'],'verifier_sha256':v['verifier_sha256'],
      'source_pdf_hashes':{s['name']:s['sha256'] for s in sources},
      'preserved_old_hashes':preserved,'new_point_counts':0,'new_pure_transports':0,
      'all_started_calculations_ended':True,'base_field':'Q','other_fibres_or_field_extensions_identified':False,
      'full_goal_complete':False}
    wr('data/fixed_pure_group/final_audit.json',json.dumps(audit,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in audit.items() if k not in ['preserved_old_hashes','source_pdf_hashes']},ensure_ascii=False))
if __name__=='__main__':main()

