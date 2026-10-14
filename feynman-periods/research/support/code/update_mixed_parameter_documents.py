"""Save the actual rank-ten mixed connection checkpoint after independent replay."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import json,hashlib,re
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/banana_mixed_parameter'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rd(name):return (ROOT/name).read_text(encoding='utf-8')
def wr(name,s):(ROOT/name).write_text(s,encoding='utf-8')
def append(name,s):
    text=rd(name)
    if s.strip() not in text:wr(name,text.rstrip()+'\n\n'+s.strip()+'\n')
def main():
    replay=json.loads((OUT/'replay.json').read_text())
    assert replay['curvature_entries']==600 and replay['direct_torus_IBP_identities']==16
    assert replay['verifier_sha256']==sha(ROOT/'code/verify_banana_mixed_parameter.py')
    assert replay['flatness_sha256']==sha(OUT/'flatness.json')
    assert replay['input_pure_replay_sha256']==sha(ROOT/'data/banana_four_parameter/replay.json')
    prior=json.loads(rd('data/banana_four_parameter/final_audit.json'))
    old=dict(prior['old_hashes_preserved'])
    for name,h in old.items():assert sha(ROOT/'data/banana_pure_ray'/name)==h
    old_proofs=dict(prior['new_proofs_sha256'])
    old_proofs['RELATIVE_TOP_COMPARISON.md']=prior['old_relative_top_proof_sha256']
    for name,h in old_proofs.items():assert sha(ROOT/name)==h
    old_pure={}
    flat=json.loads(rd('data/banana_four_parameter/flatness.json'))
    for name in ['A1.json','A2.json','A3.json','A4.json','K.json','flatness.json','replay.json','determinant.json','final_audit.json']:
        old_pure[name]=sha(ROOT/'data/banana_four_parameter'/name)
    for name,h in replay['matrix_sha256'].items():assert sha(OUT/name)==h
    proof=rd('MIXED_PARAMETER_CONNECTION.md')
    assert '\x08' not in proof and '\x0c' not in proof
    for name,h in replay['matrix_sha256'].items():assert h in proof
    stamp=datetime.now(timezone.utc).astimezone(timezone(timedelta(hours=8))).strftime('%Y-%m-%d %H:%M')
    checkpoint=f"""### 最新批次：实际秩十多参数混合对象与完整平坦连接

{stamp}（Asia/Shanghai）。MIXED_PARAMETER_CONNECTION.md为本批证明，证书在data/banana_mixed_parameter。重建命令内exit0，全曲率生成44030与独立复核45359均exit0，没有待观察计算。没有重启旧纯运输；P6/P7/P8、原顶层证明和全四参数纯证书哈希全部保持。

1. 从实际边界一动机前三个质量坐标拉回，得到E_mass秩十对象，未扭转权重分级6、3、1、Hodge维数1、5、9、10。环境商为五个分裂Tate方向；正欧氏物理链因实共轭进入此子对象。固定m=(1,2,4,8)纤维仍以原秩八活动对象为最小，不能把秩十说成该点最小。

2. 取x_i=mu_i/s、D_i=x_i partial_i、Omega=s omega、V_i=D_iOmega、beta_i=−x_i alpha_i N Omega。实际相对Cartan公式给D_i e=beta_i，sum beta_i=−Omega。边界留数常数使C_ij=D_j beta_i落在纯F1，C_ij=C_ji且sum_jC_ij=−V_i。

3. 十六个环面分部积分恒等式精确给C_ii=−x_i(Omega+sumV)。其余C_ij=aOmega+sum bV的系数全部是仿射多项式，完整式见证明§4。Hodge横截性及纯配对给完全对称三张量，四个未知全异指标项的系数行列式为−2，继而平坦性唯一固定a。由此证明这些矩阵就是实际质量提升的连接，并非只验证任意形式扩张平坦。

4. 框架(L6,beta1,beta2,beta3,e)的四个10×10有理矩阵已保存；六百曲率条目清分母逐系数为零。独立算法从紧凑公式重建，核对16个IBP、96个耦合、144个纯换基条目和全部600曲率，负对照拒绝。flatness SHA {sha(OUT/'flatness.json')}；replay SHA {sha(OUT/'replay.json')}。

5. 顶层nabla e=−Omega dlogx4+sum beta_i dlog(x_i/x4)，保持P/−P及nu=e_F−e_B的旧Kummer符号。三个单位函数域赋值矩阵为I3。沿固定质量射线，beta_i'=V_i/s及e'=Omega/s精确恢复旧系统；beta_c=3beta1+2beta2+beta3给既有秩八稳定框架。

6. 下一项先审阅本稿实际拉回与唯一性证明。可直接复用已保存纯多项式分别积分三个beta_i，构造秩九下层/秩十相对比较，但须另证三个边界Betti列及填充；不能将P8补零。质量方向/边界粘合的三种饱和类区分、固定纯R_H/R_N、有界一般图算法、独立物理基准和全文整合仍保留。完整多参数混合连接不再缺失；秩十绝对比较及其全群尚未计算。单一heartbeat继续，整体目标未完成。

"""
    s=rd('CONTINUATION.md')
    if '### 最新批次：实际秩十多参数混合对象' not in s:
        s=re.sub(r'更新于 .*?（Asia/Shanghai）',f'更新于 {stamp}（Asia/Shanghai）',s,count=1)
        marker='### 最新批次：全四参数纯平坦性、绝对杯积延伸与几何群';assert marker in s
        s=s.replace(marker,checkpoint+'### 前一批次：全四参数纯平坦性、绝对杯积延伸与几何群',1)
        start=s.index('下一批优先（以开头最新批次为准）：');end=s.index('已读取参考：',start)
        s=s[:start]+"""下一批优先（以开头最新批次为准）：

1. 审阅MIXED_PARAMETER_CONNECTION中的实际一动机拉回、Hodge提升、环面IBP和对称三张量唯一性；实际秩十多参数连接已通过600条目独立精确复核，勿重复提取或重建。旧证书和纯几何归一化保持不变。
2. 如推进秩十绝对比较，可把旧纯多项式的积分行拆为三个beta_i'=V_i/s。先明确P1/P2/P3对偶边界环面，使用Tate_basis_left_inverse及I−ZR排除环境H2(T)，再作实共轭填充，顶层Kummer项应为(2pi i)^2 log(mu_i/mu4)。第六球仍须完整混合T0增量。复用旧q周期，不重复纯ODE。现有P8本身完整，不能误称缺失。
3. 三种偶指数二补全须实际质量方向闭路或边界粘合区分。全族O(Q6)、连接平坦性和更高射线精度不能选出实际半循环。四个旧质量矩阵仍为候选。
4. 固定纯R_H/R_N、G1有界一般图算法、F1另设物理基准和P2全文继续。可以适时转向这些尚未执行的总目标，避免只扩张矩阵规模。当前十维系统不等于所有纤维最小，亦未证明十维混合群。
5. runtime/python-flint为工作区局部依赖，不改旧环境。reconstruction/flatness的false字段仅记录生成时状态，终态读取replay和final_audit。所有本批会话已exit0，无待观察项目进程。

"""+s[end:]
        wr('CONTINUATION.md',s)
    s=rd('TASKS.md')
    if '| H5e |' not in s:
        marker='| Q0 |'
        s=s.replace(marker,'| H5e | 全四参数实际混合连接 | H4b、H5c、H5d | 三个Kummer方向的实际对象、完整有理矩阵和全平坦性 | MIXED_PARAMETER_CONNECTION完成秩十拉回及仿射混合耦合，600条目独立精确复核；秩十绝对比较另待证 |\n'+marker,1)
    a=s.index('## 当前理论批次');b=s.index('## 决策门槛',a)
    s=s[:a]+"""## 当前理论批次

1. H5e已完成：实际秩十质量子对象、四个混合矩阵和600个曲率条目见MIXED_PARAMETER_CONNECTION。先审阅实际拉回与唯一性证明；后续比较可复用旧纯多项式，不重复原运输。
2. H4e：完整抽象T与实际Q6指数二已经确定。质量方向闭路或边界粘合仍须区分三种饱和类；旧质量矩阵只是候选。
3. M1：指定秩八对象的13维幺幂核、完整P8和实际coaction框架已证；固定纯R_H/R_N仍待独立识别。秩十多参数连接不自动识别这些群。
4. G1/F1：执行有界一般图分类与范围明确的独立物理基准。banana指定积分不能替代全部基准。
5. P1/P2：将定理、来源与证书整理成连贯全文；下方各批次小节保留历史，以本节及CONTINUATION开头为准。

"""+s[b:];wr('TASKS.md',s)
    append('TASKS.md',"""## 秩十实际混合连接后的任务边界

H5e通过：质量一动机子对象拉回给秩十环境，完整混合连接与Hodge过滤已明确；全部600曲率由紧凑公式独立重放。旧射线上秩八活动子对象及其P8没有改变。下一比较任务若执行，是三个独立beta的秩九下层和秩十环境P10，并非再次计算旧P8。泛单位秩三不等于指定纤维最小秩十，三种整格补全和固定纯群仍另待证。""")
    append('CLAIMS.md',"""## 实际四参数混合连接

| 主张 | 证据 | 范围 |
|---|---|---|
| 全质量秩十实际子对象存在 | 实际边界一动机前三环面坐标拉回，商Q(1)^5 | 物理扭转后；未扭转权重6、3、1，F维数1、5、9、10 |
| 三个质量提升的完整连接 | 16个环面IBP、留数常数、横截性、对称三张量det−2唯一性；仿射系数明确 | 实际几何类的连接，不只形式矩阵 |
| 全四参数混合平坦性 | 四个10×10矩阵全部600曲率条目清分母恒零，独立紧凑公式重放 | Q(x1,x2,x3,x4)，不是抽样 |
| 顶层Kummer连接与旧符号 | nabla e=−Omega dlogx4+sum beta_i dlog(xi/x4)；保持P/−P及nu=e_F−e_B | 三个泛单位独立，不改变指定纤维log2 |
| 秩十系统与旧秩八系统关系 | 固定质量beta_i'=Vi/s、e'=Omega/s；K,beta_c,e稳定 | 原P8完整且哈希未变；不称指定纤维最小秩十 |
| 秩十绝对P10及其全混合群 | 尚未计算或识别 | 不能由平坦性或旧P8补零推出 |
""")
    s=rd('PAPER.md')
    s=s.replace('The fixed pure motivic quotient, marked integral saturation and actual multiparameter mixed connection remain separate problems.',
      'Pullback along the three mass-ratio coordinates constructs an actual rank-ten relative object over the parameter family. Its four mixed connection matrices have explicit affine lower couplings, uniquely determined by torus integration by parts, Hodge transversality and the pure pairing; all six hundred curvature entries are independently verified. The special rank-eight realization remains its active subobject on the prescribed mass ray. The fixed pure motivic quotient, marked integral saturation and full rank-ten absolute comparison remain separate problems.')
    wr('PAPER.md',s)
    append('PAPER.md',r"""## 可纳入正文的实际混合连接定理

**Theorem.** Over the labelled smooth banana parameter family, the pullback of the three-coordinate mass Kummer one-motive defines an actual relative object of rank ten. Before the physical Tate twist its weight-graded ranks are (6,3,1), and its Hodge filtration dimensions are (1,5,9,10). Set \(x_i=\mu_i/s\), \(D_i=x_i\partial_i\), \(\Omega=s\omega\), \(V_i=D_i\Omega\), and \(\beta_i=-x_i\alpha_iN\Omega\). In the frame \((\Omega,V_1,\ldots,V_4,\partial_1^2\Omega,\beta_1,\beta_2,\beta_3,e)\), the complete connection consists of the previously certified pure block and explicit affine coefficients \(C_{ij}=D_j\beta_i\) given in MIXED_PARAMETER_CONNECTION, equations (8)–(9). Its top row is
\[
 \nabla e=-\Omega\,d\log x_4+\sum_{i=1}^3\beta_i\,d\log(x_i/x_4).
\]
All six integrability identities hold over the full rational function field.

**Proof.** Pull back the actual boundary quotient along its first three torus coordinates. The quotient by this rank-ten subobject is split Tate, and the actual F³ class and its mass derivatives lie in the pullback. Their constant boundary residues and Griffiths transversality put \(C_{ij}\) in the pure F¹ space. Cartan compatibility makes \(C_{ij}\) symmetric; torus integration by parts gives \(C_{ii}=-x_i(\Omega+\sum_r V_r)\). Pairing the remaining flatness equations with \(\Omega\) produces a completely symmetric cubic tensor. Its four undetermined all-distinct entries are uniquely fixed by a matrix of determinant −2. The pure connection then uniquely determines the holomorphic coefficients. The resulting affine formulas and the complete 10×10 matrices are verified independently through six hundred coefficientwise curvature identities.

On the ray \(m=(1,2,4,8)\), the subframe \(K,\beta_c,e\) with \(\beta_c=3\beta_1+2\beta_2+\beta_3\) recovers the already constructed rank-eight active object. The new family object neither changes that fibre's minimal rank nor provides its rank-ten absolute comparison by formal padding. The pure block is attributed to the published unequal-mass system; the present proof supplies the relative mixed identification and explicit lower couplings.""")
    append('README.md',"""## 当前结果：实际四参数混合系统

[MIXED_PARAMETER_CONNECTION](MIXED_PARAMETER_CONNECTION.md)构造实际秩十质量子对象，给出三个质量提升的仿射连接公式和完整Hodge框架。四个10×10矩阵的600个曲率条目、16个环面分部积分恒等式均已独立精确复核；数据终态为data/banana_mixed_parameter/replay.json。指定质量射线仍使用原最小秩八活动对象，其P8不变。下一步继续质量方向的饱和标记、固定纯群、一般图、独立物理基准和全文；若扩展绝对比较，需另行标记秩十对象的三个边界列。""")
    append('SOURCE_AUDIT.md',"""## 实际多参数混合连接的推导范围

MIXED_PARAMETER_CONNECTION直接使用前批已认证的公开纯矩阵与配对，未再下载、提取或重跑纯系统。新混合耦合从实际质量形式、相对Cartan、十六个取留数前的环面IBP及对称三张量唯一性导出；相对和一动机态射沿用已经核对的局部化构造。独立验证器从紧凑仿射公式起算600个曲率，和恢复程序的三张量路线不同，但共同使用FLINT有理多项式底层。函数域的实际拉回与其几何族实现已区分，未假设一般motivic t-structure，亦未由微分方程推断固定动机群。独立推导不自动构成文献新颖性声明。""")
    m=json.loads(rd('data/input_manifest.json'))
    m.update({
      'comparison_status':'Specified-fibre P8 and physical readout certified; actual rank-ten four-parameter mixed connection and Hodge frame proved. Rank-ten absolute comparison, marked saturation and fixed pure groups remain open.',
      'four_parameter_mixed_connection_certified':True,
      'four_parameter_mixed_connection_rank':10,
      'four_parameter_mixed_connection_certificate':'banana_mixed_parameter/flatness.json',
      'four_parameter_mixed_connection_replay':'banana_mixed_parameter/replay.json',
      'four_parameter_mixed_connection_proof':'../MIXED_PARAMETER_CONNECTION.md',
      'four_parameter_mixed_geometry_identification':'Actual mass Kummer pullback; IBP, Hodge and symmetric cubic tensor uniqueness.',
      'four_parameter_mixed_curvature_entries':600,
      'four_parameter_mixed_hodge_dimensions_untwisted':[1,5,9,10],
      'rank_ten_absolute_comparison_certified':False,
      'fixed_fibre_rank_eight_object_extended_to_all_parameters':False,
      'fixed_ray_rank_eight_embedding_in_rank_ten_verified':True})
    wr('data/input_manifest.json',json.dumps(m,ensure_ascii=False,indent=2)+'\n')
    audit={'status':'actual_rank_ten_mixed_connection_batch_complete',
      'timestamp':stamp,'proof_sha256':sha(ROOT/'MIXED_PARAMETER_CONNECTION.md'),
      'flatness_sha256':sha(OUT/'flatness.json'),'replay_sha256':sha(OUT/'replay.json'),
      'old_ray_hashes_preserved':old,'old_proof_hashes_preserved':old_proofs,
      'old_four_parameter_hashes_preserved':old_pure,
      'matrix_sha256':replay['matrix_sha256'],
      'all_started_calculations_ended':True,'successful_sessions':[44030,45359],
      'new_pure_transports_started':0,'absolute_rank10_comparison_certified':False,
      'full_goal_complete':False}
    (OUT/'final_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in audit.items() if k not in ['old_ray_hashes_preserved','old_four_parameter_hashes_preserved']},ensure_ascii=False),flush=True)
if __name__=='__main__':main()

