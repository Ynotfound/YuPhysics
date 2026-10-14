"""Publish the four-parameter checkpoint after completed independent replay."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import json,hashlib,re
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/banana_four_parameter'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rd(name):return (ROOT/name).read_text(encoding='utf-8')
def wr(name,s):(ROOT/name).write_text(s,encoding='utf-8')
def append(name,s):
    text=rd(name)
    if s.strip() not in text:wr(name,text.rstrip()+'\n\n'+s.strip()+'\n')
def main():
    replay=json.loads((OUT/'replay.json').read_text())
    assert replay['status']=='four_parameter_source_and_all_identities_independently_replayed'
    assert replay['verifier_sha256']==sha(ROOT/'code/verify_banana_four_parameter.py')
    assert replay['certificate_sha256']==sha(OUT/'flatness.json')
    prior=json.loads((ROOT/'data/banana_pure_ray/relative_top/final_audit.json').read_text())
    assert sha(ROOT/'RELATIVE_TOP_COMPARISON.md')==prior['proof_sha256']
    old=dict(prior['old_hashes_preserved'])
    old.update({
      'lower_comparison/period_matrix.json':'264ddd733392c8b4a08afefebc88a29ac3a29d956d88698aff1b20467149f939',
      'relative_top/period_matrix.json':'0eeca4619c39f5a3c3156dd1d31beeb43585e838bb824d8565037fb14c46c5b2',
      'relative_top/physical_readout.json':'8cb313ee98adccc4811b116f5c90d8068473ce862f52eebd647883204e2d33e3',
      'relative_top/replay.json':'ca66d4beb668ac4c9f6a9febe15385466ff843fb33109350ed2f4bae54c10307'})
    for name,h in old.items():assert sha(ROOT/'data/banana_pure_ray'/name)==h
    stamp=datetime.now(timezone.utc).astimezone(timezone(timedelta(hours=8))).strftime('%Y-%m-%d %H:%M')
    checkpoint=f"""### 最新批次：全四参数纯平坦性、绝对杯积延伸与几何群

{stamp}（Asia/Shanghai）。FOUR_PARAMETER_CONNECTION.md与RELATIVE_TOP_AUDIT.md为本批证明。完整系数证书在data/banana_four_parameter。提取的前次尝试仅因根式语法未支持而退出，已复用四矩阵并完成平方根覆盖下降；最终提取命令exit0。全恒等式生成24629、独立原文件复核19196均exit0；行列式核对命令内exit0。没有待观察项目计算，没有重跑任何纯运输。

1. 审阅RELATIVE_TOP_COMPARISON的相对微分、Cartan符号、膜零初值、边界环境投影和实填充。阈值支撑避让与有序T³交数相容。独立比较C_E逆的(6,8)条目和C_M的(8,6)条目，均给I=−2nu_w。原证明、P8与读出证书哈希未变。

2. 从四份原始IBP矩阵恢复完整Q(x1,x2,x3,x4)纯6×6系统，ε=0的120个ISP耦合精确为零。完整144个连接条目和36个配对条目均保存。K_raw(1,6)的根式在覆盖上逐系数证明下降，不作数值选支。

3. 清分母D次数18、E次数5后，216个曲率条目与144个配对条目所有有理系数均为零。独立验证器重新用Mathematica语法解析原文件、手写指数偏导与外积矩阵乘法，重放全部360项；180个条目与旧质量射线数据逐项相等，负对照被拒绝。flatness SHA {sha(OUT/'flatness.json')}；replay SHA {sha(OUT/'replay.json')}。

4. 紧凑四阶C满足detC=−27Delta，Delta是16个带符号平方根线性因子的积。完整迹恒等式trA=dlogq5−dlogx1−sumdlogxi−3/2 dlogDelta已精确核对；极化框架变为−sumdlogxi−1/2 dlogDelta。因此多参数闭路detT=(−1)^windDelta。

5. 实际几何杯积与公开配对都在同一全参数系统中平坦，且在旧基点完整矩阵已经由Morse定标相等，故C_geom=−2K^-1延伸到连通光滑族。只延伸纯配对尺度，不提供其他纤维的绝对周期初值。由射线群包含及全族极化上界，四参数纯几何Zariski闭包同为O(Q6)；不识别固定纯R_H/R_N或离散整群。

6. 三个质量比单位x1/x4、x2/x4、x3/x4在函数域模常数独立，赋值直接证明。因此不能将目标纤维的单一log2、秩八活动对象机械铺成全参数秩八系统。多参数混合连接须保留实际环境或另证稳定子商。

7. 下一项优先实际多参数混合相容性与质量方向/边界粘合的饱和标记；固定纯群、有界一般图算法、另设物理基准和主稿全文仍继续。全四参数“纯”平坦性与纯绝对尺度不再列为缺失。原P6/P7/P8、固定regulator非零、最小秩八和13维幺幂核均保持。旧证书未改，单一heartbeat继续。

"""
    s=rd('CONTINUATION.md')
    if '### 最新批次：全四参数纯平坦性' not in s:
        s=re.sub(r'更新于 .*?（Asia/Shanghai）',f'更新于 {stamp}（Asia/Shanghai）',s,count=1)
        marker='### 最新批次：完整秩八相对比较、物理链标记与物理积分';assert marker in s
        s=s.replace(marker,checkpoint+'### 前一批次：完整秩八相对比较、物理链标记与物理积分',1)
        start=s.index('下一批优先（以开头最新批次为准）：');end=s.index('已读取参考：',start)
        s=s[:start]+"""下一批优先（以开头最新批次为准）：

1. 使用已认证全四参数纯系统处理质量方向的几何标记。三种偶指数二補全需要实际质量方向闭路或边界粘合证据；群闭包O(Q6)和更高精度射线运输均不能选出补全。先核对旧质量候选路径和框架，勿把候选四矩阵当输入事实。
2. 构造实际多参数混合连接及其稳定子商，保留三个独立质量比单位。固定秩八对象只在指定质量关系下成立，不能只加一行beta_c就宣称全族混合平坦性。
3. 固定纯R_H/R_N仍需独立证据；全四参数几何O(Q6)不是这一问题的答案。可审阅已有有限域证书对额外Hodge自同态的约束，但不得以单值泛性特化。
4. 推进有界一般图算法、另设物理基准、文献定位与完整论文。新FOUR_PARAMETER_CONNECTION和RELATIVE_TOP_AUDIT已整合到PAPER，不重复原P6/P7/P8或纯运输。
5. 新后端在runtime/python-flint，本机局部轮包已保存；旧Python环境不变。新系数文件extraction/flatness内的“待复核”字段是生成时状态，当前应读取replay与final_audit。

"""+s[end:]
        wr('CONTINUATION.md',s)
    s=rd('TASKS.md')
    s=s.replace('裸留数Morse定标给射线代数trace为公开配对的−2倍；全四参数平坦性仍待独立核对',
        '全四参数216个曲率及144个配对条目已独立认证；裸留数几何trace=公开配对的−2倍已延伸至连通光滑族；混合多参数另列')
    a=s.index('## 当前理论批次');b=s.index('## 决策门槛',a)
    s=s[:a]+"""## 当前理论批次

1. H5d：纯四参数系统、全平坦性和配对归一化已完成，见FOUR_PARAMETER_CONNECTION。下一项为包含三个质量Kummer方向的实际混合连接，不能固定单一beta_c。
2. H4e：完整抽象T与实际Q6指数二已经确定。需要质量方向闭路或边界粘合识别实际三种饱和类之一；旧质量矩阵仍为候选。
3. M1：固定活动对象13维幺幂核、完整P8与实际物理coaction框架已证；固定纯R_H/R_N仍待独立识别。
4. G1/F1：完成有界一般图分类与范围明确的独立物理基准。banana指定积分不能替代全部基准。
5. P1/P2：把现有定理和证书整理成连贯全文，审阅一手来源和新颖性范围。下方各批次小节保留历史，以本节和CONTINUATION开头为准。

"""+s[b:]
    wr('TASKS.md',s)
    append('TASKS.md',"""## 全四参数纯认证后的任务边界

H5d的纯连接验收已通过：180个原始条目、120个零耦合、216个曲率、144个配对条目及180个旧射线限制等式均独立复核。纯绝对尺度已跨参数延伸，四参数几何群闭包O(Q6)由包含关系证明。剩余多参数任务指实际混合对象、质量方向整格标记及特定路径认证，不再指纯系统未知。RELATIVE_TOP_AUDIT完成前批膜、共轭填充、物理对偶及coaction索引审阅。""")
    append('CLAIMS.md',"""## 全四参数纯连接与归一化

| 主张 | 证据 | 适用范围 |
|---|---|---|
| 全四参数纯系统平坦 | 清分母后216个多项式条目恒零，原文件独立解析与复核 | 全Q(x1,x2,x3,x4)，不是有限参数抽样 |
| 全方向公开配对平坦 | 144个多项式条目恒零；180个旧射线限制等式核对 | 保留第六项因子二；根式项已证明有理下降 |
| 几何杯积等于公开配对的−2倍 | 全方向平坦性加已有完整基点Morse定标 | 连通光滑纯族；不推广混合周期初值 |
| 全参数det单值字符 | detT=(-1)^windDelta，精确detK和连接迹证明 | q5是微分框架退化，不产生新符号 |
| 四参数纯几何群闭包O(Q6) | 已证射线闭包包含与全族杯积上界 | 非离散整群，非固定纤维绝对纯群 |
| 固定秩八对象自动给出全参数秩八混合系统 | 不成立；三个质量比单位在函数域独立 | 全参数混合连接须明确其对象及稳定子商 |
| 新顶层几何与物理标记复核 | RELATIVE_TOP_AUDIT逐项检查相对符号、膜/填充、交数和C_E逆的读出 | 原P8证明与证书保持哈希 |""")
    s=rd('PAPER.md')
    s=s.replace('The fixed pure motivic quotient, marked integral saturation, and full multiparameter verification remain separate problems.',
      'Coefficientwise identities certify the complete four-parameter pure connection and its pairing, extending the absolute cup normalization across the connected smooth family. The geometric monodromy closure of this family is O(Q6), by the already certified ray subgroup and polarization. The fixed pure motivic quotient, marked integral saturation and actual multiparameter mixed connection remain separate problems.')
    wr('PAPER.md',s)
    append('PAPER.md',r"""## 可纳入正文的四参数纯比较定理

**Theorem.** Let W=(Ω_x,∂_1Ω_x,...,∂_4Ω_x,∂_1²Ω_x) be the bare-residue derivative frame, and let K be the inverse published pairing after the explicit last-row factor two. The four rational matrices A_i satisfy all six integrability identities and all four pairing identities over Q(x_1,x_2,x_3,x_4). The geometric algebraic cup product is −2K^{-1} throughout the connected smooth family. Its monodromy determinant is (−1)^{wind Δ}; its geometric Zariski closure, in the certified six-cycle frame, is O(Q_6).

**Proof.** Independent coefficient calculations establish 216 curvature and 144 pairing identities after clearing nonzero denominators. A second parser verifies all original entries, the vanishing ISP couplings and agreement with the earlier mass-ray certificates. The compact inverse metric yields det K=x_1²(x_1x_2x_3x_4)²Δ³/(48q_5²). Taking the logarithmic derivative gives the determinant character; q_5 disappears from the trace in the intrinsic polarized frame. The actual cup product and the published pairing are flat and agree up to the certified factor −2 at the old basepoint, hence on the connected family. Finally the ray subgroup is already Zariski dense in O(Q_6), while the family preserves the same cup product. See FOUR_PARAMETER_CONNECTION for conventions, exact data and the scope of extension.

These pure matrices originate in the published unequal-mass system and are not claimed as a new differential equation. Their explicit comparison and verification do not identify the fixed-fibre motivic group. The three independent generic Kummer units also prevent extending the special rank-eight mixed object by merely retaining its single log2 direction.""")
    append('README.md',"""## 当前结果：完整四参数纯系统

[FOUR_PARAMETER_CONNECTION](FOUR_PARAMETER_CONNECTION.md)完成四参数纯连接的216个曲率和144个配对恒等式，独立解析原文件、重放全部系数并验证与旧射线的180个等式。纯几何杯积的−2归一化已沿连通光滑族延伸；全族纯几何群闭包为O(Q6)。[RELATIVE_TOP_AUDIT](RELATIVE_TOP_AUDIT.md)完成完整P8的几何与符号审阅。下一步是实际多参数混合对象、饱和整格标记、固定纯群与其余总目标；所有旧运输及P6/P7/P8证书保持不变。

新增独立复核入口为code/verify_banana_four_parameter.py，使用本工作区局部runtime/python-flint及保存轮包；不会改变旧程序环境。data/banana_four_parameter/replay.json是复核终态。""")
    append('SOURCE_AUDIT.md',"""## 全四参数纯系统核对

重新核对Duhr–Maggio–Porkert–Semper–Stawinski arXiv:2507.23061v1原始HTML及已保存PDF/一手Bonn数据。公开矩阵和配对归作者所有；新工作为完整源条目比较、全平坦性可重放证书及与本项目裸留数的几何归一化。独立解析180条目，检查120个ISP零耦合、360个恒等式、180个射线限制式。K_raw(1,6)的平方根乘积以代数覆盖证明下降，不依赖复数数值支。python-flint仅作确定性QQ多项式运算，其0.8.0本机轮包哈希保存在证书中。原先已认证的纯比较和本次平坦性一起支持全族trace尺度延伸；不据此推广固定混合扩张或动机群。""")
    m=json.loads(rd('data/input_manifest.json'))
    m.update({
      'comparison_status':'Full P8 and physical readout certified at the specified fibre; complete four-parameter pure connection verified; fixed pure group and actual mixed parameter system remain open.',
      'banana_pairing_normalization':'Published pairing times -2 equals geometric algebraic cup on the connected smooth four-parameter family, by full flatness and fixed Morse calibration.',
      'four_parameter_pure_flatness_certified':True,
      'four_parameter_pure_pairing_certified':True,
      'four_parameter_pure_certificate':'banana_four_parameter/flatness.json',
      'four_parameter_pure_replay':'banana_four_parameter/replay.json',
      'four_parameter_pure_proof':'../FOUR_PARAMETER_CONNECTION.md',
      'four_parameter_determinant_certificate':'banana_four_parameter/determinant.json',
      'absolute_pure_cup_scale_on_connected_family_proved':True,
      'four_parameter_geometric_monodromy_Zariski_closure':'O(Q6)',
      'four_parameter_mixed_connection_certified':False,
      'relative_top_geometry_review':'../RELATIVE_TOP_AUDIT.md',
      'fixed_fibre_rank_eight_object_extended_to_all_parameters':False})
    wr('data/input_manifest.json',json.dumps(m,ensure_ascii=False,indent=2)+'\n')
    audit={'status':'four_parameter_pure_batch_complete_and_old_certificates_unchanged','timestamp':stamp,
       'old_hashes_preserved':old,'old_relative_top_proof_sha256':prior['proof_sha256'],
       'new_proofs_sha256':{n:sha(ROOT/n) for n in ['FOUR_PARAMETER_CONNECTION.md','RELATIVE_TOP_AUDIT.md']},
       'certificate_sha256':sha(OUT/'flatness.json'),'replay_sha256':sha(OUT/'replay.json'),
       'determinant_sha256':sha(OUT/'determinant.json'),
       'all_started_calculations_ended':True,'successful_sessions':[24629,19196],
       'new_pure_transports_started':0,'full_goal_complete':False}
    (OUT/'final_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(audit,ensure_ascii=False),flush=True)
if __name__=='__main__':main()
