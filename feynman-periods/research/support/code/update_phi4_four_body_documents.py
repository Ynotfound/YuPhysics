"""Checkpoint the new four-body work; preserve every prior proof/certificate hash."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import json,hashlib,re,math,tomllib
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/phi4_four_body'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rd(p):return (ROOT/p).read_text(encoding='utf-8')
def wr(p,s):(ROOT/p).write_text(s,encoding='utf-8')
def append(p,s):
    old=rd(p)
    if s.strip() not in old:wr(p,old.rstrip()+'\n\n'+s.strip()+'\n')
def main():
    prod=json.loads(rd('data/phi4_four_body/production.json'))
    checks=json.loads(rd('data/phi4_four_body/checks.json'))
    bounds=json.loads(rd('data/phi4_four_body/bounds.json'))
    replay=json.loads(rd('data/phi4_four_body/bounds_replay.json'))
    assert len(prod['rows'])==4 and all(len(r['levels'])==3 for r in prod['rows'])
    assert replay['bounds_sha256']==sha(OUT/'bounds.json')
    for f,p in [('compute_phi4_four_body.py',prod),('verify_phi4_four_body.py',checks),
                ('bound_phi4_four_body.py',bounds),('verify_phi4_four_body_bounds.py',replay)]:
        assert p['script_sha256']==sha(ROOT/'code'/f)
    prior=json.loads(rd('data/phi4_benchmark/final_audit.json'));old=dict(prior['old_hashes_preserved'])
    for path,k in [('PHI4_BENCHMARK.md','proof_sha256'),('data/phi4_benchmark/graphs.json','graphs_sha256'),
                   ('data/phi4_benchmark/benchmark.json','benchmark_sha256'),('data/phi4_benchmark/elastic_total.json','elastic_total_sha256'),
                   ('external/phi4_benchmark/manifest.json','sources_sha256')]:
        old[path]=prior[k]
    old['data/phi4_benchmark/final_audit.json']=sha(ROOT/'data/phi4_benchmark/final_audit.json')
    for p,h in old.items():assert sha(ROOT/p)==h,(p,'prior hash changed')
    elastic=json.loads(rd('data/phi4_benchmark/elastic_total.json'))['rows'][-1]
    assert elastic['s_over_m2']==20
    coeff=elastic['sigma_coefficients_inside_1_over_32pi_s']
    four=prod['rows'][0]['levels'][-1]
    inc=dict(status='inclusive_through_g4_numerical_estimate',s_over_m2=20,
             convention='1/(640*pi*m^2) times [g^2+C3*g^3+C4*g^4]',
             C3=coeff[1],C4_elastic=coeff[2],
             C4_four_body=32*math.pi*20*four['sigma_over_g4'],
             C4_inclusive=float(coeff[2])+32*math.pi*20*four['sigma_over_g4'],
             four_body_replicate_standard_error_in_C4=32*math.pi*20*four['replicate_standard_error'],
             certified_narrow_interval=False,elastic_sha256=sha(ROOT/'data/phi4_benchmark/elastic_total.json'),
             production_sha256=sha(OUT/'production.json'))
    wr('data/phi4_four_body/inclusive.json',json.dumps(inc,indent=2)+'\n')
    stamp=datetime.now(timezone(timedelta(hours=8))).strftime('%Y-%m-%d %H:%M')
    latest=f"""### 最新批次：四体产生、阈值律、严格粗界与包容截面

{stamp}（Asia/Shanghai）。PHI4_FOUR_BODY.md新增F1b实算；data/phi4_four_body保存证据。本批所有计算已结束（主生成会话81865已exit0），独立复核与文档命令均exit0，无待观察项目计算。旧PHI4_BENCHMARK及所有K3/动机/格/G1a证书保持哈希，原始资料未改。

1. 审阅旧森林减除：Schwinger雅可比T³t、U=T²t(1+tY)、F=T³t(Dx+t+t²Y)恢复两个对数的正确符号；子图及整体有限常数与同一个对称点耦合一致。PDG相空间转换为包含(2pi)^4的dPhi，递归dq²/(2pi)，入射通量2s beta，无入射2!。
2. 2→4的十张无序3|3树图已独立Wick半边注入穷尽，每图权重1；M24=−g² sum1/(Q_A²−m²)，全末态除4!。四项分母≥8m²，六项严格类空且分母<−m²；给出全物理域无内部极点与收敛证明。
3. 两对和逐次3+1两种实际四体相空间均执行，各8种子、2^14/17/19点，最终每种分解4194304事件。s/m²=20的m²sigma24/g⁴约5.75125323e−11和5.75125308e−11；32约1.14052712e−9和1.14052993e−9。独立20项振幅、全部24末态置换、二维体积、壳与守恒均核对。标准误为随机重复估计，不是严格包围，不认证打印数字。
4. 严格解析阈值：epsilon=sqrt(s)/m−4，S6→−1/m²，Phi4~m⁴sqrt2 epsilon^(7/2)/(13440pi⁴)，sigma24~g⁴/m² sqrt(2/3)epsilon^(7/2)/(5160960pi⁴)。两相空间归一化另以无质量体积s²/(24576pi⁵)精确核对。
5. s/m²=20全域混合分母[−10,−2]、S6≤−1/10，结合128×128正权Darboux包围得到严格Phi4∈(7.58e−8,7.63e−8)，m²sigma24/g⁴∈(8.8e−13,6.2e−10)。180位Arb与220位独立端点因子算法重放通过；这是严格粗界和正性证书，不能认证5.75125的小数精度。已有本地flint环境未改。
6. 20点现已加四体得到包容截面：sigma_incl≈[g²+0.0180991989932455g³+0.00014250842g⁴]/(640pi m²)。旧elastic_total仍仅弹性，新的inclusive.json才含四体。该合并仍是普通数值；32点本批未算弹性。有限耦合转换不改变四体g⁴首项，故同阶包容方案匹配成立。
7. 三圈前向g⁴光学检验明确写出二体(|A1|²−2ReA2)/2!与四体|S6|²/4!全右端，未独立计算三圈左端或全部切割；不把定义等式当额外验证。F1c窄误差认证和若要求的三圈全切割仍保留。
8. production SHA {sha(OUT/'production.json')}；bounds SHA {sha(OUT/'bounds.json')}；bounds_replay SHA {sha(OUT/'bounds_replay.json')}；proof SHA {sha(ROOT/'PHI4_FOUR_BODY.md')}。不重复已有四体采样、旧弹性图或K3运输。完整目标未完成。

"""
    s=rd('CONTINUATION.md')
    marker='### 最新批次：四维phi4两圈弹性振幅与独立物理检验已执行'
    if marker in s:
        s=s.replace(marker,latest+marker.replace('最新','前一'),1)
    else:
        assert '### 最新批次：四体产生、阈值律、严格粗界与包容截面' in s
    s=re.sub(r'更新于 .*?（Asia/Shanghai）',f'更新于 {stamp}（Asia/Shanghai）',s,count=1)
    a=s.index('下一批优先（以开头最新批次为准）：');z=s.index('已读取参考：',a)
    nextsteps="""下一批优先（以开头最新批次为准）：

1. 审阅PHI4_FOUR_BODY的无极点、阈值首项和正权Darboux界，保留普通数值与严格粗界的区别。四体生成已exit0，不因旧heartbeat文本重复生成。20点包容数值已合并，32点仅四体。
2. 可转向尚未实施的G1b：恢复原稿负边端点参数的精确全割算法，输出反例割或正性证书并独立有限重放；不能把全割coNP结果套给dual-critical。G1c须输入实际运动学及明确图族/次数界，再做割正规化、亏格/PF/不可约性。
3. F1c优先选一个弹性基准点给严格窄误差，现四体粗界只证明正性和范围。若要求光学到总截面g⁴，则必须独立完成三圈前向全切割；现有振幅g³和右端组合不替代。
4. P1/P2需将全部证明包整合成连贯全文，清理历史待证语句、列清来源。不要只累积更多矩阵或附节。G1a及F1a/F1b已有实际结果但总目标未完成。
5. 旧K3、实际相对对象、P6/P7/P8、非零性、固定Q群与饱和整基证据保留，不重复ODE/点数。旧theta/四质量矩阵及可选P10另列，不能替代剩余主线。
6. 全部本批项目进程已结束。单一heartbeat续跑；get_goal旧blocked不代表研究停滞，不能缩小目标或标为完成。

"""
    s=s[:a]+nextsteps+s[z:];wr('CONTINUATION.md',s)
    s=rd('TASKS.md');lines=[]
    for line in s.splitlines():
        if line.startswith('| F1 |'):
            line='| F1 | 独立massive phi4两圈物理基准 | A1、P1 | 全通道、反项、LSZ、交叉、匹配阶光学及方案转换 | 弹性及开通四体通道已执行，20点包容截面已合并；窄数值认证与若要求的三圈全切割仍待完成 |'
        if line.startswith('| F1b |'):
            line='| F1b | 四体产生与包容截面 | F1a可并行 | 十图树振幅、4!、四体相空间、阶数合并 | PHI4_FOUR_BODY完成两分解实际积分、精确阈值及20点合并；十图权重独立Wick复核 |'
        if line.startswith('| F1c |'):
            line='| F1c | 物理数值误差认证 | F1a、F1b | 选定基准点严格求积误差、独立重放 | 新四体严格粗界/正性通过；弹性及显示小数的窄误差仍待证 |'
        lines.append(line)
    s='\n'.join(lines)+'\n';a=s.index('## 当前理论批次');z=s.index('## 决策门槛',a)
    s=s[:a]+"""## 当前理论批次

1. F1b已执行十张2→4树图、两种四体相空间、相干干涉与20点包容截面；阈值首项有解析证明。
2. F1c已有严格粗界与独立重放，但窄数值误差认证未完成；高精度/随机重复不能替代。
3. 旧森林减除与2pi/同种因子已审阅。光学至振幅g³已有检验；总截面g⁴所需三圈前向全切割另列。
4. 下一理论主线可转G1b带权割算法、G1c运动学割几何/PF及P1/P2全文整合。不要重复已完成G1a/物理采样/K3运输。
5. 全部旧证据哈希保持，完整目标尚未完成。

"""+s[z:];wr('TASKS.md',s)
    s=rd('CLAIMS.md')
    s=s.replace('| 四体开通后的包容截面 | 未计算2→4相空间 | 不能用已有弹性总积分替代 |',
                '| 四体开通后的包容截面 | 后续PHI4_FOUR_BODY已补齐十图及相空间 | 20点已与弹性合并；普通数值与严格粗界区分 |')
    wr('CLAIMS.md',s)
    append('CLAIMS.md',"""## 四体产生的补齐与误差层级

| 主张 | 证据 | 适用范围 |
|---|---|---|
| 十个无序3\|3树图，每图权重1 | 独立Wick半边注入穷尽 | 四维单实标量phi4，六个外腿有标签 |
| 物理域无内部传播子极点 | 四个三未来向组合及六个严格类空组合 | 等正质量、真实2→4运动学 |
| 完整四体相空间与干涉 | 两对/逐次两分解、20项独立标量振幅与24个置换 | s/m²=20、32普通数值；不称区间 |
| 阈值epsilon^(7/2)系数 | 精确参数因式分解与两个Beta权重积分 | epsilon→0+，不宣称有限epsilon误差 |
| 四体严格正性及粗区间 | 180/220位独立正权Darboux界及全域振幅界 | s/m²=20，m²sigma24/g⁴∈(8.8e−13,6.2e−10) |
| 20点包容截面到g⁴ | 新四体项与旧弹性项合并 | 系数约0.00014250842/(640pi m²)，非窄区间 |
| 窄数值精度/三圈前向全切割 | 仍需独立证据 | 粗界与光学右端书写均不替代 |
""")
    append('PAPER.md',r"""## Four-particle completion of the scalar benchmark

With all six external momenta outgoing, the tree-level production amplitude is
\[
 \mathcal M_{24}^{(0)}=-g^2\sum_{|A|=3,\ 0\in A}(Q_A^2-m^2)^{-1}.
\]
The ten unordered partitions have unit weights, as verified by independent Wick injections. Four propagators have denominator at least 8m²; the other six are strictly spacelike. Thus the massive physical phase-space integral has no internal propagator singularity. The contribution to the inclusive cross section is
\[
 \sigma_{24}^{(4)}=\frac{g^4}{2s\sqrt{1-4m^2/s}\,4!}
        \int d\Phi_4\left|\sum_A(Q_A^2-m^2)^{-1}\right|^2 .
\]
PHI4_FOUR_BODY gives two independent pair/sequential phase-space constructions, a proof of the threshold law
\[
 \sigma_{24}^{(4)}\sim\frac{g^4}{m^2}
 \frac{\sqrt{2/3}}{5160960\pi^4}(\sqrt{s}/m-4)^{7/2},
\]
and both numerical crosschecks and separately certified coarse bounds. At s=20m² the production coefficient is estimated as 5.75125×10⁻¹¹/m². The rigorous bound (8.8×10⁻¹³,6.2×10⁻¹⁰)/m² proves positivity but does not certify those displayed digits.

Adding production to the earlier elastic result gives the numerical inclusive estimate
\[
 \sigma_{\rm incl}(20m^2)\simeq
 [g^2+0.0180991989932455g^3+0.00014250842g^4]/(640\pi m^2)+O(g^5).
\]
The leading production term is unchanged through g⁴ under the finite coupling conversion already established for elastic scattering. A full independent optical check at this order requires the three-loop forward cut sum; the two-loop amplitude check alone does not supply it.
""")
    append('SOURCE_AUDIT.md',"""### 四体相空间与严格算术的来源核对

PHI4_FOUR_BODY沿用已保存的PDG Kinematics原页49.12–13、49.18，明确把(2pi)^4吸收进dPhi后递归因子为dq²/(2pi)。十图、无极点、阈值及严格正性界在本项目内独立推导，不作为文献新颖性主张。https://pdg.lbl.gov/2025/reviews/rpp2025-rev-kinematics.pdf

Arb球区间语义核对官方文档 https://python-flint.readthedocs.io/en/latest/arb.html ，本次实际使用既存本地python-flint 0.8.0轮包，未修改原运行环境。
""")
    outputs=['production.json','checks.json','bounds.json','bounds_replay.json','inclusive.json']
    audit=dict(status='four_body_benchmark_and_rigorous_coarse_bound_completed',timestamp=stamp,
               proof_sha256=sha(ROOT/'PHI4_FOUR_BODY.md'),data_sha256={f:sha(OUT/f) for f in outputs},
               old_hashes_preserved=old,script_sha256=sha(__file__),
               flint_wheel_sha256=sha(ROOT/'runtime/wheels/python_flint-0.8.0-cp311-cp311-win_amd64.whl'),
               numeric_precision_certified=False,coarse_positive_interval_certified=True,
               independent_three_loop_optical_check=False,full_goal_complete=False)
    wr('data/phi4_four_body/final_audit.json',json.dumps(audit,indent=2)+'\n')
    # Prepare (do not write) the existing app automation update.
    config=tomllib.loads(Path('E:/CodexData/.codex/automations/automation/automation.toml').read_text(encoding='utf-8'))
    prompt=config['prompt'];a=prompt.index('下一步先审阅新森林减除');z=prompt.index('状态未改变且无可操作信息时',a)
    newblock=f"""最新PHI4_FOUR_BODY.md已审阅森林与相空间归一化并完成四体产生。十个3|3树图独立Wick权重全1；四类时/六类空分母给全物理域无极点。两对与逐次3+1两种相空间在s/m²=20、32均执行，各8种子×2^19事件；20的m²sigma24/g⁴约5.75125e−11，32约1.1405e−9。随机标准误不是严格包围。精确阈值首项为g⁴/m² sqrt(2/3)epsilon^(7/2)/(5160960pi⁴)。s20全域分母界结合128²正权Darboux积分，严格Phi4∈(7.58e−8,7.63e−8)、m²sigma24/g⁴∈(8.8e−13,6.2e−10)，180位生成/220位独立重放通过；粗界不能认证展示小数。20点当前包容截面已合并为[g²+0.0180991989932455g³+0.00014250842g⁴]/(640pi m²)，不要再说20只有弹性；旧elastic_total文件仍只含弹性，新inclusive.json才含四体。production SHA{sha(OUT/'production.json')}；bounds SHA{sha(OUT/'bounds.json')}；replay SHA{sha(OUT/'bounds_replay.json')}；proof SHA{sha(ROOT/'PHI4_FOUR_BODY.md')}。主生成会话81865已exit0，全部独立/文档命令也结束，无待观察项目计算。旧证据哈希保留，不重复采样、旧42图或K3运输。

下一步先审阅新无极点/阈值/正权包围证明，再重点推进G1b负边端点参数精确割算法与独立有限证书、G1c实际运动学割几何/PF及P1/P2全文整合。F1c窄数值认证仍需选点推进，当前四体粗界不等于小数精度。三圈前向g⁴全切割仍未独立执行，光学右端的正确组合不能充当左端验证。G1a及F1a/F1b已有结果，不把组合支撑等同动机分类；旧theta/四质量矩阵和可选P10不能代替未完成主线。按证明与证书范围写结论，完整总目标未完成。

"""
    prompt=prompt[:a]+newblock+prompt[z:]
    payload=dict(id=config['id'],mode='update',kind=config['kind'],name=config['name'],prompt=prompt,
                 status=config['status'],rrule=config['rrule'],targetThreadId=config['target_thread_id'])
    if 'notification_policy' in config:payload['notificationPolicy']=config['notification_policy']
    wr('data/phi4_four_body/automation_update.json',json.dumps(payload,ensure_ascii=False,indent=2))
    print(json.dumps({k:v for k,v in audit.items() if k!='old_hashes_preserved'},indent=2))
    print('old proofs/certificates preserved:',len(old))
if __name__=='__main__':main()

