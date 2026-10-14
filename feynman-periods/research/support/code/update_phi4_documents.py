"""Save the explicit elastic phi4 benchmark and preserve older certificates."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import json,hashlib,re
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/phi4_benchmark'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rd(p):return (ROOT/p).read_text(encoding='utf-8')
def wr(p,s):(ROOT/p).write_text(s,encoding='utf-8')
def append(p,s):
    old=rd(p)
    if s.strip() not in old:wr(p,old.rstrip()+'\n\n'+s.strip()+'\n')
def main():
    b=json.loads(rd('data/phi4_benchmark/benchmark.json'));t=json.loads(rd('data/phi4_benchmark/elastic_total.json'))
    assert b['graphs_sha256']==sha(OUT/'graphs.json')
    assert b['verifier_sha256']==sha(ROOT/'code/verify_phi4_benchmark.py')
    assert t['benchmark_sha256']==sha(OUT/'benchmark.json')
    assert t['script_sha256']==sha(ROOT/'code/compute_phi4_elastic_total.py')
    prior=json.loads(rd('data/graph_support/final_audit.json'));old=dict(prior['old_hashes_preserved'])
    for path,k in [('BOUNDED_GRAPH_SUPPORT.md','proof_sha256'),('data/graph_support/certificate.json','certificate_sha256'),('data/graph_support/replay.json','replay_sha256')]:old[path]=prior[k]
    old['data/graph_support/final_audit.json']=sha(ROOT/'data/graph_support/final_audit.json')
    for path,digest in old.items():assert sha(ROOT/path)==digest,(path,'changed old evidence')
    stamp=datetime.now(timezone(timedelta(hours=8))).strftime('%Y-%m-%d %H:%M')
    latest=f"""### 最新批次：四维phi4两圈弹性振幅与独立物理检验已执行

{stamp}（Asia/Shanghai）。PHI4_BENCHMARK.md给出实际计算而非基准提案；data/phi4_benchmark保存图与数值。旧K3、动机、格及G1a所有证明/证书保持哈希。本批全部命令正常结束，无待观察进程。

1. 审阅G1a极点/桥极限与奇偶压缩的范围，原证明未改。恢复原综合稿第39页F1要求，并核对PDG第2/3/6/7页相空间、通量和同种末态约定；本文dPhi包含(2pi)^4，与PDG书写约定不同。
2. 固定四维单实标量、正极点质量m、单位外腿留数、M_tree=−g；交叉对称在壳解析点s=t=u=4m²/3定义g。所有普通四价图按外腿标签穷尽：树1张，一圈7张，两圈42张。独立Wick半边配对复核全部对称因子；两圈8505个剩余配对。非局部1PI部分为三双泡(各1/4)、六wineglass(各1/2)；三带tadpole的1PI与质量反项抵消，三十外腿可约图按LSZ/单位留数处理。
3. 完整三通道振幅M=−g+g²A1+g³A2，A1=sum b/(32pi²)，A2=sum(−b²−4J)/(1024pi⁴)。B为标准泡，b=B−B(4/3)；J给一次减除色散积分及独立三参数Schwinger森林减除。后者在x=−1、0、3用24/48/96/160阶Gauss重放，最终与色散值差约10^−9。两表示等同性来自图参数/切割证明，不用数值吻合作替代。
4. wineglass F=(c+d)[m²((a+b)²+(a+b)(c+d)+cd)−sab]精确核对，给正常阈值4m²。角平均A=3−(1−beta²)atanh(beta)²/beta²−2atanh(beta)/beta；二圈吸收部−beta(ReB+2A−3B0)/(512pi³)已与完整树/一圈干涉核对。奇态由Z2禁止；四体树振幅g²，其对光学定理首次贡献g⁴，不在二圈振幅g³检验内。无规范Ward恒等式。
5. 微分截面含g⁴的|A1|²−2ReA2，完整交叉对称；阈值耦合g_th=g−a g²−b g³的有限转换和可观测量不变性已用一般符号精确核对。a≈0.0038708444461173248，b≈−0.00001420121720122013。
6. 全角总弹性截面也已计算，交换角积分/色散积分的独立式(11)重放四个能量5/8/12/20，24点嵌套差最大约1.3e−23。s/m²=8时sigma=[g²+0.00356759925535171g³−0.0000461288386699841g⁴]/(256pi m²)。5/8/12低于四体阈值，在此阶等于包容总截面；20处仅弹性。
7. 数值均为高精度求积与相容性检查，不是严格区间认证，不能把48位输出当作48位已证精度。benchmark.json的旧not_done含总积分，但已由elastic_total.json补齐。图/代数检查为精确；源页和证明范围已明确。
8. graphs SHA {sha(OUT/'graphs.json')}；benchmark SHA {sha(OUT/'benchmark.json')}；elastic_total SHA {sha(OUT/'elastic_total.json')}。F1弹性部分已有结果，但四体开通后的2→4、数值区间和必要的更高阶光学验证不冒充完成，完整总目标仍保留。

"""
    s=rd('CONTINUATION.md');marker='### 最新批次：G1a能量支撑压缩、极点证明审阅及独立有限重放'
    if '### 最新批次：四维phi4两圈弹性振幅' not in s:
        assert marker in s
        s=s.replace(marker,latest+'### 前一批次：G1a能量支撑压缩、极点证明审阅及独立有限重放',1)
        s=re.sub(r'更新于 .*?（Asia/Shanghai）',f'更新于 {stamp}（Asia/Shanghai）',s,count=1)
    a=s.index('下一批优先（以开头最新批次为准）：');z=s.index('已读取参考：',a)
    s=s[:a]+"""下一批优先（以开头最新批次为准）：

1. 审阅PHI4_BENCHMARK的森林减除、泡/wineglass系数、总体耦合条件及相空间的2pi/2!约定；独立三参数与色散数值已交叉，勿重复普通图枚举或旧K3运输。
2. F1下一步为s>16m²的2→4树振幅和四体相空间。应有十个3|3分割图，M6=−g² sum 1/(Q_A²−m²)，全末态除4!；须实际核对组合及归一化。可先生成有限相空间验证，再做严格积分包围。不能把s=20现有弹性数值当包容总截面。
3. 弹性数值虽已独立核对仍需误差认证；先选一两个基准点做严格包围，不声称所有打印小数可信。全截面g4对应三圈前向吸收部；如要同阶光学检验须另给全切割，原g3检验不替代。
4. G1剩余仍含带权割端点参数算法和带实际运动学的割几何/PF；G1a已完成。P1/P2须将全部证明整合正文、明确来源和适用域，避免只是累积附节。
5. P6/P7/P8、固定非零/最小性、固定Q活动群和饱和整基等旧证据保持。四个旧质量矩阵、theta和可选P10独立列，不重跑已完成运输。
6. 全部本批命令已结束，无待观察项目计算。单一heartbeat续跑；get_goal旧blocked不代表研究停滞，不能改成完成或缩小总目标。

"""+s[z:];wr('CONTINUATION.md',s)
    s=rd('TASKS.md');lines=[]
    for line in s.splitlines():
        if line.startswith('| F1 |'):
            line='| F1 | 独立massive phi4两圈物理基准 | A1、P1 | 全通道、反项、LSZ、交叉、匹配阶数光学定理与方案转换 | PHI4_BENCHMARK已执行弹性振幅/微分及总截面；开通四体通道与严格数值误差仍待完成 |'
        lines.append(line)
        if line.startswith('| F1 |'):
            lines.extend([
              '| F1a | 两圈弹性振幅与可观测量 | F1规范 | 全图精确权重、收敛积分、幺正性/方案匹配 | 已计算与独立交叉；数值不是区间证书 |',
              '| F1b | 四体产生与包容截面 | F1a可并行 | 十图树振幅、4!、四体相空间、与弹性项按阶组合 | 待执行 |',
              '| F1c | 物理数值误差认证 | F1a、F1b | 选定基准点的严格求积误差、独立重放 | 待执行；高工作精度不替代 |'])
    s='\n'.join(lines)+'\n';a=s.index('## 当前理论批次');z=s.index('## 决策门槛',a)
    s=s[:a]+"""## 当前理论批次

1. F1a已落地：四维实标量、极点质量、单位留数、对称点耦合；42个两圈图及权重精确Wick复核，实际三通道两圈振幅与弹性微分/总截面计算。PHI4_BENCHMARK给参数/色散独立表示。
2. 两圈光学检验到g³，有限耦合转换对截面到g⁴的代数不变性已证。指定能区4<s/m²<16只有弹性态；不把该范围扩至四体开通后。
3. 下一批优先审阅森林减除和归一化，执行F1b四体产生，再补F1c严格求积；打印48位数不是48位证书。
4. G1a已完成，G1b带权割及G1c割几何/PF仍继续。P1/P2整合全文与引用；不把组合分类等同动机分类。
5. 原K3/活动动机/整格证据全保持，无新ODE。可选P10及旧theta/质量矩阵不替代未完成主线。

"""+s[z:];wr('TASKS.md',s)
    append('CLAIMS.md',"""## 独立phi4物理基准：本批弹性范围

| 主张 | 证据 | 适用范围 |
|---|---|---|
| 普通四价图穷尽及对称因子 | 图自同构生成，与独立Wick半边配对复核 | 固定四外腿标签，V=1/2/3，共1/7/42图 |
| 两圈完整弹性振幅 | 三通道双泡/wineglass，二点OS/LSZ，子图及总体森林减除 | 单实标量、等正质量、对称点定义g |
| 独立实际wineglass求值 | Schwinger三参数减除与一次减除色散两式 | 三亚阈值点普通数值吻合约1e−9；非区间证书 |
| 交叉与光学定理 | 显式对称和全部二粒子树/一圈干涉 | 光学到振幅g³；非三圈g⁴前向验证 |
| 弹性截面到g⁴ | 包含一圈平方与树/两圈干涉；全角积分另独立换序 | s/m²=5/8/12是本阶总截面；20只给弹性 |
| 有限方案转换 | 对称点和弹性阈值耦合，符号精确核对可观测量 | 相同极点质量和留数，不宣称MS匹配 |
| 四体开通后的包容截面 | 未计算2→4相空间 | 不能用已有弹性总积分替代 |
| 高精度输出与误差 | 60位工作精度、独立表示/积分顺序对照 | 不等于严格区间，不保证全部48位打印数字 |
""")
    append('PAPER.md',r"""## 独立物理基准：重整化四维phi4弹性散射

Fix the pole mass and unit external residue and define the quartic coupling at the crossing-symmetric on-shell analytic point s=t=u=4m²/3. With \(B(x)=-\int_0^1\log(1-xv(1-v)-i0)\,dv\), \(b=B-B(4/3)\), and
\[
 J(x)=(x-4/3)\int_4^\infty
 \frac{\sqrt{1-4/r}\,[A(r)-B(4/3)]}{(r-4/3)(r-x-i0)}\,dr,
 \qquad A(r)=\tfrac12\int_{-1}^1B(-(r-4)(1-z)/2)\,dz,
\]
the two-loop amplitude is
\[
 \mathcal M=-g+\frac{g^2}{32\pi^2}\sum_{\xi=s,t,u}b(\xi/m^2)
 +\frac{g^3}{1024\pi^4}\sum_{\xi=s,t,u}[-b(\xi/m^2)^2-4J(\xi/m^2)]
 +O(g^4).
\]
PHI4_BENCHMARK derives a separately convergent Schwinger-parameter representation for J with both its subdivergence and overall subtraction explicit. An exhaustive labelled graph inventory and independent Wick pairing recover the three double-bubble and six wineglass contributions and account for tadpoles, local insertions and LSZ terms. The result is crossing-symmetric and satisfies the complete two-particle unitarity relation at orders g² and g³. A finite conversion to a threshold-defined coupling preserves the differential cross section through g⁴.

The physical cross section includes both the one-loop square and the tree/two-loop interference. For s/m²=8, numerical evaluation gives
\[
 \sigma_{\rm el}\simeq
 \frac{g^2+0.00356759925535171g^3-0.0000461288386699841g^4}
      {256\pi m^2}.
\]
Independent parameter and reordered spectral quadratures provide checks; these numerical values have not yet received rigorous interval error bounds. Below the four-particle threshold the elastic result is the inclusive cross section at this order. Above that threshold the tree-level 2→4 contribution must still be included. This benchmark is independent of the banana period and is not asserted as a new two-loop scattering result.
""")
    append('README.md',"""## 当前物理结果：两圈phi4弹性基准

[PHI4_BENCHMARK](PHI4_BENCHMARK.md)已从原稿提案推进到实际弹性计算：42张两圈连通图的权重经独立Wick配对复核，完整通道/反项/外腿规范明确，给出两圈振幅、微分和总弹性截面、匹配阶数光学检验及有限方案转换。wineglass参数与色散积分独立交叉。数值不是区间证书；四体开通后的2→4及误差认证仍待推进。""")
    append('SOURCE_AUDIT.md',"""## phi4独立物理基准的规范

直接查看原综合稿第39页，确认两圈振幅、截面到四阶、一圈平方、四体开通、方案转换和外腿均为原定范围。另查看PDG 2025 Kinematics第2/3/6/7页的外态、相空间、通量及两体运动学。本文dPhi包含(2pi)^4，而PDG把其外置；同种末态另除2!，四体将除4!。原PDF和哈希在external/phi4_benchmark/manifest.json。

两圈图与森林减除由本项目直接推导并计算，未把搜索摘要当作振幅公式依据。曾尝试打开Tong分章PDF但页面读取失败，未引用其正文。普通高精度求积只作相容性检查，不能升级为严格包围；两种表示数值吻合也不代替解析比较证明。本例是标准phi4问题的独立重放，不作新颖性主张。""")
    m=json.loads(rd('data/input_manifest.json'));m.update({
      'physical_benchmark_proof':'../PHI4_BENCHMARK.md','physical_benchmark_graphs':'phi4_benchmark/graphs.json',
      'physical_benchmark_values':'phi4_benchmark/benchmark.json','physical_total_elastic':'phi4_benchmark/elastic_total.json',
      'physical_benchmark_status':'Equal-mass D4 real phi4 elastic amplitude through g3 and cross sections through g4 executed, with explicit subtractions and independent numerical checks.',
      'physical_benchmark_numerical_interval_certified':False,
      'physical_four_body_production_computed':False,'independent_physical_benchmark_complete':False,
      'physical_elastic_part_computed':True,'physical_optical_check_order':'g^3 amplitude; no independent g^4 forward three-loop audit'})
    wr('data/input_manifest.json',json.dumps(m,ensure_ascii=False,indent=2)+'\n')
    audit={'status':'two_loop_phi4_elastic_benchmark_executed','timestamp':stamp,
      'proof_sha256':sha(ROOT/'PHI4_BENCHMARK.md'),'graphs_sha256':sha(OUT/'graphs.json'),
      'benchmark_sha256':sha(OUT/'benchmark.json'),'elastic_total_sha256':sha(OUT/'elastic_total.json'),
      'sources_sha256':sha(ROOT/'external/phi4_benchmark/manifest.json'),'old_hashes_preserved':old,
      'all_started_calculations_ended':True,'new_ODE_transports':0,
      'graph_weights_exactly_replayed':True,'numerical_values_interval_certified':False,
      'inelastic_four_body_complete':False,'full_F1_complete':False,'full_goal_complete':False}
    wr('data/phi4_benchmark/final_audit.json',json.dumps(audit,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in audit.items() if k!='old_hashes_preserved'},ensure_ascii=False))
if __name__=='__main__':main()

