"""Integrate bounded graph support while preserving all earlier research certificates."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import json,hashlib,re
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/graph_support'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rd(p):return (ROOT/p).read_text(encoding='utf-8')
def wr(p,s):(ROOT/p).write_text(s,encoding='utf-8')
def append(p,s):
    old=rd(p)
    if s.strip() not in old:wr(p,old.rstrip()+'\n\n'+s.strip()+'\n')
def main():
    cert=json.loads(rd('data/graph_support/certificate.json'))
    replay=json.loads(rd('data/graph_support/replay.json'))
    assert replay['status']=='independent_permutation_and_exact_kernel_replay_passed'
    assert replay['certificate_sha256']==sha(OUT/'certificate.json')
    assert replay['verifier_sha256']==sha(ROOT/'code/verify_graph_support.py')
    oldaudit=json.loads(rd('data/marked_saturation/final_audit.json'))
    old=dict(oldaudit['old_hashes_preserved'])
    for path,key in [('MARKED_SATURATION.md','proof_sha256'),
       ('data/marked_saturation/certificate.json','certificate_sha256'),
       ('data/marked_saturation/periods.json','periods_sha256'),
       ('data/marked_saturation/replay.json','replay_sha256')]:
        old[path]=oldaudit[key]
    old['data/marked_saturation/final_audit.json']=sha(ROOT/'data/marked_saturation/final_audit.json')
    for path,h in old.items():assert sha(ROOT/path)==h,(path,'old input changed')
    stamp=datetime.now(timezone(timedelta(hours=8))).strftime('%Y-%m-%d %H:%M')
    latest=f"""### 最新批次：G1a能量支撑压缩、极点证明审阅及独立有限重放

{stamp}（Asia/Shanghai）。新证明BOUNDED_GRAPH_SUPPORT.md与data/graph_support已保存。当前批次完成一般图任务的组合模块，不把完整G1、物理基准或总目标标为完成。无新点数、ODE、纯周期或混合比较计算；旧证明和全部证书哈希保持。

1. 先审阅MARKED_SATURATION的十一类穷尽、A_K[2]判别式型和实际L⊂T输入，未发现需修改之处。实际饱和整基结论继续有效，不重复运输。
2. 从原稿energy_character_classification.pdf第4、5、8、11页恢复单位标量能量核、字符定义与dual-critical判据，并查看原页。图论一手Király–Kisfaludi-Bak arXiv1410.1653的好排序、奇割递归Theorem3和平行边对不变性Proposition4已核对。
3. 极点因子化与无极点衰减证明已逐项写清。补充不连通侧伪极点消去的有理桥极限证明：以c−1条根eta的桥连接分量，乘(2eta)^(c−1)，用桥能量固定的乘积分解严格取零极限，不依赖分布极限。
4. 平行边按保留数奇偶压缩。q个自环强制进入S，b个非空边束只需2^b个奇偶图；每个接受图恢复2^(m−q−b)个带标签字符。给出完整按删边数的生成多项式，确定性子集算法和2^O(n²)poly(m,n)位复杂度上界。默认n≤12且2^(n+b)≤2^22，超限明确拒绝；不把预算拒绝当成空支撑。
5. 全部1099个五顶点以内简单图、432个三顶点含重边/自环图已用独立排序穷举核对。七顶点细分K4的512个奇偶图另以全排列检查；能量核512次根翻号各直接求和5040个排序，再按定义做字符投影，复现59项=9个单根+50个三根。生成树另枚举得50。模非零证明指定三外点能量切片上的非零；恒零由理论承担。百边banana压缩给2^99个字符。
6. 原综合稿第42–43页coNP归约只适用于有理带符号权图的全割严格正性。不能推广为dual-critical判定或一般动机分类的复杂度定理。2014一手文献中的开放状态仅作历史引用，不冒称当前未决状态。
7. certificate SHA {sha(OUT/'certificate.json')}；replay SHA {sha(OUT/'replay.json')}。生成和独立验证均exit0，独立验证约7秒，无待观察计算。原稿59项是重现，不宣称新结果或59维周期空间。

"""
    s=rd('CONTINUATION.md')
    marker='### 最新批次：判别式二次型唯一标记饱和格与整基绝对比较'
    if '### 最新批次：G1a能量支撑压缩' not in s:
        assert marker in s
        s=re.sub(r'更新于 .*?（Asia/Shanghai）',f'更新于 {stamp}（Asia/Shanghai）',s,count=1)
        s=s.replace(marker,latest+'### 前一批次：判别式二次型唯一标记饱和格与整基绝对比较',1)
    a=s.index('下一批优先（以开头最新批次为准）：');b=s.index('已读取参考：',a)
    s=s[:a]+"""下一批优先（以开头最新批次为准）：

1. 审阅BOUNDED_GRAPH_SUPPORT的极点证明、桥极限、平行边计数与资源契约。G1a已经完成，勿把它重新写为尚未执行，勿反复重跑已有有限枚举。
2. F1现为适合的下一主线：恢复原稿massive phi4两圈2→2的设定，明确维数、方案和振幅归一化，逐步执行独立物理基准。banana指定积分不代替此项；实标量理论不虚构规范Ward恒等式。先做有可验收输出的有限步骤。
3. G1剩余层保持：带权割安全性可按原稿负边端点参数实现精确证书；割几何需外部运动学和割方程，边表不足；亏格/Picard–Fuchs/不可约性尚未自动化。不要把字符支撑当作动机分类，也不要从指数输出推出一般判定困难。
4. P1/P2必须从累积证明包整合正文、来源归属和最新状态，清理旧段落的待证语句。尚无嵌套对象证明的可见性过滤保留条件。原稿字符支撑定理的文献优先权另需系统检索。
5. 所有K3点数、全闭路、Q6、21迹、P6/P7/P8、固定非零/最小性、固定Q群、实际饱和整基均保留。旧theta/四质量矩阵与可选P10不作为重复启动运输的理由；纯整基不自动认证混合整基。
6. 本批所有命令结束，无待观察项目计算。保持单一heartbeat；旧get_goal blocked字段不代表当前无进展，不新建目标或缩小完整目标。

"""+s[b:]
    # Replace a stale global status that otherwise contradicts the top checkpoint.
    s=s.replace('全部有限闭路、21迹及指定六循环有理基的绝对P6已认证；指定饱和半循环、完整混合比较和固定纯群识别继续；指定regulator非零已在新批次完成。一般图算法、明确物理基准、完整论文与证书均保留在总目标中。',
       '指定纤维的全闭路、21迹、实际Q6、纯/活动混合比较、固定非零和纯/活动群、实际饱和整基均已有证明与证书。G1a组合支撑已执行；一般图剩余几何层、独立物理基准和全文整合仍在总目标中。')
    wr('CONTINUATION.md',s)
    s=rd('TASKS.md');lines=[]
    for line in s.splitlines():
        if line.startswith('| G1 |'):
            parts=line.split('|');parts[-2]=' G1a已完成并独立重放；带权割、几何/PF/不可约性仍待各自有界实现 ';line='|'.join(parts)
        lines.append(line)
        if line.startswith('| G1 |'):
            lines.extend([
             '| G1a | 能量字符支撑压缩 | A1、A2 | 泛独立根/守恒能量、边束奇偶、复杂度、独立有限重放 | BOUNDED_GRAPH_SUPPORT完成；1099简单图、432多重图及59项直接能量证书通过 |',
             '| G1b | 带权割安全性 | G1a可并行 | 有理权重、负边端点参数、精确割证书与界 | 原稿已有理论，本批未实施 |',
             '| G1c | 有界割几何与微分方程 | A1、A2 | 带运动学的割方程、正规化、亏格/PF/不可约性及终止界 | 未实现，不能从组合支撑自动推导 |'])
    s='\n'.join(lines)+'\n';a=s.index('## 当前理论批次');b=s.index('## 决策门槛',a)
    s=s[:a]+"""## 当前理论批次

1. G1a：完成原稿支撑定理审阅，补充有理桥极限，证明平行边压缩并执行独立精确认证。泛支撑不等于特化支撑，59项不是59维周期空间。
2. G1剩余：带权割安全性已有原稿算法，下一层可实现端点参数证书；边表到割几何须补运动学/割方程，PF与不可约性另定范围。
3. F1：下一重点为原massive phi4两圈2→2独立物理基准，先恢复约定和有限验收。banana单点积分不替代；无规范Ward约束可虚构。
4. P1/P2：整合完整正文及文献归属。组合附节已经可审阅，整篇论文仍未完成；不给未证可见性嵌套结构冠以定理。
5. 原K3及指定活动动机证明包保持，不重复点数/运输/饱和筛选。旧theta/四质量矩阵及可选P10另列，不替代当前未完成主线。

"""+s[b:];wr('TASKS.md',s)
    s=rd('CLAIMS.md').replace(
      '| energy character 支撑等价于 G 减 S 的 dual critical 性 | 综合稿有定义、定理与附录证明，未在本轮逐行审核 | 后续专项核查一般性、残数与无极点衰减论证 |',
      '| energy character 支撑等价于 G 减 S 的 dual critical 性 | BOUNDED_GRAPH_SUPPORT已审阅残数/无极点衰减，补充有理桥极限 | 泛独立根与守恒顶点能量、单位标量传播子；不保证每个物理特化非零 |')
    wr('CLAIMS.md',s)
    append('CLAIMS.md',"""## 有界一般图的组合模块G1a

| 主张 | 证据 | 范围 |
|---|---|---|
| 平行边奇偶图完整决定泛字符支撑 | 好排序的入度奇偶不变及原支撑定理 | 每个自环必须属于删边字符 |
| 精确压缩计数及度数多项式 | 每束固定删边奇偶有2^(r−1)个子集，二项式筛选 | 按字符计数一次，不按好排序重复计数 |
| 参数n的确定性FPT | b≤n(n−1)/2，枚举2^b奇偶图，每图子集DP | 位复杂度2^O(n²)poly(m,n)，不是一般多项式或树宽FPT声明 |
| 有界实现与独立核对 | 1099简单图、432重边/自环图；排列与删除枚举 | 全部有限指定覆盖，不拿测试代替任意图证明 |
| 三圈细分K4恰59项 | 原稿结果；512次根翻号×5040排序的精确有限域直接核验证 | 9个单根、50个三根；指定三外点能量切片非零 |
| coNP复杂度的适用范围 | 原稿附录C第42–43页MAXCUT归约 | 全割带符号权严格正性，非dual-critical或动机分类 |
| 一般图几何/PF/动机分类 | 尚未执行 | 不由G1a推得，也不称不可判定 |
""")
    append('PAPER.md',r"""## 组合附节：可认证的泛能量支撑压缩

**Theorem.** Let G be a connected labelled multigraph with n vertices, m edges, q loops, and b nonempty bundles of parallel nonloop edges of multiplicities r_j. For unit scalar propagators and independent roots and conserving vertex energies, the character support is represented exactly by the parity graphs H_p which are dual-critical. Every character contains all loops and satisfies |S∩B_j|=r_j−p_j modulo two. Consequently
\[
 |\mathcal V(G)|=2^{m-q-b}\#\{p:H_p\text{ is dual-critical}\},
 \qquad
 \sum_{S\in\mathcal V(G)}z^{|S|}
 =z^q\sum_{p:H_p\ {\rm dual\mbox{-}critical}}
  \prod_j\frac{(1+z)^{r_j}+(-1)^{r_j-p_j}(1-z)^{r_j}}2.
\]
A deterministic subset recurrence constructs the compressed certificate in \(2^{O(n^2)}\operatorname{poly}(m,n)\) bit operations. Listing all individual supported characters is not part of this complexity assertion.

**Proof.** The generic energy support criterion is proved by residues of the time-ordered kernel, the odd-cut recursion, and decay in a generic external-energy direction. BOUNDED_GRAPH_SUPPORT supplies all signs and a rational bridge regularization of disconnected-side cancellation. In a fixed vertex ordering, deleting a pair of parallel edges preserves every predecessor-degree parity. Thus the criterion depends only on H_p. Each prescribed deletion parity in a nonempty bundle is satisfied by exactly half its subsets; multiplication gives the count and binomial parity filters give the polynomial. Removing the last vertex gives the exact subset recurrence and the complexity bound, since b≤n(n−1)/2.

The graph-theoretic notion and odd-cut equivalence are due to the prior dual-critical literature; see Király–Kisfaludi-Bak, arXiv:1410.1653, Theorem3 and Proposition4. The present certificate independently checks all1099 labelled simple graphs of at most five vertices and432 bounded multigraphs. It reproduces the original subdivided-K4 example with59 supported characters by evaluating512 root reflections, each through all5040 time orderings, over an exact finite field. These counts do not classify motives or count master integrals. Physical specializations and cut geometry require separate data and proofs.
""")
    append('README.md',"""## 当前进展：一般图组合支撑已可执行

[BOUNDED_GRAPH_SUPPORT](BOUNDED_GRAPH_SUPPORT.md)审阅完整泛能量字符判据，补充伪极点消去证明，并用平行边奇偶压缩实现有界精确算法。1099个简单图、432个多重图和原三圈59字符例子已独立重放。百边banana的2^99个支撑无需逐项枚举。G1a完成；带权割、割几何/PF、独立物理基准和完整论文仍须推进。证据保存在data/graph_support，旧K3与动机证明/数值证书全部保留。""")
    append('SOURCE_AUDIT.md',"""## G1a的来源与复杂度纠偏

已查看energy_character_classification.pdf第4、5、8、11页：核的规范、残数奇偶因子、动态规划和59项例子。原稿定理的本批审阅写于BOUNDED_GRAPH_SUPPORT，新增桥极限证明把不连通侧相消写成有理恒等式。59项及所用有限域参数均归属于原稿，当前工作是独立恢复与重放，不声称新颖性。

一手文献Király–Kisfaludi-Bak, Notes on dual-critical graphs, arXiv1410.1653，已查看第1–3页的定义、Theorem3和Proposition4；外部文件及原稿哈希见external/graph_algorithms/manifest.json。平行边对不变性是既有图论事实，压缩支撑及度数公式由其与原稿支撑定理推导。该文2014年摘要的复杂度未决状态仅按当时描述引用。

另查看feynman_periods_expanded.pdf第42–43页附录C，coNP归约只涉及有理带符号割的全正性，不能外推为dual-critical或动机分类困难。尚未做字符支撑定理完整优先权检索；不使用“首次”表述。所有外部关键公式直接核对原页，不仅依赖抽取文本。""")
    m=json.loads(rd('data/input_manifest.json'));m.update({
      'graph_support_status':'G1a generic unit-scalar energy support compressed, proved and independently replayed',
      'graph_support_proof':'../BOUNDED_GRAPH_SUPPORT.md','graph_support_certificate':'graph_support/certificate.json',
      'graph_support_replay':'graph_support/replay.json',
      'graph_support_coverage':{'simple_graphs':1099,'multigraphs':432,'example_characters':59},
      'graph_support_default_limits':{'vertices':12,'total_parity_subset_states':2**22},
      'graph_geometric_classification_complete':False,'independent_physical_benchmark_complete':False,
      'comparison_status':'Specified-fibre pure/active mixed comparisons, fixed-Q groups and marked saturation proved; G1a energy support certified. Graph geometric layers, independent physical benchmark and full manuscript remain.'})
    wr('data/input_manifest.json',json.dumps(m,ensure_ascii=False,indent=2)+'\n')
    audit={'status':'bounded_graph_support_module_complete','timestamp':stamp,
      'proof_sha256':sha(ROOT/'BOUNDED_GRAPH_SUPPORT.md'),
      'certificate_sha256':sha(OUT/'certificate.json'),'replay_sha256':sha(OUT/'replay.json'),
      'old_hashes_preserved':old,
      'signed_cut_source_sha256':sha(Path('E:/document/feynman_period/feynman_periods_expanded.pdf')),
      'signed_cut_pages_viewed':[42,43],'new_ODE_transports':0,'all_started_calculations_ended':True,
      'G1a_complete':True,'full_G1_complete':False,'F1_complete':False,'full_goal_complete':False}
    wr('data/graph_support/final_audit.json',json.dumps(audit,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in audit.items() if k!='old_hashes_preserved'},ensure_ascii=False))
if __name__=='__main__':main()

