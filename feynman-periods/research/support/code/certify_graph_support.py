"""Generate bounded graph certificates and an exact modular energy witness."""
from pathlib import Path
from itertools import combinations,product
import json,hashlib
from graph_support import compile_support,dual_critical_dp,contains,BudgetExceeded
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/graph_support'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name,x):(OUT/name).write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8')
def connected(n,edges):
    reached={0}
    while True:
        old=len(reached)
        for u,v in edges:
            if u in reached or v in reached:reached.update([u,v])
        if len(reached)==old:return len(reached)==n

def modular_energy(n,edges,roots,energies,prime,imaginary):
    full=(1<<n)-1;F=[0]*(1<<n);F[0]=1
    for A in range(1,full+1):
        total=0
        for v in range(n):
            if A>>v&1:total+=F[A^(1<<v)]
        if A==full:F[A]=total%prime;continue
        den=sum(roots[j] for j,(u,v) in enumerate(edges) if ((A>>u)^(A>>v))&1)
        den+=imaginary*sum(energies[v] for v in range(n) if A>>v&1)
        den%=prime
        if not den:raise ArithmeticError('singular evaluation')
        F[A]=total*pow(den,-1,prime)%prime
    pre=1
    for w in roots:pre=pre*(2*w)%prime
    return F[-1]*pow(pre,-1,prime)%prime

def walsh(values,p):
    a=list(values);step=1
    while step<len(a):
        for start in range(0,len(a),2*step):
            for j in range(step):
                u,v=a[start+j],a[start+j+step]
                a[start+j]=(u+v)%p;a[start+j+step]=(u-v)%p
        step*=2
    inv=pow(len(a),-1,p)
    return [x*inv%p for x in a]

def main():
    OUT.mkdir(exist_ok=True)
    simple=[]
    for n in range(1,6):
        pairs=list(combinations(range(n),2))
        for mask in range(1<<len(pairs)):
            edges=[e for j,e in enumerate(pairs) if mask>>j&1]
            simple.append({'n':n,'edge_mask':mask,**dual_critical_dp(n,edges)})
    save('simple_graphs.json',simple)
    multigraph=[]
    for mult in product(range(4),repeat=3):
        edges=[e for e,r in zip([(0,1),(0,2),(1,2)],mult) for _ in range(r)]
        if not connected(3,edges):continue
        for loops in range(8):
            es=edges+[(v,v) for v in range(3) if loops>>v&1]
            multigraph.append(compile_support(3,es))
    save('multigraphs.json',multigraph)
    examples={
      'single_vertex_two_loops':(1,[(0,0)]*2),
      'two_double_bubbles':(3,[(0,1)]*2+[(1,2)]*2),
      'banana4':(2,[(0,1)]*4),
      'banana100':(2,[(0,1)]*100),
      'K4':(4,list(combinations(range(4),2))),
      'subdivided_K4':(7,[(0,1),(0,2),(0,3),(1,4),(4,2),(2,5),(5,3),(3,6),(6,1)])}
    ex={name:compile_support(n,es) for name,(n,es) in examples.items()}
    save('examples.json',ex)
    p=1000000009;ii=430477711;assert ii*ii%p==p-1
    n,es=examples['subdivided_K4'];energies=[0,0,0,0,1,2,-3];roots=[2**j for j in range(len(es))]
    values=[]
    for signs in range(1<<len(es)):
        rr=[-w if signs>>j&1 else w for j,w in enumerate(roots)]
        values.append(modular_energy(n,es,rr,energies,p,ii))
    projected=walsh(values,p);support=[j for j,a in enumerate(projected) if a]
    expected=[mask for mask in range(1<<len(es)) if contains(ex['subdivided_K4'],[j for j in range(len(es)) if mask>>j&1])]
    assert support==expected and len(support)==59
    histogram={str(k):sum(mask.bit_count()==k for mask in support) for k in range(10)}
    assert histogram['1']==9 and histogram['3']==50
    save('modular_energy.json',{'prime':p,'imaginary_unit':ii,'n':n,'edges':es,
      'root_values':roots,'vertex_energies':energies,'signed_kernel_values':values,
      'projected_values':projected,'nonzero_masks':support,'degree_histogram':histogram,
      'scope':'Exact nonzero witnesses on specified three-external-vertex slice; modular zeros alone do not prove identities.'})
    budget_rejected=False
    try:compile_support(7,es,state_budget=16)
    except BudgetExceeded:budget_rejected=True
    assert budget_rejected
    save('certificate.json',{'status':'bounded_compressed_energy_support_generated',
       'generator_sha256':sha(__file__),'engine_sha256':sha(ROOT/'code/graph_support.py'),
       'data_sha256':{name:sha(OUT/name) for name in ['simple_graphs.json','multigraphs.json','examples.json','modular_energy.json']},
       'simple_graphs':len(simple),'multigraphs':len(multigraph),
       'subdivided_K4_support':59,'banana100_support':ex['banana100']['generic_support_count'],
       'budget_overrun_rejected':True,'source_manifest_sha256':sha(ROOT/'external/graph_algorithms/manifest.json'),
       'geometry_and_full_motive_classification':False})
    print(json.dumps({'simple_graphs':len(simple),'multigraphs':len(multigraph),'energy_characters':59,'banana100':ex['banana100']['generic_support_count']}))
if __name__=='__main__':main()

