"""Independent permutation and direct time-order verification; imports no generator."""
from pathlib import Path
from itertools import combinations,permutations
from collections import Counter
from math import comb,isqrt
import hashlib,json,time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/graph_support'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(n):return json.loads((OUT/n).read_text(encoding='utf-8'))
def dc_order(vertices,edges):
    if any(u==v for u,v in edges):return None
    for order in permutations(vertices):
        earlier=set();ok=True
        for j,v in enumerate(order):
            if j and sum((u==v and w in earlier) or (w==v and u in earlier) for u,w in edges)%2!=1:
                ok=False;break
            earlier.add(v)
        if ok:return list(order)
    return None
def verify_dp(n,edges,c,full_subsets=True):
    actual=dc_order(list(range(n)),edges)
    assert (actual is not None)==c['dual_critical']
    if actual is None:assert c['good_order'] is None
    else:
        order=c['good_order'];assert sorted(order)==list(range(n))
        assert all(sum((u==v and w in order[:j]) or (w==v and u in order[:j]) for u,w in edges)%2 for j,v in enumerate(order) if j)
    if full_subsets:
        bits=0
        for mask in range(1,1<<n):
            vs=[v for v in range(n) if mask>>v&1]
            es=[(u,v) for u,v in edges if u in vs and v in vs]
            if dc_order(vs,es) is not None:bits|=1<<mask
        assert format(bits,'x')==c['reachable_subsets_hex']
def compressed_member(c,S):
    if any(j not in S for j in c['forced_loop_ids']):return False
    p=0
    for j,b in enumerate(c['bundles']):
        if sum(k not in S for k in b['edge_ids'])%2:p+=2**j
    return c['cases'][p]['dual_critical']
def degree_formula(c):
    out=Counter()
    for row in c['cases']:
        if not row['dual_critical']:continue
        poly={len(c['forced_loop_ids']):1}
        for b,par in zip(c['bundles'],row['deleted_parities']):
            new=Counter()
            for a,value in poly.items():
                for k in range(par,b['multiplicity']+1,2):
                    new[a+k]+=value*comb(b['multiplicity'],k)
            poly=new
        out.update(poly)
    return dict(out)
def verify_compiled(c,enumerate_deletions=True,full_subsets=True):
    n=c['n'];es=c['edges'];loops=[j for j,(u,v) in enumerate(es) if u==v]
    assert loops==c['forced_loop_ids']
    expected={}
    for j,(u,v) in enumerate(es):
        if u!=v:expected.setdefault(tuple(sorted((u,v))),[]).append(j)
    assert [(tuple(x['endpoints']),x['edge_ids']) for x in c['bundles']]==sorted(expected.items())
    assert len(c['cases'])==2**len(expected)
    for p,row in enumerate(c['cases']):
        assert row['retained_parity_mask']==p
        active=[b['endpoints'] for j,b in enumerate(c['bundles']) if p>>j&1]
        verify_dp(n,active,row,full_subsets)
        assert row['deleted_parities']==[(b['multiplicity']-((p>>j)&1))%2 for j,b in enumerate(c['bundles'])]
    polynomial=degree_formula(c)
    assert sum(polynomial.values())==int(c['generic_support_count'])
    if enumerate_deletions:
        nonloop=[j for j in range(len(es)) if j not in loops];hist=Counter()
        for sm in range(1<<len(nonloop)):
            S=set(loops)|{j for k,j in enumerate(nonloop) if sm>>k&1}
            retained=[e for j,e in enumerate(es) if j not in S]
            oracle=dc_order(list(range(n)),retained) is not None
            assert oracle==compressed_member(c,S)
            if oracle:hist[len(S)]+=1
            if loops:assert not compressed_member(c,S-{loops[0]})
        assert dict(hist)==polynomial
    return polynomial
def main():
    started=time.monotonic();cert=read('certificate.json')
    assert sha(ROOT/'code/graph_support.py')==cert['engine_sha256']
    assert sha(ROOT/'code/certify_graph_support.py')==cert['generator_sha256']
    for name,digest in cert['data_sha256'].items():assert sha(OUT/name)==digest
    assert sha(ROOT/'external/graph_algorithms/manifest.json')==cert['source_manifest_sha256']
    simple=read('simple_graphs.json')
    for c in simple:
        pairs=list(combinations(range(c['n']),2))
        es=[e for j,e in enumerate(pairs) if c['edge_mask']>>j&1]
        verify_dp(c['n'],es,c)
    print('simple graphs independent permutation replay: '+str(len(simple)),flush=True)
    multi=read('multigraphs.json')
    for c in multi:verify_compiled(c)
    print('multigraph deletion and degree-polynomial replay: '+str(len(multi)),flush=True)
    ex=read('examples.json');histograms={}
    for name,c in ex.items():
        histograms[name]=verify_compiled(c,enumerate_deletions=name!='banana100',full_subsets=name!='subdivided_K4')
    assert ex['banana100']['generic_support_count']==str(2**99)
    assert histograms['single_vertex_two_loops']=={2:1}
    assert histograms['two_double_bubbles']=={2:4}
    # The subset recurrence is checked on all 1099 simple graphs; the larger
    # seven-vertex example is checked by all full permutations for every parity mask.
    d=read('modular_energy.json');p=d['prime'];ii=d['imaginary_unit'];n=d['n'];edges=d['edges']
    assert all(p%r for r in range(2,isqrt(p)+1))
    assert ii*ii%p==p-1 and sum(d['vertex_energies'])==0
    chains=[]
    for order in permutations(range(n)):
        mask=0;prefix=[]
        for v in order[:-1]:
            mask|=1<<v;prefix.append(mask)
        chains.append(prefix)
    values=[];full=(1<<n)-1
    for signs in range(1<<len(edges)):
        roots=[(-w if signs>>j&1 else w)%p for j,w in enumerate(d['root_values'])]
        inv=[0]*(full+1)
        for mask in range(1,full):
            cut=sum(roots[j] for j,(u,v) in enumerate(edges) if (bool(mask>>u&1)!=bool(mask>>v&1)))
            energy=sum(d['vertex_energies'][v] for v in range(n) if mask>>v&1)
            den=(cut+ii*energy)%p;assert den
            inv[mask]=pow(den,p-2,p)
        total=0
        for chain in chains:
            term=1
            for mask in chain:term=term*inv[mask]%p
            total+=term
        pref=1
        for w in roots:pref=pref*2*w%p
        values.append(total%p*pow(pref,p-2,p)%p)
    assert values==d['signed_kernel_values']
    count=len(values);projected=[]
    for mask in range(count):
        z=sum((-value if (mask&signs).bit_count()%2 else value) for signs,value in enumerate(values))
        projected.append(z*pow(count,p-2,p)%p)
    assert projected==d['projected_values']
    support=[j for j,v in enumerate(projected) if v]
    assert support==d['nonzero_masks']
    assert all(compressed_member(ex['subdivided_K4'],{j for j in range(len(edges)) if mask>>j&1})==(mask in support) for mask in range(count))
    assert Counter(mask.bit_count() for mask in support)=={1:9,3:50}
    trees=0
    for es in combinations(edges,n-1):
        seen={0}
        while True:
            old=len(seen)
            for u,v in es:
                if u in seen or v in seen:seen.update((u,v))
            if old==len(seen):break
        trees+=len(seen)==n
    assert trees==50
    tampered=list(projected);tampered[support[0]]=0
    assert tampered!=d['projected_values']
    result={'status':'independent_permutation_and_exact_kernel_replay_passed',
      'simple_graphs':len(simple),'multigraphs':len(multi),
      'seven_vertex_parity_graphs_checked':512,'signed_energy_evaluations':512,
      'orders_per_energy_evaluation':len(chains),'nonzero_characters':len(support),
      'degree_histogram':{'1':9,'3':50},'spanning_trees':trees,
      'banana100_support':str(2**99),'explicit_examples_degree_counts':histograms,
      'negative_controls':['omitted_required_loop','tampered_projected_value'],
      'certificate_sha256':sha(OUT/'certificate.json'),'verifier_sha256':sha(__file__),
      'runtime_seconds':round(time.monotonic()-started,3),
      'claims_excluded':['general_cut_curve_geometry','Picard_Fuchs_generation','motivic_classification']}
    (OUT/'replay.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['status','simple_graphs','multigraphs','nonzero_characters','runtime_seconds']}))
if __name__=='__main__':main()

