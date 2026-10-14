"""Exact connected and 1PI phi^4 four-point inventory through two loops."""
from pathlib import Path
from itertools import product,permutations
from math import factorial
from fractions import Fraction
from collections import Counter
import json,hashlib
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/phi4_benchmark'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canon(V,external,edges):
    return min((tuple(p[v] for v in external),tuple(sorted(tuple(sorted((p[a],p[b]))) for a,b in edges))) for p in permutations(range(V)))
def connected(V,edges):
    parent=list(range(V))
    def find(x):
        while x!=parent[x]:x=parent[x]
        return x
    for u,v in edges:parent[find(u)]=find(v)
    return len({find(v) for v in range(V)})==1
def internal_graphs(V,degrees):
    slots=[(u,v) for u in range(V) for v in range(u,V)]
    def rec(j,left,es):
        if j==len(slots):
            if not any(left):yield list(es)
            return
        u,v=slots[j];lim=left[u]//2 if u==v else min(left[u],left[v])
        for k in range(lim+1):
            deg=list(left);deg[u]-=k;deg[v]-=k
            yield from rec(j+1,deg,es+[(u,v)]*k)
    yield from rec(0,degrees,[])
def main():
    OUT.mkdir(exist_ok=True);allrows=[];summary={}
    for V in (1,2,3):
        seen={}
        for external in product(range(V),repeat=4):
            deg=[4-external.count(v) for v in range(V)]
            for edges in internal_graphs(V,deg):
                if not connected(V,edges):continue
                key=canon(V,external,edges)
                if key in seen:continue
                ext,es=key;mult=Counter(es)
                av=sum(tuple(p[v] for v in ext)==ext and tuple(sorted(tuple(sorted((p[a],p[b]))) for a,b in es))==es for p in permutations(range(V)))
                denominator=av
                for (u,v),k in mult.items():denominator*=factorial(k)*(2**k if u==v else 1)
                onepi=all(connected(V,es[:j]+es[j+1:]) for j in range(len(es)))
                loops=sum(u==v for u,v in es)
                edeg=sorted(Counter(ext).get(v,0) for v in range(V))
                topology='tree' if V==1 else ('bubble' if V==2 and onepi and not loops else
                    ('double_bubble' if V==3 and onepi and not loops and edeg==[0,2,2] else
                     ('wineglass' if V==3 and onepi and not loops else ('tadpole_1PI' if onepi else 'external_reducible'))))
                row={'vertices':V,'loops':V-1,'external_vertices':list(ext),'internal_edges':[list(e) for e in es],
                     'symmetry_weight':str(Fraction(1,denominator)),'one_particle_irreducible':onepi,
                     'selfloops':loops,'topology':topology}
                seen[key]=row
        rows=sorted(seen.values(),key=lambda x:(x['external_vertices'],x['internal_edges']))
        allrows.extend(rows)
        counts=Counter(x['topology'] for x in rows)
        weights={k:str(sum((Fraction(x['symmetry_weight']) for x in rows if x['topology']==k),Fraction())) for k in counts}
        summary[str(V)]={'labelled_connected_graphs':len(rows),'topology_counts':dict(counts),'weight_sums':weights}
    data={'convention':'External labels 0,1,2,3 fixed; L_int=-g phi^4/4!; weights before vertex/propagator phases.',
          'graphs':allrows,'summary':summary,'generator_sha256':sha(__file__),
          'counterterm_scope':'Ordinary quartic graphs only; local insertions and LSZ treated in proof.'}
    (OUT/'graphs.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary))
if __name__=='__main__':main()

