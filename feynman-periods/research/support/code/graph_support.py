"""Exact compressed generic-energy support for finite labelled connected multigraphs.

Input: n vertices 0..n-1 and labelled edges in input order; loops allowed.
Output: parity-bundle certificate, never a geometric or motivic label.
"""
from itertools import combinations
import json,argparse
class BudgetExceeded(ValueError):pass

def normalize(n,edges):
    if type(n)!=int or n<1:raise ValueError('n must be a positive integer')
    pairs={};loops=[]
    for j,edge in enumerate(edges):
        if len(edge)!=2 or any(type(x)!=int or x<0 or x>=n for x in edge):
            raise ValueError('invalid edge endpoints')
        u,v=sorted(edge)
        if u==v:loops.append(j)
        else:pairs.setdefault((u,v),[]).append(j)
    reached={0}
    while True:
        old=len(reached)
        for u,v in pairs:
            if u in reached or v in reached:reached.update((u,v))
        if len(reached)==old:break
    if len(reached)!=n:raise ValueError('connected input graph required')
    bundles=[{'endpoints':list(p),'edge_ids':ids,'multiplicity':len(ids)} for p,ids in sorted(pairs.items())]
    return bundles,loops

def dual_critical_dp(n,pairs):
    adj=[0]*n
    for u,v in pairs:
        if u==v:raise ValueError('loopless parity graph required')
        adj[u]^=1<<v;adj[v]^=1<<u
    good=[False]*(1<<n);last=[-1]*(1<<n)
    for A in range(1,1<<n):
        if A&(A-1)==0:
            good[A]=True;last[A]=A.bit_length()-1;continue
        rest=A
        while rest:
            bit=rest&-rest;rest-=bit;v=bit.bit_length()-1;B=A^bit
            if good[B] and (adj[v]&B).bit_count()%2:
                good[A]=True;last[A]=v;break
    order=[];A=(1<<n)-1
    if good[A]:
        while A:
            v=last[A];order.append(v);A^=1<<v
        order.reverse()
    packed=bytearray((len(good)+7)//8)
    for A,x in enumerate(good):
        if x:packed[A//8]|=1<<(A%8)
    return {'dual_critical':good[-1],
            'reachable_subsets_hex':format(int.from_bytes(packed,'little'),'x'),
            'good_order':order if good[-1] else None}

def compile_support(n,edges,*,vertex_limit=12,state_budget=1<<22):
    # No partial output or silent truncation when an explicit budget is exceeded.
    if type(n)!=int or n<1:raise ValueError('n must be a positive integer')
    if n>vertex_limit:raise BudgetExceeded('vertex bound exceeded')
    bundles,loops=normalize(n,edges);b=len(bundles)
    states=1<<(n+b)
    if states>state_budget:raise BudgetExceeded('parity/subset state budget exceeded')
    cases=[]
    for p in range(1<<b):
        active=[x['endpoints'] for j,x in enumerate(bundles) if (p>>j)&1]
        row=dual_critical_dp(n,active)
        row['retained_parity_mask']=p
        row['deleted_parities']=[(x['multiplicity']-((p>>j)&1))%2 for j,x in enumerate(bundles)]
        cases.append(row)
    hits=sum(x['dual_critical'] for x in cases);m0=len(edges)-len(loops)
    return {'schema':'energy-support-parity-v1','n':n,'edges':[list(e) for e in edges],
       'bundles':bundles,'forced_loop_ids':loops,'cases':cases,
       'positive_parity_cases':hits,'generic_support_count':str(hits*(1<<(m0-b))),
       'degree_polynomial':'z^loops sum_good_p product_b ((1+z)^r+(-1)^(r-p)(1-z)^r)/2',
       'state_budget':state_budget,'states_allocated_over_all_cases':states,
       'scope':'Independent energy roots; independent conserving vertex energies; unit scalar propagators.',
       'geometric_classification':'not_evaluated'}

def contains(cert,deleted_ids):
    ids=list(deleted_ids)
    if len(set(ids))!=len(ids) or any(type(i)!=int or i<0 or i>=len(cert['edges']) for i in ids):
        raise ValueError('deleted edge labels must be a subset')
    S=set(ids)
    if not set(cert['forced_loop_ids']).issubset(S):return False
    mask=sum(((sum(i not in S for i in b['edge_ids'])%2)<<j) for j,b in enumerate(cert['bundles']))
    return cert['cases'][mask]['dual_critical']

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input');p.add_argument('output');p.add_argument('--state-budget',type=int,default=1<<22)
    p.add_argument('--vertex-limit',type=int,default=12);a=p.parse_args()
    with open(a.input,encoding='utf-8') as f:d=json.load(f)
    result=compile_support(d['n'],d['edges'],vertex_limit=a.vertex_limit,state_budget=a.state_budget)
    with open(a.output,'w',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
if __name__=='__main__':main()

