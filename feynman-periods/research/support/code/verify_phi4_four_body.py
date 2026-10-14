"""Independent combinatorics/amplitude checks and deterministic phase-volume checks."""
from pathlib import Path
from itertools import combinations,permutations
from fractions import Fraction
import hashlib,json,math
import numpy as np
from scipy.special import roots_legendre
from compute_phi4_four_body import pair_events,sequential_events,amplitude
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/phi4_four_body'
def scalar_amplitude(s,ks):
    p=math.sqrt(s/4-1);E=math.sqrt(s)/2
    momenta=[np.array([-E,0,0,-p]),np.array([-E,0,0,p])]+list(ks)
    terms=[]
    for A in combinations(range(6),3):
        q=sum((momenta[i] for i in A),start=np.zeros(4))
        terms.append(1/(q[0]*q[0]-sum(q[j]*q[j] for j in (1,2,3))-1))
    return sum(terms)/2

def volume(s,N,kind):
    x,w=roots_legendre(N);z=(x+1)/2;w=w/2
    u=(1-np.cos(np.pi*z))/2;jac=np.pi*np.sin(np.pi*z)/2
    U,V=np.meshgrid(u,u,indexing='ij');W=np.outer(w*jac,w*jac)
    E=np.sqrt(s)
    if kind=='pair':
        a=4+((E-2)**2-4)*U;db=(E-np.sqrt(a))**2-4;b=4+db*V
        lr=np.sqrt(np.maximum(0,(s-a-b)**2-4*a*b))
        f=((E-2)**2-4)*db*lr/s*np.sqrt(1-4/a)*np.sqrt(1-4/b)
    else:
        a=9+((E-1)**2-9)*U;db=(np.sqrt(a)-1)**2-4;b=4+db*V
        f=((E-1)**2-9)*db*np.sqrt(np.maximum(0,(s-a-1)**2-4*a))/s
        f*=np.sqrt(np.maximum(0,(a-b-1)**2-4*b))/a*np.sqrt(1-4/b)
    return float(np.sum(W*f)/(2048*np.pi**5))

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    # Enumerate labelled Wick injections, not graph automorphisms.
    counts={}
    for perm in permutations(range(8),6):
        verts=[x//4 for x in perm]
        if verts.count(0)!=3:continue
        A=tuple(i for i,v in enumerate(verts) if v==0);B=tuple(i for i,v in enumerate(verts) if v==1)
        key=min(A,B);counts[key]=counts.get(key,0)+1
    weights={str(k):str(Fraction(v,2*24**2)) for k,v in counts.items()}
    assert len(weights)==10 and set(weights.values())=={'1'}
    threshold=[Fraction(1,8) if sum(i<2 for i in A) in (0,2) else Fraction(-1,4) for A in counts]
    assert sum(threshold)==-1
    rng=np.random.default_rng(88617);errs=[];sym=[]
    for factory in [pair_events,sequential_events]:
        ks,w=factory(20.,rng.random((64,8)));ss,d=amplitude(20.,ks)
        errs.append(max(abs(scalar_amplitude(20.,k)-v) for k,v in zip(ks,ss)))
        for perm in permutations(range(4)):
            vv=[scalar_amplitude(20.,k[list(perm)]) for k in ks[:4]]
            sym.append(float(np.max(np.abs(vv-ss[:4]))))
        assert np.all(np.sum(d>0,axis=1)==4) and np.all(np.sum(d<0,axis=1)==6)
        assert np.all(d[d>0]>=8-1e-11) and np.all(d[d<0]<=-2+1e-11) and np.all(d[d<0]>=-10-1e-11)
    assert max(errs+sym)<2e-13
    rows=[dict(s=s,N=N,pair=volume(s,N,'pair'),sequential=volume(s,N,'sequential')) for s in [20.,32.] for N in [24,48,96]]
    threshold_rows=[]
    for eps in [0.1,0.01,0.001]:
        s=(4+eps)**2;v=volume(s,96,'pair')
        leading=math.sqrt(2)*eps**3.5/(13440*math.pi**4)
        threshold_rows.append(dict(epsilon=eps,phase_volume=v,leading_ratio=v/leading))
    out=dict(status='exact_Wick_and_independent_numerical_checks',weights=weights,
             threshold_sum=str(sum(threshold)),amplitude_difference=max(errs),
             permutation_difference=max(sym),volumes=rows,threshold=threshold_rows,
             script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             generator_sha256=hashlib.sha256((ROOT/'code/compute_phi4_four_body.py').read_bytes()).hexdigest(),
             numerical_checks_are_intervals=False)
    (OUT/'checks.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()

