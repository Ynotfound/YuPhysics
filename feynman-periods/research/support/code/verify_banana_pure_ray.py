"""Exact local structure and an independent holomorphic-series check of the ray system."""
from pathlib import Path
from itertools import product
from math import factorial
import json,hashlib
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
folder=ROOT/'data/banana_pure_ray'
path=folder/'system.json'
data=json.loads(path.read_text(encoding='utf-8'))
s=sp.Symbol('ss');t=sp.Symbol('t');z=sp.Symbol('z')
A=sp.Matrix(data['matrix']);thresholds=[1,9,25,49,81,121,169,225]
records=[]
for a in thresholds:
    R=A.applyfunc(lambda e:sp.cancel((s-a)*e).subs(s,a))
    assert R.rank()==1 and sp.expand(R.charpoly(z).as_expr()-z**5*(z+sp.Rational(3,2)))==0
    assert R*R==-sp.Rational(3,2)*R
    records.append({'s':a,'rank':1,'eigenvalues':['0 (multiplicity 5)','-3/2']})
q=next(sp.sympify(p) for p,n in data['common_denominator_factors'] if sp.degree(sp.sympify(p),s)==5)
assert sp.gcd(q,sp.diff(q,s))==1
intervals=sp.Poly(q,s).intervals(eps=sp.Rational(1,10**6))
assert len(intervals)==3 and all(n==1 for ab,n in intervals)
assert sp.Poly(q,s).count_roots(-sp.oo,sp.oo)==3
metric=json.loads((ROOT/'data/banana_metric.json').read_text(encoding='utf-8'))
x=sp.symbols('m12 m22 m32 m42')
N=sp.Matrix(metric['inverse_C_numerator']);d=sp.sympify(metric['inverse_C_denominator'])
Y11=3*N[0,0]/(x[0]**2*d)
ray=dict(zip(x,[a/s for a in [1,4,16,64]]))
Y11ray=-12*s**4*q/sp.prod(s-a for a in thresholds)
assert sp.cancel(Y11.subs(ray)-Y11ray)==0
print('node_residues_apparent_roots_and_pairing_factor_passed',flush=True)

# Independent scalar period: constant term of (M*N)^n gives multinomial squares.
MAX=8;mu=[1,4,16,64]
H=[sp.S.Zero for i in range(6)]
for ns in product(range(MAX+1),repeat=4):
    n=sum(ns)
    if n>MAX:continue
    coef=sp.Rational(factorial(n)**2,sp.prod(factorial(k)**2 for k in ns))*sp.prod(m**k for m,k in zip(mu,ns))
    H[0]+=coef*t**(n+1)
    for i in range(4):
        if ns[i]:H[i+1]+=coef*ns[i]/mu[i]*t**n
    if ns[0]>=2:H[5]+=coef*ns[0]*(ns[0]-1)/mu[0]**2*t**(n-1)
den=sp.lcm([sp.denom(e) for e in A])
P=(den*A).applyfunc(sp.cancel)
degree=max([sp.degree(den,s)]+[sp.degree(v,s) for v in P if v!=0])
L=sp.Poly(sp.expand(t**degree*den.subs(s,1/t)),t)
Pt=[[sp.Poly(sp.expand(t**degree*P[i,j].subs(s,1/t)),t) for j in range(6)] for i in range(6)]
valuation=lambda p:min(mon[0] for mon,coef in p.terms()) if not p.is_zero else 10**6
tail=[MAX+2]+[MAX+1]*4+[MAX]
bounds=[]
for i in range(6):
    residual=sp.Poly(-t**2*L.as_expr()*sp.diff(H[i],t)-sum(Pt[i][j].as_expr()*H[j] for j in range(6)),t)
    bound=min([valuation(L)+tail[i]+1]+[valuation(Pt[i][j])+tail[j] for j in range(6)])
    assert all(coef==0 or mon[0]>=bound for mon,coef in residual.terms()),(i,bound,residual)
    bounds.append(bound-1)
out={'status':'pure_ray_local_structure_and_holomorphic_jet_verified',
     'system_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
     'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'node_residues':records,'apparent_polynomial':str(q),
     'apparent_real_root_intervals':[[str(v) for v in ab]+[int(n)] for ab,n in intervals],
     'number_of_nonreal_apparent_roots':2,'complex_root_rectangles_certified':False,
     'Yx11_along_ray':str(Y11ray),'independent_period_total_degree':MAX,
     'cleared_residual_verified_through_t_degrees':bounds,
     'all_six_series_residuals_verified':True,'global_transport_certified':False,
     'global_four_variable_flatness_certified':False}
(folder/'verification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':out['status'],'series_residual_degrees':bounds,'apparent_real_roots':3,'nonreal_roots':2}),flush=True)
