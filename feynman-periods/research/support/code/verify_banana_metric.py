"""Recover the ten rational couplings in the published Frobenius normalization.

Checks: immutable-source hash, compact inverse formula, discriminant identity,
exact inversion and independent multivariate Frobenius jets. This does not
calibrate the pairing against a geometrically marked integral Betti basis.
"""
from pathlib import Path
from itertools import product, combinations
from fractions import Fraction
from math import factorial
import hashlib,json
import sympy as sp
from sympy.parsing.mathematica import parse_mathematica
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'external/bonn-banana-2025'
source=SOURCE/'ZCYInverse.m'
manifest=json.loads((SOURCE/'manifest.json').read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(source)==next(f['sha256'] for f in manifest['files'] if f['name']==source.name)
Zinv=sp.Matrix(parse_mathematica(source.read_text(encoding='utf-8')))
x=sp.symbols('m12 m22 m32 m42')
assert Zinv.shape==(6,6) and Zinv.free_symbols==set(x)
assert all(Zinv[i,5]==Zinv[5,i]==0 for i in range(1,6))
C=sp.zeros(4)
for i in range(4):
    other=[j for j in range(4) if j!=i]
    C[i,i]=2*((x[i]-1)**2+3*(sum(x[j]**2 for j in other)
                 -2*sum(x[j]*x[k] for j,k in combinations(other,2))))
for i,j in combinations(range(4),2):
    k,l=[v for v in range(4) if v not in (i,j)]
    C[i,j]=C[j,i]=3*(x[i]+x[j]-x[k]-x[l])**2-4*x[i]*x[j]-12*x[k]*x[l]-2*(x[i]+x[j])+6*(x[k]+x[l])-1
C=C.applyfunc(sp.expand)
B=sp.diag(*x)*C*sp.diag(*x)/3
assert (B-Zinv[1:5,1:5]).applyfunc(sp.expand)==sp.zeros(4)
print('source_block_and_compact_formula_passed',flush=True)

detC=sp.expand(C.det(method='domain-ge'))
Delta=sp.Poly(-detC/27,*x)
# Independent sign-product definition of the threshold discriminant.
q=sp.symbols('q1:5')
norm=sp.Poly(1,*q)
for signs in product((-1,1),repeat=4):
    norm*=sp.Poly(1+sum(e*v for e,v in zip(signs,q)),*q)
assert all(all(e%2==0 for e in mon) for mon,coef in norm.terms())
delta_from_signs=sp.Poly.from_dict({tuple(e//2 for e in mon):coef for mon,coef in norm.terms()},x)
assert Delta==delta_from_signs
inverse_numerator,inverse_denominator=DomainMatrix.from_Matrix(C).inv_den()
N=inverse_numerator.to_Matrix()
d=inverse_numerator.domain.to_sympy(inverse_denominator)
assert (C*N-d*sp.eye(4)).applyfunc(sp.expand)==sp.zeros(4)
print('discriminant_and_exact_inverse_passed',flush=True)

# Independent Frobenius construction: no published matrix in the series input.
ORDER=4
zero=(0,0,0,0)
indices=sorted((a for a in product(range(ORDER+1),repeat=4) if sum(a)<=ORDER),key=lambda a:(sum(a),a))
one={zero:Fraction(1)}
def add(*ss):
    out={}
    for s in ss:
        for a,c in s.items():out[a]=out.get(a,Fraction(0))+c
    return {a:c for a,c in out.items() if c}
def scale(a,c):return {m:v*c for m,v in a.items() if v*c}
def mul(a,b):
    out={}
    for i,c in a.items():
        for j,dj in b.items():
            k=tuple(i[t]+j[t] for t in range(4))
            if sum(k)<=ORDER:out[k]=out.get(k,Fraction(0))+c*dj
    return {i:c for i,c in out.items() if c}
def theta(a,i):return {m:v*m[i] for m,v in a.items() if m[i]}
harm=[sum((Fraction(1,j) for j in range(1,n+1)),Fraction(0)) for n in range(ORDER+1)]
F={a:Fraction(factorial(sum(a))**2,sp.prod(factorial(k)**2 for k in a)) for a in indices}
G=[{a:2*c*(harm[sum(a)]-harm[a[i]]) for a,c in F.items() if harm[sum(a)]!=harm[a[i]]} for i in range(4)]
Fminus=add(F,scale(one,-1));power=one;Finverse=one
for n in range(1,ORDER+1):
    power=mul(power,Fminus);Finverse=add(Finverse,scale(power,(-1)**n))
assert mul(F,Finverse)==one
regular=[mul(g,Finverse) for g in G]
E=[[add(one if a==i else {},theta(regular[i],a)) for i in range(4)] for a in range(4)]
F2=mul(F,F)
Ytheta=[[mul(F2,scale(add(*(mul(E[a][i],E[b][j]) for i in range(4) for j in range(4) if i!=j)),-1)) for b in range(4)] for a in range(4)]
def polynomial_series(expr):
    return {a:Fraction(int(c.p),int(c.q)) for a,c in sp.Poly(expr,*x).terms() if sum(a)<=ORDER}
Btheta=[[polynomial_series(C[a,b]/3) for b in range(4)] for a in range(4)]
for a in range(4):
    for b in range(4):
        assert add(*(mul(Btheta[a][c],Ytheta[c][b]) for c in range(4)))==(one if a==b else {})
print('independent_frobenius_series_degree_4_passed',flush=True)

s=sp.Symbol('s',nonzero=True)
mu=sp.symbols('mu1:5',nonzero=True)
threshold=sp.prod(s-t for t in [1,9,25,49,81,121,169,225])
ray=dict(zip(x,[a/s for a in [1,4,16,64]]))
assert sp.cancel(Delta.as_expr().subs(ray)-threshold/s**8)==0
point=dict(zip(x,[sp.Rational(-1,64),sp.Rational(-1,16),sp.Rational(-1,4),sp.Integer(-1)]))
metric_x=3*sp.diag(*[1/v for v in x])*N*sp.diag(*[1/v for v in x])/d
metric_point=metric_x.subs(point)
J=sp.diag(*[point[v] for v in x]);J[3,:]=sp.Matrix([[-point[v]/(-64) for v in x]])
metric_physical=J*metric_point*J.T/4096
det_physical=sp.factor(metric_physical.det())
assert det_physical==-sp.Rational(3,1)/(4096*threshold.subs(s,-64))
def strings(M):return [[str(v) for v in row] for row in M.tolist()]
out={'status':'published_normalization_rational_metric_verified',
     'source_dataset':'https://doi.org/10.60507/FK2/J4OPZY',
     'source_file_sha256':sha(source),'code_sha256':sha(Path(__file__)),
     'variables':[str(v) for v in x],
     'inverse_metric_formula':'Y_x^(-1)=diag(x) C diag(x)/3',
     'C':strings(C),'inverse_C_numerator':strings(N),'inverse_C_denominator':str(d),
     'discriminant':str(Delta.as_expr()),'det_C_equals_minus_27_Delta':True,
     'sign_product_discriminant_verified':True,'frobenius_jet_total_degree':ORDER,
     'frobenius_jet_monomials':len(indices),'all_16_series_identities_verified':True,
     'physical_change':'Yhat=s^(-2) J Y_x J^t; J_ai=x_i delta_ai for a<=3, J_4i=-x_i/s',
     'physical_determinant':'-3/(s^10 Delta(mu/s))',
     'physical_metric_at_m1248_s_minus64':strings(metric_physical),
     'physical_determinant_at_m1248_s_minus64':str(det_physical),
     'absolute_geometric_cup_scale_certified':False,
     'geometrically_marked_integral_basis_certified':False,
     'fixed_basepoint_regulator_nonzero_certified':False}
(ROOT/'data/banana_metric.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':out['status'],'independent_series_monomials':len(indices),'basepoint_determinant':str(det_physical),'absolute_cup_scale_certified':False}),flush=True)
