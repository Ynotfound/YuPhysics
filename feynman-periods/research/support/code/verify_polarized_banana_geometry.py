"""Independent metric/Levi-Civita check and exact regular-singular analysis."""
from pathlib import Path
import sympy as sp
import json,hashlib

ROOT=Path(__file__).resolve().parents[1]
folder=ROOT/'data/banana_pure_ray'
path=folder/'polarized_system.json'
data=json.loads(path.read_text(encoding='utf-8'))
B=sp.Matrix(data['matrix']);KB=sp.Matrix(data['inverse_pairing'])
s=sp.Symbol('ss');z=sp.Symbol('z');m1,m2,m3,m4=mu=sp.symbols('u1 u2 u3 u4')
met=json.loads((ROOT/'data/banana_metric.json').read_text(encoding='utf-8'))
x=sp.symbols('m12 m22 m32 m42')
C=sp.Matrix(met['C']).subs(dict(zip(x,[m/s for m in mu])))
L=sp.Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[-1,-1,-1,-s]])
Ki=(s**2*L.T*C*L/3).applyfunc(sp.cancel)
ray=dict(zip(mu,[1,4,16,64]))
K=Ki.subs(ray);Y=K.inv().applyfunc(sp.cancel)
assert (KB[1:5,1:5]-K).applyfunc(sp.cancel)==sp.zeros(4)
derivatives=[lambda e,m=m:m*sp.diff(e,m) for m in mu[:3]]+[lambda e:sp.diff(e,s)]
dK=[Ki.applyfunc(D).subs(ray) for D in derivatives]
dY=[(-Y*d*Y).applyfunc(sp.cancel) for d in dK]
for b in range(4):
    for c in range(4):
        g=sp.cancel(sum(K[c,d]*(dY[3][b,d]+dY[b][3,d]-dY[d][3,b]) for d in range(4))/2)
        assert sp.cancel(B[b+1,c+1]-g)==0,(b,c)
assert all(sp.cancel(B[b+1,5]+Y[3,b])==0 for b in range(4))
print('metric_and_all_16_Levi_Civita_entries_verified',flush=True)
nodes=[1,9,25,49,81,121,169,225]
residues={}
for a in [0]+nodes:
    R=B.applyfunc(lambda e:sp.cancel((s-a)*e).subs(s,a))
    assert all(e.is_Rational for e in R)
    if a:
        assert R.rank()==1 and R*R==-R/2
    residues[str(a)]=R
R0=residues['0'];E=R0**2
assert E*E==E and E.rank()==1 and R0*E==-E
# Resonant shearing at zero: erase the integral exponent -1.
J=sp.eye(6)-E+s*E;Ji=sp.eye(6)-E+E/s
Jconn=(J.diff(s)*Ji+J*B*Ji).applyfunc(sp.cancel)
den=sp.lcm([sp.denom(v) for v in Jconn]);factors=sp.factor_list(den)[1]
assert all(n==1 for f,n in factors)
N0=Jconn.applyfunc(lambda e:sp.cancel(s*e).subs(s,0))
assert N0.rank()==2 and N0**2==sp.zeros(6)
print('zero_resonance_removed_rank2_squarezero_residue',flush=True)
# At infinity use its own chart, not the zero shearing.
S=sp.diag(s,s,s,s,s*s,1/s)
Cinf=(S.diff(s)*S.inv()+S*B*S.inv()).applyfunc(sp.cancel)
for e in Cinf:
    if e:assert sp.degree(sp.numer(e),s)-sp.degree(sp.denom(e),s)<=-1
def leading(e):
    if not e or sp.degree(sp.numer(e),s)-sp.degree(sp.denom(e),s)<-1:return sp.S.Zero
    return sp.LC(sp.Poly(sp.numer(e),s))/sp.LC(sp.Poly(sp.denom(e),s))
Ninf=-Cinf.applyfunc(leading)
assert (Ninf.rank(),(Ninf**2).rank(),(Ninf**3).rank())==(2,1,0)
mat=lambda M:[[str(e) for e in row] for row in M.tolist()]
out={'status':'geometric_metric_Levi_Civita_and_all_local_types_verified',
     'polarized_system_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
     'metric_sha256':hashlib.sha256((ROOT/'data/banana_metric.json').read_bytes()).hexdigest(),
     'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'all_16_Levi_Civita_entries_verified':True,'physical_inverse_metric_matched':True,
     'finite_residues':{a:mat(R) for a,R in residues.items()},
     'node_minimal_polynomial':'z*(z+1/2)',
     'zero_projector':mat(E),'zero_shearing':mat(J),'zero_sheared_connection':mat(Jconn),
     'zero_sheared_denominator_factors':[[str(f),int(n)] for f,n in factors],
     'zero_nilpotent_residue':mat(N0),'zero_nilpotent_ranks':[2,0],
     'infinity_shearing_exponents':[1,1,1,1,2,-1],
     'infinity_nilpotent_residue':mat(Ninf),'infinity_nilpotent_ranks':[2,1,0],
     'four_parameter_flatness_verified':False,'global_transport_certified':False,
     'absolute_period_comparison_certified':False}
(folder/'polarized_geometry_check.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':out['status'],'zero_ranks':[2,0],'infinity_ranks':[2,1,0]}),flush=True)
