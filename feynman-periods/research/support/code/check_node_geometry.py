"""Exact Morse Hessians at the eight interior banana critical points."""
from pathlib import Path
import sympy as s
import itertools,json,hashlib

ROOT=Path(__file__).resolve().parents[1]
x=s.symbols('x1:4');r=s.symbols('r1:5');S=sum(r)
f=(sum(r[i]*x[i] for i in range(3))+r[3])*(sum(r[i]/x[i] for i in range(3))+r[3])
at={v:1 for v in x}
gradient=s.Matrix([s.diff(f,v).subs(at) for v in x])
assert gradient==s.zeros(3,1)
H=s.hessian(f,x).subs(at)
expected=2*(S*s.diag(*r[:3])-s.Matrix(r[:3])*s.Matrix(1,3,r[:3]))
assert H==expected
assert s.factor(H.det()-8*S**2*s.prod(r))==0
records=[]
for signs in itertools.product([-1,1],repeat=3):
    signed=[e*m for e,m in zip(signs,[1,2,4])]+[8]
    total=sum(signed);a=total**2
    point=[str(s.Rational(8*e,m)) for e,m in zip(signs,[1,2,4])]+['1']
    det=int(H.det().subs(dict(zip(r,signed))))
    assert total>0 and det!=0
    records.append({'node':a,'signed_masses':signed,'torus_projective_point':point,
                    'critical_value':a,'logarithmic_hessian_determinant':det})
records.sort(key=lambda d:d['node'])
assert [d['node'] for d in records]==[1,9,25,49,81,121,169,225]
out={'status':'eight_unique_transverse_interior_Morse_nodes_checked',
     'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'Hessian_formula':'2*(S*diag(r1,r2,r3)-r*r^t)',
     'Hessian_determinant':'8*S^2*r1*r2*r3*r4',
     'transverse_parameter_derivative':-1,'nodes':records,
     'geometric_topology_not_certified_by_this_program':
         'Picard-Lefschetz integral vanishing sphere, self-intersection -2, and orthogonality to fixed boundary classes are proved in the companion text.'}
p=ROOT/'data/banana_pure_ray/node_branches/morse_nodes.json'
p.parent.mkdir(exist_ok=True)
p.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out))
