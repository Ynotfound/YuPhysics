"""Exact rational compatibility on the whole mass ray, plus complete root boxes.

No numerical transport is performed. The inverse pairing check avoids taking
the inverse of a large rational matrix. Root isolation uses exact QQ arithmetic.
"""
from pathlib import Path
import json,hashlib
import sympy as sp
from sympy.parsing.mathematica import parse_mathematica

ROOT=Path(__file__).resolve().parents[1]
folder=ROOT/'data/banana_pure_ray'
source=ROOT/'external/bonn-banana-2025/ZCYInverse.m'
system=folder/'system.json'
s=sp.Symbol('ss')
x=sp.symbols('m12 m22 m32 m42')
A=sp.Matrix(json.loads(system.read_text(encoding='utf-8'))['matrix'])
K=sp.Matrix(parse_mathematica(source.read_text(encoding='utf-8')))
K=K.subs(dict(zip(x,[sp.Rational(m)/s for m in [1,4,16,64]]))).applyfunc(sp.cancel)
Ti=sp.diag(1,1,1,1,1,sp.Rational(1,2))
KH=(s**2*Ti*K*Ti).applyfunc(sp.cancel)
residual=(KH.diff(s)+A.T*KH+KH*A).applyfunc(sp.cancel)
assert residual==sp.zeros(6)
assert KH==KH.T
assert KH.subs(s,-64).det()!=0
print('global_rational_pairing_compatibility_passed',flush=True)
data=json.loads(system.read_text(encoding='utf-8'))
q=next(sp.sympify(p) for p,n in data['common_denominator_factors'] if sp.degree(sp.sympify(p),s)==5)
poly=sp.Poly(q,s,domain=sp.QQ)
real,complexes=poly.intervals(all=True,eps=sp.Rational(1,10**6))
assert len(real)==3 and len(complexes)==2
boxes=[]
for (lo,hi),multiplicity in real:
    assert multiplicity==1 and poly.count_roots(lo,hi)==1
    boxes.append({'kind':'real','re':[str(lo),str(hi)],'im':['0','0'],'count':1})
for (lo,hi),multiplicity in complexes:
    assert multiplicity==1 and poly.count_roots(lo,hi)==1
    boxes.append({'kind':'nonreal','re':[str(sp.re(lo)),str(sp.re(hi))],
                  'im':[str(sp.im(lo)),str(sp.im(hi))],'count':1})
for i,b in enumerate(boxes):
    for c in boxes[:i]:
        disjoint=any(sp.Rational(b[axis][1])<sp.Rational(c[axis][0]) or
                     sp.Rational(c[axis][1])<sp.Rational(b[axis][0]) for axis in ['re','im'])
        assert disjoint
    for rational_pole in [0,1,9,25,49,81,121,169,225]:
        assert not(sp.Rational(b['re'][0])<=rational_pole<=sp.Rational(b['re'][1]) and
                   sp.Rational(b['im'][0])<=0<=sp.Rational(b['im'][1]))
out={'status':'ray_pairing_identity_and_all_finite_pole_boxes_verified',
     'system_sha256':hashlib.sha256(system.read_bytes()).hexdigest(),
     'source_pairing_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
     'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'sympy_version':sp.__version__,
     'inverse_pairing_basis':'H=s^(-1)*(phi,partial_xi phi,partial_x1^2 phi)',
     'inverse_pairing':[[str(v) for v in row] for row in KH.tolist()],
     'identity':'d_s K_H + A_s^t K_H + K_H A_s = 0',
     'identity_scope':'all rational functions on the mass ray, 36 entries',
     'q5':str(q),'q5_boxes':boxes,'q5_root_count':5,
     'rational_finite_poles':[0,1,9,25,49,81,121,169,225],
     'root_isolation_method':'SymPy exact rational real/complex isolation and exact rectangle root counts',
     'path_clearance_certified':False,'four_parameter_flatness_certified':False,
     'global_transport_certified':False,'absolute_cup_scale_certified':False}
(folder/'pairing_and_roots.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':out['status'],'boxes':boxes}),flush=True)
