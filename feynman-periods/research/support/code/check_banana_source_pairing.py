"""Exact pointwise compatibility of source connection and full source pairing.

This checks a finite algebraic identity at the specified point, not global
flatness and not numerical continuation. It also resolves the last-row factor.
"""
from pathlib import Path
import json,hashlib
import sympy as sp
from sympy.parsing.mathematica import parse_mathematica

ROOT=Path(__file__).resolve().parents[1]
p=ROOT/'external/bonn-banana-2025/ZCYInverse.m'
K=sp.Matrix(parse_mathematica(p.read_text(encoding='utf-8')))
x=sp.symbols('m12 m22 m32 m42')
point=dict(zip(x,[sp.Rational(-1,64),sp.Rational(-1,16),sp.Rational(-1,4),sp.Integer(-1)]))
K0=K.subs(point).applyfunc(sp.simplify)
assert K0==K0.T and all(v.is_Rational for v in K0)
Z0=K0.inv()
data=json.loads((ROOT/'data/banana_pure_ray/A_x1.json').read_text(encoding='utf-8'))
A=sp.Matrix(data['matrix']).subs(sp.Symbol('ss'),-64)
dK=K.diff(x[0]).subs(point).applyfunc(sp.simplify)
assert all(v.is_Rational for v in dK)
dZ=-Z0*dK*Z0
matches=[]
for factor in [sp.Rational(1),sp.Rational(2),sp.Rational(1,2)]:
    T=sp.diag(1,1,1,1,1,factor)
    Z=T*Z0*T;derivative=T*dZ*T
    if derivative==A*Z+Z*A.T:matches.append(str(factor))
assert len(matches)==1,matches
out={'status':'pointwise_source_connection_pairing_compatibility_passed',
     'point_x':['-1/64','-1/16','-1/4','-1'],'direction':'partial_x1',
     'last_row_rescaling_for_W_basis':matches[0],
     'identity':'partial_x1 Z_W = A_x1 Z_W + Z_W A_x1^t',
     'source_file_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
     'connection_source_sha256':data['source_sha256'],
     'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'global_flatness_certified':False,'absolute_cup_scale_certified':False}
(ROOT/'data/banana_source_pairing_check.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out),flush=True)
