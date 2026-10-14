"""Replay a saved rational step witness without regenerating its Taylor series."""
from pathlib import Path
import json,hashlib
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
folder=ROOT/'data/banana_pure_ray'
p=folder/'first_trunk_step.json'
d=json.loads(p.read_text(encoding='utf-8'))
system=folder/'polarized_system.json'
assert d['input_system_sha256']==hashlib.sha256(system.read_bytes()).hexdigest()
assert d['code_sha256']==hashlib.sha256((ROOT/'code/certify_banana_ordinary_step.py').read_bytes()).hexdigest()
s,z=sp.symbols('ss z');rat=sp.Rational
assert d['start']==['-64','0'] and d['end']==['-64','1'] and d['initial_matrix']=='I_6'
D=sp.diag(*map(rat,d['D']));B=sp.Matrix(json.loads(system.read_text())['matrix'])
q=sum(rat(v)*z**k for k,v in enumerate(d['denominator_coefficients_ascending']))
N=sum((sp.Matrix(c)*z**k for k,c in enumerate(d['numerator_matrix_coefficients_ascending'])),sp.zeros(6))
P=sum((sp.Matrix(c)*z**k for k,c in enumerate(d['witness_coefficients_ascending'])),sp.zeros(6))
assert P.subs(z,0)==sp.eye(6)
assert (N-q*(D*B*D.inv()).subs(s,z-64)).applyfunc(sp.cancel)==sp.zeros(6)
assert sp.expand(q-sp.prod(z-64-a for a in [0,1,9,25,49,81,121,169,225]))==0
r=rat(d['bounds']['radius']);assert r==1
lower=sp.prod(64+a-r for a in [0,1,9,25,49,81,121,169,225])
assert lower>0
def norm_coeff(M):
    return max(sum(sum(abs(coef)*r**mon[0] for mon,coef in sp.Poly(M[i,j],z).terms()) for j in range(6)) for i in range(6))
residual=(q*P.diff(z)-N*P).applyfunc(sp.expand)
for e in residual:
    pol=sp.Poly(e,z)
    assert all(coef==0 or mon[0]>=d['exact_residual_zero_below_degree'] for mon,coef in pol.terms())
M=norm_coeff(N)/lower;res=norm_coeff(residual)/lower
error=r*3**sp.ceiling(M*r)*res
assert lower==rat(d['bounds']['denominator_lower'])
assert M==rat(d['bounds']['matrix_norm_upper'])
assert res==rat(d['bounds']['residual_norm_upper'])
assert error==rat(d['bounds']['error_norm_upper'])
endpoint=P.subs(z,sp.I).applyfunc(sp.expand)
assert endpoint==sp.Matrix(d['endpoint_center_real'])+sp.I*sp.Matrix(d['endpoint_center_imag'])
assert 0<error<=sp.Rational(2)**d['error_bound_dyadic_exponent']
out={'status':'independent_polynomial_residual_replay_passed',
     'certificate_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
     'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'Taylor_coefficient_generator_called':False,
     'source_rational_matrix_matched':True,'bounds_recomputed_exactly':True,
     'endpoint_real_imag_radius':d['uniform_entry_real_and_imag_error_bound'],
     'global_loop_certified':False}
(folder/'first_trunk_step_replay.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out),flush=True)
