"""Change to the geometric polarized frame, removing derivative-chart poles.

The inverse metric constructs the unique isotropic complement algebraically.
All operations are in QQ(s); no numerical evidence is used for cancellation.
"""
from pathlib import Path
import json,hashlib,time
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
folder=ROOT/'data/banana_pure_ray'
s=sp.Symbol('ss')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
system=folder/'system.json';pairing=folder/'pairing_and_roots.json'
A=sp.Matrix(json.loads(system.read_text(encoding='utf-8'))['matrix'])
KH=sp.Matrix(json.loads(pairing.read_text(encoding='utf-8'))['inverse_pairing'])
mu=[1,4,16,64]
simp=lambda M:M.applyfunc(sp.cancel)
R=sp.zeros(6)
R[0,0]=1
for i in range(1,4):R[i,i]=mu[i-1]/s
R[4,0]=-1/s
for i in range(1,5):R[4,i]=-mu[i-1]/s**2
R[5,5]=1
Ri=R.inv()
# The dual of omega in the temporary basis (omega,v1,v2,v3,v4,H5).
dual=simp((Ri[:,0].T*KH)*Ri)
c=dual.copy();c[0,0]/=2
P=R.copy();P[5,:]=simp(c*R)
Ti=sp.eye(6)
for j in range(5):Ti[5,j]=-c[0,j]/c[0,5]
Ti[5,5]=1/c[0,5]
Pi=simp(Ri*Ti)
assert simp(P*Pi)==sp.eye(6)
print('polarized_gauge_and_inverse_constructed',flush=True)
KB=simp(Pi.T*KH*Pi)
assert KB[0,0]==KB[5,5]==0 and KB[0,5]==KB[5,0]==1
assert all(KB[i,j]==0 for i in [0,5] for j in range(1,5))
assert KB==KB.T
rows=[]
for i in range(6):
    row=simp((P[i,:].diff(s)+P[i,:]*A)*Pi)
    rows.append(list(row))
    print(f'connection_row_{i}_completed',flush=True)
B=sp.Matrix(rows)
den=sp.lcm([sp.denom(e) for e in B])
factors=sp.factor_list(den)[1]
q=next(sp.sympify(p) for p,n in json.loads(system.read_text(encoding='utf-8'))['common_denominator_factors'] if sp.degree(sp.sympify(p),s)==5)
assert sp.gcd(den,q)==1
assert B[0,:]==sp.eye(6)[4,:]
assert simp(KB.diff(s)+B.T*KB+KB*B)==sp.zeros(6)
out={'status':'polarized_ray_gauge_exactly_verified_and_apparent_q5_removed',
     'basis':'(omega,mu1*d_mu1 omega,mu2*d_mu2 omega,mu3*d_mu3 omega,d_s omega,w)',
     'w_conditions':'<omega,w>=1, <v_a,w>=0, <w,w>=0 in published normalization',
     'old_system_sha256':sha(system),'old_pairing_sha256':sha(pairing),
     'code_sha256':sha(Path(__file__)),
     'gauge_P':[[str(v) for v in row] for row in P.tolist()],
     'inverse_gauge':[[str(v) for v in row] for row in Pi.tolist()],
     'matrix':[[str(v) for v in row] for row in B.tolist()],
     'inverse_pairing':[[str(v) for v in row] for row in KB.tolist()],
     'common_denominator_factors':[[str(f),int(n)] for f,n in factors],
     'apparent_q5_removed':True,'polarization_identity_verified':True,
     'absolute_cup_scale_certified':False,'global_transport_certified':False}
(folder/'polarized_system.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':out['status'],'common_denominator_factors':out['common_denominator_factors']}),flush=True)
