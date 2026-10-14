"""Check connection reconstruction against an explicit flat null-cone lift.

This is an independent sign/factor example, NOT the banana metric Y.
"""
from pathlib import Path
import hashlib
import json
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
x = sp.symbols('x1:5', nonzero=True)
n = 4
ambient = sp.zeros(6)
ambient[0,5] = ambient[5,0] = 1
ambient[1:5,1:5] = sp.eye(4)
omega = x[0]*sp.Matrix([1, *x, -sum(t*t for t in x)/2])
v = [omega.diff(t) for t in x]

def pair(a,b):
    return sp.simplify((a.T*ambient*b)[0])

Y = sp.Matrix(n,n,lambda a,b: pair(v[a],v[b]))
assert Y == x[0]**2*sp.eye(n)
Yi = Y.inv()
tau = sp.zeros(6,1)
for a in range(n):
    for b in range(n):
        tau -= Yi[a,b]*v[b].diff(x[a])/n
assert pair(omega,tau) == 1
coeffs = Yi*sp.Matrix([pair(z,tau) for z in v])
w = tau - sp.Matrix.hstack(*v)*coeffs
w = sp.simplify(w - pair(w,w)*omega/2)
assert pair(omega,w) == 1 and pair(w,w) == 0
assert all(pair(z,w) == 0 for z in v)

Gamma = [[[sp.simplify(sum(Yi[c,d]*(sp.diff(Y[b,d],x[a])
           +sp.diff(Y[a,d],x[b])-sp.diff(Y[a,b],x[d]))/2 for d in range(n)))
           for c in range(n)] for b in range(n)] for a in range(n)]
def curvature(d,c,a,b):
    return sp.simplify(sp.diff(Gamma[b][c][d],x[a])-sp.diff(Gamma[a][c][d],x[b])
        +sum(Gamma[b][c][e]*Gamma[a][e][d]-Gamma[a][c][e]*Gamma[b][e][d] for e in range(n)))

Ric = sp.Matrix(n,n,lambda c,b: sp.simplify(sum(curvature(a,c,a,b) for a in range(n))))
scalar = sp.simplify(sp.trace(Yi*Ric))
C = sp.simplify(-(Ric - scalar*Y/6)/2)
assert C == sp.diag(-sp.Rational(3,2)/x[0]**2,
                    sp.Rational(1,2)/x[0]**2,
                    sp.Rational(1,2)/x[0]**2,
                    sp.Rational(1,2)/x[0]**2)
for a in range(n):
    for b in range(n):
        rhs = sum((Gamma[a][b][c]*v[c] for c in range(n)), sp.zeros(6,1))
        rhs += C[a,b]*omega-Y[a,b]*w
        assert sp.simplify(v[b].diff(x[a])-rhs) == sp.zeros(6,1)
    rhs = -sum((Yi[b,c]*C[a,c]*v[b] for b in range(n) for c in range(n)), sp.zeros(6,1))
    assert sp.simplify(w.diff(x[a])-rhs) == sp.zeros(6,1)

output = {'status':'connection_reconstruction_sign_and_factor_check_passed',
          'model':'omega=x1*(1,x1,x2,x3,x4,-sum(x_i^2)/2) in flat six-dimensional quadric space',
          'metric':'Y=x1^2 I4',
          'C_diagonal':[str(C[i,i]) for i in range(n)],
          'frame_gram_verified':True,'all_twenty_vector_derivatives_verified':True,
          'banana_metric_computed':False,'absolute_periods_computed':False,
          'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(ROOT/'data/polarized_connection_check.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':output['status'],'banana_metric_computed':False}))
