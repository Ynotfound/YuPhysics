"""Certify the full first trunk segment -64 -> -64+i by a rational residual.

All polynomial witnesses and all bounds use exact rational arithmetic.
No floating-point value is used as an error certificate. The current driver
handles this one real-centred step; it is not yet a complete path integrator.
"""
from pathlib import Path
from fractions import Fraction as F
from math import ceil
import json,hashlib,time
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
folder=ROOT/'data/banana_pure_ray'

def zero(d):return [[F(0) for j in range(d)] for i in range(d)]
def eye(d):return [[F(i==j) for j in range(d)] for i in range(d)]
def add(A,B):return [[a+b for a,b in zip(u,v)] for u,v in zip(A,B)]
def scale(a,A):return [[a*x for x in row] for row in A]
def mul(A,B):
    d=len(A)
    return [[sum((A[i][k]*B[k][j] for k in range(d)),F(0)) for j in range(d)] for i in range(d)]
def coeff_norm(matrices,r):
    d=len(matrices[0])
    return max(sum(sum(abs(M[i][j])*r**k for k,M in enumerate(matrices)) for j in range(d)) for i in range(d))

def witness(q,N,order,r,q_lower):
    d=len(N[0]);A=[]
    for n in range(order):
        v=N[n] if n<len(N) else zero(d)
        for j in range(1,min(n,len(q)-1)+1):v=add(v,scale(-q[j],A[n-j]))
        A.append(scale(1/q[0],v))
    P=[eye(d)]
    for n in range(order):
        v=zero(d)
        for k in range(n+1):v=add(v,mul(A[k],P[n-k]))
        P.append(scale(F(1,n+1),v))
    residual=[zero(d) for _ in range(max(len(q)+len(P)-2,len(N)+len(P)-1))]
    for j,a in enumerate(q):
        for k in range(1,len(P)):residual[j+k-1]=add(residual[j+k-1],scale(a*k,P[k]))
    for j,M in enumerate(N):
        for k,Pc in enumerate(P):residual[j+k]=add(residual[j+k],scale(-1,mul(M,Pc)))
    assert all(R==zero(d) for R in residual[:order])
    M=coeff_norm(N,r)/q_lower
    residual_upper=coeff_norm(residual,r)/q_lower
    exponential_upper=F(3)**ceil(M*r)
    error=r*exponential_upper*residual_upper
    return P,residual,{'radius':r,'denominator_lower':q_lower,'matrix_norm_upper':M,
                       'residual_norm_upper':residual_upper,'exp_upper':exponential_upper,
                       'error_norm_upper':error}

def at_imag(P,h):
    real=zero(len(P[0]));imag=zero(len(P[0]))
    for k,M in enumerate(P):
        if k%2:imag=add(imag,scale((-1)**((k-1)//2)*h**k,M))
        else:real=add(real,scale((-1)**(k//2)*h**k,M))
    return real,imag

def self_test():
    # Exact solution y=1/(1-z), tested at z=i/2, not an implementation-shaped fixture.
    P,_,bound=witness([F(1),F(-1)],[[[F(1)]]],4,F(1,2),F(1,2))
    re,im=at_imag(P,F(1,2));error=bound['error_norm_upper']
    assert abs(re[0][0]-F(4,5))<=error and abs(im[0][0]-F(2,5))<=error
    assert error>0

def main():
    self_test()
    source=folder/'polarized_system.json'
    data=json.loads(source.read_text(encoding='utf-8'))
    s,z=sp.symbols('ss z')
    B=sp.Matrix(data['matrix'])
    D=sp.diag(-64,-64,-64,-64,4096,sp.Rational(-1,64))
    A=D*B*D.inv()
    poles=[0,1,9,25,49,81,121,169,225]
    q=sp.prod(s-a for a in poles)
    N=(q*A).applyfunc(sp.cancel)
    assert all(sp.denom(e).is_Integer for e in N)
    qp=sp.Poly(q.subs(s,z-64),z)
    np=[[sp.Poly(sp.expand(N[i,j].subs(s,z-64)),z) for j in range(6)] for i in range(6)]
    frac=lambda v:F(int(sp.numer(v)),int(sp.denom(v)))
    qc=[frac(qp.nth(k)) for k in range(qp.degree()+1)]
    degree=max(p.degree() for row in np for p in row if not p.is_zero)
    nc=[[[frac(np[i][j].nth(k)) for j in range(6)] for i in range(6)] for k in range(degree+1)]
    radius=F(1);lower=F(1)
    for a in poles:lower*=abs(-64-a)-radius
    assert lower>0
    order=20
    started=time.time()
    P,residual,bound=witness(qc,nc,order,radius,lower)
    print(f'exact_polynomial_witness_seconds_{time.time()-started:.2f}',flush=True)
    real,imag=at_imag(P,F(1))
    error=bound['error_norm_upper']
    assert 0<error<F(1,10**20)
    k=error.denominator.bit_length()-error.numerator.bit_length()
    dyadic=lambda k:F(1,2**k) if k>=0 else F(2**(-k))
    while dyadic(k)<error:k-=1
    while dyadic(k+1)>=error:k+=1
    enc=lambda M:[[str(v) for v in row] for row in M]
    out={'status':'first_complete_trunk_segment_certified_by_exact_residual',
         'input_system_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
         'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'basis':'constant D times polarized frame','D':['-64','-64','-64','-64','4096','-1/64'],
         'start':['-64','0'],'end':['-64','1'],'initial_matrix':'I_6','order':order,
         'denominator_coefficients_ascending':list(map(str,qc)),
         'numerator_matrix_coefficients_ascending':[enc(M) for M in nc],
         'witness_coefficients_ascending':[enc(M) for M in P],
         'exact_residual_zero_below_degree':order,
         'bounds':{a:str(b) for a,b in bound.items()},
         'endpoint_center_real':enc(real),'endpoint_center_imag':enc(imag),
         'uniform_entry_real_and_imag_error_bound':f'2^(-{k})',
         'error_bound_dyadic_exponent':-k,
         'error_includes_all_truncation':True,'floating_point_used_in_certificate':False,
         'global_loop_transport_certified':False,'absolute_initial_periods_certified':False,
         'independent_exact_solution_self_test_passed':True}
    (folder/'first_trunk_step.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':out['status'],'error_at_most':out['uniform_entry_real_and_imag_error_bound'],
                      'matrix_norm_upper_approx':float(bound['matrix_norm_upper']),
                      'approximation_is_display_only':True}),flush=True)

if __name__=='__main__':main()
