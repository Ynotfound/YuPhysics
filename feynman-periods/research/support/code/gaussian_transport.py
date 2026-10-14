"""Exact arithmetic for replayable complex ODE residual certificates.

Complex numbers are pairs of integers/rationals. No floating point enters
the certificate checks. Matrices are row-major lists of such pairs.
"""
from fractions import Fraction as F
from math import comb, lcm, ceil
from pathlib import Path
import hashlib, json
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / 'data/banana_pure_ray'
POLES = [0, 1, 9, 25, 49, 81, 121, 169, 225]
D_VALUES = ['-64', '-64', '-64', '-64', '4096', '-1/64']
Z = (0, 0)

def ga(a,b): return (a[0]+b[0], a[1]+b[1])
def gn(a): return (-a[0],-a[1])
def gm(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gs(a,b): return (a*b[0],a*b[1])
def ab(a): return abs(a[0])+abs(a[1])
def gp(a,n):
    v=(1,0)
    for _ in range(n): v=gm(v,a)
    return v
def mz(d=6): return [Z for _ in range(d*d)]
def mi(scale=1,d=6): return [(scale if i==j else 0,0) for i in range(d) for j in range(d)]
def ma(A,B): return [ga(a,b) for a,b in zip(A,B)]
def ms(a,A): return [gm(a,b) for b in A]
def mm(A,B,d=6):
    C=mz(d)
    for i in range(d):
        for k in range(d):
            a=A[d*i+k]
            if a==Z: continue
            for j in range(d):
                b=B[d*k+j]
                if b!=Z: C[d*i+j]=ga(C[d*i+j],gm(a,b))
    return C
def mnorm(A,d=6): return max(sum(ab(A[d*i+j]) for j in range(d)) for i in range(d))
def cnorm(A,d=6): return max(sum(ab(C[d*i+j]) for C in A for j in range(d)) for i in range(d))
def encode(A): return [[str(a),str(b)] for a,b in A]
def decode(A): return [(int(a),int(b)) for a,b in A]
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dyadic_upper(x):
    x=F(x)
    if x==0: return F(0)
    assert x>0
    k=x.numerator.bit_length()-x.denominator.bit_length()
    v=F(2)**k
    return v if v>=x else 2*v
def dyadic_round(x,bits):
    t=F(x)*(1<<bits)
    # Deterministic nearest rounding; exact error is measured separately.
    return (2*t.numerator+t.denominator)//(2*t.denominator)

def source_polynomials():
    """Read the actual preferred connection, not a numerical approximation."""
    s=sp.Symbol('ss')
    source=FOLDER/'polarized_system.json'
    B=sp.Matrix(json.loads(source.read_text(encoding='utf-8'))['matrix'])
    D=sp.diag(*map(sp.Rational,D_VALUES))
    q=sp.prod(s-a for a in POLES)
    N=(q*D*B*D.inv()).applyfunc(sp.cancel)
    polys=[sp.Poly(e,s) for e in N]
    cv=lambda v:F(int(sp.numer(v)),int(sp.denom(v)))
    qco=[cv(sp.Poly(q,s).nth(k)) for k in range(10)]
    deg=max(p.degree() for p in polys if not p.is_zero)
    nco=[[cv(p.nth(k)) for p in polys] for k in range(deg+1)]
    return qco,nco

def shift(poly,c,h):
    # p(c+h*t), allowing arbitrary Gaussian rational c and h.
    m=len(poly)
    cp=[gp(c,j) for j in range(m)];hp=[gp(h,j) for j in range(m)]
    return [gs(1,sumg(gs(poly[k]*comb(k,j),gm(cp[k-j],hp[j]))
                    for k in range(j,m))) for j in range(m)]
def sumg(values):
    out=Z
    for v in values: out=ga(out,v)
    return out

def local_polynomials(source,c,h):
    """dt equation: q(t)Y_t=N(t)Y, q,N Gaussian integral."""
    qraw,nraw=source
    q=shift(qraw,c,h)
    columns=[shift([M[j] for M in nraw],c,h) for j in range(36)]
    N=[[gm(h,columns[j][k]) for j in range(36)] for k in range(len(nraw))]
    scale=1
    for pair in q+[v for M in N for v in M]:
        for v in pair: scale=lcm(scale,F(v).denominator)
    cast=lambda a:tuple(int(F(v)*scale) for v in a)
    return list(map(cast,q)),[[cast(v) for v in M] for M in N],scale

def disk_lower(c,h,scale):
    # l1(h) >= |h|; linfty(c-a) <= |c-a|. All inequalities rational.
    radius=ab(h)
    distances=[max(abs(c[0]-a),abs(c[1]))-radius for a in POLES]
    assert all(v>0 for v in distances), 'disk may meet a pole'
    value=F(scale)
    for v in distances: value*=v
    return value,distances

def residual(q,N,P):
    R=[mz() for _ in range(max(len(q)+len(P)-2,len(N)+len(P)-1))]
    for j,a in enumerate(q):
        for k in range(1,len(P)):
            R[j+k-1]=ma(R[j+k-1],ms(gs(k,a),P[k]))
    for j,A in enumerate(N):
        for k,B in enumerate(P):
            R[j+k]=ma(R[j+k],ms((-1,0),mm(A,B)))
    return R

def bounds(q,N,P,bits,lower):
    assert P[0]==mi(1<<bits)
    M=F(cnorm(N))/lower
    defect=F(cnorm(residual(q,N,P)),1<<bits)/lower
    exponential=F(3)**ceil(M)
    error=exponential*defect
    endpoint=[sumg(C[j] for C in P) for j in range(36)]
    return endpoint,{'denominator_lower':lower,'matrix_norm_upper':M,
                      'residual_norm_upper':defect,'exp_upper':exponential,
                      'error_norm_upper':error}

def compose(L,le,C,ce,bits):
    """L,C have common dyadic scale. Radii use usual complex infinity norm."""
    S=1<<bits;product=mm(L,C)
    rounded=[(dyadic_round(F(v[0],S*S),bits),dyadic_round(F(v[1],S*S),bits)) for v in product]
    rd=max(sum(F(ab(ga(product[6*i+j],gn(gs(S,rounded[6*i+j])))),S*S)
               for j in range(6)) for i in range(6))
    err=(F(mnorm(L),S)+le)*ce+le*F(mnorm(C),S)+rd
    return rounded,dyadic_upper(err),rd

def import_first_step(bits):
    p=FOLDER/'first_trunk_step.json';d=json.loads(p.read_text(encoding='utf-8'))
    assert d['input_system_sha256']==sha(FOLDER/'polarized_system.json')
    C=[(F(a),F(b)) for row,icol in zip(d['endpoint_center_real'],d['endpoint_center_imag']) for a,b in zip(row,icol)]
    S=1<<bits;rounded=[(dyadic_round(a,bits),dyadic_round(b,bits)) for a,b in C]
    rd=max(sum(ab(ga(C[6*i+j],gn(gs(F(1,S),rounded[6*i+j])))) for j in range(6)) for i in range(6))
    return rounded,dyadic_upper(F(2)**d['error_bound_dyadic_exponent']+rd)
