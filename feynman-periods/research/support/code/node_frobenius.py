"""Exact local node model, Frobenius vectors and reflection error bounds."""
from gaussian_transport import *
from refined_transport import radius_up,save_json
import sys
sys.set_int_max_str_digits(100000)

def gd(x,y):
    n=y[0]*y[0]+y[1]*y[1]
    assert n!=0
    return (F(x[0]*y[0]+x[1]*y[1])/n,F(x[1]*y[0]-x[0]*y[1])/n)
def mv(A,v): return [sumg(gm(A[6*i+j],v[j]) for j in range(6)) for i in range(6)]
def va(a,b): return [ga(x,y) for x,y in zip(a,b)]
def vscale(a,v): return [gm(a,x) for x in v]
def vfdecode(v): return [(F(a),F(b)) for a,b in v]
def vnorm(v): return max(map(ab,v))
def vcoeff_norm(vs): return max(sum(ab(v[i]) for v in vs) for i in range(6))
def rat(v): return F(int(sp.numer(v)),int(sp.denom(v)))
def real_matrix(M): return [(rat(v),F(0)) for v in M]

def node_model(source,a):
    assert a in POLES[1:]
    s=sp.Symbol('ss');qraw,nraw=source
    qrest=sp.Poly(sp.prod(s-b for b in POLES if b!=a),s)
    qr=[rat(qrest.nth(k)) for k in range(qrest.degree()+1)]
    den=sum(c*F(a)**k for k,c in enumerate(qr))
    R=sp.Matrix(6,6,[sp.Rational(sum(M[j]*F(a)**k for k,M in enumerate(nraw)))/sp.Rational(den) for j in range(36)])
    assert R.rank()==1 and R*R==-R/2
    vectors=R.columnspace()+R.nullspace();assert len(vectors)==6
    Q=sp.Matrix.hstack(*(v/max(abs(x) for x in v) for v in vectors))
    Qi=Q.inv();Lambda=sp.diag(sp.Rational(-1,2),0,0,0,0,0)
    assert Qi*R*Q==Lambda
    # Synthetic division of N(s)-R qrest(s) by s-a, entry by entry.
    polys=[]
    for j in range(36):
        f=[M[j]-(rat(R[j])*qr[k] if k<len(qr) else 0) for k,M in enumerate(nraw)]
        out=[F(0)]*(len(f)-1);out[-1]=f[-1]
        for k in range(len(f)-2,0,-1): out[k-1]=f[k]+a*out[k]
        assert f[0]+a*out[0]==0
        polys.append(out)
    transformed=[]
    for k in range(len(polys[0])):
        coefficient=Qi*sp.Matrix(6,6,[sp.Rational(p[k]) for p in polys])*Q
        transformed.append([rat(v) for v in coefficient])
    return Q,Qi,R,(qr,transformed)

def node_disk_lower(a,h,scale):
    ds=[F(abs(a-b))-ab(h) for b in POLES if b!=a]
    assert min(ds)>0
    lower=F(scale)
    for x in ds: lower*=x
    return lower,ds

D0=[F(0)]+[F(1,2)]*5

def exact_coefficients(q,N,order):
    assert q[0][1]==0 and q[0][0]!=0
    P=[[(F(i==0),F(0)) for i in range(6)]]
    for n in range(1,order+1):
        v=[Z]*6
        for j in range(min(n-1,len(N)-1)+1): v=va(v,mv(N[j],P[n-1-j]))
        for j in range(1,min(n,len(q)-1)+1):
            v=va(v,[gn(gm(q[j],gs(F(n-j)-D0[i],P[n-j][i]))) for i in range(6)])
        P.append([gs(1/(F(q[0][0])*(F(n)-D0[i])),v[i]) for i in range(6)])
    return P

def fuchs_residual(q,N,P):
    out=[[Z]*6 for _ in range(max(len(q)+len(P)-1,len(N)+len(P)))]
    for j,qj in enumerate(q):
        for n,v in enumerate(P):
            out[j+n]=va(out[j+n],[gm(qj,gs(F(n)-D0[i],v[i])) for i in range(6)])
    for j,Nj in enumerate(N):
        for n,v in enumerate(P): out[j+n+1]=va(out[j+n+1],vscale((-1,0),mv(Nj,v)))
    return out

def fuchs_bounds(q,N,P,lower):
    m=len(P);res=fuchs_residual(q,N,P)
    assert all(v==[Z]*6 for v in res[:m])
    M=F(cnorm(N))/lower;delta=vcoeff_norm(res[m:])/lower
    kappa=M/(F(m)+F(1,2));assert kappa<1
    error=delta/((F(m)-F(1,2))*(1-kappa))
    return {'vanishing_degree':m,'denominator_lower':lower,'M':M,'delta':delta,'kappa':kappa,'vector_error_Q_basis':error}

def pairing_at(a,h):
    s=sp.Symbol('ss');d=json.loads((FOLDER/'polarized_system.json').read_text())
    K=sp.Matrix(d['inverse_pairing']);D=sp.diag(*map(sp.Rational,D_VALUES))
    point=sp.Rational(a)+sp.Rational(h[0])+sp.I*sp.Rational(h[1])
    val=(D.inv()*K*D.inv()).subs(s,point).applyfunc(sp.expand)
    return [(rat(sp.re(v)),rat(sp.im(v))) for v in val]

def reflection_from_vector(P,Q,error,K,bits):
    S=1<<bits
    pc=[sumg(v[i] for v in P) for i in range(6)]
    exact_center=mv(real_matrix(Q),pc)
    Vint=[(dyadic_round(a,bits),dyadic_round(b,bits)) for a,b in exact_center]
    V=[gs(F(1,S),v) for v in Vint]
    rho=vnorm([ga(a,gn(b)) for a,b in zip(exact_center,V)])
    ev=radius_up(F(mnorm(real_matrix(Q)))*error+rho,bits)
    cov=[sumg(gm(V[i],K[6*i+j]) for i in range(6)) for j in range(6)]
    dc=sumg(gm(cov[i],V[i]) for i in range(6))
    B=[gm(V[i],cov[j]) for i in range(6) for j in range(6)]
    k=mnorm(K);u=vnorm(V);Delta=6*k*ev*(2*u+ev)
    d0=max(abs(dc[0]),abs(dc[1]));low=d0-Delta
    assert low>0, 'pairing denominator not separated from zero'
    error_ref=2*(Delta/low+mnorm(B)*Delta/(d0*low))
    exact_ref=[ga(x,gs(-2,gd(y,dc))) for x,y in zip(mi(),B)]
    Tint=[(dyadic_round(a,bits),dyadic_round(b,bits)) for a,b in exact_ref]
    rhor=mnorm([ga(v,gn(gs(F(1,S),w))) for v,w in zip(exact_ref,Tint)])
    er=radius_up(error_ref+rhor,bits)
    return {'bits':bits,'vector_center_gaussian_integer':encode(Vint),'vector_rounding_error':str(rho),
            'vector_error_upper':str(ev),'pairing_center':encode([dc])[0],
            'pairing_error_upper':str(Delta),'pairing_modulus_lower':str(low),
            'reflection_center_gaussian_integer':encode(Tint),'reflection_rounding_error':str(rhor),
            'reflection_error_upper':str(er)}
