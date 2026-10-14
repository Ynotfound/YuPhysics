"""Nilpotent Frobenius gauge after the exact resonant shearing at zero."""
from node_frobenius import *
from certify_horizontal_landmarks import inverse


def zero_model():
    s=sp.Symbol('ss');p=FOLDER/'polarized_geometry_check.json'
    data=json.loads(p.read_text());B0=sp.Matrix(data['zero_sheared_connection'])
    J=sp.Matrix(data['zero_shearing']);N0=sp.Matrix(data['zero_nilpotent_residue'])
    Q,N= N0.jordan_form()
    # A common scale on each Jordan chain leaves the unit superdiagonal intact.
    scales=[];k=0
    while k<6:
        size=2 if k<5 and N[k,k+1]==1 else 1
        scale=max(abs(v) for v in Q[:,k:k+size])
        scales.extend([1/scale]*size);k+=size
    Q=Q*sp.diag(*scales);Qi=Q.inv()
    assert Qi*N0*Q==N and N*N==sp.zeros(6) and N.rank()==2
    assert mnorm(real_matrix(N))==1
    A=(Qi*B0*Q).applyfunc(sp.cancel)
    q=sp.prod(s-a for a in POLES)
    num=(q*A).applyfunc(sp.cancel)
    pol=[sp.Poly(v,s) for v in num];degree=max(v.degree() for v in pol if not v.is_zero)
    source=([rat(sp.Poly(q,s).nth(k)) for k in range(10)],
            [[rat(v.nth(k)) for v in pol] for k in range(degree+1)])
    qrest=sp.cancel(q/s);regnum=((num-N*qrest)/s).applyfunc(sp.cancel)
    pol=[sp.Poly(v,s) for v in regnum];degree=max(v.degree() for v in pol if not v.is_zero)
    regular=([rat(sp.Poly(qrest,s).nth(k)) for k in range(9)],
             [[rat(v.nth(k)) for v in pol] for k in range(degree+1)])
    D=sp.diag(*map(sp.Rational,D_VALUES))
    to_D=(D*J.inv()*Q).applyfunc(sp.cancel)
    return Q,Qi,N,source,regular,to_D


def comm(N,X): return ma(mm(N,X),ms((-1,0),mm(X,N)))


def inv_euler(n,N,X):
    C=comm(N,X);CC=comm(N,C)
    return ma(ma(ms((F(1,n),0),X),ms((F(1,n*n),0),C)),ms((F(1,n**3),0),CC))


def exact_gauge_coefficients(q,H,N,order):
    assert q[0][1]==0 and q[0][0]!=0
    P=[mi()]
    for n in range(1,order+1):
        rhs=mz()
        for j in range(min(n-1,len(H)-1)+1): rhs=ma(rhs,mm(H[j],P[n-1-j]))
        for j in range(1,min(n,len(q)-1)+1):
            term=ma(ms((n-j,0),P[n-j]),ms((-1,0),comm(N,P[n-j])))
            rhs=ma(rhs,ms(gn(q[j]),term))
        P.append(ms((F(1,q[0][0]),0),inv_euler(n,N,rhs)))
    return P


def gauge_residual(q,H,N,P):
    R=[mz() for _ in range(max(len(q)+len(P)-1,len(H)+len(P)))]
    for j,qj in enumerate(q):
        for n,M in enumerate(P):
            R[j+n]=ma(R[j+n],ms(qj,ma(ms((n,0),M),ms((-1,0),comm(N,M)))))
    for j,Hj in enumerate(H):
        for n,M in enumerate(P): R[j+n+1]=ma(R[j+n+1],ms((-1,0),mm(Hj,M)))
    return R


def Jbound(m): return F(1,m)+F(2,m*m)+F(2,m**3)


def gauge_bounds(q,H,N,P,lower):
    m=len(P);R=gauge_residual(q,H,N,P)
    assert all(M==mz() for M in R[:m])
    M=F(cnorm(H))/lower;delta=F(cnorm(R[m:]))/lower
    kappa=M*Jbound(m+1);assert kappa<1
    error=delta*Jbound(m)/(1-kappa)
    return {'vanishing_degree':m,'denominator_lower':lower,'regular_norm_upper':M,
            'weighted_defect_upper':delta,'kappa':kappa,'gauge_error_upper':error}


def endpoint_gauge(P,error,bits):
    S=1<<bits;exact=[sumg(M[j] for M in P) for j in range(36)]
    C=[(dyadic_round(x,bits),dyadic_round(y,bits)) for x,y in exact]
    rd=mnorm([ga(v,gn(gs(F(1,S),w))) for v,w in zip(exact,C)])
    return C,radius_up(error+rd,bits),rd
