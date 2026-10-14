"""Independent known-solution tests and exact geometry checks for transport."""
from gaussian_transport import *
from certify_high_trunk import proposal
from verify_high_trunk import check_residual
import sympy as sp

def known_solutions():
    s,t=sp.symbols('s t');c=2+3*sp.I;h=sp.Rational(1,8)-sp.I/16
    Fm=sp.eye(6);Fm[:2,:2]=sp.Matrix([[1,s],[s,1+s*s]])
    A=Fm.diff(s)*Fm.inv()
    exact=(Fm.subs(s,c+h*t)*Fm.subs(s,c).inv()).applyfunc(sp.expand)
    Nt=(h*A.subs(s,c+h*t)).applyfunc(sp.expand)
    vals=[sp.Poly(v,t) for v in Nt];degree=max(p.degree() for p in vals if not p.is_zero)
    pairs=[];scale=1
    for k in range(degree+1):
        row=[]
        for p in vals:
            re,im=sp.expand(p.nth(k)).as_real_imag();pair=(F(re),F(im))
            row.append(pair)
            scale=lcm(scale,pair[0].denominator,pair[1].denominator)
        pairs.append(row)
    N=[[tuple(int(scale*v) for v in p) for p in row] for row in pairs];q=[(scale,0)]
    bits=240;P=proposal(q,N,4,bits)
    endpoint,b=bounds(q,N,P,bits,F(scale))
    assert b['error_norm_upper']==0
    expected=exact.subs(t,1)
    for v,known in zip(endpoint,expected):
        re,im=sp.expand(known).as_real_imag()
        assert (F(v[0],1<<bits),F(v[1],1<<bits))==(F(re),F(im))
    # A different exact solution with a genuine pole and nonzero truncation:
    # y'=-h/(i+h*t)y, h=(2+i)/8; y(1)=8i/(2+9i).
    q=[(0,8),(2,1)];N=[mz()];N[0][0]=(-2,-1)
    P=proposal(q,N,16,bits)
    endpoint,b=bounds(q,N,P,bits,F(5))
    re=F(72,85);im=F(16,85)
    dre=F(endpoint[0][0],1<<bits)-re;dim=F(endpoint[0][1],1<<bits)-im
    assert 0<b['error_norm_upper']<F(1,1000)
    assert dre*dre+dim*dim<=b['error_norm_upper']**2
    assert b==check_residual(q,N,P,bits,F(5))
    rejected=False
    try: disk_lower((F(1),F(0)),(F(1,4),F(0)),1)
    except AssertionError: rejected=True
    assert rejected
    return {'complex_noncommuting_polynomial_solution':'exactly reproduced',
            'complex_rational_pole_solution':'exact error within residual bound',
            'pole_intersecting_disk':'rejected'}

def main():
    tests=known_solutions()
    system=FOLDER/'polarized_system.json';d=json.loads(system.read_text(encoding='utf-8'))
    s=sp.Symbol('ss');B=sp.Matrix(d['matrix']);K=sp.Matrix(d['inverse_pairing'])
    p8=sp.prod(s-a for a in POLES[1:])
    assert sp.cancel(sp.trace(B)+1/s+sp.diff(p8,s)/(2*p8))==0
    assert sp.factor(K.det()-s*s*p8/3)==0
    p=FOLDER/'high_trunk_transport.json';tr=json.loads(p.read_text(encoding='utf-8'))
    rep=json.loads((FOLDER/'high_trunk_replay.json').read_text(encoding='utf-8'))
    assert rep['certificate_sha256']==sha(p)
    bits=tr['bits'];S=1<<bits;eps=F(tr['error_upper'])
    centre=decode(tr['endpoint_center_gaussian_integer'])
    U=sp.Matrix(6,6,[sp.Rational(a,S)+sp.I*sp.Rational(b,S) for a,b in centre])
    D=sp.diag(*map(sp.Rational,D_VALUES));KD=D.inv()*K*D.inv()
    end=sp.Rational(tr['end'][0])+sp.I*sp.Rational(tr['end'][1])
    K0=KD.subs(s,-64);Ke=KD.subs(s,end).applyfunc(sp.expand)
    R=(U.T*Ke*U-K0).applyfunc(sp.expand)
    norm=lambda M:max(sum(F(abs(sp.re(M[i,j])))+F(abs(sp.im(M[i,j]))) for j in range(6)) for i in range(6))
    # Transpose infinity norms can be up to dimension times the original norm.
    allowed=6*eps*norm(Ke)*(2*F(mnorm(centre),S)+eps)
    residual_norm=norm(R)
    assert residual_norm<=allowed
    invp=FOLDER/'high_trunk_inverse.json';iv=json.loads(invp.read_text(encoding='utf-8'))
    assert iv['transport_sha256']==sha(p)
    V=decode(iv['inverse_center_gaussian_integer'])
    # Independent exact Neumann residual, using explicit scalar sums.
    rr=[]
    for i in range(6):
        for j in range(6):
            val=sumg(gm(V[6*i+k],centre[6*k+j]) for k in range(6))
            rr.append(ga((S*S if i==j else 0,0),gn(val)))
    residual_inverse=F(mnorm(rr),S*S);vn=F(mnorm(V),S)
    eta=residual_inverse+vn*eps
    assert str(residual_inverse)==iv['residual_upper'] and str(eta)==iv['eta_upper']
    assert eta<1 and F(iv['inverse_error_upper'])>=eta*vn/(1-eta)
    assert residual_inverse+vn*F(1)>=1  # this deliberately huge input ball is refused
    out={'status':'known_solutions_inverse_and_polarization_checks_passed',
         'code_sha256':sha(__file__),'transport_sha256':sha(p),'inverse_sha256':sha(invp),
         'known_solution_tests':tests,'inverse_recomputed_without_numerical_inversion':True,
         'excessive_inverse_uncertainty_rejected':True,
         'trace_connection':'-1/s-P8_prime/(2*P8)',
         'det_inverse_pairing':'s^2*P8(s)/3',
         'exact_transport_determinant':'det(U)^2 = b^2*P8(b)/(s^2*P8(s)); branch fixed by U(b)=I',
         'pairing_residual_norm':str(residual_norm),'pairing_allowed_norm':str(allowed),
         'pairing_interval_contains_zero':True,
         'normalization':'published pairing, absolute geometric scale not certified',
         'global_loop_certified':False,'absolute_periods_certified':False}
    (FOLDER/'high_trunk_geometry_check.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k not in ['pairing_residual_norm','pairing_allowed_norm']}),flush=True)

if __name__=='__main__': main()
