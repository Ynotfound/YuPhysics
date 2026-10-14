"""Independent scalar-residual replay for the resonant zero gauge."""
from node_frobenius import *
from verify_high_trunk import conv,padded_add,horner


def check_source(d):
    s=sp.Symbol('ss');geo=json.loads((FOLDER/'polarized_geometry_check.json').read_text())
    original=json.loads((FOLDER/'polarized_system.json').read_text())
    B=sp.Matrix(original['matrix']);B0=sp.Matrix(geo['zero_sheared_connection']);J=sp.Matrix(geo['zero_shearing'])
    Ji=J.inv();assert (B0-J.diff(s)*Ji-J*B*Ji).applyfunc(sp.cancel)==sp.zeros(6)
    Q=sp.Matrix(d['Q']);Qi=sp.Matrix(d['Q_inverse']);N0=sp.Matrix(d['nilpotent'])
    assert Q*Qi==sp.eye(6) and N0*N0==sp.zeros(6) and N0.rank()==2
    assert mnorm(real_matrix(N0))==1
    assert Qi*sp.Matrix(geo['zero_nilpotent_residue'])*Q==N0
    D=sp.diag(*map(sp.Rational,D_VALUES))
    assert (sp.Matrix(d['to_D'])-D*Ji*Q).applyfunc(sp.cancel)==sp.zeros(6)
    qglobal=sp.prod(s-a for a in POLES)
    actual=(qglobal*Qi*B0*Q).applyfunc(sp.cancel)
    assert all(sp.Poly(v,s).degree()<=12 for v in actual if v!=0)
    qr=list(map(F,d['full_source_q']));nr=[list(map(F,M)) for M in d['full_source_N']]
    assert len(qr)==10 and len(nr)<=13
    assert sum(sp.Rational(v)*s**k for k,v in enumerate(qr))==sp.expand(qglobal)
    for j in range(36):
        assert sp.expand(actual[j]-sum(sp.Rational(M[j])*s**k for k,M in enumerate(nr)))==0
    q=decode(d['q_coefficients_gaussian_integer']);H=[decode(M) for M in d['H_coefficients_gaussian_integer']]
    assert len(q)<=9 and len(H)<=12
    scale=int(d['common_q_N_scale']);h=tuple(map(F,d['h']));assert scale>0 and h==(0,F(1,4))
    Nm=real_matrix(N0)
    for t in range(1,23):
        point=gs(t,h);qg=horner([(v,0) for v in qr],point)
        val=[horner([(M[j],0) for M in nr],point) for j in range(36)]
        qt=horner(q,(t,0));ht=[horner([M[j] for M in H],(t,0)) for j in range(36)]
        expected=(scale,0)
        for a in POLES[1:]: expected=gm(expected,ga(point,(-a,0)))
        assert qt==expected
        assert ms(gs(t,qg),ht)==ms(qt,ma(ms(gs(t,h),val),ms(gn(qg),Nm)))
    lower=F(scale)
    for a in POLES[1:]: lower*=F(a)-ab(h)
    return q,H,Nm,lower


def check_residual(q,H,N,P,lower):
    assert P[0]==mi();m=len(P);totals=[]
    for i in range(6):
        for j in range(6):
            euler=[gs(n,M[6*i+j]) for n,M in enumerate(P)]
            for k in range(6):
                for n,M in enumerate(P):
                    euler[n]=ga(euler[n],gn(gm(N[6*i+k],M[6*k+j])))
                    euler[n]=ga(euler[n],gm(M[6*i+k],N[6*k+j]))
            R=conv(q,euler)
            for k in range(6):
                R=padded_add(R,[Z]+conv([A[6*i+k] for A in H],[M[6*k+j] for M in P]),-1)
            assert all(v==Z for v in R[:m])
            totals.append(sum(ab(v) for v in R[m:]))
    M=F(max(sum(ab(A[6*i+j]) for A in H for j in range(6)) for i in range(6)))/lower
    delta=max(sum(totals[6*i:6*i+6]) for i in range(6))/lower
    J=lambda n:F(1,n)+F(2,n*n)+F(2,n*n*n)
    kappa=M*J(m+1);assert kappa<1
    return {'vanishing_degree':m,'denominator_lower':lower,'regular_norm_upper':M,
            'weighted_defect_upper':delta,'kappa':kappa,'gauge_error_upper':delta*J(m)/(1-kappa)}


def main():
    p=FOLDER/'zero_monodromy/frobenius.json';d=json.loads(p.read_text())
    assert d['system_sha256']==sha(FOLDER/'polarized_system.json')
    assert d['geometry_sha256']==sha(FOLDER/'polarized_geometry_check.json')
    assert d['code_sha256']==sha(ROOT/'code/certify_zero_frobenius.py')
    assert d['arithmetic_sha256']==sha(ROOT/'code/zero_frobenius.py')
    q,H,N,lower=check_source(d);print('zero_source_and_resonant_shearing_replayed',flush=True)
    P=[vfdecode(M) for M in d['gauge_coefficients_exact']]
    assert len(P)==d['order']+1
    bd=check_residual(q,H,N,P,lower);assert {k:str(v) for k,v in bd.items()}==d['bounds']
    bits=d['bits'];S=1<<bits;C=decode(d['endpoint_center_gaussian_integer'])
    exact=[horner([M[j] for M in P],(1,0)) for j in range(36)]
    rd=mnorm([ga(v,gn(gs(F(1,S),w))) for v,w in zip(exact,C)])
    raw=bd['gauge_error_upper']+rd;e=F(d['endpoint_error_upper'])
    assert rd==F(d['endpoint_rounding_error']) and raw<=e<raw+F(1,S)
    bad=[list(M) for M in P];bad[1][0]=ga(bad[1][0],(1,0))
    rejected=False
    try: check_residual(q,H,N,bad,lower)
    except AssertionError: rejected=True
    assert rejected
    out={'status':'zero_nilpotent_gauge_independently_replayed','certificate_sha256':sha(p),
         'verifier_sha256':sha(__file__),'coefficient_generator_called':False,
         'floating_point_used':False,'resonant_shearing_checked_against_original':True,
         'low_order_tampering_rejected':rejected,'gauge_error_upper':str(e),
         'gauge_error_dyadic':str(dyadic_upper(e)),'based_monodromy_certified':False}
    save_json(p.with_name('frobenius_replay.json'),out);print(json.dumps(out),flush=True)


if __name__=='__main__': main()
