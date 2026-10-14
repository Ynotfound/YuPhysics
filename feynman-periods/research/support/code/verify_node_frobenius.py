"""Replay actual local node certificates without generating Frobenius series."""
from node_frobenius import *
from verify_high_trunk import conv,padded_add,horner

def source_check(d,source):
    qraw,nraw=source;a=d['node'];h=tuple(map(F,d['h']));scale=int(d['common_q_N_scale'])
    Q=sp.Matrix(d['Q']);Qi=sp.Matrix(d['Q_inverse']);R=sp.Matrix(d['residue_in_D_basis'])
    Lambda=sp.diag(sp.Rational(-1,2),0,0,0,0,0)
    assert Q*Qi==sp.eye(6) and Qi*R*Q==Lambda and R.rank()==1
    qa=F(1)
    for b in POLES:
        if b!=a: qa*=a-b
    for j in range(36): assert rat(R[j])==sum(M[j]*F(a)**k for k,M in enumerate(nraw))/qa
    q=decode(d['q_coefficients_gaussian_integer']);N=[decode(M) for M in d['N_coefficients_gaussian_integer']]
    assert len(q)<=9 and len(N)<=11
    qm=real_matrix(Q);qim=real_matrix(Qi);lm=real_matrix(Lambda)
    # Polynomial identity of degree <=20 after clearing the global denominator.
    # 21 exact sample values suffice; no Taylor or shift routine is used.
    for t in range(1,22):
        point=ga((F(a),F(0)),gs(t,h))
        qglobal=horner([(v,0) for v in qraw],point)
        nsource=[horner([(M[j],0) for M in nraw],point) for j in range(36)]
        qt=horner(q,(t,0));nt=[horner([M[j] for M in N],(t,0)) for j in range(36)]
        expected=(F(scale),F(0))
        for b in POLES:
            if b!=a: expected=gm(expected,ga(point,(-b,0)))
        assert qt==expected
        rhs=ms(qt,ma(ms(gs(t,h),mm(mm(qim,nsource),qm)),ms(gn(qglobal),lm)))
        lhs=ms(gs(t,qglobal),nt)
        assert lhs==rhs
    return q,N,Q

def residual_check(q,N,P,lower):
    m=len(P);assert P[0]==[(F(i==0),F(0)) for i in range(6)]
    coeff_bounds=[]
    for i in range(6):
        shifted=[gs(F(n)-D0[i],v[i]) for n,v in enumerate(P)]
        r=conv(q,shifted)
        for j in range(6): r=padded_add(r,[Z]+conv([M[6*i+j] for M in N],[v[j] for v in P]),-1)
        assert all(x==Z for x in r[:m]), 'low-order Fuchs residual is nonzero'
        coeff_bounds.append(sum(ab(v) for v in r[m:]))
    M=F(max(sum(ab(A[6*i+j]) for A in N for j in range(6)) for i in range(6)))/lower
    delta=max(coeff_bounds)/lower;kappa=M/(F(m)+F(1,2));assert kappa<1
    error=delta/((F(m)-F(1,2))*(1-kappa))
    return {'vanishing_degree':m,'denominator_lower':lower,'M':M,'delta':delta,'kappa':kappa,'vector_error_Q_basis':error}

def reflection_check(d,P,Q,error):
    ref=d['reflection'];bits=ref['bits'];S=1<<bits
    # Evaluate the witness and change basis with exact SymPy matrices,
    # independently of the generator's row-major matrix-vector operation.
    pc=sp.Matrix([sum(sp.Rational(v[i][0])+sp.I*sp.Rational(v[i][1]) for v in P) for i in range(6)])
    vc=(Q*pc).applyfunc(sp.expand)
    vexact=[(rat(sp.re(v)),rat(sp.im(v))) for v in vc]
    V=[gs(F(1,S),v) for v in decode(ref['vector_center_gaussian_integer'])]
    rho=max(ab(ga(a,gn(b))) for a,b in zip(vexact,V))
    assert str(rho)==ref['vector_rounding_error']
    raw=F(mnorm(real_matrix(Q)))*error+rho;ev=F(ref['vector_error_upper'])
    assert raw<=ev<raw+F(1,S)
    s=sp.Symbol('ss');data=json.loads((FOLDER/'polarized_system.json').read_text())
    K0=sp.Matrix(data['inverse_pairing']);D=sp.diag(*map(sp.Rational,D_VALUES))
    h=list(map(sp.Rational,d['h']));point=d['node']+h[0]+sp.I*h[1]
    Ks=(D.inv()*K0*D.inv()).subs(s,point).applyfunc(sp.expand)
    K=[(rat(sp.re(v)),rat(sp.im(v))) for v in Ks]
    assert encode(K)==d['pairing_at_endpoint_exact']
    cov=[sumg(gm(V[i],K[6*i+j]) for i in range(6)) for j in range(6)]
    dc=sumg(gm(cov[i],V[i]) for i in range(6));B=[gm(V[i],cov[j]) for i in range(6) for j in range(6)]
    Delta=6*mnorm(K)*ev*(2*vnorm(V)+ev);d0=max(abs(dc[0]),abs(dc[1]));low=d0-Delta
    assert low>0 and encode([dc])[0]==ref['pairing_center']
    assert str(Delta)==ref['pairing_error_upper'] and str(low)==ref['pairing_modulus_lower']
    exact=[ga(x,gs(-2,gd(y,dc))) for x,y in zip(mi(),B)]
    saved=[gs(F(1,S),v) for v in decode(ref['reflection_center_gaussian_integer'])]
    rd=mnorm([ga(v,gn(w)) for v,w in zip(exact,saved)])
    assert str(rd)==ref['reflection_rounding_error']
    raw=2*(Delta/low+mnorm(B)*Delta/(d0*low))+rd
    er=F(ref['reflection_error_upper']);assert raw<=er<raw+F(1,S)
    assert ab(ga(sumg(saved[6*i+i] for i in range(6)),(-4,0)))<=6*er
    return er

def main():
    dest=FOLDER/'node_frobenius';source=source_polynomials();items=[];negative=False
    for a in POLES[1:]:
        p=dest/f'node_{a}.json';d=json.loads(p.read_text())
        assert d['system_sha256']==sha(FOLDER/'polarized_system.json')
        assert d['code_sha256']==sha(ROOT/'code/certify_node_frobenius.py')
        assert d['arithmetic_sha256']==sha(ROOT/'code/node_frobenius.py')
        assert d['source_arithmetic_sha256']==sha(ROOT/'code/gaussian_transport.py')
        assert d['node']==a and d['endpoint']==[str(a),d['h'][1]] and F(d['h'][0])==0
        q,N,Q=source_check(d,source)
        h=tuple(map(F,d['h']));ds=[F(abs(a-b))-ab(h) for b in POLES if b!=a]
        assert min(ds)>0 and list(map(str,ds))==d['pole_distance_lower']
        lower=F(d['common_q_N_scale'])
        for v in ds: lower*=v
        P=[vfdecode(v) for v in d['frobenius_vector_coefficients_exact']]
        assert len(P)==d['order']+1
        bd=residual_check(q,N,P,lower)
        assert {k:str(v) for k,v in bd.items()}==d['bounds']
        er=reflection_check(d,P,Q,bd['vector_error_Q_basis'])
        items.append({'node':a,'file':p.name,'sha256':sha(p),'order':d['order'],'endpoint':d['endpoint'],
                      'local_reflection_error_upper':str(er),'dyadic_upper':str(dyadic_upper(er))})
        if a==1:
            bad=[list(v) for v in P];bad[1][0]=ga(bad[1][0],(1,0))
            try: residual_check(q,N,bad,lower)
            except AssertionError: negative=True
            assert negative
        print(json.dumps({'node_replayed':a,'order':d['order'],'dyadic_error':str(dyadic_upper(er))}),flush=True)
    out={'status':'eight_actual_node_local_reflections_independently_replayed',
         'verifier_sha256':sha(__file__),'source_system_sha256':sha(FOLDER/'polarized_system.json'),
         'certificates':items,'Frobenius_coefficient_generator_called':False,'floating_point_used':False,
         'source_check':'21 exact evaluations of the degree at most 20 identity',
         'residual_check':'independent scalar convolutions','low_order_tampering_rejected':negative,
         'all_local_seeds':['a','1/4'],'uniform_reflection_error_upper':str(max(F(x['dyadic_upper']) for x in items)),
         'common_basepoint_node_monodromies_certified':False,'zero_monodromy_certified':False,
         'Q6_traces_certified':False,'absolute_amplitudes_certified':False}
    save_json(dest/'replay.json',out)
    print(json.dumps({k:v for k,v in out.items() if k!='certificates'}),flush=True)

if __name__=='__main__': main()
