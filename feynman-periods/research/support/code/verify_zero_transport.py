"""Independent scalar residual, frame, inverse and Machin-bound replay."""
from node_frobenius import *
from verify_high_trunk import check_source,check_residual
from verify_horizontal_landmarks import inverse_check
from verify_node_branches import check_action
from verify_node_frobenius import reflection_check


def check_product(L,le,C,ce,new,ne,bits):
    S=1<<bits
    exact=[sumg(gm(L[6*i+k],C[6*k+j]) for k in range(6)) for i in range(6) for j in range(6)]
    rd=max(sum(F(ab(ga(exact[6*i+j],gn(gs(S,new[6*i+j])))),S*S) for j in range(6)) for i in range(6))
    raw=(F(mnorm(L),S)+le)*ce+le*F(mnorm(C),S)+rd
    assert raw<=ne<raw+F(1,S)
    return rd


def check_pi(d):
    n=d['terms'];assert n==100
    intervals=[]
    # Successive-term ratio recurrence, independent of the generator power sums.
    for q in [5,239]:
        term=F(1,q);total=F(0)
        for k in range(n):
            total+=term
            term*=-F(2*k+1,(2*k+3)*q*q)
        intervals.append(sorted([total,total+term]))
    low=16*intervals[0][0]-4*intervals[1][1]
    high=16*intervals[0][1]-4*intervals[1][0]
    assert F(d['lower'])==low and F(d['upper'])==high
    assert F(d['center'])==(low+high)/2 and F(d['radius'])==(high-low)/2
    # Exact tangent identity; the angle is in (0, pi/2).
    t=F(1,5)
    t=2*t/(1-t*t);t=2*t/(1-t*t)
    assert (t-F(1,239))/(1+t/F(239))==1
    return (low+high)/2,(high-low)/2


def main():
    dest=FOLDER/'zero_monodromy';tp=dest/'transport.json';d=json.loads(tp.read_text())
    seedp=dest/'frobenius.json';seed=json.loads(seedp.read_text())
    assert json.loads((dest/'frobenius_replay.json').read_text())['certificate_sha256']==sha(seedp)
    planp=dest/'plan.json';plan=json.loads(planp.read_text())
    assert d['seed_sha256']==plan['seed_sha256']==sha(seedp)
    assert d['plan_sha256']==sha(planp)
    assert d['code_sha256']==plan['code_sha256']==sha(ROOT/'code/certify_zero_transport.py')
    assert d['system_sha256']==plan['system_sha256']==sha(FOLDER/'polarized_system.json')
    assert d['status']=='zero_based_monodromy_and_sixth_reflection_certified'
    assert len(d['segments'])==len(plan['steps'])==64
    bits=d['bits'];S=1<<bits;assert bits==seed['bits']
    C=decode(seed['endpoint_center_gaussian_integer']);e=F(seed['endpoint_error_upper'])
    source=source_polynomials()
    zsource=(list(map(F,seed['full_source_q'])),[list(map(F,A)) for A in seed['full_source_N']])
    toD=sp.Matrix(seed['to_D']);s=sp.Symbol('ss');frame='zero_Q';last=seed['endpoint']
    assert plan['start']==last==['0','1/4'] and plan['switch_to_D_height']=='2'
    checked=[];negative=False
    for ix,(entry,st) in enumerate(zip(d['segments'],plan['steps'])):
        p=dest/entry['file'];assert sha(p)==entry['sha256'];v=json.loads(p.read_text())
        assert v['index']==ix+1 and v['bits']==bits
        assert v['start']==st['start']==last and v['end']==st['end']==entry['end']
        assert v['system_sha256']==d['system_sha256'] and v['plan_sha256']==sha(planp)
        c=tuple(map(F,last));h=tuple(map(F,st['step']))
        assert ga(c,h)==tuple(map(F,v['end']))
        switched=frame=='zero_Q' and c[1]>=2
        assert switched==v['switched_before_step']
        if switched:
            evaluated=toD.subs(s,sp.Rational(c[0])+sp.I*sp.Rational(c[1])).applyfunc(sp.expand)
            T=[(rat(sp.re(x)),rat(sp.im(x))) for x in evaluated]
            assert T==vfdecode(v['switch_matrix_exact'])
            exact=[sumg(gm(T[6*i+k],gs(F(1,S),C[6*k+j])) for k in range(6)) for i in range(6) for j in range(6)]
            new=decode(v['input_center_gaussian_integer']);ne=F(v['input_error_upper'])
            rd=mnorm([ga(x,gn(gs(F(1,S),y))) for x,y in zip(exact,new)])
            raw=mnorm(T)*e+rd
            assert raw<=ne<raw+F(1,S) and rd==F(v['switch_rounding_error'])
            C,e,frame=new,ne,'D'
        else:
            assert decode(v['input_center_gaussian_integer'])==C and F(v['input_error_upper'])==e
            assert v['switch_matrix_exact'] is None and F(v['switch_rounding_error'])==0
        assert v['frame']==frame
        q=decode(v['q_coefficients_gaussian_integer']);N=[decode(A) for A in v['N_coefficients_gaussian_integer']]
        P=[decode(A) for A in v['polynomial_coefficients_gaussian_integer']]
        assert len(P)==v['order']+1
        scale=int(v['common_q_N_scale']);assert scale>0
        check_source(q,N,scale,c,h,zsource if frame=='zero_Q' else source)
        lower,ds=disk_lower(c,h,scale)
        assert list(map(str,ds))==st['pole_distance_lower']
        assert lower==scale*F(st['monic_denominator_lower'])
        bd=check_residual(q,N,P,bits,lower)
        assert {k:str(x) for k,x in bd.items()}==v['bounds']
        le=F(v['local_matrix_error_upper'])
        assert bd['error_norm_upper']<=le<bd['error_norm_upper']+F(1,S)
        L=[sumg(A[j] for A in P) for j in range(36)]
        assert L==decode(v['transition_center_gaussian_integer'])
        new=decode(v['output_center_gaussian_integer']);ne=F(v['output_error_upper'])
        rd=check_product(L,le,C,e,new,ne,bits)
        assert rd==F(v['matrix_rounding_error'])
        if ix==0:
            try: check_product(L,le,C,e,new,F(0),bits)
            except AssertionError: negative=True
            assert negative
        C,e,last=new,ne,v['end'];checked.append(entry['sha256'])
        if (ix+1)%16==0: print(json.dumps({'zero_steps_replayed':ix+1}),flush=True)
    assert C==decode(d['center_gaussian_integer']) and e==F(d['error_upper'])
    assert last==d['end']==plan['end']==['0','64'] and frame==d['frame']=='D'
    hp=FOLDER/'horizontal_trunk/landmarks.json';hor=json.loads(hp.read_text())
    assert d['landmarks_sha256']==sha(hp)
    assert json.loads(hp.with_name('landmarks_replay.json').read_text())['landmarks_sha256']==sha(hp)
    item=next(x for x in hor['landmarks'] if x['pole']==0)
    H=[gs(1<<(bits-hor['bits']),v) for v in decode(item['high_inverse']['center_gaussian_integer'])]
    he=F(item['high_inverse']['error_upper'])
    W=decode(d['based_frame_center_gaussian_integer']);we=F(d['based_frame_error_upper'])
    assert d['based_point']==['-64','64']
    assert check_product(H,he,C,e,W,we,bits)==F(d['rebasing_rounding_error'])
    V=decode(d['based_frame_inverse']['center_gaussian_integer'])
    ve=inverse_check(W,we,d['based_frame_inverse'],bits)
    Nint=[gs(S,v) for v in real_matrix(sp.Matrix(seed['nilpotent']))]
    # Right multiplication by the 0/1 Jordan matrix is exact on dyadic centres.
    WN=[gs(F(1,S),sumg(gm(W[6*i+k],Nint[6*k+j]) for k in range(6))) for i in range(6) for j in range(6)]
    assert all(x.denominator==y.denominator==1 for x,y in WN)
    raw=F(mnorm(Nint),S)*we
    wne=F(ceil(raw*S),S)
    assert check_product(W,we,Nint,F(0),WN,wne,bits)==F(d['nilpotent_first_rounding_error'])==0
    A=decode(d['conjugated_nilpotent_center_gaussian_integer']);ae=F(d['conjugated_nilpotent_error_upper'])
    assert check_product(WN,wne,V,ve,A,ae,bits)==F(d['nilpotent_second_rounding_error'])
    pc,pe=check_pi(d['pi_certificate'])
    Ac=[gs(F(1,S),v) for v in A]
    exact=[ga((int(i==j),0),gm((0,2*pc),Ac[6*i+j])) for i in range(6) for j in range(6)]
    T=decode(d['monodromy_center_gaussian_integer']);te=F(d['monodromy_error_upper'])
    rd=mnorm([ga(x,gn(gs(F(1,S),y))) for x,y in zip(exact,T)])
    raw=2*((abs(pc)+pe)*ae+pe*mnorm(Ac))+rd
    assert raw<=te<raw+F(1,S) and rd==F(d['monodromy_rounding_error'])
    n1p=FOLDER/'node_branches/node_1/transport.json';n1=json.loads(n1p.read_text())
    assert d['node_one_sha256']==sha(n1p)
    assert json.loads(n1p.with_name('replay.json').read_text())['transport_sha256']==sha(n1p)
    v1=[gs(1<<(bits-n1['bits']),v) for v in decode(n1['based_center_gaussian_integer'])]
    e1=F(n1['based_error_upper']);v6=decode(d['sixth_vector_center_gaussian_integer']);e6=F(d['sixth_vector_error_upper'])
    assert check_action(T,te,v1,e1,v6,e6,bits)==F(d['sixth_vector_rounding_error'])
    fake={'node':-64,'h':['0','64'],'reflection':d['sixth_reflection'],'pairing_at_endpoint_exact':d['based_pairing_exact']}
    er=reflection_check(fake,[[gs(F(1,S),v) for v in v6]],sp.eye(6),e6)
    out={'status':'zero_transport_monodromy_and_sixth_reflection_independently_replayed',
         'transport_sha256':sha(tp),'verifier_sha256':sha(__file__),'segments':64,'segment_sha256':checked,
         'monodromy_error_upper':str(te),'monodromy_error_dyadic':str(dyadic_upper(te)),
         'sixth_reflection_error_upper':str(er),'sixth_reflection_error_dyadic':str(dyadic_upper(er)),
         'pi_bound_replayed':True,'residual_algorithm':'scalar convolution',
         'coefficient_generator_called':False,'floating_point_used':False,
         'zero_radius_tampering_rejected':negative,'absolute_periods_certified':False}
    save_json(dest/'replay.json',out)
    print(json.dumps({k:v for k,v in out.items() if k!='segment_sha256'}),flush=True)


if __name__=='__main__': main()
