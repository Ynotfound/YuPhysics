"""Independent exact replay of the horizontal route and both normalizations."""
from gaussian_transport import *
from verify_high_trunk import check_source,check_residual,horner

def propagated_check(L,le,C,ce,new,ne,rounding,bits):
    S=1<<bits
    product=[sumg(gm(L[6*i+k],C[6*k+j]) for k in range(6)) for i in range(6) for j in range(6)]
    rd=max(sum(F(ab(ga(product[6*i+j],gn(gs(S,new[6*i+j])))),S*S) for j in range(6)) for i in range(6))
    assert rd==rounding
    raw=(F(mnorm(L),S)+le)*ce+le*F(mnorm(C),S)+rd
    assert raw<=ne<raw+F(1,S)

def main():
    dest=FOLDER/'horizontal_trunk';p=dest/'transport.json';d=json.loads(p.read_text())
    plan=json.loads((dest/'plan.json').read_text());seedp=FOLDER/'high_trunk_transport.json'
    seed=json.loads(seedp.read_text());rep=json.loads((FOLDER/'high_trunk_replay.json').read_text())
    assert d['plan_sha256']==sha(dest/'plan.json')
    assert d['seed_sha256']==plan['seed_sha256']==rep['certificate_sha256']==sha(seedp)
    assert d['source_system_sha256']==sha(FOLDER/'polarized_system.json')
    assert d['code_sha256']==sha(ROOT/'code/certify_horizontal_trunk.py')
    for name,digest in d['dependency_sha256'].items(): assert digest==sha(ROOT/'code'/name)
    bits=d['bits'];S=1<<bits;source=source_polynomials()
    states={'physical':(decode(seed['endpoint_center_gaussian_integer']),F(seed['error_upper'])),
            'high':(mi(S),F(0))}
    last=seed['end'];seen=[]
    for ix,item in enumerate(d['segments']):
        f=dest/item['file'];assert sha(f)==item['sha256']
        v=json.loads(f.read_text());step=plan['steps'][ix]
        assert v['start']==last==step['start'] and v['end']==step['end']
        assert v['bits']==bits and v['step']==step['step'] and v['landmark_pole']==step['landmark_pole']
        assert v['system_sha256']==d['source_system_sha256'] and v['plan_sha256']==d['plan_sha256']
        q=decode(v['q_coefficients_gaussian_integer']);N=[decode(M) for M in v['N_coefficients_gaussian_integer']]
        P=[decode(M) for M in v['polynomial_coefficients_gaussian_integer']]
        c=tuple(map(F,v['start']));h=tuple(map(F,v['step']));scale=int(v['common_q_N_scale'])
        assert ga(c,h)==tuple(map(F,v['end']))
        check_source(q,N,scale,c,h,source)
        lower,distances=disk_lower(c,h,scale)
        assert list(map(str,distances))==step['pole_distance_lower_on_disk']
        assert lower==scale*F(step['monic_denominator_lower'])
        bd=check_residual(q,N,P,bits,lower)
        assert {k:str(x) for k,x in bd.items()}==v['bounds']
        le=F(v['local_error_upper']);assert bd['error_norm_upper']<=le<bd['error_norm_upper']+F(1,S)
        L=[horner([M[j] for M in P],(1,0)) for j in range(36)]
        assert L==decode(v['endpoint_center_gaussian_integer'])
        for label,(C,ce) in list(states.items()):
            new=decode(v[label+'_center_gaussian_integer']);ne=F(v[label+'_error_upper'])
            propagated_check(L,le,C,ce,new,ne,F(v[label+'_rounding_error']),bits)
            states[label]=(new,ne)
            assert item[label+'_error_upper']==str(ne)
        if step['landmark_pole'] is not None: seen.append(dict(item,pole=step['landmark_pole']))
        last=v['end']
        if (ix+1)%16==0: print(json.dumps({'replayed':ix+1,'end':last}),flush=True)
    assert d['status']=='horizontal_trunk_to_all_nine_landmarks_certified'
    assert len(d['segments'])==len(plan['steps'])==73 and seen==d['landmarks']
    assert [v['pole'] for v in seen]==POLES and d['end']==last
    for label,(C,e) in states.items():
        assert decode(d[label+'_center_gaussian_integer'])==C and F(d[label+'_error_upper'])==e
    out={'status':'all_horizontal_certificates_independently_replayed','transport_sha256':sha(p),
         'verifier_sha256':sha(__file__),'source_checker_sha256':sha(ROOT/'code/verify_high_trunk.py'),
         'steps':73,'landmark_poles':POLES,'end':last,'both_initial_errors_and_rounding_checked':True,
         'source_matching':'exact evaluation and degree bounds','residuals':'scalar polynomial convolution',
         'numerical_proposal_called':False,'floating_arithmetic_used':False,
         'high_error_upper':str(states['high'][1]),'physical_sequential_error_upper':str(states['physical'][1]),
         'closed_loops_certified':False,'absolute_periods_certified':False}
    (dest/'replay.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if 'error_upper' not in k}),flush=True)

if __name__=='__main__': main()
