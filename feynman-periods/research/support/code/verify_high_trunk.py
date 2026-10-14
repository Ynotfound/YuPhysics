"""Independent replay: no proposal generator, no Taylor recurrence, no mpmath.

Source matching uses exact polynomial evaluations and degree bounds, whereas
the generator uses binomial coefficient shifts. Residuals are recomputed
entry-by-entry by scalar convolution rather than matrix coefficient products.
"""
from gaussian_transport import *
from math import ceil
import argparse

def horner(p,t):
    v=Z
    for a in reversed(p): v=ga(gm(v,t),a)
    return v
def conv(a,b):
    c=[Z]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x==Z: continue
        for j,y in enumerate(b):
            if y!=Z: c[i+j]=ga(c[i+j],gm(x,y))
    return c
def padded_add(a,b,sign=1):
    return [ga(a[k] if k<len(a) else Z,gs(sign,b[k]) if k<len(b) else Z)
            for k in range(max(len(a),len(b)))]

def check_source(q,N,scale,c,h,source):
    qraw,nraw=source
    assert len(q)<=len(qraw) and len(N)<=len(nraw)
    # Each component difference is a polynomial of degree <= degree_limit.
    # A degree+1 exact evaluation grid proves coefficient equality.
    for t in range(max(len(qraw),len(nraw))):
        s=ga(c,gs(t,h))
        assert horner(q,(t,0))==gs(scale,horner([(v,0) for v in qraw],s))
        for j in range(36):
            lhs=horner([A[j] for A in N],(t,0))
            rhs=gs(scale,gm(h,horner([(A[j],0) for A in nraw],s)))
            assert lhs==rhs

def check_residual(q,N,P,bits,lower):
    assert P[0]==mi(1<<bits)
    component_bounds=[]
    for i in range(6):
        for j in range(6):
            derivative=[gs(k,P[k][6*i+j]) for k in range(1,len(P))]
            r=conv(q,derivative)
            for k in range(6):
                r=padded_add(r,conv([A[6*i+k] for A in N],[C[6*k+j] for C in P]),-1)
            component_bounds.append(sum(ab(v) for v in r))
    residual_bound=max(sum(component_bounds[6*i:6*i+6]) for i in range(6))
    M=F(max(sum(ab(A[6*i+j]) for A in N for j in range(6)) for i in range(6)))/lower
    defect=F(residual_bound,1<<bits)/lower
    exponential=F(3)**ceil(M)
    return {'denominator_lower':lower,'matrix_norm_upper':M,'residual_norm_upper':defect,
            'exp_upper':exponential,'error_norm_upper':exponential*defect}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--reject-tampering',action='store_true');args=parser.parse_args()
    p=FOLDER/'high_trunk_transport.json';d=json.loads(p.read_text(encoding='utf-8'))
    assert d['source_system_sha256']==sha(FOLDER/'polarized_system.json')
    assert d['plan_sha256']==sha(FOLDER/'high_trunk_plan.json')
    assert d['seed_sha256']==sha(FOLDER/'first_trunk_step.json')
    assert d['code_sha256']==sha(ROOT/'code/certify_high_trunk.py')
    assert d['arithmetic_code_sha256']==sha(ROOT/'code/gaussian_transport.py')
    seed_replay=json.loads((FOLDER/'first_trunk_step_replay.json').read_text(encoding='utf-8'))
    assert seed_replay['certificate_sha256']==d['seed_sha256']
    assert seed_replay['status']=='independent_polynomial_residual_replay_passed'
    plan=json.loads((FOLDER/'high_trunk_plan.json').read_text(encoding='utf-8'))
    bits=d['bits'];S=1<<bits;source=source_polynomials()
    C,e=import_first_step(bits);last=['-64','1'];checked=[]
    for ix,item in enumerate(d['segments']):
        f=FOLDER/item['file'];assert sha(f)==item['sha256']
        v=json.loads(f.read_text(encoding='utf-8'));step=plan['steps'][ix]
        assert v['start']==last==step['start'] and v['end']==step['end']
        assert v['bits']==bits and v['system_sha256']==d['source_system_sha256']
        assert v['plan_sha256']==d['plan_sha256']
        q=decode(v['q_coefficients_gaussian_integer']);N=[decode(A) for A in v['N_coefficients_gaussian_integer']]
        P=[decode(A) for A in v['polynomial_coefficients_gaussian_integer']]
        scale=int(v['common_q_N_scale']);c=tuple(map(F,v['start']));h=tuple(map(F,step['step']))
        assert ga(c,h)==tuple(map(F,v['end'])) and scale>0
        assert v['order']==len(P)-1
        check_source(q,N,scale,c,h,source)
        lower,distances=disk_lower(c,h,scale)
        assert list(map(str,distances))==step['pole_distance_lower_on_disk']
        b=check_residual(q,N,P,bits,lower)
        assert {k:str(x) for k,x in b.items()}==v['bounds']
        le=F(v['local_error_upper']);assert le>=b['error_norm_upper']>=0
        L=[horner([A[j] for A in P],(1,0)) for j in range(36)]
        assert L==decode(v['endpoint_center_gaussian_integer'])
        new=decode(v['propagated_center_gaussian_integer'])
        # Independent product and bound, allowing any saved rounded centre.
        product=[sumg(gm(L[6*i+k],C[6*k+j]) for k in range(6)) for i in range(6) for j in range(6)]
        rounding=max(sum(F(ab(ga(product[6*i+j],gn(gs(S,new[6*i+j])))),S*S) for j in range(6)) for i in range(6))
        assert rounding==F(v['composition_rounding_error'])
        bound=(F(mnorm(L),S)+le)*e+le*F(mnorm(C),S)+rounding
        new_error=F(v['propagated_error_upper']);assert new_error>=bound
        assert item['end']==v['end'] and F(item['propagated_error_upper'])==new_error
        C,e,last=new,new_error,v['end'];checked.append(item['sha256'])
        print(json.dumps({'replayed_segment':ix+1,'end':last}),flush=True)
    assert d['end']==last and F(d['error_upper'])==e
    assert decode(d['endpoint_center_gaussian_integer'])==C
    negative=False
    if args.reject_tampering:
        tampered=[list(A) for A in P];tampered[1][0]=ga(tampered[1][0],(S,0))
        bad=check_residual(q,N,tampered,bits,lower)
        assert bad['error_norm_upper']>le
        try:
            qbad=list(q);qbad[0]=ga(qbad[0],(1,0));check_source(qbad,N,scale,c,h,source)
        except AssertionError: negative=True
        assert negative
    out={'status':'all_saved_complex_steps_independently_replayed','certificate_sha256':sha(p),
         'verifier_sha256':sha(__file__),'arithmetic_sha256':sha(ROOT/'code/gaussian_transport.py'),
         'segments':len(checked),'segment_sha256':checked,'source_matching':'exact evaluation + degree bounds',
         'residual_algorithm':'independent scalar polynomial convolutions','Taylor_generator_called':False,
         'floating_point_used':False,'end':last,'error_upper':str(e),'tampering_tests_passed':negative,
         'global_loop_certified':False,'absolute_periods_certified':False}
    (FOLDER/'high_trunk_replay.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out),flush=True)

if __name__=='__main__': main()
