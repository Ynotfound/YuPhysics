"""Replay saved branch matrices, vector propagation and rebased reflections."""
from node_frobenius import *
from verify_high_trunk import check_residual,check_source,horner
from verify_node_frobenius import reflection_check
import argparse


def check_frame_source(q,N,scale,c,h,source,Q,Qi,frame):
    if frame=='D':
        return check_source(q,N,scale,c,h,source)
    assert frame=='Q'
    qraw,nraw=source;assert len(q)<=10 and len(N)<=12
    qm=real_matrix(Q);qim=real_matrix(Qi)
    for t in range(12):
        s=ga(c,gs(t,h))
        assert horner(q,(t,0))==gs(scale,horner([(v,0) for v in qraw],s))
        val=[horner([(M[j],0) for M in nraw],s) for j in range(36)]
        wanted=ms(gs(scale,h),mm(mm(qim,val),qm))
        got=[horner([M[j] for M in N],(t,0)) for j in range(36)]
        assert got==wanted


def check_action(A,ae,C,ce,new,new_error,bits):
    S=1<<bits
    exact=[sumg(gm(A[6*i+j],C[j]) for j in range(6)) for i in range(6)]
    rd=max(F(ab(ga(v,gn(gs(S,w)))),S*S) for v,w in zip(exact,new))
    raw=(F(mnorm(A),S)+ae)*ce+ae*F(vnorm(C),S)+rd
    assert raw<=new_error<raw+F(1,S)
    return rd


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--nodes',type=int,nargs='+',default=[1,9,25,49,81])
    args=parser.parse_args();source=source_polynomials()
    planp=FOLDER/'node_frobenius/branch_plan.json';plan=json.loads(planp.read_text())
    seedrep=json.loads((FOLDER/'node_frobenius/replay.json').read_text())
    hrp=FOLDER/'horizontal_trunk/landmarks.json';hor=json.loads(hrp.read_text())
    assert json.loads(hrp.with_name('landmarks_replay.json').read_text())['landmarks_sha256']==sha(hrp)
    for a in args.nodes:
        dest=FOLDER/'node_branches'/f'node_{a}';tp=dest/'transport.json';d=json.loads(tp.read_text())
        seedp=FOLDER/f'node_frobenius/node_{a}.json';seed=json.loads(seedp.read_text())
        assert next(c for c in seedrep['certificates'] if c['node']==a)['sha256']==sha(seedp)
        assert d['node']==a and d['system_sha256']==sha(FOLDER/'polarized_system.json')
        assert d['seed_sha256']==sha(seedp) and d['plan_sha256']==sha(planp)
        assert d['landmarks_sha256']==sha(hrp)
        assert d['code_sha256']==sha(ROOT/'code/certify_node_branches.py')
        assert d['arithmetic_sha256']==sha(ROOT/'code/node_branch_transport.py')
        assert d['source_arithmetic_sha256']==sha(ROOT/'code/gaussian_transport.py')
        assert d['status']=='node_branch_and_common_basepoint_reflection_certified'
        bits=d['bits'];S=1<<bits
        Q=sp.Matrix(seed['Q']);Qi=sp.Matrix(seed['Q_inverse'])
        assert Q*Qi==sp.eye(6)
        coeff=[vfdecode(v) for v in seed['frobenius_vector_coefficients_exact']]
        exact=[sumg(v[i] for v in coeff) for i in range(6)]
        C=[(dyadic_round(x,bits),dyadic_round(y,bits)) for x,y in exact]
        rd=max(ab(ga(v,gn(gs(F(1,S),w)))) for v,w in zip(exact,C))
        assert F(d['initial_rounding_error'])==rd
        raw=F(seed['bounds']['vector_error_Q_basis'])+rd;e=F(d['initial_error_upper'])
        assert raw<=e<raw+F(1,S)
        branch=next(b for b in plan['branches'] if b['node']==a)
        steps=list(reversed(branch['steps_top_to_bottom']))
        frame='Q';last=seed['endpoint'];checked=[];negative=False
        assert len(d['segments'])==64
        for ix,(ref,st) in enumerate(zip(d['segments'],steps)):
            p=dest/ref['file'];assert sha(p)==ref['sha256'];v=json.loads(p.read_text())
            assert v['node']==a and v['index']==ix+1 and v['bits']==bits
            assert v['start']==last==st['end'] and v['end']==st['start']==ref['end']
            assert v['seed_sha256']==sha(seedp) and v['plan_sha256']==sha(planp)
            assert v['system_sha256']==d['system_sha256']
            c=tuple(map(F,v['start']));h=tuple(map(F,st['up']['step']))
            assert ga(c,h)==tuple(map(F,v['end']))
            switched=frame=='Q' and c[1]>=F(d['switch_height'])
            assert switched==v['switched_before_step']
            if switched:
                exact=[sumg(gs(rat(Q[i,j]),C[j]) for j in range(6)) for i in range(6)]
                new=decode(v['input_center_gaussian_integer'])
                rd=max(F(ab(ga(x,gn(y))),S) for x,y in zip(exact,new))
                raw=mnorm(real_matrix(Q))*e+rd;ne=F(v['input_error_upper'])
                assert rd==F(v['switch_rounding_error']) and raw<=ne<raw+F(1,S)
                C,e,frame=new,ne,'D'
            else:
                assert decode(v['input_center_gaussian_integer'])==C
                assert F(v['input_error_upper'])==e and F(v['switch_rounding_error'])==0
            assert v['frame']==frame
            q=decode(v['q_coefficients_gaussian_integer'])
            N=[decode(M) for M in v['N_coefficients_gaussian_integer']]
            P=[decode(M) for M in v['polynomial_coefficients_gaussian_integer']]
            assert len(P)==v['order']+1
            scale=int(v['common_q_N_scale']);assert scale>0
            check_frame_source(q,N,scale,c,h,source,Q,Qi,frame)
            lower,ds=disk_lower(c,h,scale)
            assert list(map(str,ds))==st['up']['pole_distance_lower']
            bd=check_residual(q,N,P,bits,lower)
            assert {k:str(x) for k,x in bd.items()}==v['bounds']
            le=F(v['local_matrix_error_upper'])
            assert bd['error_norm_upper']<=le<bd['error_norm_upper']+F(1,S)
            L=[sumg(M[j] for M in P) for j in range(36)]
            assert L==decode(v['transition_center_gaussian_integer'])
            new=decode(v['output_center_gaussian_integer']);ne=F(v['output_error_upper'])
            rd=check_action(L,le,C,e,new,ne,bits)
            assert rd==F(v['vector_rounding_error'])
            if ix==0:
                try: check_action(L,le,C,e,new,F(0),bits)
                except AssertionError: negative=True
                assert negative
            C,e,last=new,ne,v['end'];checked.append(ref['sha256'])
            if (ix+1)%16==0: print(json.dumps({'node_replayed':a,'segments':ix+1}),flush=True)
        assert last==d['end']==[str(a),'64'] and frame==d['frame']=='D'
        assert C==decode(d['center_gaussian_integer']) and e==F(d['error_upper'])
        point=next(p for p in hor['landmarks'] if p['pole']==a)
        V=decode(point['high_inverse']['center_gaussian_integer']);ve=F(point['high_inverse']['error_upper'])
        based=decode(d['based_center_gaussian_integer']);be=F(d['based_error_upper'])
        rd=check_action(V,ve,C,e,based,be,bits)
        assert rd==F(d['rebasing_rounding_error']) and d['based_point']==['-64','64']
        fake={'node':-64,'h':['0','64'],'reflection':d['based_reflection'],
              'pairing_at_endpoint_exact':d['based_pairing_exact']}
        er=reflection_check(fake,[[gs(F(1,S),v) for v in based]],sp.eye(6),be)
        out={'status':'node_branch_and_based_reflection_independently_replayed','node':a,
             'transport_sha256':sha(tp),'verifier_sha256':sha(__file__),'segments':64,
             'segment_sha256':checked,'based_point':d['based_point'],
             'based_reflection_error_upper':str(er),'based_reflection_error_dyadic':str(dyadic_upper(er)),
             'source_check':'12 exact points and degree bounds, with frame conjugation',
             'residual_check':'independent scalar convolutions',
             'coefficient_generator_called':False,'floating_point_used':False,
             'zero_radius_tampering_rejected':negative,'zero_monodromy_certified':False,
             'Q6_21_traces_certified':False,'absolute_periods_certified':False}
        save_json(dest/'replay.json',out)
        print(json.dumps({k:v for k,v in out.items() if k!='segment_sha256'}),flush=True)


if __name__=='__main__': main()
