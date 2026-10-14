"""Generate complex-centred dyadic proposals, then certify EXACT residuals.

mpmath only proposes polynomial coefficients; it supplies no trusted bound.
Saved segments are reusable checkpoints and are never silently overwritten.
"""
from gaussian_transport import *
import mpmath as mp
import argparse,time

def proposal(q,N,order,bits):
    mp.mp.dps=100
    cv=lambda z:mp.mpc(z[0],z[1])
    qm=list(map(cv,q));nm=[mp.matrix(6,6) for _ in N]
    for A,M in zip(nm,N):
        for j,v in enumerate(M): A[j//6,j%6]=cv(v)
    P=[mp.eye(6)]
    for n in range(order):
        V=mp.zeros(6)
        for j in range(min(n,len(N)-1)+1): V+=nm[j]*P[n-j]
        for j in range(1,min(n,len(q)-1)+1): V-=qm[j]*(n-j+1)*P[n-j+1]
        P.append(V/(qm[0]*(n+1)))
    S=1<<bits
    return [[(int(mp.nint(mp.re(v)*S)),int(mp.nint(mp.im(v)*S))) for v in M] for M in P]

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--steps',type=int,default=8)
    args=parser.parse_args()
    planp=FOLDER/'high_trunk_plan.json';plan=json.loads(planp.read_text(encoding='utf-8'))
    assert plan['system_sha256']==sha(FOLDER/'polarized_system.json')
    assert 0<args.steps<=len(plan['steps'])
    source=source_polynomials();bits=240;order=64
    dest=FOLDER/'high_trunk';dest.mkdir(exist_ok=True)
    C,e=import_first_step(bits)
    summary={'status':'partial_high_trunk','source_system_sha256':sha(FOLDER/'polarized_system.json'),
             'plan_sha256':sha(planp),'seed_sha256':sha(FOLDER/'first_trunk_step.json'),
             'code_sha256':sha(__file__),'arithmetic_code_sha256':sha(Path(__file__).with_name('gaussian_transport.py')),
             'bits':bits,'start':['-64','0'],'seed_end':['-64','1'],'initial_matrix':'I_6',
             'basis':'constant D times polarized frame','D':D_VALUES,'segments':[],
             'global_loop_transport_certified':False,'absolute_periods_certified':False}
    for idx,step in enumerate(plan['steps'][:args.steps]):
        started=time.time();path=dest/f'step_{idx+1:03d}.json'
        c=tuple(map(F,step['start']));h=tuple(map(F,step['step']))
        q,N,scale=local_polynomials(source,c,h);lower,ds=disk_lower(c,h,scale)
        assert [str(v) for v in ds]==step['pole_distance_lower_on_disk']
        if path.exists():
            old=json.loads(path.read_text(encoding='utf-8'))
            assert old['system_sha256']==summary['source_system_sha256'] and old['plan_sha256']==sha(planp)
            assert old['start']==step['start'] and old['end']==step['end'] and old['bits']==bits
            P=[decode(M) for M in old['polynomial_coefficients_gaussian_integer']]
        else: P=proposal(q,N,order,bits)
        endpoint,bound=bounds(q,N,P,bits,lower)
        le=dyadic_upper(bound['error_norm_upper'])
        assert le<F(2)**-120, 'proposal needs refinement; no transport result written'
        C,e,rd=compose(endpoint,le,C,e,bits)
        entry={'status':'exact_complex_residual_step_certified','index':idx+1,
               'system_sha256':summary['source_system_sha256'],'plan_sha256':sha(planp),
               'start':step['start'],'end':step['end'],'normalized_coordinate':'s=start+step*t, endpoint t=1',
               'bits':bits,'order':len(P)-1,'common_q_N_scale':str(scale),
               'q_coefficients_gaussian_integer':encode(q),'N_coefficients_gaussian_integer':[encode(M) for M in N],
               'polynomial_coefficients_gaussian_integer':[encode(M) for M in P],
               'bounds':{k:str(v) for k,v in bound.items()},'local_error_upper':str(le),
               'endpoint_center_gaussian_integer':encode(endpoint),
               'propagated_center_gaussian_integer':encode(C),'propagated_error_upper':str(e),
               'composition_rounding_error':str(rd),'floating_proposal':'mpmath 100 decimal digits, untrusted',
               'certificate_arithmetic':'Gaussian integers and rational inequalities only',
               'rigorous_source':'exact residual, denominator distance, Gronwall and matrix-product error bound'}
        if not path.exists(): path.write_text(json.dumps(entry,indent=2)+'\n',encoding='utf-8')
        else: assert old==entry, 'existing certificate differs; investigate rather than overwrite'
        summary['segments'].append({'file':f'high_trunk/{path.name}','sha256':sha(path),'end':step['end'],
                                    'local_error_upper':str(le),'propagated_error_upper':str(e)})
        summary.update(end=step['end'],endpoint_center_gaussian_integer=encode(C),error_upper=str(e))
        (FOLDER/'high_trunk_transport.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
        print(json.dumps({'step':idx+1,'end':step['end'],'seconds':round(time.time()-started,2),
                          'matrix_norm_upper_display':float(bound['matrix_norm_upper']),
                          'local_error_display':float(le),'propagated_error_display':float(e)}),flush=True)
    summary['status']='requested_high_trunk_segments_certified'
    (FOLDER/'high_trunk_transport.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__': main()
