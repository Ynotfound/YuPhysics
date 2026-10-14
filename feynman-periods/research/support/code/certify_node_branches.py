"""Certify ascending node vectors, then rebase using saved horizontal inverses."""
from node_branch_transport import *
import argparse,time


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--nodes',type=int,nargs='+',default=[1,9,25,49,81])
    parser.add_argument('--steps',type=int,default=64)
    args=parser.parse_args();assert 1<=args.steps<=64
    planp=FOLDER/'node_frobenius/branch_plan.json';plan=json.loads(planp.read_text())
    rp=planp.with_name('branch_plan_replay.json')
    assert json.loads(rp.read_text())['plan_sha256']==sha(planp)
    hrp=FOLDER/'horizontal_trunk/landmarks.json';hor=json.loads(hrp.read_text())
    source=source_polynomials();bits=240;S=1<<bits;order=64
    for a in args.nodes:
        started=time.time();seedp=FOLDER/f'node_frobenius/node_{a}.json'
        seed=json.loads(seedp.read_text());Q=sp.Matrix(seed['Q']);Qi=sp.Matrix(seed['Q_inverse'])
        qsource=transformed_source(source,Q,Qi)
        branch=next(b for b in plan['branches'] if b['node']==a)
        assert branch['local_sha256']==sha(seedp)
        dest=FOLDER/'node_branches'/f'node_{a}';dest.mkdir(parents=True,exist_ok=True)
        C,e,initial_rounding=initial_Q_vector(seed,bits);current_frame='Q'
        summary={'status':'partial_node_branch_vector','node':a,'bits':bits,
                 'system_sha256':sha(FOLDER/'polarized_system.json'),
                 'seed_sha256':sha(seedp),'plan_sha256':sha(planp),
                 'landmarks_sha256':sha(hrp),'code_sha256':sha(__file__),
                 'arithmetic_sha256':sha(ROOT/'code/node_branch_transport.py'),
                 'source_arithmetic_sha256':sha(ROOT/'code/gaussian_transport.py'),
                 'start':seed['endpoint'],'initial_frame':'Q',
                 'switch_height':str(SWITCH_HEIGHT[a]),
                 'initial_rounding_error':str(initial_rounding),
                 'initial_error_upper':str(e),'segments':[],
                 'absolute_periods_certified':False}
        for ix,st in enumerate(list(reversed(branch['steps_top_to_bottom']))[:args.steps]):
            tick=time.time();up=st['up'];c=tuple(map(F,up['center']));h=tuple(map(F,up['step']))
            switched=False;switch_rounding=F(0)
            if current_frame=='Q' and c[1]>=SWITCH_HEIGHT[a]:
                C,e,switch_rounding=switch_vector(C,e,Q,bits);current_frame='D';switched=True
            src=qsource if current_frame=='Q' else source
            q,N,scale=local_polynomials(src,c,h);lower,ds=disk_lower(c,h,scale)
            assert list(map(str,ds))==up['pole_distance_lower']
            path=dest/f'step_{ix+1:03d}.json'
            if path.exists():
                old=json.loads(path.read_text());assert old['system_sha256']==summary['system_sha256']
                assert old['frame']==current_frame and old['plan_sha256']==sha(planp)
                P=[decode(M) for M in old['polynomial_coefficients_gaussian_integer']]
            else:
                P=proposal(q,N,order,bits)
            L,bd=bounds(q,N,P,bits,lower);le=radius_up(bd['error_norm_upper'],bits)
            assert le<F(2)**-135, f'node {a}, step {ix+1}: refine proposal before saving'
            new,enew,rd=propagate_vector(L,le,C,e,bits)
            entry={'status':'node_vector_ordinary_step_certified','node':a,'index':ix+1,
                   'system_sha256':summary['system_sha256'],'plan_sha256':sha(planp),
                   'seed_sha256':sha(seedp),'start':st['end'],'end':st['start'],
                   'frame':current_frame,'switched_before_step':switched,
                   'switch_rounding_error':str(switch_rounding),
                   'input_center_gaussian_integer':encode(C),'input_error_upper':str(e),
                   'bits':bits,'order':len(P)-1,'common_q_N_scale':str(scale),
                   'q_coefficients_gaussian_integer':encode(q),
                   'N_coefficients_gaussian_integer':[encode(M) for M in N],
                   'polynomial_coefficients_gaussian_integer':[encode(M) for M in P],
                   'bounds':{k:str(v) for k,v in bd.items()},
                   'local_matrix_error_upper':str(le),'transition_center_gaussian_integer':encode(L),
                   'output_center_gaussian_integer':encode(new),'output_error_upper':str(enew),
                   'vector_rounding_error':str(rd),'proposal_precision':'mpmath 100 decimal digits, untrusted'}
            if path.exists(): assert old==entry, 'saved step differs; do not overwrite'
            else: save_json(path,entry)
            C,e=new,enew
            summary['segments'].append({'file':path.name,'sha256':sha(path),'end':entry['end']})
            summary.update(end=entry['end'],frame=current_frame,center_gaussian_integer=encode(C),error_upper=str(e))
            save_json(dest/'transport.json',summary)
            if (ix+1)%8==0 or ix==0 or switched:
                print(json.dumps({'node':a,'step':ix+1,'height':entry['end'][1],
                                  'frame':current_frame,'M_display':float(bd['matrix_norm_upper']),
                                  'vector_error_display':float(e),'seconds':round(time.time()-tick,2)}),flush=True)
        if args.steps==64:
            assert current_frame=='D' and summary['end']==[str(a),'64']
            point=next(p for p in hor['landmarks'] if p['pole']==a)
            V=decode(point['high_inverse']['center_gaussian_integer']);ve=F(point['high_inverse']['error_upper'])
            based,be,rd=propagate_vector(V,ve,C,e,bits)
            K=pairing_at(-64,(F(0),F(64)))
            ref=reflection_at_vector([gs(F(1,S),v) for v in based],be,K,bits)
            summary.update(status='node_branch_and_common_basepoint_reflection_certified',
                           based_point=['-64','64'],based_center_gaussian_integer=encode(based),
                           based_error_upper=str(be),rebasing_rounding_error=str(rd),
                           based_pairing_exact=encode(K),based_reflection=ref)
            save_json(dest/'transport.json',summary)
            print(json.dumps({'node':a,'status':summary['status'],
                              'based_reflection_error_dyadic':str(dyadic_upper(F(ref['reflection_error_upper']))),
                              'total_seconds':round(time.time()-started,2)}),flush=True)


if __name__=='__main__': main()
