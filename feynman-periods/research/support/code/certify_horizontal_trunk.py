"""Certify one physical route, simultaneously normalized at two basepoints.

The old (-64)-normalization includes its certified seed error. The high
basepoint normalization starts exactly at I and is useful for word traces.
"""
from refined_transport import *
from certify_high_trunk import proposal
import time

def main():
    dest=FOLDER/'horizontal_trunk';planp=dest/'plan.json';plan=json.loads(planp.read_text())
    seedp=FOLDER/'high_trunk_transport.json';seed=json.loads(seedp.read_text())
    assert plan['system_sha256']==seed['source_system_sha256']==sha(FOLDER/'polarized_system.json')
    assert plan['seed_sha256']==sha(seedp)
    bits=seed['bits'];assert bits==240
    source=source_polynomials();order=64
    physical=decode(seed['endpoint_center_gaussian_integer']);ep=F(seed['error_upper'])
    high=mi(1<<bits);eh=F(0)
    deps=['gaussian_transport.py','refined_transport.py','certify_high_trunk.py']
    summary={'status':'partial_horizontal_transport','source_system_sha256':plan['system_sha256'],
             'plan_sha256':sha(planp),'seed_sha256':sha(seedp),'code_sha256':sha(__file__),
             'dependency_sha256':{n:sha(ROOT/'code'/n) for n in deps},'bits':bits,
             'physical_basepoint':['-64','0'],'high_basepoint':plan['start'],'D':D_VALUES,
             'radius_policy':'ceil(radius*2^240)/2^240','segments':[],'landmarks':[],
             'closed_loop_certified':False,'absolute_periods_certified':False}
    for ix,step in enumerate(plan['steps'],1):
        begin=time.time();p=dest/f'step_{ix:03d}.json'
        c=tuple(map(F,step['start']));h=tuple(map(F,step['step']))
        q,N,scale=local_polynomials(source,c,h);lower,ds=disk_lower(c,h,scale)
        assert list(map(str,ds))==step['pole_distance_lower_on_disk']
        if p.exists():
            old=json.loads(p.read_text());assert old['plan_sha256']==sha(planp)
            assert old['system_sha256']==plan['system_sha256'] and old['bits']==bits
            P=[decode(M) for M in old['polynomial_coefficients_gaussian_integer']]
        else: P=proposal(q,N,order,bits)
        L,bd=bounds(q,N,P,bits,lower);le=radius_up(bd['error_norm_upper'],bits)
        assert le<F(2)**-120, f'local residual too large at step {ix}; no certificate saved'
        physical,ep,rp=compose_fine(L,le,physical,ep,bits)
        high,eh,rh=compose_fine(L,le,high,eh,bits)
        out={'status':'complex_residual_horizontal_step_certified','index':ix,
             'system_sha256':plan['system_sha256'],'plan_sha256':sha(planp),
             'start':step['start'],'end':step['end'],'step':step['step'],'bits':bits,'order':len(P)-1,
             'common_q_N_scale':str(scale),'q_coefficients_gaussian_integer':encode(q),
             'N_coefficients_gaussian_integer':[encode(M) for M in N],
             'polynomial_coefficients_gaussian_integer':[encode(M) for M in P],
             'bounds':{k:str(v) for k,v in bd.items()},'local_error_upper':str(le),
             'endpoint_center_gaussian_integer':encode(L),
             'physical_center_gaussian_integer':encode(physical),'physical_error_upper':str(ep),
             'physical_rounding_error':str(rp),
             'high_center_gaussian_integer':encode(high),'high_error_upper':str(eh),
             'high_rounding_error':str(rh),'landmark_pole':step['landmark_pole'],
             'float_role':'untrusted proposal only; all accepting bounds exact'}
        if p.exists(): assert old==out, 'saved step differs; investigate rather than overwrite'
        else: save_json(p,out)
        entry={'file':p.name,'sha256':sha(p),'end':step['end'],'physical_error_upper':str(ep),'high_error_upper':str(eh)}
        summary['segments'].append(entry)
        if step['landmark_pole'] is not None: summary['landmarks'].append(dict(entry,pole=step['landmark_pole']))
        summary.update(end=step['end'],physical_center_gaussian_integer=encode(physical),physical_error_upper=str(ep),
                       high_center_gaussian_integer=encode(high),high_error_upper=str(eh))
        save_json(dest/'transport.json',summary)
        if ix==1 or ix%8==0 or step['landmark_pole'] is not None:
            print(json.dumps({'step':ix,'end':step['end'],'seconds':round(time.time()-begin,2),
                              'M_display':float(bd['matrix_norm_upper']),'physical_error_display':float(ep),
                              'high_error_display':float(eh)}),flush=True)
    summary['status']='horizontal_trunk_to_all_nine_landmarks_certified';save_json(dest/'transport.json',summary)

if __name__=='__main__': main()
