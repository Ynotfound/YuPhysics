"""Actual banana node coefficients, exact weighted tails and local reflections."""
from node_frobenius import *
import argparse,time

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--nodes',type=int,nargs='+',default=[1,9,25,49,81,121,169,225])
    args=parser.parse_args();source=source_polynomials();dest=FOLDER/'node_frobenius';dest.mkdir(exist_ok=True)
    for a in args.nodes:
        p=dest/f'node_{a}.json'
        if p.exists():
            old=json.loads(p.read_text());assert old['system_sha256']==sha(FOLDER/'polarized_system.json')
            print(json.dumps({'node':a,'status':'existing_certificate_retained'}),flush=True);continue
        started=time.time();Q,Qi,R,local=node_model(source,a)
        # Make the contraction genuinely effective in the chosen exact eigenbasis.
        height=F(1,4)
        for attempt in range(20):
            h=(F(0),height);q,N,scale=local_polynomials(local,(F(a),F(0)),h)
            lower,ds=node_disk_lower(a,h,scale);M=F(cnorm(N))/lower
            if M<=32: break
            height/=2
        else: raise AssertionError('no suitable local disk found within bounded search')
        for order in [64,96,128]:
            P=exact_coefficients(q,N,order);bd=fuchs_bounds(q,N,P,lower)
            K=pairing_at(a,h);ref=reflection_from_vector(P,Q,bd['vector_error_Q_basis'],K,240)
            if F(ref['reflection_error_upper'])<=F(2)**-120: break
        else: raise AssertionError('local reflection needs further precision; no result saved')
        matrix_encode=lambda M:[[str(v) for v in row] for row in M.tolist()]
        out={'status':'actual_node_frobenius_tail_and_local_reflection_certified','node':a,
             'system_sha256':sha(FOLDER/'polarized_system.json'),'code_sha256':sha(__file__),
             'arithmetic_sha256':sha(ROOT/'code/node_frobenius.py'),
             'source_arithmetic_sha256':sha(ROOT/'code/gaussian_transport.py'),
             'D':D_VALUES,'Q':matrix_encode(Q),'Q_inverse':matrix_encode(Qi),'residue_in_D_basis':matrix_encode(R),
             'h':list(map(str,h)),'endpoint':[str(a),str(height)],'coordinate':'s=a+h*t, endpoint t=1',
             'order':order,'common_q_N_scale':str(scale),'pole_distance_lower':list(map(str,ds)),
             'q_coefficients_gaussian_integer':encode(q),'N_coefficients_gaussian_integer':[encode(M) for M in N],
             'frobenius_vector_coefficients_exact':[encode(v) for v in P],
             'bounds':{k:str(v) for k,v in bd.items()},'pairing_at_endpoint_exact':encode(K),'reflection':ref,
             'coefficient_arithmetic':'exact Gaussian rational; no numerical Taylor proposal',
             'local_CCW_monodromy_reflection_certified':True,'based_high_monodromy_certified':False,
             'absolute_geometric_amplitude_certified':False,'zero_monodromy_certified_by_this_file':False}
        save_json(p,out)
        print(json.dumps({'node':a,'height':str(height),'order':order,'M_display':float(M),
                          'kappa_display':float(bd['kappa']),'reflection_error_dyadic':str(dyadic_upper(F(ref['reflection_error_upper']))),
                          'seconds':round(time.time()-started,2)}),flush=True)

if __name__=='__main__': main()
