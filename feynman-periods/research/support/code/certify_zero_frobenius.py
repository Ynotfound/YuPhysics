"""Exact Frobenius gauge with a weighted nilpotent-adjoint tail."""
from zero_frobenius import *
import time


def main():
    dest=FOLDER/'zero_monodromy';dest.mkdir(exist_ok=True);target=dest/'frobenius.json'
    if target.exists():
        print('existing zero Frobenius certificate retained');return
    start=time.time();Q,Qi,N,source,regular,toD=zero_model();h=(F(0),F(1,4))
    q,H,scale=local_polynomials(regular,Z,h);lower,ds=node_disk_lower(0,h,scale)
    print(json.dumps({'stage':'exact_nilpotent_gauge_coefficients','order':128,
                      'regular_norm_display':float(F(cnorm(H))/lower)}),flush=True)
    P=exact_gauge_coefficients(q,H,real_matrix(N),128)
    bd=gauge_bounds(q,H,real_matrix(N),P,lower)
    C,e,rd=endpoint_gauge(P,bd['gauge_error_upper'],320)
    assert e<F(2)**-210
    mat=lambda M:[[str(v) for v in row] for row in M.tolist()]
    out={'status':'zero_nilpotent_Frobenius_gauge_certified','endpoint':['0','1/4'],'h':['0','1/4'],
         'system_sha256':sha(FOLDER/'polarized_system.json'),
         'geometry_sha256':sha(FOLDER/'polarized_geometry_check.json'),
         'code_sha256':sha(__file__),'arithmetic_sha256':sha(ROOT/'code/zero_frobenius.py'),
         'Q':mat(Q),'Q_inverse':mat(Qi),'nilpotent':mat(N),'to_D':mat(toD),
         'full_source_q':[str(v) for v in source[0]],
         'full_source_N':[[str(v) for v in M] for M in source[1]],
         'order':128,'bits':320,'common_q_N_scale':str(scale),
         'q_coefficients_gaussian_integer':encode(q),'H_coefficients_gaussian_integer':[encode(M) for M in H],
         'gauge_coefficients_exact':[encode(M) for M in P],'bounds':{k:str(v) for k,v in bd.items()},
         'endpoint_center_gaussian_integer':encode(C),'endpoint_error_upper':str(e),
         'endpoint_rounding_error':str(rd),
         'norm_convention':'usual complex matrix infinity norm; l1 Gaussian modulus is only an upper bound',
         'local_monodromy_formula':'G(1)*(I+2*pi*i*N)*G(1)^-1',
         'based_zero_monodromy_certified':False,'absolute_periods_certified':False}
    save_json(target,out)
    print(json.dumps({'status':out['status'],'gauge_error_dyadic':str(dyadic_upper(e)),
                      'kappa_display':float(bd['kappa']),'seconds':round(time.time()-start,2)}),flush=True)


if __name__=='__main__': main()
