"""Certify a dyadic approximate inverse using an exact Neumann residual."""
from gaussian_transport import *
import mpmath as mp

def verify(d,cert):
    S=1<<d['bits'];C=decode(d['endpoint_center_gaussian_integer']);e=F(d['error_upper'])
    V=decode(cert['inverse_center_gaussian_integer'])
    E=ma(mi(S*S),ms((-1,0),mm(V,C)))
    rn=F(mnorm(E),S*S);vn=F(mnorm(V),S)
    eta=rn+vn*e
    assert 0<=eta<1
    bound=eta*vn/(1-eta)
    assert cert['residual_upper']==str(rn) and cert['eta_upper']==str(eta)
    assert F(cert['inverse_error_upper'])>=bound
    return eta,bound

def main():
    p=FOLDER/'high_trunk_transport.json';d=json.loads(p.read_text(encoding='utf-8'))
    replay=json.loads((FOLDER/'high_trunk_replay.json').read_text(encoding='utf-8'))
    assert replay['certificate_sha256']==sha(p)
    S=1<<d['bits'];C=decode(d['endpoint_center_gaussian_integer']);mp.mp.dps=100
    A=mp.matrix([[mp.mpc(*C[6*i+j])/S for j in range(6)] for i in range(6)])
    proposed=A**-1
    V=[(int(mp.nint(mp.re(v)*S)),int(mp.nint(mp.im(v)*S))) for v in proposed]
    E=ma(mi(S*S),ms((-1,0),mm(V,C)))
    rn=F(mnorm(E),S*S);vn=F(mnorm(V),S);eta=rn+vn*F(d['error_upper'])
    assert eta<1
    bound=eta*vn/(1-eta)
    out={'status':'endpoint_true_transport_inverse_certified','transport_sha256':sha(p),
         'code_sha256':sha(__file__),'end':d['end'],'bits':d['bits'],
         'inverse_center_gaussian_integer':encode(V),'residual_upper':str(rn),
         'eta_upper':str(eta),'inverse_error_upper':str(dyadic_upper(bound)),
         'proof':'||I-VU|| <= ||I-VC||+||V|| e = eta < 1; ||U^-1-V|| <= eta ||V||/(1-eta)',
         'interpretation':'reverse open path from endpoint to -64; not absolute comparison',
         'floating_point_in_bounds':False,'global_loop_certified':False}
    verify(d,out)
    (FOLDER/'high_trunk_inverse.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':out['status'],'error_upper':out['inverse_error_upper'],
                      'error_display':float(F(out['inverse_error_upper'])),'eta_display':float(eta)}),flush=True)

if __name__=='__main__': main()
