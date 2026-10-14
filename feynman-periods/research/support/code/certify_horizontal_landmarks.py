"""Recombine the common seed once; certify inverses at nine upper landmarks."""
from refined_transport import *
import mpmath as mp

def inverse(C,e,bits):
    S=1<<bits;mp.mp.dps=100
    M=mp.matrix([[mp.mpc(*C[6*i+j])/S for j in range(6)] for i in range(6)])
    V=[(int(mp.nint(mp.re(x)*S)),int(mp.nint(mp.im(x)*S))) for x in M**-1]
    residual=ma(mi(S*S),ms((-1,0),mm(V,C)))
    rn=F(mnorm(residual),S*S);vn=F(mnorm(V),S);eta=rn+vn*e
    assert eta<1
    return {'center_gaussian_integer':encode(V),'residual_upper':str(rn),'eta_upper':str(eta),
            'error_upper':str(radius_up(eta*vn/(1-eta),bits))}

def main():
    dest=FOLDER/'horizontal_trunk';tp=dest/'transport.json';d=json.loads(tp.read_text())
    rp=dest/'replay.json';rep=json.loads(rp.read_text());assert rep['transport_sha256']==sha(tp)
    seedp=FOLDER/'high_trunk_transport.json';seed=json.loads(seedp.read_text())
    assert d['seed_sha256']==sha(seedp)
    bits=d['bits'];SC=decode(seed['endpoint_center_gaussian_integer']);se=F(seed['error_upper'])
    items=[]
    for point in d['landmarks']:
        f=dest/point['file'];assert sha(f)==point['sha256'];step=json.loads(f.read_text())
        H=decode(step['high_center_gaussian_integer']);he=F(step['high_error_upper'])
        C,e,rounding=compose_fine(H,he,SC,se,bits)
        ivh=inverse(H,he,bits);ivb=inverse(C,e,bits)
        v={'pole':point['pole'],'point':point['end'],'segment_file':point['file'],'segment_sha256':sha(f),
           'high_center_gaussian_integer':encode(H),'high_error_upper':str(he),'high_inverse':ivh,
           'physical_factored_center_gaussian_integer':encode(C),'physical_factored_error_upper':str(e),
           'physical_factored_rounding_error':str(rounding),'physical_inverse':ivb,
           'physical_sequential_error_upper':step['physical_error_upper']}
        items.append(v)
        print(json.dumps({'pole':v['pole'],'high_error_display':float(he),'physical_factored_error_display':float(e),
                          'physical_sequential_error_display':float(F(step['physical_error_upper'])),
                          'high_inverse_error_display':float(F(ivh['error_upper']))}),flush=True)
    out={'status':'nine_upper_landmark_transports_and_inverses_certified','code_sha256':sha(__file__),
         'transport_sha256':sha(tp),'transport_replay_sha256':sha(rp),'seed_sha256':sha(seedp),
         'bits':bits,'physical_basepoint':['-64','0'],'high_basepoint':['-64','64'],'landmarks':items,
         'physical_error_method':'one exact recombination of the shared seed, retaining its original error',
         'trace_rebasing':'simultaneous conjugation only; no geometric absolute initial periods',
         'closed_loops_certified':False,'absolute_periods_certified':False}
    save_json(dest/'landmarks.json',out)

if __name__=='__main__': main()
