"""Replay factored seed bounds and all 18 inverse enclosures exactly."""
from gaussian_transport import *
from verify_horizontal_trunk import propagated_check

def inverse_check(C,e,cert,bits):
    S=1<<bits;V=decode(cert['center_gaussian_integer'])
    r=[ga((S*S if i==j else 0,0),gn(sumg(gm(V[6*i+k],C[6*k+j]) for k in range(6))))
       for i in range(6) for j in range(6)]
    rn=F(mnorm(r),S*S);vn=F(mnorm(V),S);eta=rn+vn*e
    assert cert['residual_upper']==str(rn) and cert['eta_upper']==str(eta)
    assert 0<=eta<1
    err=F(cert['error_upper']);raw=eta*vn/(1-eta)
    assert raw<=err<raw+F(1,S)
    return err

def main():
    dest=FOLDER/'horizontal_trunk';lp=dest/'landmarks.json';d=json.loads(lp.read_text())
    tp=dest/'transport.json';tr=json.loads(tp.read_text())
    rp=dest/'replay.json';r=json.loads(rp.read_text())
    assert r['transport_sha256']==d['transport_sha256']==sha(tp)
    assert d['transport_replay_sha256']==sha(rp)
    assert d['code_sha256']==sha(ROOT/'code/certify_horizontal_landmarks.py')
    seedp=FOLDER/'high_trunk_transport.json';seed=json.loads(seedp.read_text())
    assert d['seed_sha256']==tr['seed_sha256']==sha(seedp)
    bits=d['bits'];SC=decode(seed['endpoint_center_gaussian_integer']);se=F(seed['error_upper'])
    assert [v['pole'] for v in d['landmarks']]==POLES
    maxima={'high_forward':F(0),'high_inverse':F(0),'physical_factored_forward':F(0),'physical_inverse':F(0)}
    for v,item in zip(d['landmarks'],tr['landmarks']):
        f=dest/v['segment_file'];assert v['segment_sha256']==item['sha256']==sha(f)
        step=json.loads(f.read_text());H=decode(v['high_center_gaussian_integer']);he=F(v['high_error_upper'])
        assert H==decode(step['high_center_gaussian_integer']) and he==F(step['high_error_upper'])
        C=decode(v['physical_factored_center_gaussian_integer']);ce=F(v['physical_factored_error_upper'])
        propagated_check(H,he,SC,se,C,ce,F(v['physical_factored_rounding_error']),bits)
        ih=inverse_check(H,he,v['high_inverse'],bits);ib=inverse_check(C,ce,v['physical_inverse'],bits)
        for key,err in [('high_forward',he),('high_inverse',ih),('physical_factored_forward',ce),('physical_inverse',ib)]:
            maxima[key]=max(maxima[key],err)
    out={'status':'all_nine_landmark_factorizations_and_18_inverses_replayed',
         'landmarks_sha256':sha(lp),'verifier_sha256':sha(__file__),'arithmetic':'exact rational, no numerical inversion',
         'uniform_dyadic_upper_bounds':{k:str(dyadic_upper(v)) for k,v in maxima.items()},
         'exact_maxima':{k:str(v) for k,v in maxima.items()},'closed_loops_certified':False,
         'absolute_periods_certified':False}
    (dest/'landmarks_replay.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k!='exact_maxima'}),flush=True)

if __name__=='__main__': main()
