"""Independent cyclic scalar trace sums, telescoping errors and Gram recovery."""
from node_frobenius import *
from itertools import product,combinations
from math import floor,isqrt


def main():
    p=FOLDER/'q6_traces.json';d=json.loads(p.read_text())
    assert d['code_sha256']==sha(ROOT/'code/certify_q6_traces.py')
    assert d['enclosure_code_sha256']==sha(ROOT/'code/certify_five_node_traces.py')
    assert d['system_sha256']==sha(FOLDER/'polarized_system.json')
    assert [v['node'] for v in d['sources']]==[1,9,25,49,81,0]
    matrices=[];radii=[]
    for source in d['sources']:
        a=source['node'];folder=FOLDER/'node_branches'/f'node_{a}' if a else FOLDER/'zero_monodromy'
        tp=folder/'transport.json';rp=folder/'replay.json'
        assert sha(tp)==source['transport_sha256'] and sha(rp)==source['replay_sha256']
        tr=json.loads(tp.read_text());assert tr['based_point']==d['basepoint']==['-64','64']
        assert json.loads(rp.read_text())['transport_sha256']==sha(tp)
        ref=tr['based_reflection'] if a else tr['sixth_reflection'];S=1<<ref['bits']
        T=[gs(F(1,S),v) for v in decode(ref['reflection_center_gaussian_integer'])]
        matrices.append(ma(T,ms((-1,0),mi())));radii.append(F(ref['reflection_error_upper']))
    assert [v['indices'] for v in d['pair_traces']]==[list(ix) for ix in combinations(range(1,7),2)]
    assert [v['indices'] for v in d['triple_traces']]==[[1,2,4],[1,2,6],[1,3,4],[1,3,6],[1,4,6],[1,5,6]]
    values={};max_error=F(0);checks=[]
    for v in d['pair_traces']+d['triple_traces']:
        ix=[i-1 for i in v['indices']];k=len(ix);center=Z
        for indices in product(range(6),repeat=k):
            term=(F(1),F(0))
            for j in range(k): term=gm(term,matrices[ix[j]][6*indices[j]+indices[(j+1)%k]])
            center=ga(center,term)
        error=F(0)
        for j in range(k):
            term=radii[ix[j]]
            for l in range(k):
                if l<j: term*=mnorm(matrices[ix[l]])+radii[ix[l]]
                if l>j: term*=mnorm(matrices[ix[l]])
            error+=6*term
        assert encode([center])[0]==v['center'] and error==F(v['error_upper'])
        re=[center[0]-error,center[0]+error];im=[center[1]-error,center[1]+error]
        assert list(map(str,re))==v['real_interval'] and list(map(str,im))==v['imaginary_interval']
        lo,hi=ceil(re[0]),floor(re[1]);assert lo==hi and im[0]<=0<=im[1]
        assert lo==v['unique_integer_if_integral']
        values[tuple(v['indices'])]=lo;max_error=max(max_error,error)
        checks.append({'indices':v['indices'],'integer':lo,'error_dyadic':str(dyadic_upper(error))})
    # Orient cycles 2..5 relative to cycle 1. Cycle 6=T0 cycle1 has no extra sign freedom.
    # Since (T0-I)^2=0 and T0 preserves the pairing, <delta1,T0 delta1>=-2.
    gram=sp.diag(*([-2]*6));signs=[-1,-1,1,-1,-1]
    for j in range(1,6):
        q=isqrt(values[(1,j+1)]);assert q*q==values[(1,j+1)] and q>0
        gram[0,j]=gram[j,0]=signs[j-1]*q
    assert gram[0,5]==-2
    for i in range(1,6):
        for j in range(i+1,6):
            square=values[(i+1,j+1)]
            if square==0: gram[i,j]=gram[j,i]=0
            else:
                val=F(values[(1,i+1,j+1)],int(gram[0,i]*gram[0,j]))
                assert val.denominator==1 and val*val==square
                gram[i,j]=gram[j,i]=int(val)
    qp=ROOT/'data/recovered_lattice.json'
    assert sha(qp)==d['candidate_sha256']
    candidate=sp.Matrix(json.loads(qp.read_text())['Q6'])
    assert gram==candidate and gram.det()==192 and gram.rank()==6
    out={'status':'all_21_trace_enclosures_independently_replayed',
         'trace_certificate_sha256':sha(p),'verifier_sha256':sha(__file__),
         'trace_center_method':'explicit cyclic scalar sums','trace_error_method':'telescoping product norm',
         'checks':checks,'maximum_trace_error_upper':str(max_error),
         'maximum_trace_error_dyadic':str(dyadic_upper(max_error)),
         'geometric_Gram_up_to_specified_orientations':[[int(x) for x in row] for row in gram.tolist()],
         'Gram_determinant':str(gram.det()),'Gram_rank':gram.rank(),
         'geometric_inputs':'Integral vanishing spheres of square -2, polarization, and actual T0 with (T0-I)^2=0.',
         'arithmetic_alone_proves_geometric_inputs':False,
         'index_in_T_given_proved_det_T_48':2,
         'marked_saturation_coset_certified':False,'absolute_comparison_certified':False}
    save_json(p.with_name('q6_traces_replay.json'),out)
    print(json.dumps({k:v for k,v in out.items() if k not in ['maximum_trace_error_upper','checks']}))


if __name__=='__main__': main()
