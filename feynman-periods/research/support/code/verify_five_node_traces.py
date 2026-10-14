"""Independent cyclic-index evaluation and exact trace enclosure replay."""
from node_frobenius import *
from itertools import product
from math import floor,isqrt


def main():
    p=FOLDER/'node_branches/five_node_traces.json';d=json.loads(p.read_text())
    assert d['code_sha256']==sha(ROOT/'code/certify_five_node_traces.py')
    matrices=[];radii=[]
    for source in d['sources']:
        a=source['node'];folder=FOLDER/'node_branches'/f'node_{a}'
        assert sha(folder/'transport.json')==source['transport_sha256']
        assert sha(folder/'replay.json')==source['replay_sha256']
        transport=json.loads((folder/'transport.json').read_text())
        assert transport['based_point']==d['basepoint']==['-64','64']
        ref=transport['based_reflection'];scale=1<<ref['bits']
        M=[gs(F(1,scale),v) for v in decode(ref['reflection_center_gaussian_integer'])]
        matrices.append(ma(M,ms((-1,0),mi())));radii.append(F(ref['reflection_error_upper']))
    maximum_error=F(0);values={}
    for entry in d['pair_traces']+d['triple_traces']:
        ix=[i-1 for i in entry['indices']];k=len(ix);assert k in (2,3)
        center=Z
        for indices in product(range(6),repeat=k):
            term=(F(1),F(0))
            for j in range(k):
                term=gm(term,matrices[ix[j]][6*indices[j]+indices[(j+1)%k]])
            center=ga(center,term)
        assert encode([center])[0]==entry['center']
        # Telescoping the product, independent of subtracting two products.
        e=F(0)
        for j in range(k):
            term=radii[ix[j]]
            for l in range(k):
                if l<j: term*=mnorm(matrices[ix[l]])+radii[ix[l]]
                if l>j: term*=mnorm(matrices[ix[l]])
            e+=term
        e*=6;assert e==F(entry['error_upper'])
        re=[center[0]-e,center[0]+e];im=[center[1]-e,center[1]+e]
        assert list(map(str,re))==entry['real_interval']
        assert list(map(str,im))==entry['imaginary_interval']
        lo,hi=ceil(re[0]),floor(re[1]);unique=lo if lo==hi and im[0]<=0<=im[1] else None
        assert unique==entry['unique_integer_if_integral']
        assert unique is not None
        values[tuple(entry['indices'])]=unique
        maximum_error=max(maximum_error,e)
    # Reconstruct the first five cycle pairings up to independent orientations.
    # Use first-cycle edges to choose orientations, then the two triangles.
    gram=sp.diag(*([-2]*5));signs=[-1,-1,1,-1]
    for j in range(1,5):
        square=values[(1,j+1)];q=isqrt(square);assert q*q==square and q>0
        gram[0,j]=gram[j,0]=signs[j-1]*q
    for i in range(1,5):
        for j in range(i+1,5):
            square=values[(i+1,j+1)]
            if square==0: gram[i,j]=gram[j,i]=0
            else:
                v=F(values[(1,i+1,j+1)],int(gram[0,i]*gram[0,j]))
                assert v.denominator==1 and v*v==square
                gram[i,j]=gram[j,i]=int(v)
    candidate=json.loads((ROOT/'data/recovered_lattice.json').read_text())
    assert d['candidate_sha256']==sha(ROOT/'data/recovered_lattice.json')
    assert gram==sp.Matrix(candidate['Q6'])[:5,:5]
    assert d['all_twelve_enclosures_isolate_an_integer'] and d['all_candidate_targets_match']
    out={'status':'ten_pair_and_two_triple_trace_enclosures_independently_replayed',
         'trace_certificate_sha256':sha(p),'verifier_sha256':sha(__file__),
         'trace_center_method':'explicit cyclic scalar sums','trace_error_method':'telescoping product norm',
         'maximum_trace_error_upper':str(maximum_error),'maximum_trace_error_dyadic':str(dyadic_upper(maximum_error)),
         'first_five_Gram_up_to_orientation':[[int(v) for v in row] for row in gram.tolist()],
         'first_five_Gram_determinant':str(gram.det()),'first_five_Gram_rank':gram.rank(),
         'geometry_required':'Integral vanishing spheres of square -2 and comparison with this pure K3 connection.',
         'Q6_21_traces_certified':False,'zero_monodromy_certified':False,'absolute_periods_certified':False}
    save_json(p.with_name('five_node_traces_replay.json'),out)
    print(json.dumps(out))


if __name__=='__main__': main()
