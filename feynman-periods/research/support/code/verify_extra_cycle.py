"""Independent scalar trace and half-integer coordinate replay."""
from node_frobenius import *
from verify_node_branches import check_action
from math import floor,isqrt
import argparse


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--node',type=int,required=True);a=ap.parse_args().node
    p=FOLDER/f'marked_lattice/cycle_{a}.json';d=json.loads(p.read_text())
    bp=p.with_name('marked_zero.json');b=json.loads(bp.read_text())
    np=FOLDER/f'node_branches/node_{a}/transport.json';n=json.loads(np.read_text())
    assert d['node']==a and d['code_sha256']==sha(ROOT/'code/certify_extra_cycle.py')
    assert d['marked_basis_sha256']==sha(bp) and d['node_transport_sha256']==sha(np)
    assert json.loads(bp.with_name('marked_zero_replay.json').read_text())['certificate_sha256']==sha(bp)
    assert json.loads(np.with_name('replay.json').read_text())['transport_sha256']==sha(np)
    bits=d['bits'];S=1<<bits;ref=n['based_reflection'];re=F(ref['reflection_error_upper'])
    R=ma([gs(F(1,1<<ref['bits']),v) for v in decode(ref['reflection_center_gaussian_integer'])],ms((-1,0),mi()))
    values=[]
    for i,(src,t) in enumerate(zip(d['reference_sources'],d['trace_pair_enclosures'])):
        node=src['node'];path=FOLDER/f'node_branches/node_{node}/transport.json' if node else FOLDER/'zero_monodromy/transport.json'
        assert sha(path)==src['sha256'] and t['indices']==[i+1,7]
        tr=json.loads(path.read_text());rf=tr['based_reflection'] if node else tr['sixth_reflection']
        A=ma([gs(F(1,1<<rf['bits']),v) for v in decode(rf['reflection_center_gaussian_integer'])],ms((-1,0),mi()))
        ae=F(rf['reflection_error_upper'])
        c=sumg(gm(A[6*j+k],R[6*k+j]) for j in range(6) for k in range(6))
        e=6*((mnorm(A)+ae)*re+ae*mnorm(R))
        assert encode([c])[0]==t['center'] and e==F(t['error_upper'])
        rl=[c[0]-e,c[0]+e];il=[c[1]-e,c[1]+e]
        assert list(map(str,rl))==t['real_interval'] and list(map(str,il))==t['imaginary_interval']
        lo,hi=ceil(rl[0]),floor(rl[1]);assert lo==hi==t['unique_integer_if_integral'] and il[0]<=0<=il[1]
        assert lo>=0 and isqrt(lo)**2==lo;values.append(lo)
    anchor=d['anchor_cycle']-1;q=d['anchor_pairing_positive'];assert q>0 and q*q==values[anchor]
    P=decode(b['basis_center_gaussian_integer']);col=[P[6*i+anchor] for i in range(6)];ce=F(b['column_error_upper'][anchor])
    raw=decode(d['action_center_gaussian_integer']);ae=F(d['action_error_upper'])
    assert check_action([gs(S,v) for v in R],re,col,ce,raw,ae,bits)==F(d['action_rounding_error'])
    v=decode(d['scaled_vector_center_gaussian_integer']);ve=F(d['scaled_vector_error_upper'])
    rd=max(F(ab(ga(gs(F(1,q),x),gn(y))),S) for x,y in zip(raw,v))
    assert rd==F(d['scaling_rounding_error']) and ae/q+rd<=ve<ae/q+rd+F(1,S)
    Vi=decode(b['basis_inverse']['center_gaussian_integer']);ie=F(b['basis_inverse']['error_upper'])
    c=decode(d['coordinate_center_gaussian_integer']);e=F(d['coordinate_error_upper'])
    assert check_action(Vi,ie,v,ve,c,e,bits)==F(d['coordinate_rounding_error'])
    values2=[]
    for z,t in zip(c,d['half_integer_intervals']):
        z=gs(F(1,S),z);rl=[2*(z[0]-e),2*(z[0]+e)];il=[2*(z[1]-e),2*(z[1]+e)]
        assert list(map(str,rl))==t['twice_real'] and list(map(str,il))==t['twice_imag']
        lo,hi=ceil(rl[0]),floor(rl[1]);assert lo==hi==t['unique_integer'] and il[0]<=0<=il[1]
        values2.append(sp.Rational(lo,2))
    x=sp.Matrix(values2);Q=sp.Matrix(json.loads((FOLDER/'q6_traces_replay.json').read_text())['geometric_Gram_up_to_specified_orientations'])
    assert list(map(str,x))==d['geometric_coordinates'] and (x.T*Q*x)[0]==-2
    pairings=Q*x;assert list(map(str,pairings))==d['pairings_with_six_cycles']
    assert [v*v for v in pairings]==values and pairings[anchor]==q
    coset=[int(2*v)%2 for v in x];outside=any(coset)
    assert coset==d['twice_coordinate_parity'] and outside==d['outside_L']
    assert d['identified_saturation_coset']==(coset if outside else None)
    out={'status':'additional_geometric_cycle_independently_replayed','node':a,'certificate_sha256':sha(p),
         'verifier_sha256':sha(__file__),'geometric_coordinates':d['geometric_coordinates'],
         'coordinate_error_dyadic':str(dyadic_upper(e)),'outside_L':outside,
         'identified_saturation_coset':coset if outside else None,
         'floating_point_used':False,'absolute_periods_certified':False}
    save_json(p.with_name(f'cycle_{a}_replay.json'),out);print(json.dumps(out))


if __name__=='__main__':main()
