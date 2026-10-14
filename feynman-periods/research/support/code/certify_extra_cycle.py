"""Recover another geometric vanishing sphere in the certified Q6 frame."""
from node_branch_transport import propagate_vector
from node_frobenius import *
from certify_five_node_traces import enclosure
from math import floor,isqrt
import argparse


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--node',type=int,required=True);args=ap.parse_args();a=args.node
    dest=FOLDER/'marked_lattice';bp=dest/'marked_zero.json';b=json.loads(bp.read_text())
    assert json.loads(bp.with_name('marked_zero_replay.json').read_text())['certificate_sha256']==sha(bp)
    np=FOLDER/f'node_branches/node_{a}/transport.json';n=json.loads(np.read_text())
    assert json.loads(np.with_name('replay.json').read_text())['transport_sha256']==sha(np)
    bits=b['bits'];S=1<<bits;ref=n['based_reflection'];qrep=json.loads((FOLDER/'q6_traces_replay.json').read_text())
    R=ma([gs(F(1,1<<ref['bits']),v) for v in decode(ref['reflection_center_gaussian_integer'])],ms((-1,0),mi()))
    re=F(ref['reflection_error_upper']);matrices=[];radii=[];sources=[]
    for node in [1,9,25,49,81,0]:
        p=FOLDER/f'node_branches/node_{node}/transport.json' if node else FOLDER/'zero_monodromy/transport.json'
        d=json.loads(p.read_text());rf=d['based_reflection'] if node else d['sixth_reflection']
        matrices.append(ma([gs(F(1,1<<rf['bits']),v) for v in decode(rf['reflection_center_gaussian_integer'])],ms((-1,0),mi())))
        radii.append(F(rf['reflection_error_upper']));sources.append({'node':node,'sha256':sha(p)})
    matrices.append(R);radii.append(re)
    traces=[enclosure((i,6),matrices,radii) for i in range(6)]
    squares=[v['unique_integer_if_integral'] for v in traces]
    assert all(v is not None and v>=0 and isqrt(v)**2==v for v in squares)
    anchor=next(i for i,v in enumerate(squares) if v)
    q=isqrt(squares[anchor]) # Choose the new sphere's orientation through this positive pairing.
    P=decode(b['basis_center_gaussian_integer']);col=[P[6*i+anchor] for i in range(6)]
    ce=F(b['column_error_upper'][anchor])
    Rint=[gs(S,v) for v in R];assert all(x.denominator==y.denominator==1 for x,y in Rint)
    acted,ae,ard=propagate_vector(Rint,re,col,ce,bits)
    exact=[gs(F(1,q),v) for v in acted]
    v=[(dyadic_round(F(x,S),bits),dyadic_round(F(y,S),bits)) for x,y in exact]
    rd=max(F(ab(ga(x,gn(y))),S) for x,y in zip(exact,v));ve=radius_up(ae/q+rd,bits)
    Vi=decode(b['basis_inverse']['center_gaussian_integer']);ie=F(b['basis_inverse']['error_upper'])
    coords,e,crd=propagate_vector(Vi,ie,v,ve,bits)
    exactcoords=[];intervals=[]
    for c in coords:
        z=gs(F(1,S),c);rl=[2*(z[0]-e),2*(z[0]+e)];il=[2*(z[1]-e),2*(z[1]+e)]
        lo,hi=ceil(rl[0]),floor(rl[1]);assert lo==hi and il[0]<=0<=il[1]
        exactcoords.append(sp.Rational(lo,2))
        intervals.append({'twice_real':list(map(str,rl)),'twice_imag':list(map(str,il)),'unique_integer':lo})
    Q=sp.Matrix(qrep['geometric_Gram_up_to_specified_orientations']);x=sp.Matrix(exactcoords);pairings=Q*x
    assert (x.T*Q*x)[0]==-2 and pairings[anchor]==q
    assert [v*v for v in pairings]==squares
    coset=[int(2*v)%2 for v in x];is_outside=any(coset)
    if is_outside:assert coset in b['actual_T0_fixed_classes']
    out={'status':'additional_geometric_cycle_coordinates_certified','node':a,'bits':bits,'code_sha256':sha(__file__),
         'marked_basis_sha256':sha(bp),'node_transport_sha256':sha(np),'reference_sources':sources,
         'trace_pair_enclosures':traces,'anchor_cycle':anchor+1,'anchor_pairing_positive':q,
         'action_center_gaussian_integer':encode(acted),'action_error_upper':str(ae),'action_rounding_error':str(ard),
         'scaled_vector_center_gaussian_integer':encode(v),'scaled_vector_error_upper':str(ve),'scaling_rounding_error':str(rd),
         'coordinate_center_gaussian_integer':encode(coords),'coordinate_error_upper':str(e),'coordinate_rounding_error':str(crd),
         'half_integer_intervals':intervals,'geometric_coordinates':list(map(str,x)),
         'pairings_with_six_cycles':list(map(str,pairings)),'square':-2,
         'twice_coordinate_parity':coset,'outside_L':is_outside,
         'identified_saturation_coset':coset if is_outside else None,'absolute_periods_certified':False}
    save_json(dest/f'cycle_{a}.json',out)
    print(json.dumps({'node':a,'coordinates':out['geometric_coordinates'],'pairings':out['pairings_with_six_cycles'],
                      'coordinate_error_dyadic':str(dyadic_upper(e)),'outside_L':is_outside,'coset':coset}))


if __name__=='__main__':main()
