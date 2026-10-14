"""Exact trace enclosures for the five rebased nodal reflections."""
from node_frobenius import *
from itertools import combinations
from math import floor


def enclosure(indices,matrices,radii):
    product=mi();n0=F(1);n1=F(1)
    for i in indices:
        A=matrices[i];norm=mnorm(A)
        product=mm(product,A);n0*=norm;n1*=norm+radii[i]
    center=sumg(product[6*i+i] for i in range(6));error=6*(n1-n0)
    real=[center[0]-error,center[0]+error];imag=[center[1]-error,center[1]+error]
    lo,hi=ceil(real[0]),floor(real[1])
    unique=lo if lo==hi and imag[0]<=0<=imag[1] else None
    return {'indices':[i+1 for i in indices],'center':encode([center])[0],
            'error_upper':str(error),'real_interval':list(map(str,real)),
            'imaginary_interval':list(map(str,imag)),'unique_integer_if_integral':unique}


def main():
    nodes=[1,9,25,49,81];matrices=[];radii=[];sources=[];point=None
    for a in nodes:
        dest=FOLDER/'node_branches'/f'node_{a}';p=dest/'transport.json'
        d=json.loads(p.read_text());r=json.loads((dest/'replay.json').read_text())
        assert r['transport_sha256']==sha(p)
        if point is None: point=d['based_point']
        assert point==d['based_point']==['-64','64']
        ref=d['based_reflection'];S=1<<ref['bits']
        M=[gs(F(1,S),v) for v in decode(ref['reflection_center_gaussian_integer'])]
        matrices.append(ma(M,ms((-1,0),mi())));radii.append(F(ref['reflection_error_upper']))
        sources.append({'node':a,'transport_sha256':sha(p),'replay_sha256':sha(dest/'replay.json')})
    pairs=[enclosure(i,matrices,radii) for i in combinations(range(5),2)]
    triples=[enclosure(i,matrices,radii) for i in [(0,1,3),(0,2,3)]]
    qp=ROOT/'data/recovered_lattice.json';Q=json.loads(qp.read_text())['Q6']
    for t in pairs:
        i,j=[v-1 for v in t['indices']];t['candidate_target']=Q[i][j]**2
        t['candidate_matches_unique_integer']=t['unique_integer_if_integral']==t['candidate_target']
    for t in triples:
        i,j,k=[v-1 for v in t['indices']];t['candidate_target']=Q[i][j]*Q[j][k]*Q[k][i]
        t['candidate_matches_unique_integer']=t['unique_integer_if_integral']==t['candidate_target']
    allcert=all(v['unique_integer_if_integral'] is not None for v in pairs+triples)
    out={'status':'five_node_trace_enclosures_computed','nodes':nodes,'basepoint':point,
         'system_sha256':sha(FOLDER/'polarized_system.json'),'code_sha256':sha(__file__),
         'sources':sources,'candidate_sha256':sha(qp),'pair_traces':pairs,'triple_traces':triples,
         'all_twelve_enclosures_isolate_an_integer':allcert,
         'all_candidate_targets_match':all(v['candidate_matches_unique_integer'] for v in pairs+triples),
         'integrality_basis':'Geometric K3 monodromies preserve the integral transcendental lattice; all products of T_i-I have integer trace.',
         'integrality_verified_by_this_arithmetic_program':False,
         'Q6_21_traces_certified':False,'zero_monodromy_certified':False,
         'absolute_comparison_certified':False}
    save_json(FOLDER/'node_branches/five_node_traces.json',out)
    print(json.dumps({'status':out['status'],'all_unique':allcert,
                      'all_candidates_match':out['all_candidate_targets_match'],
                      'pairs':[(t['indices'],t['unique_integer_if_integral'],float(F(t['error_upper']))) for t in pairs],
                      'triples':[(t['indices'],t['unique_integer_if_integral'],float(F(t['error_upper']))) for t in triples]}))


if __name__=='__main__': main()
