"""All 21 based trace enclosures, with no absolute-period claim."""
from certify_five_node_traces import enclosure
from node_frobenius import *
from itertools import combinations


def main():
    sources=[];matrices=[];radii=[]
    for a in [1,9,25,49,81,0]:
        folder=FOLDER/'node_branches'/f'node_{a}' if a else FOLDER/'zero_monodromy'
        p=folder/'transport.json';rp=folder/'replay.json'
        d=json.loads(p.read_text());r=json.loads(rp.read_text())
        assert r['transport_sha256']==sha(p) and d['based_point']==['-64','64']
        ref=d['based_reflection'] if a else d['sixth_reflection']
        S=1<<ref['bits']
        T=[gs(F(1,S),v) for v in decode(ref['reflection_center_gaussian_integer'])]
        matrices.append(ma(T,ms((-1,0),mi())));radii.append(F(ref['reflection_error_upper']))
        sources.append({'cycle':len(sources)+1,'node':a,'transport_sha256':sha(p),'replay_sha256':sha(rp)})
    pairs=[enclosure(ix,matrices,radii) for ix in combinations(range(6),2)]
    triples=[enclosure(ix,matrices,radii) for ix in [(0,1,3),(0,1,5),(0,2,3),(0,2,5),(0,3,5),(0,4,5)]]
    qp=ROOT/'data/recovered_lattice.json';Q=json.loads(qp.read_text())['Q6']
    for t in pairs+triples:
        ix=[i-1 for i in t['indices']];target=1
        for j in range(len(ix)): target*=Q[ix[j]][ix[(j+1)%len(ix)]]
        t['candidate_target']=target
        t['candidate_matches_unique_integer']=t['unique_integer_if_integral']==target
    out={'status':'all_21_based_trace_enclosures_computed','basepoint':['-64','64'],
         'system_sha256':sha(FOLDER/'polarized_system.json'),'code_sha256':sha(__file__),
         'enclosure_code_sha256':sha(ROOT/'code/certify_five_node_traces.py'),
         'sources':sources,'candidate_sha256':sha(qp),'pair_traces':pairs,'triple_traces':triples,
         'all_21_enclosures_isolate_an_integer':all(v['unique_integer_if_integral'] is not None for v in pairs+triples),
         'all_candidate_targets_match':all(v['candidate_matches_unique_integer'] for v in pairs+triples),
         'integrality_basis':'Integral K3 Picard-Lefschetz reflections and the actual zero monodromy.',
         'integrality_verified_by_this_arithmetic_program':False,'absolute_comparison_certified':False}
    save_json(FOLDER/'q6_traces.json',out)
    print(json.dumps({k:out[k] for k in ['status','all_21_enclosures_isolate_an_integer','all_candidate_targets_match']}))
    print(json.dumps({'pairs':[(v['indices'],v['unique_integer_if_integral']) for v in pairs],
                      'triples':[(v['indices'],v['unique_integer_if_integral']) for v in triples]}))


if __name__=='__main__': main()
