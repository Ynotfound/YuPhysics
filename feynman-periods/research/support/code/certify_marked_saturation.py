"""Select actual index-two saturation by the two-torsion discriminant quadratic form."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import Counter
import json,hashlib,sympy as sp
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/marked_saturation'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def rows(M):return [[int(x) if x.q==1 else str(x) for x in row] for row in M.tolist()]
def norm(M):return max(sum(abs(x) for x in M.row(i)) for i in range(M.rows))
def q2(G):
    data=[]
    for w in product(range(2),repeat=G.rows):
        W=sp.Matrix(w)
        if any(x%2 for x in G*W):continue
        q=(W.T*G*W)[0]/4
        data.append({'numerator':list(w),'q_mod_2':str(q%2)})
    return data
def periods(P,W):
    b=P['bits'];scale=2
    old=[tuple(map(int,c)) for c in P['physical_B_center_gaussian_integer']]
    new=[]
    for i in range(6):
        for j in range(6):
            c=[sum(old[6*i+k][a]*(2*W[k,j]) for k in range(6)) for a in range(2)]
            assert all(x.q==1 for x in c)
            new.append(list(map(lambda x:str(int(x)),c)))
    return {'bits':b+1,'center_gaussian_integer':new,
        'error_upper':str(F(P['physical_B_error_upper'])*F(str(norm(W)))),
        'right_basis_matrix':rows(W),'right_basis_norm':str(norm(W)),
        'rounding_error':'0','basepoint':-64,'de_rham_basis':'bare-residue rational B frame'}
def main():
    OUT.mkdir(exist_ok=True)
    paths=['data/banana_pure_ray/q6_traces_replay.json','data/hessian_lattice.json',
      'LATTICE_SATURATION.md','Q6_GEOMETRIC_CERTIFICATE.md','data/recovered_lattice.json',
      'data/banana_pure_ray/marked_lattice/absolute_periods.json','data/fixed_pure_group/final_audit.json']
    inputs={p:sha(ROOT/p) for p in paths}
    qdata=read(paths[0]);Q=sp.Matrix(qdata['geometric_Gram_up_to_specified_orientations'])
    H=sp.Matrix(read(paths[1])['transcendental_gram'])
    assert Q.det()==192 and H.det()==48
    target=q2(H);assert len(target)==16 and all(sp.Rational(t['q_mod_2']).q==1 for t in target)
    candidates=[]
    for v in product(range(2),repeat=6):
        if not any(v):continue
        V=sp.Matrix(v)
        if any(x%2 for x in Q*V) or (V.T*Q*V)[0]%8:continue
        B=sp.eye(6);B[:,v.index(1)]=V/2;G=B.T*Q*B
        assert all(x.q==1 for x in G) and all(G[i,i]%2==0 for i in range(6))
        assert G.det()==48
        tors=q2(G);bad=[t for t in tors if sp.Rational(t['q_mod_2']).q!=1]
        candidates.append({'numerator':list(v),'basis':rows(B),'gram':rows(G),
          'two_torsion':tors,'histogram':dict(Counter(t['q_mod_2'] for t in tors)),
          'nonintegral_count':len(bad),'first_obstruction':bad[0] if bad else None})
    good=[c for c in candidates if c['nonintegral_count']==0]
    assert len(candidates)==11 and len(good)==1 and good[0]['numerator']==[1,1,1,1,0,0]
    B=sp.Matrix(good[0]['basis']);G=sp.Matrix(good[0]['gram'])
    R=sp.Matrix(read('data/recovered_lattice.json')['standard_isometry'])
    assert R.det()==1 and R.T*G*R==H
    P=read('data/banana_pure_ray/marked_lattice/absolute_periods.json')
    transformed={'status':'absolute_pure_periods_in_complete_integral_betti_bases',
      'source_sha256':sha(ROOT/paths[-2]),'saturated_basis':periods(P,B),
      'standard_basis':periods(P,B*R),'source_error_unchanged':P['physical_B_error_upper']}
    # Fix source lookup explicitly: last-but-one path above is the absolute comparison.
    assert transformed['source_sha256']==sha(ROOT/'data/banana_pure_ray/marked_lattice/absolute_periods.json')
    (OUT/'periods.json').write_text(json.dumps(transformed,indent=2)+'\n')
    cert={'status':'unique_index_two_saturation_from_geometric_discriminant_form',
      'code_sha256':sha(__file__),'inputs_sha256':inputs,'Q6':rows(Q),'target_gram':rows(H),
      'target_two_torsion':target,'index':2,'even_candidates':candidates,
      'selected_half_numerator':good[0]['numerator'],'saturated_gram':rows(G),
      'standard_isometry':rows(R),'standard_isometry_determinant':1,
      'periods_sha256':sha(OUT/'periods.json'),
      'geometric_inputs':'Actual Q6 inclusion and full integral T from LATTICE_SATURATION; not determinant alone.',
      'mass_loop_matrices_certified':False,'new_ODE_transports':0}
    (OUT/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    for key in ['saturated_basis','standard_basis']:
        print(key,'error',float(F(transformed[key]['error_upper'])))
    print('selected',good[0]['numerator'],'certificate',sha(OUT/'certificate.json'),flush=True)
if __name__=='__main__':main()

