"""Independent integer/Fraction replay; no SymPy and no generator imports."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product,permutations
from collections import Counter
import hashlib,json
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/marked_saturation'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def mat(a):return [[F(x) for x in r] for r in a]
def trans(a):return [list(c) for c in zip(*a)]
def mul(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def eye(n):return [[F(i==j) for j in range(n)] for i in range(n)]
def det(a):
    n=len(a);ans=F(0)
    for perm in permutations(range(n)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n))
        term=F(sign)
        for i,j in enumerate(perm):term*=a[i][j]
        ans+=term
    return ans
def norm(a):return max(sum(abs(x) for x in row) for row in a)
def quad(a,v):return sum(v[i]*a[i][j]*v[j] for i in range(len(v)) for j in range(len(v)))
def torsion(a):
    out=[]
    for v in product((0,1),repeat=len(a)):
        dual=[sum(row[j]*F(v[j],2) for j in range(len(v))) for row in a]
        if any(x.denominator!=1 for x in dual):continue
        out.append({'numerator':list(v),'q_mod_2':str(quad(a,v)/4%2)})
    return out
def main():
    c=json.loads((OUT/'certificate.json').read_text())
    assert c['code_sha256']==sha(ROOT/'code/certify_marked_saturation.py')
    for name,h in c['inputs_sha256'].items():assert sha(ROOT/name)==h
    Q=mat(c['Q6']);H=mat(c['target_gram'])
    assert Q==mat(json.loads((ROOT/'data/banana_pure_ray/q6_traces_replay.json').read_text())['geometric_Gram_up_to_specified_orientations'])
    assert H==mat(json.loads((ROOT/'data/hessian_lattice.json').read_text())['transcendental_gram'])
    assert det(Q)==192 and det(H)==48
    tt=torsion(H);assert tt==c['target_two_torsion'] and len(tt)==16
    assert all(F(x['q_mod_2']).denominator==1 for x in tt)
    actual=[]
    for v in product((0,1),repeat=6):
        if any(v) and quad(Q,v)%8==0 and all(sum(row[j]*v[j] for j in range(6))%2==0 for row in Q):
            actual.append(list(v))
    assert actual==[x['numerator'] for x in c['even_candidates']] and len(actual)==11
    survivors=[];obstructions=[]
    for row in c['even_candidates']:
        v=row['numerator'];B=eye(6);pivot=v.index(1)
        for i in range(6):B[i][pivot]=F(v[i],2)
        G=mul(trans(B),mul(Q,B))
        assert B==mat(row['basis']) and G==mat(row['gram'])
        assert abs(det(B))==F(1,2) and det(G)==48
        assert all(x.denominator==1 for rr in G for x in rr) and all(G[i][i]%2==0 for i in range(6))
        tt=torsion(G);bad=[x for x in tt if F(x['q_mod_2']).denominator!=1]
        assert tt==row['two_torsion'] and dict(Counter(x['q_mod_2'] for x in tt))==row['histogram']
        assert len(bad)==row['nonintegral_count'] and (bad[0] if bad else None)==row['first_obstruction']
        if not bad:survivors.append(v)
        else:obstructions.append({'candidate':v,'witness':bad[0]})
    assert survivors==[c['selected_half_numerator']]==[[1,1,1,1,0,0]]
    selected=next(x for x in c['even_candidates'] if x['numerator']==survivors[0])
    B=mat(selected['basis']);G=mat(selected['gram']);R=mat(c['standard_isometry'])
    assert G==mat(c['saturated_gram']) and det(R)==1 and mul(trans(R),mul(G,R))==H
    assert c['periods_sha256']==sha(OUT/'periods.json')
    p=json.loads((ROOT/'data/banana_pure_ray/marked_lattice/absolute_periods.json').read_text())
    dst=json.loads((OUT/'periods.json').read_text())
    assert dst['source_sha256']==sha(ROOT/'data/banana_pure_ray/marked_lattice/absolute_periods.json')
    old=[list(map(int,a)) for a in p['physical_B_center_gaussian_integer']]
    checks=0
    for label,W,cap in [('saturated_basis',B,84),('standard_basis',mul(B,R),81)]:
        d=dst[label];assert d['bits']==p['bits']+1 and mat(d['right_basis_matrix'])==W
        raw=[[sum(old[6*i+k][axis]*2*W[k][j] for k in range(6)) for axis in range(2)]
             for i in range(6) for j in range(6)]
        assert all(v.denominator==1 for entry in raw for v in entry)
        assert [[int(v) for v in entry] for entry in raw]==[list(map(int,a)) for a in d['center_gaussian_integer']]
        assert F(d['error_upper'])==norm(W)*F(p['physical_B_error_upper'])<F(1,2**cap)
        assert F(d['right_basis_norm'])==norm(W) and d['rounding_error']=='0'
        checks+=36
    # A wrong half-cycle is rejected by an actual nonintegral discriminant value.
    wrong=next(x for x in c['even_candidates'] if x['numerator']==[1,0,0,0,0,1])
    assert wrong['first_obstruction'] and F(wrong['first_obstruction']['q_mod_2']).denominator==2
    output={'status':'marked_saturation_independently_replayed','certificate_sha256':sha(OUT/'certificate.json'),
      'periods_sha256':sha(OUT/'periods.json'),'verifier_sha256':sha(__file__),
      'arithmetic':'integers and Fraction; independent bilinear sums and Leibniz determinants',
      'even_classes_exhausted':11,'excluded_by_discriminant_quadratic_form':10,
      'selected_half_numerator':survivors[0],'explicit_exclusion_witnesses':obstructions,
      'standard_isometry_replayed':True,'transformed_period_entries_replayed':checks,
      'saturated_period_error_less_than':'2^-84','standard_period_error_less_than':'2^-81',
      'wrong_half_cycle_negative_control_rejected':True,'new_ODE_transports':0,
      'geometric_input_scope':'Uses proved actual inclusion, evenness, index two and full abstract T; no new claim about mass loops.'}
    (OUT/'replay.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k!='explicit_exclusion_witnesses'}))
if __name__=='__main__':main()

