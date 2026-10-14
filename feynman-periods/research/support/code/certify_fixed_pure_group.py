"""Exact arithmetic witnesses; geometric arguments are in FIXED_PURE_GROUP.md."""
from pathlib import Path
import sympy as s,json,hashlib,math
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/fixed_pure_group';OUT.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    source=ROOT/'data/picard_point_counts.json'
    old=json.loads(source.read_text())
    sample=next(v for v in old['samples'] if v['masses']==[1,2,4,8] and v['s']==-64)
    assert sample['p']==11 and sample['no_residual_root_is_p_times_root_of_unity']
    X=s.Symbol('X');P=s.Poly.from_list(sample['residual_polynomial_coefficients'],X)
    expected=X**6+10*X**5+55*X**4-484*X**3+6655*X**2+146410*X+1771561
    assert P.as_expr()==expected
    mod5=s.Poly(P,X,modulus=5);xx=s.Poly(X,X,modulus=5)
    def power(e):
        a=s.Poly(1,X,modulus=5);b=xx
        while e:
            if e%2:a=(a*b).rem(mod5)
            b=(b*b).rem(mod5);e//=2
        return a
    assert s.gcd(mod5,power(25)-xx).degree()==0
    assert s.gcd(mod5,power(125)-xx).degree()==0
    assert (power(15625)-xx).rem(mod5).is_zero
    mod11=s.Poly(P,X,modulus=11)
    assert mod11==s.Poly(X**5*(X-1),X,modulus=11)
    assert int(P.diff().eval(1))%11==1
    mod29=s.Poly(P,X,modulus=29)
    q4=X**4-13*X**3+X**2-7*X-4
    assert mod29==s.Poly((X-1)*(X-5)*q4,X,modulus=29)
    assert s.gcd(mod29,mod29.diff()).degree()==0
    hensel={'root_mod29':1,'P_at_root_mod29':int(P.eval(1))%29,'derivative_mod29':int(P.diff().eval(1))%29,
      'minus_three_euler_mod29':pow(26,14,29)}
    assert hensel['P_at_root_mod29']==0 and hensel['derivative_mod29']!=0
    assert hensel['minus_three_euler_mod29']==28
    Q=s.Matrix([[-2,-2,-2,4,-2,-2],[-2,-2,0,2,0,-4],[-2,0,-2,2,0,-4],
       [4,2,2,-2,0,8],[-2,0,0,0,-2,-4],[-2,-4,-4,8,-4,-2]])
    assert Q.det()==192
    assert s.Rational(Q.det(),3)==8**2
    sysfile=ROOT/'data/banana_pure_ray/polarized_system.json'
    sys=json.loads(sysfile.read_text());ss=s.Symbol('ss')
    K=s.Matrix(sys['inverse_pairing']).subs(ss,-64)
    dK=K.det();period_square=s.cancel(64*Q.det()/dK)
    assert period_square==s.Rational(9,86198192929486625)
    # Discriminant class is unchanged by rational basis changes, including index two.
    squarefree=s.prod(int(p) for p,e in s.factorint(s.denom(period_square)).items() if e%2)
    squarepart=s.sqrt(s.denom(period_square)/squarefree)
    assert squarepart.is_Integer and squarepart==85
    assert s.numer(period_square)==9 and squarefree!=1
    # Independent four-parameter determinant, with exact rational gauge to bare B.
    # Prefer direct legacy K for the determinant witness; matrix is saved for replay.
    out={'status':'fixed_pure_group_arithmetic_witnesses_certified',
      'base_field':'Q','masses':[1,2,4,8],'s':-64,'input_point_counts_sha256':sha(source),
      'input_polarized_system_sha256':sha(sysfile),
      'input_absolute_comparison_sha256':sha(ROOT/'data/banana_pure_ray/marked_lattice/absolute_periods.json'),
      'code_sha256':sha(__file__),'new_point_counts_computed':0,
      'frobenius_coefficients':list(map(int,P.all_coeffs())),
      'irreducibility_mod5_rabin':True,
      'mod11':'X^5*(X-1)','unique_unit_root_mod11':1,'unit_root_derivative_mod11':1,
      'powers_have_distinct_roots_argument':'irreducible P plus unique 11-adic unit root',
      'auxiliary_prime':29,'mod29_factorization':'(X-1)(X-5)(X^4-13X^3+X^2-7X-4)',
      'hensel_and_inertness':hensel,'Q6':[[int(v) for v in row] for row in Q.tolist()],'det_Q6':192,
      'imaginary_quadratic_candidate_squarefree_d':3,
      'K_B_at_minus64':[[str(v) for v in row] for row in K.tolist()],
      'det_K_B':str(dK),'det_Tate_normalized_comparison_squared':str(period_square),
      'discriminant_squarefree_part':int(squarefree),'discriminant_square_factor':int(squarepart),
      'geometry_proof':'FIXED_PURE_GROUP.md','Hodge_endomorphism_field_conclusion_requires_proof':'Q',
      'Nori_group_conclusion_requires_proof':'O(Q6) x Gm, over Q on T(3) plus Q(1)',
      'independent_replay_completed':False}
    (OUT/'certificate.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k not in ['K_B_at_minus64','Q6']},ensure_ascii=False),flush=True)
if __name__=='__main__':main()

