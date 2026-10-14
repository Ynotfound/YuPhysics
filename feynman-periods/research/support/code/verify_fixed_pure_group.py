"""Independent integer/Fraction replay, without SymPy or the certificate generator."""
from pathlib import Path
from fractions import Fraction as F
from itertools import permutations
import json,hashlib,ast,math
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/fixed_pure_group'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0:a.pop()
    return a
def add(a,b,p):return trim([((a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0))%p for i in range(max(len(a),len(b)))])
def mul(a,b,p):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%p
    return trim(c)
def rem(a,b,p):
    a=trim([v%p for v in a]);b=trim([v%p for v in b])
    while a!=[0] and len(a)>=len(b):
        c=a[-1]*pow(b[-1],-1,p)%p;k=len(a)-len(b)
        for i,x in enumerate(b):a[i+k]=(a[i+k]-c*x)%p
        a=trim(a)
    return a
def gcd(a,b,p):
    while b!=[0]:a,b=b,rem(a,b,p)
    return [(x*pow(a[-1],-1,p))%p for x in a]
def power(n,f,p):
    a=[1];b=[0,1]
    while n:
        if n&1:a=rem(mul(a,b,p),f,p)
        b=rem(mul(b,b,p),f,p);n//=2
    return a
def evaluate(a,x,p):
    v=0
    for c in reversed(a):v=(v*x+c)%p
    return v
def rat_expr(v):
    if not isinstance(v,str):return F(v)
    def walk(n):
        if isinstance(n,ast.Constant) and type(n.value)==int:return F(n.value)
        if isinstance(n,ast.Name) and n.id=='ss':return F(-64)
        if isinstance(n,ast.UnaryOp):
            a=walk(n.operand)
            if isinstance(n.op,ast.USub):return -a
            if isinstance(n.op,ast.UAdd):return a
        if isinstance(n,ast.BinOp):
            a,b=walk(n.left),walk(n.right)
            if isinstance(n.op,ast.Add):return a+b
            if isinstance(n.op,ast.Sub):return a-b
            if isinstance(n.op,ast.Mult):return a*b
            if isinstance(n.op,ast.Div):return a/b
            if isinstance(n.op,ast.Pow) and b.denominator==1:return a**int(b)
        raise ValueError(ast.dump(n))
    return walk(ast.parse(v,mode='eval').body)
def det(A):
    ans=F(0);n=len(A)
    for pm in permutations(range(n)):
        term=F((-1)**sum(pm[i]>pm[j] for i in range(n) for j in range(i+1,n)))
        for i in range(n):term*=A[i][pm[i]]
        ans+=term
    return ans
def main():
    d=json.loads((OUT/'certificate.json').read_text())
    for key,name in [('input_point_counts_sha256','data/picard_point_counts.json'),
      ('input_polarized_system_sha256','data/banana_pure_ray/polarized_system.json'),
      ('input_absolute_comparison_sha256','data/banana_pure_ray/marked_lattice/absolute_periods.json')]:
        assert d[key]==sha(ROOT/name)
    sample=next(v for v in json.loads((ROOT/'data/picard_point_counts.json').read_text())['samples'] if v['s']==-64)
    tr=[r['smooth_count']-1-r['q']**2-16*r['q'] for r in sample['extensions']]
    cs=[F(1)]
    for k in range(1,4):cs.append(-sum(cs[k-i]*tr[i-1] for i in range(1,k+1))/k)
    coeff=cs+[121*cs[2],11**4*cs[1],F(11**6)]
    assert all(x.denominator==1 for x in coeff) and list(map(int,coeff))==d['frobenius_coefficients']
    p=list(reversed(list(map(int,coeff))));p5=[x%5 for x in p]
    assert gcd(p5,add(power(25,p5,5),[0,-1],5),5)==[1]
    assert gcd(p5,add(power(125,p5,5),[0,-1],5),5)==[1]
    assert rem(add(power(15625,p5,5),[0,-1],5),p5,5)==[0]
    assert [x%11 for x in p]==[0,0,0,0,0,10,1]
    derivative=[i*p[i] for i in range(1,len(p))]
    assert evaluate(derivative,1,11)==1
    factor=mul(mul([-1,1],[-5,1],29),[-4,-7,1,-13,1],29)
    assert factor==[x%29 for x in p]
    assert evaluate(p,1,29)==0 and evaluate(derivative,1,29)==1
    assert 26 not in {x*x%29 for x in range(29)}
    assert gcd([x%29 for x in p],[x%29 for x in derivative],29)==[1]
    Q=[[F(v) for v in r] for r in d['Q6']]
    assert det(Q)==192 and det(Q)/3==64
    sys=json.loads((ROOT/'data/banana_pure_ray/polarized_system.json').read_text())
    K=[[rat_expr(v) for v in row] for row in sys['inverse_pairing']]
    assert K==[[F(v) for v in row] for row in d['K_B_at_minus64']]
    dk=det(K);w=64*det(Q)/dk
    assert dk==F(d['det_K_B']) and w==F(d['det_Tate_normalized_comparison_squared'])
    factors=[5,13,29,37,73,89,113,233]
    for q in factors:assert all(q%k for k in range(2,math.isqrt(q)+1))
    rad=math.prod(factors)
    assert rad==d['discriminant_squarefree_part']==11930545734185
    assert w==F(9,85**2*rad) and math.isqrt(w.denominator)**2!=w.denominator
    bad=p[:];bad[0]+=1
    assert evaluate(bad,1,29)!=0
    sources=json.loads((ROOT/'external/fixed-pure-group/manifest.json').read_text())
    for r in sources:assert sha(ROOT/'external/fixed-pure-group'/r['name'])==r['sha256']
    out={'status':'fixed_pure_arithmetic_independently_replayed',
      'certificate_sha256':sha(OUT/'certificate.json'),'verifier_sha256':sha(__file__),
      'arithmetic':'Python integers, Fraction, independent AST rational evaluator and Leibniz determinant',
      'irreducibility_mod5':True,'unique_11_adic_unit_root':True,
      'simple_29_adic_eigenline_and_inert_quadratic':True,
      'det_comparison_square':str(w),'nonrational_determinant_period_proved_given_absolute_comparison':True,
      'new_point_counts_computed':0,'old_transports_repeated':0,'negative_control_rejected':True,
      'external_source_hashes_checked':True,
      'theorem_conditions_and_group_identification':'FIXED_PURE_GROUP.md'}
    (OUT/'replay.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out),flush=True)
if __name__=='__main__':main()

