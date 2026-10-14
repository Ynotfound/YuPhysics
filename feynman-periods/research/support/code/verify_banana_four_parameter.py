"""Independent parser and coefficient replay of the four-parameter certificate.

This file does not import the extractor or the certifier.  Published expressions
are parsed by SymPy's Mathematica parser; curvature derivatives are reconstructed
coefficientwise, and matrix products use column-row outer sums.
"""
from pathlib import Path
import sys,json,re,hashlib,time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'runtime/python-flint'))
from flint import fmpq,fmpq_mpoly_ctx,fmpq_poly
import sympy as sp
from sympy.parsing.mathematica import parse_mathematica
OUT=ROOT/'data/banana_four_parameter';SRC=ROOT/'external/bonn-banana-2025'
ctx=fmpq_mpoly_ctx.get(('m12','m22','m32','m42'));zero=ctx.constant(0)
xx=sp.symbols('m12 m22 m32 m42');uu=sp.symbols('u1:5',positive=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def poly(data):return ctx.from_dict({tuple(v[:4]):fmpq(v[4]) for v in data})
def split(s):
    s=s.strip();assert s[0]=='{' and s[-1]=='}'
    stack=[];out=[];start=1
    for k,c in enumerate(s[1:-1],1):
        if c in '([{':stack.append(c)
        elif c in ')]}':stack.pop()
        elif c==',' and not stack:out.append(s[start:k]);start=k+1
    assert not stack;out.append(s[start:-1]);return out
def rows(path):
    t=re.sub(r'\(\*.*?\*\)','',path.read_text(encoding='utf-8'),flags=re.S)
    return [split(s) for s in split(t)]
def from_expr(e,cover=False):
    if cover:
        e=e.subs(dict(zip(xx,[u*u for u in uu])),simultaneous=True)
    n,d=e.as_numer_denom()
    def convert(p):
        items=sp.Poly(p,*(uu if cover else xx),domain=sp.QQ).as_dict()
        if cover:
            assert all(all(k%2==0 for k in m) for m in items)
            items={tuple(k//2 for k in m):c for m,c in items.items()}
        return ctx.from_dict({m:fmpq(str(c)) for m,c in items.items()})
    return convert(n),convert(d)
def dp(p,k):
    out={}
    for e,c in p.to_dict().items():
        e=tuple(map(int,e))
        if e[k]:
            f=list(e);f[k]-=1;out[tuple(f)]=c*e[k]
    return ctx.from_dict(out)
def prod(A,B):
    C=[[zero for j in range(6)] for i in range(6)]
    for k in range(6):
        for i in range(6):
            for j in range(6):C[i][j]+=A[i][k]*B[k][j]
    return C
def transpose(A):return list(map(list,zip(*A)))
def raypoly(p):
    deg=max((sum(e) for e in p.to_dict()),default=0)
    coeff=[fmpq(0)]*(int(deg)+1)
    for e,c in p.to_dict().items():
        z=c
        for mu,k in zip((1,4,16,64),e):z*=mu**int(k)
        coeff[int(deg-sum(e))]+=z
    return fmpq_poly(coeff),int(deg)
def oldpoly(e,s):
    p=sp.Poly(e,s,domain=sp.QQ)
    return fmpq_poly([fmpq(str(c)) for c in reversed(p.all_coeffs())])
def ray_equal(n,d,text,shift=0):
    p,a=raypoly(n);q,b=raypoly(d);shift+=b-a
    ss=fmpq_poly([0,1])
    if shift>=0:p*=ss**shift
    else:q*=ss**(-shift)
    s=sp.Symbol('ss');v=sp.sympify(text);u,w=v.as_numer_denom()
    return p*oldpoly(w,s)==q*oldpoly(u,s)
def main():
    start=time.time();cert=read(OUT/'flatness.json');clear=read(OUT/'cleared_system.json')
    assert cert['cleared_system_sha256']==sha(OUT/'cleared_system.json')
    D=poly(clear['D']);E=poly(clear['E']);assert D and E
    N=[[[poly(v) for v in row] for row in A] for A in clear['N']]
    M=[[poly(v) for v in row] for row in clear['M']]
    pure=(4,5,6,7,8,14);scale=(1,1,1,1,1,2);source_hashes={}
    ray_checks=0;source_checks=0
    for k in range(4):
        path=SRC/f'Diffeqm{k+1}.m';source_hashes[path.name]=sha(path);source=rows(path)
        saved=read(OUT/f'A{k+1}.json')
        assert saved['source_sha256']==sha(path)
        assert cert['input_sha256'][f'A{k+1}.json']==sha(OUT/f'A{k+1}.json')
        old=read(ROOT/f'data/banana_pure_ray/A_x{k+1}.json')
        for i in range(6):
            for j in range(6):
                e=parse_mathematica(source[pure[i]][pure[j]].replace(r'\[Epsilon]','0'))
                n,d=from_expr(e*sp.Rational(scale[i],scale[j]))
                sn,sd=poly(saved['matrix'][i][j]['n']),poly(saved['matrix'][i][j]['d'])
                assert n*sd==sn*d and N[k][i][j]*sd==D*sn,(k,i,j,'source')
                assert ray_equal(sn,sd,old['matrix'][i][j]),(k,i,j,'ray')
                source_checks+=1;ray_checks+=1
            for j in range(9,14):
                e=parse_mathematica(source[pure[i]][j].replace(r'\[Epsilon]','0'))
                assert e==0,(k,i,j,'ISP')
            print('independent source',k+1,'row',i,round(time.time()-start,1),flush=True)
    path=SRC/'ZCYInverse.m';source_hashes[path.name]=sha(path);source=rows(path)
    saved=read(OUT/'K.json');old=read(ROOT/'data/banana_pure_ray/pairing_and_roots.json')
    assert saved['source_sha256']==sha(path) and cert['input_sha256']['K.json']==sha(OUT/'K.json')
    for i in range(6):
        for j in range(6):
            e=parse_mathematica(source[i][j])/scale[i]/scale[j]
            n,d=from_expr(e,'Sqrt' in source[i][j])
            sn,sd=poly(saved['matrix'][i][j]['n']),poly(saved['matrix'][i][j]['d'])
            assert n*sd==sn*d and M[i][j]*sd==E*sn,(i,j,'source_pairing')
            assert ray_equal(sn,sd,old['inverse_pairing'][i][j],2),(i,j,'ray_pairing')
            source_checks+=1;ray_checks+=1
    assert M==transpose(M)
    count=0
    for a in range(4):
        for b in range(a+1,4):
            AB=prod(N[a],N[b]);BA=prod(N[b],N[a])
            for i in range(6):
                for j in range(6):
                    r=D*(dp(N[a][i][j],b)-dp(N[b][i][j],a))-N[a][i][j]*dp(D,b)+N[b][i][j]*dp(D,a)+AB[i][j]-BA[i][j]
                    assert r==0,('flatness',a,b,i,j);count+=1
            print('independent curvature',a+1,b+1,round(time.time()-start,1),flush=True)
    pcount=0
    for a in range(4):
        NM=prod(transpose(N[a]),M);MN=prod(M,N[a])
        for i in range(6):
            for j in range(6):
                r=D*E*dp(M[i][j],a)-D*dp(E,a)*M[i][j]+E*NM[i][j]+E*MN[i][j]
                assert r==0,('pairing',a,i,j);pcount+=1
    assert count==216 and pcount==144 and source_checks==ray_checks==180
    # Independently require the changed A_1(0,0)+x_2 to break the identity.
    changed=[r.copy() for r in N[0]];changed[0][0]+=D*ctx.gen(1)
    r=D*(dp(changed[0][0],1)-dp(N[1][0][0],0))-changed[0][0]*dp(D,1)+N[1][0][0]*dp(D,0)+prod(changed,N[1])[0][0]-prod(N[1],changed)[0][0]
    assert r!=0
    out={'status':'four_parameter_source_and_all_identities_independently_replayed',
        'source_entries':source_checks,'ISP_zero_entries':120,'old_ray_equalities':ray_checks,
        'flatness_entries':count,'pairing_entries':pcount,
        'source_parser':'SymPy parse_mathematica; exact root-cover descent for K05',
        'derivatives':'independent exponentwise rule','products':'column-row outer accumulation',
        'negative_control_rejected':True,'source_sha256':source_hashes,
        'certificate_sha256':sha(OUT/'flatness.json'),'cleared_system_sha256':sha(OUT/'cleared_system.json'),
        'verifier_sha256':sha(__file__),'fixed_pure_group_identified':False,
        'mixed_four_parameter_flatness_certified':False}
    (OUT/'replay.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out),flush=True)
if __name__=='__main__':main()
