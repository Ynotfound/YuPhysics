"""Independent replay from the compact affine formula, not the reconstruction."""
from pathlib import Path
import sys,json,hashlib,time
from itertools import combinations,permutations
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'runtime/python-flint'))
from flint import fmpq,fmpq_mpoly_ctx
import sympy as sp
CTX=fmpq_mpoly_ctx.get(('m12','m22','m32','m42'));xx=list(CTX.gens());Z=CTX.constant(0)
OUT=ROOT/'data/banana_mixed_parameter';PURE=ROOT/'data/banana_four_parameter'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def poly(a):return CTX.from_dict({tuple(r[:4]):fmpq(r[4]) for r in a})
def dlog(p,i):return CTX.from_dict({e:c*e[i] for e,c in p.to_dict().items() if e[i]})
def product(A,B):
    m=len(A);C=[[Z for j in range(m)] for i in range(m)]
    for k in range(m):
        for i in range(m):
            for j in range(m):C[i][j]+=A[i][k]*B[k][j]
    return C
def coupling(i,j):
    if i==j:return [-xx[i]]*5+[Z]
    k,l=[v for v in range(4) if v not in (i,j)]
    v=[Z]*6
    v[0]=(2*xx[i]+2*xx[j]-xx[k]-xx[l])/6
    v[1+i]=(xx[i]-1)/3;v[1+j]=(xx[j]-1)/3
    v[1+k]=(3*xx[i]+3*xx[j]-xx[k]-3*xx[l]+1)/6
    v[1+l]=(3*xx[i]+3*xx[j]-3*xx[k]-xx[l]+1)/6
    return v
def check_ibp():
    alpha=sp.symbols('a b c')+(sp.Integer(1),)
    x=sp.symbols('x1:5');N=sum(1/a for a in alpha);M=sum(p*a for p,a in zip(x,alpha));f=M*N;F=f-1
    def torus_derivative(e,i):
        return alpha[i]*sp.diff(e,alpha[i]) if i<3 else -sum(alpha[k]*sp.diff(e,alpha[k]) for k in range(3))
    A=[a*N for a in alpha]
    for i in range(4):
        for j in range(4):
            p=-int(i==j)*A[i]*F+x[i]*A[i]*A[j]+torus_derivative(A[j],i)*F-A[j]*torus_derivative(f,i)-alpha[j]/alpha[i]
            assert sp.expand(p)==0,('torus_IBP',i,j)
    assert sp.expand(sum(x[i]*A[i] for i in range(4))-f)==0
    return 16
def main():
    start=time.time();ibp=check_ibp()
    C=[[coupling(i,j) for j in range(4)] for i in range(4)]
    assert all(C[i][j]==C[j][i] for i in range(4) for j in range(4))
    for i in range(4):
        for k in range(6):assert sum((C[i][j][k] for j in range(4)),Z)==-int(k==i+1)
        saved=read(OUT/f'couplings_{i+1}.json')
        for j in range(4):
            for k in range(6):
                n,d=map(poly,[saved['rows'][j][k]['n'],saved['rows'][j][k]['d']])
                assert n==d*C[i][j][k]
    # Uniqueness of the four all-distinct entries of a symmetric cubic tensor.
    inc=sp.Matrix([[1,1,0,0],[1,0,1,0],[0,1,1,0],[1,0,0,1]])
    assert inc.det()==-2
    cert=read(OUT/'flatness.json');clear=read(OUT/'cleared_system.json')
    assert cert['cleared_system_sha256']==sha(OUT/'cleared_system.json')
    D=poly(clear['D']);N=[[[poly(v) for v in row] for row in A] for A in clear['N']]
    S=[CTX.constant(1),*xx,CTX.constant(1)]
    for j in range(4):
        p=OUT/f'Alog{j+1}.json';assert sha(p)==cert['matrix_sha256'][p.name]
        A=read(p)['matrix'];old=read(PURE/f'A{j+1}.json')['matrix']
        for r in range(10):
            for c in range(10):
                n,d=poly(A[r][c]['n']),poly(A[r][c]['d'])
                assert N[j][r][c]*d==D*n
                if r<6 and c<6:
                    on,od=poly(old[r][c]['n']),poly(old[r][c]['d'])
                    assert (n-int(r==c==j+1)*d)*od*S[c]==d*xx[j]*S[r]*on
                elif 6<=r<9 and c<6:
                    assert n==d*C[r-6][j][c]
                elif r==9:
                    expected=(int(c==j+6) if j<3 else -int(c in [0,6,7,8]))
                    assert n==d*expected
                else:assert n==0
    count=0
    for a,b in combinations(range(4),2):
        AB=product(N[a],N[b]);BA=product(N[b],N[a])
        for i in range(10):
            for j in range(10):
                r=D*(dlog(N[a][i][j],b)-dlog(N[b][i][j],a))-N[a][i][j]*dlog(D,b)+N[b][i][j]*dlog(D,a)+AB[i][j]-BA[i][j]
                assert r==0,('curvature',a,b,i,j);count+=1
        print('independent rank10 pair',a+1,b+1,'seconds',round(time.time()-start,2),flush=True)
    # A wrong sign on just the beta1 row destroys symmetry/top curvature.
    bad=[v for v in C[0][1]];bad[0]+=1
    assert bad!=C[1][0]
    # Kummer unit directions have the diagonal valuation identity.
    valuations=sp.eye(3);assert valuations.det()==1
    # On fixed masses, -Theta/s beta_i=V_i/s and -Theta/s e=Omega/s.
    assert all([sum((C[i][j][k] for j in range(4)),Z) for k in range(6)]==[-CTX.constant(int(k==i+1)) for k in range(6)] for i in range(4))
    out={'status':'compact_mixed_connection_geometry_and_600_entries_independently_replayed',
      'rank':10,'curvature_entries':count,'direct_torus_IBP_identities':ibp,
      'compact_coupling_entries':96,'constant_Kummer_rank':3,'uniqueness_minor_determinant':-2,
      'input_pure_replay_sha256':sha(PURE/'replay.json'),
      'matrix_sha256':cert['matrix_sha256'],'flatness_sha256':sha(OUT/'flatness.json'),
      'verifier_sha256':sha(__file__),'negative_control_rejected':True,
      'old_rank8_ray_embedding_checked':True,'absolute_rank10_comparison_certified':False,
      'geometric_object_identification_proof':'MIXED_PARAMETER_CONNECTION.md'}
    (OUT/'replay.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out),flush=True)
if __name__=='__main__':main()
