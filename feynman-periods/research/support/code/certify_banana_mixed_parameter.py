"""Coefficientwise certification of the reconstructed actual rank-ten system."""
from pathlib import Path
import sys,json,hashlib,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'runtime/python-flint'))
from flint import fmpq,fmpq_mpoly_ctx
from itertools import combinations
CTX=fmpq_mpoly_ctx.get(('m12','m22','m32','m42'));Z=CTX.constant(0);I=CTX.constant(1)
OUT=ROOT/'data/banana_mixed_parameter'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def poly(v):return CTX.from_dict({tuple(e[:4]):fmpq(e[4]) for e in v})
def po(p):return [[*map(int,e),str(c)] for e,c in sorted(p.to_dict().items())]
def save(n,d):
    p=OUT/n;q=p.with_suffix('.tmp');q.write_text(json.dumps(d,separators=(',',':'))+'\n');q.replace(p)
def theta(p,k):return CTX.gen(k)*p.derivative(k)
def main():
    start=time.time();D=I;seen=set();raw=[]
    recon=json.loads((OUT/'reconstruction.json').read_text())
    for k in range(4):
        p=OUT/f'Alog{k+1}.json';assert sha(p)==recon['matrix_sha256'][p.name]
        M=[[(poly(v['n']),poly(v['d'])) for v in row] for row in json.loads(p.read_text())['matrix']]
        raw.append(M)
        for row in M:
            for n,d in row:
                key=str(d)
                if key not in seen:
                    seen.add(key);D=D*(d/D.gcd(d));D=D/D.coeffs()[0]
    N=[[[n*(D/d) for n,d in row] for row in A] for A in raw]
    save('cleared_system.json',{'D':po(D),'N':[[[po(p) for p in r] for r in A] for A in N],
       'matrix_sha256':recon['matrix_sha256'],'D_degree':int(D.total_degree())})
    done=[]
    for a,b in combinations(range(4),2):
        count=0
        for i in range(10):
            for j in range(10):
                r=D*(theta(N[b][i][j],a)-theta(N[a][i][j],b))-N[b][i][j]*theta(D,a)+N[a][i][j]*theta(D,b)
                r+=sum((N[b][i][k]*N[a][k][j]-N[a][i][k]*N[b][k][j] for k in range(10)),Z)
                assert r==0,('nonzero_rank10_curvature',a,b,i,j,str(r))
                count+=1
        done.append({'directions':[a+1,b+1],'zero_entries':count})
        save(f'curvature_{a+1}{b+1}.json',done[-1])
        print('rank10 curvature',a+1,b+1,'passed',round(time.time()-start,2),flush=True)
    out={'status':'full_four_parameter_rank_ten_mixed_flatness_certified',
      'basis':'L6,beta1,beta2,beta3,e','directions':'theta_i=x_i partial_x_i',
      'curvature_entries':600,'pairs':done,'reconstruction_sha256':sha(OUT/'reconstruction.json'),
      'matrix_sha256':recon['matrix_sha256'],'cleared_system_sha256':sha(OUT/'cleared_system.json'),
      'certifier_sha256':sha(__file__),'independent_replay_completed':False,
      'absolute_rank10_comparison_certified':False}
    save('flatness.json',out);print(json.dumps(out),flush=True)
if __name__=='__main__':main()
