"""Deterministic coefficientwise certification over QQ[x1,x2,x3,x4]."""
from pathlib import Path
import sys,json,hashlib,time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'runtime/python-flint'))
import flint
from flint import fmpq,fmpq_mpoly_ctx
OUT=ROOT/'data/banana_four_parameter'
CTX=fmpq_mpoly_ctx.get(('m12','m22','m32','m42'))
ZERO=CTX.constant(0); ONE=CTX.constant(1)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name,obj):
    p=OUT/name; q=p.with_suffix('.tmp')
    q.write_text(json.dumps(obj,ensure_ascii=False,separators=(',',':'))+'\n',encoding='utf-8');q.replace(p)
def poly(rows):return CTX.from_dict({tuple(r[:4]):fmpq(r[4]) for r in rows})
def enc(p):return [[*map(int,e),str(c)] for e,c in sorted(p.to_dict().items())]
def read(name):
    return [[(poly(v['n']),poly(v['d'])) for v in row]
            for row in json.loads((OUT/name).read_text())['matrix']]
def denominator(matrices):
    D=ONE;seen=set()
    for M in matrices:
        for row in M:
            for n,d in row:
                key=str(d)
                if key not in seen:
                    seen.add(key);D=D*(d/D.gcd(d));D=D/D.coeffs()[0]
    assert D!=0
    return D
def clear(M,D):return [[n*(D/d) for n,d in row] for row in M]
def dotrow(A,B,i,j):return sum((A[i][k]*B[k][j] for k in range(6)),ZERO)
def main():
    start=time.time();names=[f'A{k+1}.json' for k in range(4)]+['K.json']
    hashes={n:sha(OUT/n) for n in names}
    raw=[read(n) for n in names]
    D=denominator(raw[:4]);E=denominator(raw[4:])
    N=[clear(M,D) for M in raw[:4]];M=clear(raw[4],E)
    save('cleared_system.json',{'input_sha256':hashes,
        'D':enc(D),'E':enc(E),'N':[[[enc(v) for v in row] for row in mat] for mat in N],
        'M':[[enc(v) for v in row] for row in M],
        'D_degree':int(D.total_degree()),'E_degree':int(E.total_degree()),
        'D_terms':len(D),'E_terms':len(E)})
    print('cleared denominators',D.total_degree(),E.total_degree(),'seconds',time.time()-start,flush=True)
    flat=[]
    for a in range(4):
        for b in range(a+1,4):
            done=[];terms=[]
            for i in range(6):
                for j in range(6):
                    pieces=[D*(N[b][i][j].derivative(a)-N[a][i][j].derivative(b)),
                        -N[b][i][j]*D.derivative(a), N[a][i][j]*D.derivative(b),
                        dotrow(N[b],N[a],i,j),-dotrow(N[a],N[b],i,j)]
                    residual=sum(pieces,ZERO)
                    assert residual==0,('curvature',a,b,i,j,str(residual))
                    done.append([i,j]);terms.append(max(len(p) for p in pieces))
                print('curvature',a+1,b+1,'row',i,'seconds',round(time.time()-start,1),flush=True)
            item={'directions':[a+1,b+1],'zero_entries':done,'largest_summand_terms':max(terms)}
            flat.append(item);save(f'curvature_{a+1}{b+1}.json',item)
    paired=[]
    for a in range(4):
        done=[]
        NT=list(map(list,zip(*N[a])))
        for i in range(6):
            for j in range(6):
                residual=D*(E*M[i][j].derivative(a)-M[i][j]*E.derivative(a))+E*(dotrow(NT,M,i,j)+dotrow(M,N[a],i,j))
                assert residual==0,('pairing',a,i,j,str(residual))
                done.append([i,j])
            print('pairing',a+1,'row',i,'seconds',round(time.time()-start,1),flush=True)
        item={'direction':a+1,'zero_entries':done};paired.append(item);save(f'pairing_{a+1}.json',item)
    # This altered connection fails the exact curvature check.
    altered=[r.copy() for r in N[0]];altered[0][0]+=D*CTX.gen(1)
    a,b=0,1;i,j=0,0
    bad=D*(N[b][i][j].derivative(a)-altered[i][j].derivative(b))-N[b][i][j]*D.derivative(a)+altered[i][j]*D.derivative(b)+dotrow(N[b],altered,i,j)-dotrow(altered,N[b],i,j)
    assert bad!=0
    out={'status':'complete_four_parameter_pure_flatness_and_pairing_certified',
       'variables':['m12','m22','m32','m42'],'method':'exact rational polynomial coefficients; no sampling',
       'curvature_convention':'d_a A_b-d_b A_a+A_b A_a-A_a A_b',
       'flatness_entries':216,'pairing_entries':144,'flatness':flat,'pairing':paired,
       'source_extraction_sha256':sha(OUT/'extraction.json'),'input_sha256':hashes,
       'cleared_system_sha256':sha(OUT/'cleared_system.json'),'code_sha256':sha(__file__),
       'python_flint_version':flint.__version__,'wheel_sha256':sha(ROOT/'runtime/wheels/python_flint-0.8.0-cp311-cp311-win_amd64.whl'),
       'negative_control_rejected':True,'mixed_system_certified':False,
       'independent_replay_completed':False}
    save('flatness.json',out)
    print(json.dumps({k:out[k] for k in ['status','flatness_entries','pairing_entries','negative_control_rejected']}),flush=True)
if __name__=='__main__':main()
