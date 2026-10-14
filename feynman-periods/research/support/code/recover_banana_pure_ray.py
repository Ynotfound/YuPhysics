"""Extract the epsilon=0 pure block and restrict published matrices to m=1248.

Only parses mathematical expressions; never executes the downloaded code.
No claim of rigorous analytic transport or of an integral cycle marking.
"""
from pathlib import Path
import re,json,hashlib,time
import sympy as sp
from sympy.parsing.mathematica import parse_mathematica

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'external/bonn-banana-2025'
DEST=ROOT/'data/banana_pure_ray'
DEST.mkdir(exist_ok=True)
s=sp.Symbol('ss')
pure=[4,5,6,7,8,14]
isp=list(range(9,14))
mass_squares=[1,4,16,64]
manifest=json.loads((SOURCE/'manifest.json').read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def split_list(text):
    text=text.strip()
    assert text[0]=='{' and text[-1]=='}'
    out=[];start=1;depth=0
    for i,c in enumerate(text[1:-1],1):
        if c in '{[(':depth+=1
        elif c in '}])':depth-=1
        elif c==',' and depth==0:
            out.append(text[start:i].strip());start=i+1
        assert depth>=0
    assert depth==0
    out.append(text[start:-1].strip())
    return out

def restrict_parse(text):
    text=text.replace(r'\[Epsilon]','0')
    for i,mu in enumerate(mass_squares,1):
        text=re.sub(r'\bm'+str(i)+r'2\b','('+str(mu)+'/ss)',text)
    val=parse_mathematica(text)
    assert val.free_symbols<= {s}
    return sp.cancel(val)

matrices=[];input_hashes=[]
T=sp.diag(1,1,1,1,1,2)
for k in range(1,5):
    path=SOURCE/f'Diffeqm{k}.m'
    assert sha(path)==next(f['sha256'] for f in manifest['files'] if f['name']==path.name)
    input_hashes.append(sha(path))
    cache=DEST/f'A_x{k}.json'
    if cache.exists():
        data=json.loads(cache.read_text(encoding='utf-8'))
        assert data['source_sha256']==sha(path)
        A=sp.Matrix([[sp.sympify(v) for v in row] for row in data['matrix']])
        print(f'reused_verified_A_x{k}',flush=True)
    else:
        started=time.time()
        txt=re.sub(r'\(\*.*?\*\)','',path.read_text(encoding='utf-8'),flags=re.S)
        rows=[split_list(row) for row in split_list(txt)]
        assert len(rows)==15 and all(len(row)==15 for row in rows)
        for i in pure:
            for j in isp:
                assert restrict_parse(rows[i][j])==0,(k,i,j,'nonzero ISP coupling')
        A=T*sp.Matrix([[restrict_parse(rows[i][j]) for j in pure] for i in pure])*T.inv()
        assert A[0,:]==sp.eye(6)[k,:]
        if k==1:assert A[1,:]==sp.eye(6)[5,:]
        cache.write_text(json.dumps({'source_sha256':sha(path),'variables':['ss'],
            'basis':'(phi,d_x1 phi,d_x2 phi,d_x3 phi,d_x4 phi,d_x1^2 phi)',
            'ISP_couplings_vanish_at_epsilon_zero':True,
            'matrix':[[str(v) for v in row] for row in A.tolist()]},indent=2),encoding='utf-8')
        print(f'extracted_A_x{k}_seconds_{time.time()-started:.1f}',flush=True)
    matrices.append(A)

physical=-sp.eye(6)/s
for mu,A in zip(mass_squares,matrices):physical-=mu*A/s**2
physical=physical.applyfunc(sp.cancel)
den=sp.lcm([sp.denom(e) for e in physical])
factors=sp.factor_list(den)
point=physical.subs(s,-64)
assert all(v.is_Rational for v in point)
out={'status':'published_pure_ray_system_extracted',
     'masses':[1,2,4,8],'variable':'ss = s',
     'basis':'s^(-1)*(phi,d_x1 phi,d_x2 phi,d_x3 phi,d_x4 phi,d_x1^2 phi)',
     'equation':'dH/ds=A_s H',
     'matrix':[[str(v) for v in row] for row in physical.tolist()],
     'common_denominator_factors':[[str(p),int(e)] for p,e in factors[1]],
     'at_s_minus64':[[str(v) for v in row] for row in point.tolist()],
     'source_hashes':input_hashes,'code_sha256':sha(Path(__file__)),
     'global_loop_transport_certified':False,'independent_flatness_checked':False,
     'absolute_periods_certified':False,'geometric_cycle_marking_certified':False}
(DEST/'system.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':out['status'],'common_denominator_factors':out['common_denominator_factors']}),flush=True)
