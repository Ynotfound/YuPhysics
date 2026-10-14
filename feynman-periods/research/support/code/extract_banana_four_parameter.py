"""Extract all four pure matrices over QQ(x1,x2,x3,x4), with checkpoints."""
import json,time
from four_parameter_algebra import *

def main():
    manifest=json.loads((SOURCE/'manifest.json').read_text(encoding='utf-8'))
    pure=[4,5,6,7,8,14]; isp=list(range(9,14)); scale=[1,1,1,1,1,2]
    for k in range(4):
        path=SOURCE/f'Diffeqm{k+1}.m'; target=f'A{k+1}.json'
        expected=next(r['sha256'] for r in manifest['files'] if r['name']==path.name)
        assert sha(path)==expected
        if (OUT/target).exists():
            assert json.loads((OUT/target).read_text())['source_sha256']==expected
            print('reuse',target,flush=True); continue
        started=time.time(); rows=source_rows(path)
        assert len(rows)==15 and all(len(r)==15 for r in rows)
        for i in pure:
            for j in isp: assert parse(rows[i][j])==0,(k,i,j)
        A=[]
        for i in range(6):
            A.append([FF(scale[i])*parse(rows[pure[i]][pure[j]])/scale[j] for j in range(6)])
            print('extracted',k+1,'row',i,'seconds',round(time.time()-started,2),flush=True)
        assert A[0]==[FF(int(j==k+1)) for j in range(6)]
        if k==0: assert A[1]==[FF(int(j==5)) for j in range(6)]
        save(target,{'source_sha256':expected,'variables':list(NAMES),'epsilon':0,
            'basis':'W=(phi,d1 phi,d2 phi,d3 phi,d4 phi,d1^2 phi)',
            'ISP_couplings_zero':True,'matrix':[[enc(x) for x in row] for row in A],
            'entry_num_degrees':[[degree(x.numer) for x in row] for row in A],
            'entry_den_degrees':[[degree(x.denom) for x in row] for row in A]})
        print('saved',target,sha(OUT/target),flush=True)
    path=SOURCE/'ZCYInverse.m'
    expected=next(r['sha256'] for r in manifest['files'] if r['name']==path.name)
    assert sha(path)==expected
    if not (OUT/'K.json').exists():
        rows=source_rows(path); assert len(rows)==6 and all(len(r)==6 for r in rows)
        def rational_entry(text):
            if 'Sqrt' not in text: return parse(text)
            # Work on the square-root cover, then prove descent coefficientwise.
            t=re.sub(r'(m[1-4]2)\^\(3/2\)',lambda m:'ROOT'+m[1]+'^3',text)
            t=re.sub(r'Sqrt\[(m[1-4]2)\]',lambda m:'ROOT'+m[1],t)
            t=re.sub(r'\bm[1-4]2\b',lambda m:'('+m[0]+'^2)',t)
            t=t.replace('ROOT','')
            v=parse(t)
            assert all(all(k%2==0 for k in e) for p in (v.numer,v.denom) for e in p)
            descend=lambda p:RR.from_dict({tuple(k//2 for k in e):c for e,c in p.items()})
            return FF.new(descend(v.numer),descend(v.denom))
        K=[[rational_entry(rows[i][j])/scale[i]/scale[j] for j in range(6)] for i in range(6)]
        assert K==transpose(K)
        save('K.json',{'source_sha256':expected,'variables':list(NAMES),
            'meaning':'inverse published pairing in W basis: T^(-t) K_raw T^(-1)',
            'matrix':[[enc(x) for x in row] for row in K]})
        print('saved K',sha(OUT/'K.json'),flush=True)
    save('extraction.json',{'status':'four_variable_pure_matrices_and_inverse_pairing_extracted',
         'files':{p.name:sha(p) for p in [*(OUT/f'A{k+1}.json' for k in range(4)),OUT/'K.json']},
         'extractor_sha256':sha(__file__),'algebra_sha256':sha(ROOT/'code/four_parameter_algebra.py'),
         'flatness_certified':False,'mixed_four_parameter_system_certified':False})

if __name__=='__main__': main()
