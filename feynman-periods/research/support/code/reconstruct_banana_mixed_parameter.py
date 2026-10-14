"""Reconstruct the actual mass lifts by transversality and flatness.

The geometric uniqueness proof is in MIXED_PARAMETER_CONNECTION.md.
"""
from rational_flint import *
from itertools import permutations
import time
OUT=ROOT/'data/banana_mixed_parameter';PURE=ROOT/'data/banana_four_parameter'
def sgn(p):return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
def det(A):
    out=RF(0)
    for p in permutations(range(len(A))):
        v=RF(sgn(p))
        for i in range(len(A)):v=v*A[i][p[i]]
        out+=v
    return out
def minor(A,i,j):return [[x for b,x in enumerate(row) if b!=j] for a,row in enumerate(A) if a!=i]
def main():
    start=time.time();S=[RF(1),*X,RF(1)]
    A=[]
    for j in range(4):
        raw=read_matrix(PURE/f'A{j+1}.json')
        L=[[X[j]*S[r]*raw[r][c]/S[c]+int(r==c==j+1) for c in range(6)] for r in range(6)]
        assert L[0]==[RF(int(k==j+1)) for k in range(6)]
        A.append(L)
    d=json.loads((PURE/'determinant.json').read_text())
    C=[[RF(pi(v)) for v in row] for row in d['C']]
    dc=det(C)
    Y=[[3*sgn([0,1])*(-1)**(i+j)*det(minor(C,j,i))/dc for j in range(4)] for i in range(4)]
    # sgn([0,1]) is +1; the formula above is the adjugate transpose.
    assert all(dot(Y[i],[C[k][j]/3 for k in range(4)])==int(i==j) for i in range(4) for j in range(4))
    y=[sum(Y[k][l] for k in range(4)) for l in range(4)]
    total=sum(y)
    assert all(Y[i][i]==X[i]*total for i in range(4))
    R=[[(-Y[i][j]+X[i]*y[j]+X[j]*y[i]) for j in range(4)] for i in range(4)]
    distinct={(0,1,2):(R[0][1]+R[0][2]-R[0][3])/2,
      (0,1,3):(R[0][1]-R[0][2]+R[0][3])/2,
      (0,2,3):(-R[0][1]+R[0][2]+R[0][3])/2}
    distinct[1,2,3]=R[1][2]-distinct[0,1,2]
    def tensor(i,j,l):
        ids=sorted((i,j,l))
        if len(set(ids))==3:return distinct[tuple(ids)]
        repeated=next(k for k in ids if ids.count(k)>=2)
        other=next((k for k in ids if k!=repeated),repeated)
        return -X[repeated]*y[other]
    b=[[[dot([tensor(i,j,l) for l in range(4)],[C[l][k]/3 for l in range(4)]) for k in range(4)]
        for j in range(4)] for i in range(4)]
    assert all(sum(b[i][j][k] for j in range(4))==-int(i==k) for i in range(4) for k in range(4))
    assert all(b[i][i][k]==-X[i] for i in range(4) for k in range(4))
    print('cubic tensor and tangent coefficients constructed',round(time.time()-start,2),flush=True)
    couplings=[]
    for i in range(4):
        row=[]
        for j in range(4):
            l=(j+1)%4
            a=-b[i][j][l].theta(l)+b[i][l][l].theta(j)
            a-=dot(b[i][j],[A[l][r+1][l+1] for r in range(4)])
            a+=dot(b[i][l],[A[j][r+1][l+1] for r in range(4)])
            row.append([a,*b[i][j],RF(0)])
        couplings.append(row)
        save(OUT/f'couplings_{i+1}.json',{'mass_index':i+1,'rows':[[v.encode() for v in r] for r in row]})
        print('couplings',i+1,'seconds',round(time.time()-start,2),flush=True)
    assert all(couplings[i][j]==couplings[j][i] for i in range(4) for j in range(4))
    for i in range(4):
        assert couplings[i][i]==[-X[i]]*5+[RF(0)]
        assert [sum(couplings[i][j][k] for j in range(4)) for k in range(6)]==[RF(-int(k==i+1)) for k in range(6)]
    mixed=[]
    for j in range(4):
        B=[[RF(0) for _ in range(10)] for _ in range(10)]
        for r in range(6):
            for c in range(6):B[r][c]=A[j][r][c]
        for i in range(3):B[6+i][:6]=couplings[i][j]
        if j<3:B[9][6+j]=RF(1)
        else:
            B[9][0]=RF(-1)
            for k in range(6,9):B[9][k]=RF(-1)
        mixed.append(B)
        save(OUT/f'Alog{j+1}.json',{'direction':'theta_'+str(j+1),
          'basis':'(Omega,theta1 Omega,...,theta4 Omega,d1^2 Omega,beta1,beta2,beta3,e)',
          'matrix':[[v.encode() for v in row] for row in B]})
    inputs={f'A{k+1}.json':sha(PURE/f'A{k+1}.json') for k in range(4)}
    inputs['determinant.json']=sha(PURE/'determinant.json')
    save(OUT/'reconstruction.json',{'status':'geometric_mass_lifts_reconstructed_and_matrix_written',
      'input_sha256':inputs,'pure_replay_sha256':sha(PURE/'replay.json'),
      'generator_sha256':sha(__file__),'arithmetic_sha256':sha(ROOT/'code/rational_flint.py'),
      'cubic_tensor':{''.join(str(i+1) for i in key):v.encode() for key,v in distinct.items()},
      'coupling_sha256':{f'couplings_{i+1}.json':sha(OUT/f'couplings_{i+1}.json') for i in range(4)},
      'matrix_sha256':{f'Alog{j+1}.json':sha(OUT/f'Alog{j+1}.json') for j in range(4)},
      'full_flatness_certified':False,'absolute_rank10_comparison_certified':False})
    print('saved all four rank-ten matrices',round(time.time()-start,2),flush=True)
if __name__=='__main__':main()
