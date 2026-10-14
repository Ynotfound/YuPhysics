"""Fraction-only replay of the Lie witness, without SymPy rank routines."""
import json,hashlib
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
ROOT=Path(__file__).resolve().parents[1];DEST=ROOT/'data/banana_pure_ray'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
parse=lambda A:[[F(v) for v in row] for row in A]
def mul(A,B):
    return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def add(A,B,scale=1):
    return [[a+scale*b for a,b in zip(ar,br)] for ar,br in zip(A,B)]
def transpose(A): return [list(row) for row in zip(*A)]
def eye(n): return [[F(i==j) for j in range(n)] for i in range(n)]
def det(A):
    A=[row[:] for row in A];out=F(1);n=len(A)
    for j in range(n):
        pivot=next((i for i in range(j,n) if A[i][j]),None)
        if pivot is None:return F(0)
        if pivot!=j:A[j],A[pivot]=A[pivot],A[j];out=-out
        p=A[j][j];out*=p
        for i in range(j+1,n):
            c=A[i][j]/p
            for k in range(j,n):A[i][k]-=c*A[j][k]
    return out
def main():
    p=DEST/'geometric_monodromy.json';d=json.loads(p.read_text());rp=DEST/'q6_traces_replay.json'
    assert d['Q6_trace_replay_sha256']==sha(rp)
    assert d['code_sha256']==sha(ROOT/'code/certify_geometric_monodromy.py')
    Q=parse(d['Gram']);assert Q==json.loads(rp.read_text())['geometric_Gram_up_to_specified_orientations']
    I=eye(6);Z=[[F(0)]*6 for _ in range(6)]
    S=[]
    for i,encoded in enumerate(d['reflections']):
        R=[row[:] for row in I];R[i]=[x+y for x,y in zip(R[i],Q[i])]
        assert R==parse(encoded) and mul(R,R)==I
        assert mul(mul(transpose(R),Q),R)==Q and det(R)==-1
        S.append(R)
    U=mul(S[0],S[1]);assert U==parse(d['unipotent'])
    A=add(U,I,-1);A2=mul(A,A);assert mul(A2,A)==Z and A2!=Z
    X=add(A,A2,-F(1,2));assert X==parse(d['logarithm'])
    X2=mul(X,X);assert mul(X2,X)==Z and add(add(I,X),X2,F(1,2))==U
    basis=[]
    for word,encoded in zip(d['conjugation_words'],d['Lie_basis']):
        Y=X
        for i in word:
            assert 1<=i<=6
            Y=mul(mul(S[i-1],Y),S[i-1])
        assert Y==parse(encoded) and add(mul(transpose(Y),Q),mul(Q,Y))==Z
        basis.append(Y)
    assert len(basis)==15
    positions=list(combinations(range(6),2));assert [list(v) for v in positions]==d['coordinate_positions_zero_based']
    C=[[mul(Q,Y)[i][j] for Y in basis] for i,j in positions]
    assert C==parse(d['coordinate_matrix']);detC=det(C)
    assert detC!=0 and detC==F(d['coordinate_determinant'])
    bad=[row[:] for row in C]
    for row in bad:row[1]=row[0]
    assert det(bad)==0
    out={'status':'full_orthogonal_Lie_witness_independently_replayed',
         'certificate_sha256':sha(p),'verifier_sha256':sha(__file__),
         'arithmetic':'Fraction-only matrix products and Gaussian elimination',
         'nonzero_determinant':str(detC),'Lie_dimension':15,
         'duplicate_column_negative_control_rejected':True,
         'absolute_motivic_group_identified':False}
    (DEST/'geometric_monodromy_replay.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out))
if __name__=='__main__':main()
