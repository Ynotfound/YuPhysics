"""Exact Lie-algebra witness from the now geometrically certified Q6."""
from node_frobenius import *
from itertools import combinations


def main():
    rp=FOLDER/'q6_traces_replay.json';r=json.loads(rp.read_text())
    assert r['trace_certificate_sha256']==sha(FOLDER/'q6_traces.json')
    Q=sp.Matrix(r['geometric_Gram_up_to_specified_orientations']);I=sp.eye(6)
    S=[I+I[:,i]*I[:,i].T*Q for i in range(6)]
    for R in S: assert R.T*Q*R==Q and R*R==I and R.det()==-1
    U=S[0]*S[1];A=U-I;assert A**3==sp.zeros(6)
    X=A-A*A/2;assert X.T*Q+Q*X==sp.zeros(6) and X**3==sp.zeros(6)
    assert I+X+X*X/2==U
    words=[[],[3],[4],[5],[6],[3,1],[3,4],[3,6],[4,1],[4,5],[4,6],[5,1],[5,6],[6,1],[3,1,5]]
    basis=[]
    for word in words:
        Y=X
        for i in word: Y=S[i-1]*Y*S[i-1]
        basis.append(Y)
    positions=list(combinations(range(6),2))
    C=sp.Matrix([[int((Q*Y)[i,j]) for Y in basis] for i,j in positions])
    assert C.det()!=0
    enc=lambda M:[[str(x) for x in row] for row in M.tolist()]
    out={'status':'full_orthogonal_Lie_algebra_witness_computed',
         'Q6_trace_replay_sha256':sha(rp),'code_sha256':sha(__file__),
         'Gram':enc(Q),'reflections':[enc(v) for v in S],'unipotent':enc(U),'logarithm':enc(X),
         'conjugation_words':words,'Lie_basis':[enc(Y) for Y in basis],
         'coordinate_positions_zero_based':[list(v) for v in positions],
         'coordinate_matrix':enc(C),'coordinate_determinant':str(C.det()),'Lie_dimension':15,
         'Zariski_closure_conclusion':'O(Q6) over Q, hence O(6,C) after complexification',
         'proof_obligation':'Polynomial one-parameter subgroup from powers of U; SO is connected in characteristic zero.',
         'arithmetic_index_in_integral_orthogonal_group_computed':False,
         'absolute_motivic_group_identified':False}
    save_json(FOLDER/'geometric_monodromy.json',out)
    print(json.dumps({'status':out['status'],'coordinate_determinant':out['coordinate_determinant'],'dimension':15}))


if __name__=='__main__': main()
