import sympy as s,json
from pathlib import Path
from sympy.matrices.normalforms import smith_normal_form
root=Path(__file__).resolve().parents[1]
Q=s.Matrix(json.loads((root/'data/banana_pure_ray/q6_traces_replay.json').read_text())['geometric_Gram_up_to_specified_orientations'])
C=s.Matrix([[2,2,2,-3,2,3],[0,0,1,0,1,0],[0,1,0,0,1,0],[0,1,1,-1,1,0],[0,0,0,0,-1,0],[-1,-1,-1,1,-1,-2]])
M=s.Matrix(json.loads((root/'data/banana_pure_ray/marked_lattice/marked_zero.json').read_text())['actual_T0'])
for name,v in zip('ABC',[[0,1,1,1,0,1],[1,0,0,0,0,1],[1,1,1,1,0,0]]):
    B=s.eye(6);B[:,v.index(1)]=s.Matrix(v)/2
    q=B.T*Q*B
    cc=B.inv()*C*B;mm=B.inv()*M*B
    print(name,'SNF',list(smith_normal_form(q,domain=s.ZZ).diagonal()))
    for label,g,sign in [('conjugation',cc,-1),('T0',mm,1)]:
        test=(g-sign*s.eye(6))*q.inv()
        print(label,all(x.q==1 for x in test),test)

