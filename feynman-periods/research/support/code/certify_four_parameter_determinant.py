"""Exact compact metric and determinant identities on the four-parameter base."""
from certify_banana_four_parameter import *
from itertools import permutations

def sign(p):return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
def determinant(A):
    return sum((sign(p)*__import__('functools').reduce(lambda a,b:a*b,(A[i][p[i]] for i in range(len(A))),ONE)
                for p in permutations(range(len(A)))),ZERO)
def main():
    data=json.loads((OUT/'cleared_system.json').read_text())
    D=poly(data['D']);E=poly(data['E']);N=[[[poly(v) for v in row] for row in A] for A in data['N']]
    M=[[poly(v) for v in row] for row in data['M']]
    x=list(CTX.gens());C=[[ZERO for j in range(4)] for i in range(4)]
    for i in range(4):
        others=[k for k in range(4) if k!=i]
        C[i][i]=2*((x[i]-1)**2+3*(sum((x[k]**2 for k in others),ZERO)-2*sum((x[others[a]]*x[others[b]] for a in range(3) for b in range(a+1,3)),ZERO)))
        for j in range(i+1,4):
            k,l=[a for a in range(4) if a not in (i,j)]
            C[i][j]=C[j][i]=3*(x[i]+x[j]-x[k]-x[l])**2-4*x[i]*x[j]-12*x[k]*x[l]-2*(x[i]+x[j])+6*(x[k]+x[l])-1
    Delta=-determinant(C)/27
    # Check the discriminant product on the root cover, not merely at sample masses.
    root_product=ONE
    for mask in range(16):
        root_product*=1+sum((((-1)**((mask>>j)&1))*x[j] for j in range(4)),ZERO)
    roots=Delta.compose(*[v*v for v in x])
    assert roots==root_product
    q=3*E
    assert Delta.total_degree()==8 and q.total_degree()==5
    for i in range(4):
        for j in range(4):assert 3*M[i+1][j+1]==E*x[i]*C[i][j]*x[j]
    assert M[5][5]==0 and all(M[5][i]==M[i][5]==0 for i in range(1,5))
    assert 4*q*M[0][5]==E*x[0]*Delta
    for i in range(4):
        t=sum((N[i][j][j] for j in range(6)),ZERO)
        r=2*x[i]*Delta*q*t+D*(3*x[i]*q*Delta.derivative(i)+2*(1+int(i==0))*q*Delta-2*x[i]*Delta*q.derivative(i))
        assert r==0,('trace',i)
    # Independent algebraic determinant route: use the special final row/column.
    detK_num=-M[0][5]**2*determinant([r[1:5] for r in M[1:5]])
    xprod=x[0]*x[1]*x[2]*x[3]
    assert 48*q*q*detK_num==E**6*Delta**3*x[0]**2*xprod**2
    out={'status':'compact_metric_discriminant_determinant_and_trace_exactly_checked',
      'cleared_system_sha256':sha(OUT/'cleared_system.json'),'code_sha256':sha(__file__),
      'Delta':enc(Delta),'q5':enc(q),
      'C':[[enc(v) for v in row] for row in C],
      'inverse_middle_pairing':'(1/3) diag(x) C diag(x)',
      'K05':'x1 Delta/(4 q5)',
      'detK':'x1^2 (x1 x2 x3 x4)^2 Delta^3/(48 q5^2)',
      'connection_trace':'dlog(q5)-dlog(x1)-sum_i dlog(xi)-(3/2)dlog(Delta)',
      'polarized_gauge_determinant':'K05',
      'polarized_connection_trace':'-sum_i dlog(xi)-(1/2)dlog(Delta)',
      'root_cover_discriminant_product_verified':True,
      'determinant_monodromy':'(-1)^(winding(Delta))',
      'geometric_trace_scale_in_published_frame':-2,
      'scale_extension_uses_geometric_proof':'FOUR_PARAMETER_CONNECTION.md',
      'mixed_scale_certified_by_this_file':False}
    save('determinant.json',out)
    print(json.dumps({k:out[k] for k in ['status','connection_trace','polarized_connection_trace']}),flush=True)
if __name__=='__main__':main()
