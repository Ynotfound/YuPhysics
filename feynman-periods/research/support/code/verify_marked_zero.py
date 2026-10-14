"""Independent replay of relative geometric normalization and half-integrality."""
from node_frobenius import *
from verify_node_branches import check_action
from verify_zero_transport import check_product
from verify_horizontal_landmarks import inverse_check
from itertools import product
from math import floor


def main():
    p=FOLDER/'marked_lattice/marked_zero.json';d=json.loads(p.read_text())
    zp=FOLDER/'zero_monodromy/transport.json';z=json.loads(zp.read_text())
    qp=FOLDER/'q6_traces_replay.json';Q=sp.Matrix(json.loads(qp.read_text())['geometric_Gram_up_to_specified_orientations'])
    n1p=FOLDER/'node_branches/node_1/transport.json';n1=json.loads(n1p.read_text())
    assert d['code_sha256']==sha(ROOT/'code/certify_marked_zero.py')
    assert d['zero_transport_sha256']==sha(zp) and d['Q6_replay_sha256']==sha(qp) and d['node_one_sha256']==sha(n1p)
    assert json.loads(zp.with_name('replay.json').read_text())['transport_sha256']==sha(zp)
    bits=d['bits'];S=1<<bits;assert bits==z['bits']
    v1=[gs(1<<(bits-n1['bits']),v) for v in decode(n1['based_center_gaussian_integer'])]
    e1=F(n1['based_error_upper']);cols=[v1];radii=[e1]
    assert [v['node'] for v in d['normalizations']]==[9,25,49,81]
    for i,(v,source) in enumerate(zip(d['normalizations'],d['sources']),1):
        a=v['node'];src=FOLDER/f'node_branches/node_{a}/transport.json';n=json.loads(src.read_text())
        assert source['node']==a and source['transport_sha256']==sha(src)
        assert source['replay_sha256']==sha(src.with_name('replay.json'))
        ref=n['based_reflection']
        R=[gs(1<<(bits-ref['bits']),x) for x in decode(ref['reflection_center_gaussian_integer'])]
        R=ma(R,ms((-1,0),mi(S)));re=F(ref['reflection_error_upper'])
        raw=decode(v['action_center_gaussian_integer']);er=F(v['action_error_upper'])
        assert check_action(R,re,v1,e1,raw,er,bits)==F(v['action_rounding_error'])
        q=F(v['q1i']);assert q==int(Q[0,i]) and q!=0
        col=decode(v['column_center_gaussian_integer']);err=F(v['column_error_upper'])
        rd=max(F(ab(ga(gs(1/q,x),gn(y))),S) for x,y in zip(raw,col))
        bound=er/abs(q)+rd
        assert bound<=err<bound+F(1,S) and rd==F(v['scaling_rounding_error'])
        cols.append(col);radii.append(err)
    cols.append(decode(z['sixth_vector_center_gaussian_integer']));radii.append(F(z['sixth_vector_error_upper']))
    P=decode(d['basis_center_gaussian_integer']);pe=F(d['basis_error_upper'])
    assert P==[cols[j][i] for i in range(6) for j in range(6)]
    assert radii==list(map(F,d['column_error_upper'])) and pe==sum(radii)
    ve=inverse_check(P,pe,d['basis_inverse'],bits);V=decode(d['basis_inverse']['center_gaussian_integer'])
    T=decode(z['monodromy_center_gaussian_integer']);te=F(z['monodromy_error_upper'])
    TP=decode(d['T0_times_basis_center_gaussian_integer']);tpe=F(d['T0_times_basis_error_upper'])
    assert check_product(T,te,P,pe,TP,tpe,bits)==F(d['first_rounding_error'])
    M=decode(d['marked_center_gaussian_integer']);me=F(d['marked_error_upper'])
    assert check_product(V,ve,TP,tpe,M,me,bits)==F(d['second_rounding_error'])
    values=[]
    for c,entry in zip(M,d['half_integer_certificates']):
        center=gs(F(1,S),c);re=[2*(center[0]-me),2*(center[0]+me)];im=[2*(center[1]-me),2*(center[1]+me)]
        assert list(map(str,re))==entry['twice_real_interval'] and list(map(str,im))==entry['twice_imaginary_interval']
        low,high=ceil(re[0]),floor(re[1]);assert low==high==entry['unique_integer'] and im[0]<=0<=im[1]
        values.append(sp.Rational(low,2))
    actual=sp.Matrix(6,6,values)
    assert actual==sp.Matrix(d['actual_T0']) and actual.T*Q*actual==Q
    assert (actual-sp.eye(6))**2==sp.zeros(6) and actual.det()==1 and all(x.q==1 for x in actual)
    assert actual[:,0]==sp.eye(6)[:,5]
    candidates=[]
    for mask in range(1,64):
        v=[(mask>>i)&1 for i in range(6)]
        sq=sum(v[i]*int(Q[i,j])*v[j] for i in range(6) for j in range(6))
        if sq%8:continue
        fixed=all(sum(int(actual[i,j]-(i==j))*v[j] for j in range(6))%2==0 for i in range(6))
        candidates.append({'numerator':v,'fixed_by_actual_T0':fixed})
    assert sorted(candidates,key=lambda v:v['numerator'])==d['even_index_two_classes']
    fixed=sorted(v['numerator'] for v in candidates if v['fixed_by_actual_T0'])
    assert fixed==d['actual_T0_fixed_classes'] and len(candidates)==11 and len(fixed)==3
    out={'status':'actual_marked_T0_and_three_possible_saturation_cosets_independently_replayed',
         'certificate_sha256':sha(p),'verifier_sha256':sha(__file__),
         'matrix_error_dyadic':str(dyadic_upper(me)),'half_integer_isolation_entries':36,
         'actual_T0':d['actual_T0'],'remaining_even_cosets':fixed,
         'geometric_input':'[T:L]=2 and T0 integral on T imply coordinates in (1/2)Z.',
         'coefficient_generator_called':False,'floating_point_used':False,
         'specific_half_cycle_completion_certified':False,'absolute_periods_certified':False}
    save_json(p.with_name('marked_zero_replay.json'),out);print(json.dumps(out))


if __name__=='__main__':main()
