"""Recover the actual T0 in the geometric six-cycle rational basis."""
from node_branch_transport import propagate_vector
from node_frobenius import *
from certify_horizontal_landmarks import inverse
from refined_transport import compose_fine
from itertools import product
from math import floor


def main():
    dest=FOLDER/'marked_lattice';dest.mkdir(exist_ok=True)
    zp=FOLDER/'zero_monodromy/transport.json';z=json.loads(zp.read_text())
    assert json.loads(zp.with_name('replay.json').read_text())['transport_sha256']==sha(zp)
    qp=FOLDER/'q6_traces_replay.json';qr=json.loads(qp.read_text())
    assert qr['trace_certificate_sha256']==sha(FOLDER/'q6_traces.json')
    Q=sp.Matrix(qr['geometric_Gram_up_to_specified_orientations'])
    bits=z['bits'];S=1<<bits;sources=[]
    n1p=FOLDER/'node_branches/node_1/transport.json';n1=json.loads(n1p.read_text())
    v1=[gs(1<<(bits-n1['bits']),v) for v in decode(n1['based_center_gaussian_integer'])]
    e1=F(n1['based_error_upper']);cols=[v1];radii=[e1];normalizations=[]
    for a,i in [(9,1),(25,2),(49,3),(81,4)]:
        p=FOLDER/f'node_branches/node_{a}/transport.json';d=json.loads(p.read_text())
        assert json.loads(p.with_name('replay.json').read_text())['transport_sha256']==sha(p)
        ref=d['based_reflection']
        T=[gs(1<<(bits-ref['bits']),v) for v in decode(ref['reflection_center_gaussian_integer'])]
        R=ma(T,ms((-1,0),mi(S)));re=F(ref['reflection_error_upper'])
        raw,er,rd=propagate_vector(R,re,v1,e1,bits)
        factor=F(1,int(Q[0,i]));exact=[gs(factor,v) for v in raw]
        col=[(dyadic_round(F(x,S),bits),dyadic_round(F(y,S),bits)) for x,y in exact]
        srd=max(F(ab(ga(x,gn(y))),S) for x,y in zip(exact,col))
        err=radius_up(abs(factor)*er+srd,bits)
        cols.append(col);radii.append(err)
        normalizations.append({'node':a,'q1i':str(Q[0,i]),
            'action_center_gaussian_integer':encode(raw),'action_error_upper':str(er),'action_rounding_error':str(rd),
            'column_center_gaussian_integer':encode(col),'column_error_upper':str(err),'scaling_rounding_error':str(srd)})
        sources.append({'node':a,'transport_sha256':sha(p),'replay_sha256':sha(p.with_name('replay.json'))})
    cols.append(decode(z['sixth_vector_center_gaussian_integer']))
    radii.append(F(z['sixth_vector_error_upper']))
    P=[cols[j][i] for i in range(6) for j in range(6)];pe=sum(radii)
    inv=inverse(P,pe,bits);V=decode(inv['center_gaussian_integer']);ve=F(inv['error_upper'])
    T=decode(z['monodromy_center_gaussian_integer']);te=F(z['monodromy_error_upper'])
    TP,tpe,rd1=compose_fine(T,te,P,pe,bits)
    M,me,rd2=compose_fine(V,ve,TP,tpe,bits)
    entries=[];marked=[]
    for c in M:
        center=gs(F(1,S),c);lo=ceil(2*(center[0]-me));hi=floor(2*(center[0]+me))
        assert lo==hi and abs(center[1])<=me
        marked.append(F(lo,2))
        entries.append({'twice_real_interval':[str(2*(center[0]-me)),str(2*(center[0]+me))],
                        'twice_imaginary_interval':[str(2*(center[1]-me)),str(2*(center[1]+me))],
                        'unique_integer':lo})
    Mq=sp.Matrix(6,6,[sp.Rational(v) for v in marked])
    assert Mq.T*Q*Mq==Q and (Mq-sp.eye(6))**2==sp.zeros(6)
    assert all(v.q==1 for v in Mq) and Mq.det()==1
    # In this now-proved integral L automorphism, every index-two T/L class is fixed.
    candidates=[]
    for v in product(range(2),repeat=6):
        if not any(v):continue
        w=sp.Matrix(v)
        if (w.T*Q*w)[0]%8:continue
        fixed=all(int(x)%2==0 for x in (Mq-sp.eye(6))*w)
        candidates.append({'numerator':list(v),'fixed_by_actual_T0':fixed})
    old=json.loads((ROOT/'data/recovered_lattice.json').read_text())
    out={'status':'actual_T0_in_six_cycle_basis_certified',
         'bits':bits,'zero_transport_sha256':sha(zp),'Q6_replay_sha256':sha(qp),
         'node_one_sha256':sha(n1p),'code_sha256':sha(__file__),'sources':sources,
         'basis_formula':'p1=v1; pi=(Ti-I)v1/q1i for i=2..5; p6=T0 v1; all pi share one unknown amplitude',
         'normalizations':normalizations,'basis_center_gaussian_integer':encode(P),
         'column_error_upper':list(map(str,radii)),'basis_error_upper':str(pe),'basis_inverse':inv,
         'T0_times_basis_center_gaussian_integer':encode(TP),'T0_times_basis_error_upper':str(tpe),
         'first_rounding_error':str(rd1),'second_rounding_error':str(rd2),
         'marked_center_gaussian_integer':encode(M),'marked_error_upper':str(me),
         'half_integer_certificates':entries,
         'actual_T0':[[str(v) for v in row] for row in Mq.tolist()],
         'matches_old_candidate':Mq==sp.Matrix(old['T0']),
         'even_index_two_classes':candidates,
         'actual_T0_fixed_classes':[x['numerator'] for x in candidates if x['fixed_by_actual_T0']],
         'specific_half_cycle_completion_certified':False,'absolute_periods_certified':False}
    save_json(dest/'marked_zero.json',out)
    print(json.dumps({'status':out['status'],'error_dyadic':str(dyadic_upper(me)),
                      'matches_old_candidate':out['matches_old_candidate'],'T0':out['actual_T0'],
                      'fixed_classes':out['actual_T0_fixed_classes']}))


if __name__=='__main__':main()
