"""Independent scalar replay; never calls the new comparison generator."""
from verify_lower_regulator import *
from verify_zero_transport import check_pi
OUT=FOLDER/'lower_comparison'

def fixed(d,key,bits=None):
    B=d['bits'] if bits is None else bits
    return [(F(a)/(1<<B),F(b)/(1<<B)) for a,b in d[key]]
def action_check(M,me,v,ve,out,oe):
    err=max(distance(scalar_sum(M[6*i:6*i+6],v),out[i]) for i in range(6))
    check_radius((mnorm(M)+me)*ve+me*vnorm(v)+err,oe)
def main():
    op=OUT/'sixth_period.json';mp=OUT/'period_matrix.json';o=read(op);m=read(mp)
    assert o['code_sha256']==m['code_sha256']==sha(ROOT/'code/certify_lower_comparison.py')
    for cert in [o,m]:
        for f,h in cert['input_sha256'].items():assert sha(ROOT/f)==h
    for file,replay,key in [('zero_monodromy/frobenius.json','zero_monodromy/frobenius_replay.json','certificate_sha256'),
                            ('zero_monodromy/transport.json','zero_monodromy/replay.json','transport_sha256'),
                            ('marked_lattice/absolute_periods.json','marked_lattice/absolute_periods_replay.json','certificate_sha256')]:
        assert read(FOLDER/replay)[key]==sha(FOLDER/file)
    lower=read(DEST/'replay.json')
    assert lower['trunks_sha256']==sha(DEST/'trunk_integrals.json')
    assert lower['verifier_sha256']==sha(ROOT/'code/verify_lower_regulator.py')
    for a,h in lower['node_certificates'].items():assert h==sha(DEST/f'node_{a}.json')
    seed=read(FOLDER/'zero_monodromy/frobenius.json');zt=read(FOLDER/'zero_monodromy/transport.json')
    tr=read(DEST/'trunk_integrals.json');absd=read(FOLDER/'marked_lattice/absolute_periods.json')
    # Recover the coupling from the geometric lift and both saved gauges.
    source=sp.Matrix(read(FOLDER/'active_lower_system.json')['matrix'])[6,:6]
    Pi=sp.Matrix(read(FOLDER/'polarized_system.json')['inverse_gauge'])
    Di=sp.diag(*[1/sp.Rational(x) for x in D_VALUES])
    gd=(source*Pi*Di).applyfunc(sp.cancel)
    assert gd==sp.Matrix([[sp.Rational(a) for a,b in GG]])
    gqsp=(gd*sp.Matrix(seed['to_D'])).applyfunc(sp.cancel)
    assert all(not x.free_symbols for x in gqsp)
    gq=[(rat(x),F(0)) for x in gqsp];assert vfdecode(o['coupling_zero_Q'])==gq
    N=real_matrix(sp.Matrix(seed['nilpotent']));assert mm(N,N)==mz()
    coeff=[vfdecode(x) for x in seed['gauge_coefficients_exact']]
    exact=[Z]*6
    for n,A in enumerate(coeff):
        # Solve k J_k+J_k N=h g G_(k-1), then check the identity entry by entry.
        k=n+1
        rhs=[scalar_sum(gq,[A[6*i+j] for i in range(6)]) for j in range(6)]
        rhs=[(-v[1]/4,v[0]/4) for v in rhs]
        rn=[scalar_sum(rhs,[N[6*i+j] for i in range(6)]) for j in range(6)]
        row=[(rhs[j][0]/k-rn[j][0]/k**2,rhs[j][1]/k-rn[j][1]/k**2) for j in range(6)]
        for j in range(6):
            v=scalar_sum(row,[N[6*i+j] for i in range(6)])
            assert (k*row[j][0]+v[0],k*row[j][1]+v[1])==rhs[j]
        exact=[(a+c,b+d) for (a,b),(c,d) in zip(exact,row)]
    J,je=vals(o['local_primitive'])
    localraw=rnorm(gq)*F(seed['bounds']['gauge_error_upper'])*(1+mnorm(N))/4+sum(distance(x,y) for x,y in zip(exact,J))
    check_radius(localraw,je)
    rejected=False
    try:check_radius(localraw,F(0))
    except AssertionError:rejected=True
    assert rejected
    assert len(o['steps'])==len(zt['segments'])==64
    for entry,archive in zip(o['steps'],zt['segments']):
        p=FOLDER/'zero_monodromy'/entry['file'];st=read(p)
        assert entry['file']==archive['file'] and entry['sha256']==archive['sha256']==sha(p)
        R,re=calc_integral(st,gq if st['frame']=='zero_Q' else GG)
        C=fixed(st,'input_center_gaussian_integer');ce=F(st['input_error_upper'])
        out,oe=vals(entry['primitive_row']);check_row(J,je,R,re,C,ce,out,oe);J,je=out,oe
    assert (J,je)==vals(o['primitive_at_top'])
    W=fixed(zt,'based_frame_center_gaussian_integer');we=F(zt['based_frame_error_upper'])
    R,re=vals(tr['landmarks']['0']);R=[(-x,-y) for x,y in R]
    out,oe=vals(o['primitive_at_high_basepoint']);check_row(J,je,R,re,W,we,out,oe);J,je=out,oe
    nd=read(FOLDER/'node_branches/node_1/transport.json')
    v1=fixed(nd,'based_center_gaussian_integer');e1=F(nd['based_error_upper'])
    v6=fixed(zt,'sixth_vector_center_gaussian_integer');e6=F(zt['sixth_vector_error_upper'])
    dv,de=vals(o['delta_pure_at_high']);assert dv==[(a-c,b-d) for (a,b),(c,d) in zip(v6,v1)] and de==e1+e6
    Vi=fixed(zt['based_frame_inverse'],'center_gaussian_integer',zt['bits']);vie=F(zt['based_frame_inverse']['error_upper'])
    w,we=vals(o['delta_in_zero_constants']);action_check(Vi,vie,dv,de,w,we)
    inc,ie=vals(o['beta_monodromy_increment_normalized'])
    assert inc[0]==scalar_sum(J,w) and ie==rnorm(J)*we+je*(vnorm(w)+we)
    # Reject omission of the mixed monodromy increment.
    assert distance(inc[0],Z)>ie
    n1=read(DEST/'node_1.json');old,olde=vals(n1['beta_at_high_basepoint'])
    b,be=vals(o['beta_at_high_basepoint']);check_scalar(old[0],olde,J,je,w,we,b[0],be)
    inv=read(FOLDER/'high_trunk_inverse.json')
    V=fixed(inv,'inverse_center_gaussian_integer');Ve=F(inv['inverse_error_upper'])
    vb,vbe=vals(o['pure_at_physical_basepoint']);action_check(V,Ve,v6,e6,vb,vbe)
    R,re=vals(tr['vertical']);out,oe=vals(o['unscaled_beta_at_physical_basepoint'])
    check_scalar(b[0],be,R,re,vb,vbe,out[0],oe,-1);b,be=out[0],oe
    rho=vfdecode([absd['rho_center']])[0];rhoe=F(absd['rho_error_upper'])
    out,oe=vals(o['beta_period']);raw=scalar_sum([rho],[b])
    check_radius((ab(rho)+rhoe)*be+rhoe*ab(b)+distance(raw,out[0]),oe)
    assert oe<F(2)**-88
    # Pure cycle orientation: each old node is the unique vanishing sphere line
    # on the same path, hence differs by +/-1; intervals exclude the minus sign.
    P6=fixed(absd,'physical_B_center_gaussian_integer');p6e=F(absd['physical_B_error_upper'])
    Di=[F(-1,64)]*4+[F(1,4096),F(-64)]
    betas=[];errors=[];orientation=[]
    for j,a in enumerate([1,9,25,49,81]):
        d=read(DEST/f'node_{a}.json');pure,pe=vals(d['pure_at_physical_basepoint']);r,re=vals(d['rho'])
        valid=[]
        for sign in [1,-1]:
            ok=True
            for i in range(6):
                z=scalar_sum([r[0]],[pure[i]]);z=(sign*Di[i]*z[0],sign*Di[i]*z[1])
                e=abs(Di[i])*((ab(r[0])+re)*pe+re*ab(pure[i]))+p6e
                ref=P6[6*i+j]
                ok=ok and all(abs(z[k]-ref[k])<=e for k in [0,1])
            if ok:valid.append(sign)
        assert valid==[1];orientation.append({'node':a,'sign':1})
        v,e=vals(d['beta_period']);betas+=v;errors.append(e)
    # Same scaling and actual T0 image determine the sixth pure column.
    for i in range(6):
        z=scalar_sum([rho],[vb[i]]);z=(Di[i]*z[0],Di[i]*z[1])
        e=abs(Di[i])*((ab(rho)+rhoe)*vbe+rhoe*ab(vb[i]))+p6e
        assert all(abs(z[k]-P6[6*i+5][k])<=e for k in [0,1])
    betas+=out;errors.append(oe)
    pc,pe=check_pi(zt['pi_certificate'])
    t,te=vals(m['boundary_tate'])
    check_radius(4*(2*abs(pc)*pe+pe**2)+distance((-4*pc*pc,F(0)),t[0]),te)
    assert t[0][0]+te<0
    flags=read(ROOT/'data/minimality_cycle.json')['flags']
    flag=next(x for x in flags if x['edge']==['E12','F3'])
    assert flag==m['boundary_flag'] and flag['cycle_coefficient']==flag['divisor_on_first']==1
    expected=[]
    for i in range(6):expected.extend(P6[6*i:6*i+6]+[Z])
    expected+=betas+t
    assert vfdecode(m['matrix_center'])==expected
    assert list(map(F,m['beta_entry_error_upper']))==errors
    err=max(p6e,sum(errors)+te);assert F(m['matrix_infinity_error_upper'])==err<F(2)**-85
    final={'status':'sixth_mixed_period_and_lower_seven_comparison_independently_replayed',
           'sixth_certificate_sha256':sha(op),'matrix_certificate_sha256':sha(mp),'verifier_sha256':sha(__file__),
           'old_zero_steps_reused':64,'primitive_identity_checked_orders':len(coeff),
           'zero_Q_coupling_has_no_pole':True,'nonzero_mixed_monodromy_increment_checked':True,
           'negative_control_zero_tail_rejected':rejected,'cycle_orientation':orientation,
           'sixth_beta_error_less_than':'2^(-88)','matrix_infinity_error_less_than':'2^(-85)',
           'Tate_period':'(2*pi*i)^2','boundary_flag':['E12','F3'],'full_integral_basis_certified':False,
           'rank_eight_comparison_complete':False,'pure_ODE_restarted':False,'floating_point_used':False}
    save_json(OUT/'replay.json',final)
    print(json.dumps(final,ensure_ascii=False),flush=True)
if __name__=='__main__':main()

