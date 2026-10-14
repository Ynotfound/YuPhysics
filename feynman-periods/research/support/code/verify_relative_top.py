"""Scalar replay of all new top quadratures, P8 and the physical readout."""
from verify_lower_regulator import *
from verify_lower_comparison import fixed,action_check
from verify_zero_transport import check_pi
from verify_horizontal_landmarks import inverse_check
OUT=FOLDER/'relative_top'
K=[(F(-1,64),F(0))]+[Z]*5
def scaled_check(raw,e,rho,re,saved):
    out,oe=vals(saved);z=scalar_sum([rho],[raw])
    check_radius((ab(rho)+re)*e+re*ab(raw)+distance(z,out[0]),oe)
def main():
    mp=OUT/'period_matrix.json';m=read(mp);tr=read(OUT/'trunks.json')
    assert m['code_sha256']==tr['code_sha256']==sha(ROOT/'code/certify_relative_top.py')
    for f,h in m['top_period_certificate_sha256'].items():assert sha(OUT/f)==h
    lp=FOLDER/'lower_comparison/period_matrix.json';lo=read(lp)
    assert m['lower_comparison_sha256']==sha(lp)
    lr=read(lp.with_name('replay.json'));assert lr['matrix_certificate_sha256']==sha(lp)
    assert lr['verifier_sha256']==sha(ROOT/'code/verify_lower_comparison.py')
    assert vfdecode(tr['coupling_D'])==K
    # Pure frame's first form is the fixed -64 times bare residue.
    system=read(FOLDER/'polarized_system.json')
    assert sp.Matrix(system['gauge_P'])[0,:]==sp.Matrix([[1,0,0,0,0,0]])
    archive={}
    for e in read(FOLDER/'high_trunk_transport.json')['segments']:archive[e['file'].replace('\\','/')]=e['sha256']
    for e in read(FOLDER/'horizontal_trunk/transport.json')['segments']:archive['horizontal_trunk/'+e['file']]=e['sha256']
    J=[Z]*6;je=F(0);C=mi();ce=F(0)
    for ix,e in enumerate(tr['vertical_steps']):
        p=FOLDER/e['file'];d=read(p);assert sha(p)==e['sha256']
        if ix:assert archive[e['file']]==sha(p)
        else:assert read(FOLDER/'first_trunk_step_replay.json')['certificate_sha256']==sha(p)
        R,re=calc_integral(d,K,ix==0);out,oe=vals(e['row']);check_row(J,je,R,re,C,ce,out,oe);J,je=out,oe
        if ix==0:
            C=[(F(x),F(y)) for rr,ii in zip(d['endpoint_center_real'],d['endpoint_center_imag']) for x,y in zip(rr,ii)]
            ce=F(d['bounds']['error_norm_upper'])
        else:C=fixed(d,'propagated_center_gaussian_integer');ce=F(d['propagated_error_upper'])
    assert (J,je)==vals(tr['vertical'])
    J=[Z]*6;je=F(0);C=mi();ce=F(0)
    for e in tr['horizontal_steps']:
        p=FOLDER/e['file'];d=read(p);assert sha(p)==e['sha256']==archive[e['file']]
        R,re=calc_integral(d,K);out,oe=vals(e['row']);check_row(J,je,R,re,C,ce,out,oe);J,je=out,oe
        C=fixed(d,'high_center_gaussian_integer');ce=F(d['high_error_upper'])
        if d['landmark_pole'] is not None:assert (J,je)==vals(tr['landmarks'][str(d['landmark_pole'])])
    top=[];errs=[];negative=False
    oldlower=read(DEST/'replay.json')
    for a in [1,9,25,49,81]:
        o=read(OUT/f'node_{a}.json');spth=FOLDER/f'node_frobenius/node_{a}.json';seed=read(spth)
        bp=FOLDER/f'node_branches/node_{a}/transport.json';bd=read(bp)
        lowp=DEST/f'node_{a}.json';low=read(lowp)
        assert o['seed_sha256']==sha(spth) and o['branch_sha256']==sha(bp)
        assert o['lower_node_sha256']==oldlower['node_certificates'][str(a)]==sha(lowp)
        assert o['trunks_sha256']==sha(OUT/'trunks.json')
        Q=real_matrix(sp.Matrix(seed['Q']));gq=[scalar_sum(K,[Q[6*i+j] for i in range(6)]) for j in range(6)]
        coeff=[vfdecode(x) for x in seed['frobenius_vector_coefficients_exact']]
        assert scalar_sum(gq,coeff[0])==Z
        x=y=F(0)
        for n,P in enumerate(coeff):
            v=scalar_sum(gq,P);x-=v[1]/(4*(F(n)+F(1,2)));y+=v[0]/(4*(F(n)+F(1,2)))
        q,qe=vals(o['local_q']);q=q[0]
        raw=rnorm(gq)*F(seed['bounds']['vector_error_Q_basis'])/2+distance((x,y),q)
        check_radius(raw,qe)
        if a==1:
            try:check_radius(raw,F(0))
            except AssertionError:negative=True
        assert len(o['steps'])==64
        for e,arch in zip(o['steps'],bd['segments']):
            p=bp.parent/e['file'];d=read(p);assert sha(p)==e['sha256']==arch['sha256']
            R,re=calc_integral(d,gq if d['frame']=='Q' else K)
            v=fixed(d,'input_center_gaussian_integer')
            out,oe=vals(e['q']);check_scalar(q,qe,R,re,v,F(d['input_error_upper']),out[0],oe);q,qe=out[0],oe
        cv=fixed(bd,'based_center_gaussian_integer');ce=F(bd['based_error_upper'])
        R,re=vals(tr['landmarks'][str(a)]);out,oe=vals(o['q_at_high_basepoint'])
        check_scalar(q,qe,R,re,cv,ce,out[0],oe,-1);q,qe=out[0],oe
        vb,vbe=vals(low['pure_at_physical_basepoint'])
        R,re=vals(tr['vertical']);out,oe=vals(o['unscaled_q_at_physical_basepoint'])
        check_scalar(q,qe,R,re,vb,vbe,out[0],oe,-1);q,qe=out[0],oe
        rv,re=vals(low['rho']);scaled_check(q,qe,rv[0],re,o['period'])
        v,e=vals(o['period']);top+=v;errs.append(e)
        assert a<225
    assert negative
    o=read(OUT/'sixth.json')
    for f,h in o['input_sha256'].items():assert sha(ROOT/f)==h
    seed=read(FOLDER/'zero_monodromy/frobenius.json');zt=read(FOLDER/'zero_monodromy/transport.json')
    g=(sp.Matrix([[sp.Rational(x) for x,y in K]])*sp.Matrix(seed['to_D'])).applyfunc(sp.cancel)
    assert all(not x.free_symbols for x in g)
    g0=[(rat(x),F(0)) for x in g];assert g0==vfdecode(o['zero_coupling'])
    N=real_matrix(sp.Matrix(seed['nilpotent']));assert mm(N,N)==mz()
    exact=[Z]*6
    for n,A in enumerate(seed['gauge_coefficients_exact']):
        A=vfdecode(A);r=[scalar_sum(g0,[A[6*i+j] for i in range(6)]) for j in range(6)]
        rhs=[(-y/4,x/4) for x,y in r];k=n+1
        rn=[scalar_sum(rhs,[N[6*i+j] for i in range(6)]) for j in range(6)]
        v=[(rhs[j][0]/k-rn[j][0]/k**2,rhs[j][1]/k-rn[j][1]/k**2) for j in range(6)]
        for j in range(6):
            w=scalar_sum(v,[N[6*i+j] for i in range(6)])
            assert (k*v[j][0]+w[0],k*v[j][1]+w[1])==rhs[j]
        exact=[(x+u,y+v) for (x,y),(u,v) in zip(exact,v)]
    J,je=vals(o['local_primitive'])
    check_radius(rnorm(g0)*(1+mnorm(N))*F(seed['bounds']['gauge_error_upper'])/4+sum(distance(x,y) for x,y in zip(J,exact)),je)
    for e,arch in zip(o['steps'],zt['segments']):
        p=FOLDER/'zero_monodromy'/e['file'];d=read(p);assert sha(p)==e['sha256']==arch['sha256']
        R,re=calc_integral(d,g0 if d['frame']=='zero_Q' else K)
        C=fixed(d,'input_center_gaussian_integer');ce=F(d['input_error_upper'])
        out,oe=vals(e['row']);check_row(J,je,R,re,C,ce,out,oe);J,je=out,oe
    R,re=vals(tr['landmarks']['0']);R=[(-x,-y) for x,y in R]
    W=fixed(zt,'based_frame_center_gaussian_integer');we=F(zt['based_frame_error_upper'])
    out,oe=vals(o['primitive_at_high_basepoint']);check_row(J,je,R,re,W,we,out,oe);J,je=out,oe
    lc=read(FOLDER/'lower_comparison/sixth_period.json')
    assert lr['sixth_certificate_sha256']==sha(FOLDER/'lower_comparison/sixth_period.json')
    w,we=vals(lc['delta_in_zero_constants'])
    q,qe=vals(read(OUT/'node_1.json')['q_at_high_basepoint'])
    out,oe=vals(o['q_at_high_basepoint']);check_scalar(q[0],qe,J,je,w,we,out[0],oe);q,qe=out[0],oe
    vb,vbe=vals(lc['pure_at_physical_basepoint']);R,re=vals(tr['vertical'])
    out,oe=vals(o['unscaled_q_at_physical_basepoint']);check_scalar(q,qe,R,re,vb,vbe,out[0],oe,-1)
    absp=FOLDER/'marked_lattice/absolute_periods.json';absolute=read(absp)
    rv=vfdecode([absolute['rho_center']])[0];re=F(absolute['rho_error_upper'])
    scaled_check(out[0],oe,rv,re,o['period']);v,e=vals(o['period']);top+=v;errs.append(e)
    # Independent rational log and Tate bounds.
    log=m['log2_interval'];n=log['terms'];assert n==120
    power=F(1,3);partial=F(0)
    for k in range(n):partial+=2*power/(2*k+1);power/=9
    upper=partial+2*power/(2*n+1)*F(9,8)
    assert F(log['lower'])==partial and F(log['upper'])==upper
    check_pi(zt['pi_certificate']);pl=F(zt['pi_certificate']['lower']);pu=F(zt['pi_certificate']['upper'])
    for key,low,high,component in [('boundary_top_period',8*pl*pl*partial,8*pu*pu*upper,0),
                                  ('top_Tate_period',-8*pu**3,-8*pl**3,1)]:
        v,e=vals(m[key]);check_radius((high-low)/2+abs(v[0][component]-(low+high)/2),e)
        assert v[0][1-component]==0;top+=v;errs.append(e)
    boundary=read(ROOT/'data/boundary_kummer_flags.json')
    assert m['boundary_data_sha256']==sha(ROOT/'data/boundary_kummer_flags.json')
    Zm=sp.Matrix(boundary['torus_double_residues_Z']);R=sp.eye(30)[boundary['anchor_rows_zero_based'],:]
    assert R*Zm==sp.eye(3)
    idx=boundary['oriented_edges'].index(['E12','F3']);tau=sp.eye(30)[idx,:]*(sp.eye(30)-Zm*R)
    assert list(map(str,tau))==m['corrected_boundary_torus_coefficients'] and tau*Zm==sp.zeros(1,3)
    mass=sp.Matrix(boundary['mass_flag_matrix'])
    assert (tau*mass*sp.Matrix([3,2,1,0]))[0]==1
    assert (tau*mass*sp.Matrix([0,2,4,6]))[0]==-2
    P7=vfdecode(lo['matrix_center']);P8=[]
    for i in range(7):P8+=P7[7*i:7*i+7]+[Z]
    P8+=top;assert P8==vfdecode(m['matrix_center'])
    assert list(map(F,m['top_row_errors']))==errs
    err=max(F(lo['matrix_infinity_error_upper']),sum(errs));assert err==F(m['matrix_infinity_error_upper'])<F(2)**-85
    # Independent inverse residual and vector readout; floating proposal is not trusted.
    rp=OUT/'physical_readout.json';r=read(rp)
    assert r['top_matrix_sha256']==sha(mp) and r['pure_absolute_sha256']==sha(absp) and r['lower_matrix_sha256']==sha(lp)
    PD=decode(absolute['physical_D_center_gaussian_integer']);PT=[PD[6*j+i] for i in range(6) for j in range(6)]
    pe=6*F(absolute['physical_D_error_upper']);assert pe==F(r['transpose_pure_error_upper'])
    ve=inverse_check(PT,pe,r['transpose_pure_inverse'],320)
    V=fixed(r['transpose_pure_inverse'],'center_gaussian_integer',320)
    k=-(partial+upper);ke=upper-partial
    rhs,re=vals(r['rhs']);rerrors=[]
    for j in range(6):
        b=P7[42+j];q=top[j];be=F(lo['beta_entry_error_upper'][j])
        z=(q[0]-k*b[0],q[1]-k*b[1])
        rerrors.append(errs[j]+(abs(k)+ke)*be+ke*ab(b)+distance(z,rhs[j]))
    assert re==max(rerrors)
    nu,ne=vals(r['nu_balanced_D']);action_check(V,ve,rhs,re,nu,ne)
    D=[F(-64)]*4+[F(4096),F(-1,64)]
    nb=[(D[i]*nu[i][0],D[i]*nu[i][1]) for i in range(6)]+[(k,F(0))]
    nes=[abs(x)*ne for x in D]+[ke]
    assert vfdecode(r['nu_B_center'])==nb and list(map(F,r['nu_B_entry_error_upper']))==nes
    I,Ie=vals(r['physical_integral']);assert I[0]==(-2*nb[5][0],-2*nb[5][1]) and Ie==2*nes[5]
    left,right=map(F,r['physical_integral_strict_real_interval'])
    assert left<I[0][0]-Ie and I[0][0]+Ie<right and abs(I[0][1])<=Ie and Ie<F(2)**-82
    # Exact pure cup row, with the geometric -2 calibration already established.
    Di=sp.diag(*map(sp.Rational,D))
    cup=-2*(Di.T*sp.Matrix(system['inverse_pairing'])*Di).inv()
    assert cup[0,:].applyfunc(sp.cancel)==sp.Matrix([[0,0,0,0,0,-2]])
    out={'status':'full_relative_rank_eight_comparison_and_physical_readout_independently_replayed',
       'matrix_sha256':sha(mp),'readout_sha256':sha(rp),'verifier_sha256':sha(__file__),
       'P8_error_less_than':'2^(-85)','physical_integral_error_less_than':'2^(-82)',
       'physical_integral_interval':r['physical_integral_strict_real_interval'],
       'top_scalar_node_steps':320,'top_scalar_zero_steps':64,'old_pure_ODE_restarted':False,
       'boundary_projection_and_log_coefficient_checked':True,'zero_tail_negative_control_rejected':negative,
       'physical_marking_and_duality_require_geometric_proof':'RELATIVE_TOP_COMPARISON.md',
       'full_integral_basis_certified':False,'fixed_pure_motivic_group_identified':False}
    save_json(OUT/'replay.json',out);print(json.dumps(out),flush=True)
if __name__=='__main__':main()

