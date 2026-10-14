"""Actual relative log-volume periods on collapsing thimbles; reuse pure certificates."""
from certify_lower_comparison import read,unpack,fixed,action
from certify_lower_regulator import *
OUT=FOLDER/'relative_top'
K=[(F(-1,64),F(0))]+[Z]*5

def trunks():
    f=read(FOLDER/'first_trunk_step.json')
    J,je=integrated_row(f,K,True);logs=[{'file':'first_trunk_step.json','sha256':sha(FOLDER/'first_trunk_step.json'),'row':record(J,je)}]
    C=[(F(x),F(y)) for rr,ii in zip(f['endpoint_center_real'],f['endpoint_center_imag']) for x,y in zip(rr,ii)]
    ce=F(f['bounds']['error_norm_upper'])
    for p in sorted((FOLDER/'high_trunk').glob('step_*.json')):
        d=read(p);R,re=integrated_row(d,K);J,je=append_row(J,je,R,re,C,ce)
        C=fixed(d,'propagated_center_gaussian_integer');ce=F(d['propagated_error_upper'])
        logs.append({'file':str(p.relative_to(FOLDER)).replace('\\','/'),'sha256':sha(p),'row':record(J,je)})
    vertical=record(J,je);J=[Z]*6;je=F(0);C=mi();ce=F(0);landmarks={};hl=[]
    for p in sorted((FOLDER/'horizontal_trunk').glob('step_*.json')):
        d=read(p);R,re=integrated_row(d,K);J,je=append_row(J,je,R,re,C,ce)
        C=fixed(d,'high_center_gaussian_integer');ce=F(d['high_error_upper'])
        hl.append({'file':str(p.relative_to(FOLDER)).replace('\\','/'),'sha256':sha(p),'row':record(J,je)})
        if d['landmark_pole'] is not None:landmarks[str(d['landmark_pole'])]=record(J,je)
    out={'code_sha256':sha(__file__),'coupling_D':encode(K),'vertical':vertical,'landmarks':landmarks,
         'vertical_steps':logs,'horizontal_steps':hl}
    save_json(OUT/'trunks.json',out);return out

def scaled_result(q,e,rho,rhoe):
    raw=gm(rho,q);out=roundv([raw])[0]
    return record([out],radius_up((ab(rho)+rhoe)*e+rhoe*ab(q)+ab(ga(raw,gn(out))),BITS))

def node(a,tr):
    spath=FOLDER/f'node_frobenius/node_{a}.json';d=read(spath);Q=real_matrix(sp.Matrix(d['Q']))
    P=[vfdecode(v) for v in d['frobenius_vector_coefficients_exact']];gq=rowmat(K,Q);h=(F(0),F(1,4))
    assert dot(gq,P[0])==Z
    exact=gm(h,sumg(gs(1/(F(n)+F(1,2)),dot(gq,v)) for n,v in enumerate(P)))
    q=roundv([exact])[0];qe=radius_up(2*ab(h)*rn(gq)*F(d['bounds']['vector_error_Q_basis'])+ab(ga(exact,gn(q))),BITS)
    init=record([q],qe);logs=[]
    bp=FOLDER/f'node_branches/node_{a}/transport.json';bd=read(bp)
    for entry in bd['segments']:
        p=bp.parent/entry['file'];st=read(p);assert sha(p)==entry['sha256']
        R,re=integrated_row(st,gq if st['frame']=='Q' else K)
        v=fixed(st,'input_center_gaussian_integer')
        q,qe=scalar_update(q,qe,R,re,v,F(st['input_error_upper']))
        logs.append({'file':p.name,'sha256':sha(p),'q':record([q],qe)})
    cv=fixed(bd,'based_center_gaussian_integer');ce=F(bd['based_error_upper'])
    row,re=unpack(tr['landmarks'][str(a)]);q,qe=scalar_update(q,qe,row,re,cv,ce,-1);atc=record([q],qe)
    # The already-certified actual node and its Morse amplitude are unchanged.
    lowerp=DEST/f'node_{a}.json';low=read(lowerp)
    vb,vbe=unpack(low['pure_at_physical_basepoint'])
    row,re=unpack(tr['vertical']);q,qe=scalar_update(q,qe,row,re,vb,vbe,-1)
    rv,re=unpack(low['rho']);period=scaled_result(q,qe,rv[0],re)
    out={'status':'actual_relative_thimble_log_volume_certified','code_sha256':sha(__file__),'node':a,
         'seed_sha256':sha(spath),'branch_sha256':sha(bp),'lower_node_sha256':sha(lowerp),
         'trunks_sha256':sha(OUT/'trunks.json'),'local_q':init,'steps':logs,'q_at_high_basepoint':atc,
         'unscaled_q_at_physical_basepoint':record([q],qe),'period':period}
    save_json(OUT/f'node_{a}.json',out)
    print(json.dumps({'node':a,'q':[float(x) for x in vfdecode(period['center'])[0]],'error':float(F(period['error_upper']))}),flush=True)
    return out

def sixth(tr):
    seedp=FOLDER/'zero_monodromy/frobenius.json';zp=FOLDER/'zero_monodromy/transport.json'
    seed=read(seedp);zt=read(zp);n1=read(OUT/'node_1.json');toD=sp.Matrix(seed['to_D'])
    g0s=(sp.Matrix([[sp.Rational(x) for x,y in K]])*toD).applyfunc(sp.cancel)
    assert all(not x.free_symbols for x in g0s)
    g0=[(rat(x),F(0)) for x in g0s];N=real_matrix(sp.Matrix(seed['nilpotent']))
    coeff=[vfdecode(x) for x in seed['gauge_coefficients_exact']]
    raw=[Z]*6
    for n,C in enumerate(coeff):
        v=rowmat(g0,C);vN=rowmat(v,N)
        raw=va(raw,[gm((0,F(1,4)),ga(gs(F(1,n+1),x),gn(gs(F(1,(n+1)**2),y)))) for x,y in zip(v,vN)])
    J=roundv(raw);je=radius_up(rn(g0)*(1+mnorm(N))*F(seed['bounds']['gauge_error_upper'])/4+rn([ga(x,gn(y)) for x,y in zip(raw,J)]),BITS)
    local=record(J,je);logs=[]
    for e in zt['segments']:
        p=zp.parent/e['file'];st=read(p);R,re=integrated_row(st,g0 if st['frame']=='zero_Q' else K)
        C=fixed(st,'input_center_gaussian_integer');ce=F(st['input_error_upper'])
        J,je=append_row(J,je,R,re,C,ce);logs.append({'file':p.name,'sha256':sha(p),'row':record(J,je)})
    W=fixed(zt,'based_frame_center_gaussian_integer');we=F(zt['based_frame_error_upper'])
    row,re=unpack(tr['landmarks']['0']);J,je=append_row(J,je,[gn(x) for x in row],re,W,we)
    # Reuse the separately verified W^-1(v6-v1) and v6 at the physical basepoint.
    lcpath=FOLDER/'lower_comparison/sixth_period.json';lc=read(lcpath)
    w,we=unpack(lc['delta_in_zero_constants'])
    q,qe=unpack(n1['q_at_high_basepoint']);q,qe=scalar_update(q[0],qe,J,je,w,we);atc=record([q],qe)
    vb,vbe=unpack(lc['pure_at_physical_basepoint'])
    row,re=unpack(tr['vertical']);q,qe=scalar_update(q,qe,row,re,vb,vbe,-1)
    absp=FOLDER/'marked_lattice/absolute_periods.json';absolute=read(absp)
    rho=vfdecode([absolute['rho_center']])[0];rhoe=F(absolute['rho_error_upper'])
    period=scaled_result(q,qe,rho,rhoe)
    out={'status':'actual_relative_T0_thimble_log_volume_certified','code_sha256':sha(__file__),
         'input_sha256':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in [seedp,zp,lcpath,absp,OUT/'node_1.json',OUT/'trunks.json']},
         'zero_coupling':encode(g0),'local_primitive':local,'steps':logs,'primitive_at_high_basepoint':record(J,je),
         'q_at_high_basepoint':atc,'unscaled_q_at_physical_basepoint':record([q],qe),'period':period}
    save_json(OUT/'sixth.json',out)
    print(json.dumps({'cycle':6,'q':[float(x) for x in vfdecode(period['center'])[0]],'error':float(F(period['error_upper']))}),flush=True)
    return out

def assemble(ns,six):
    lowpath=FOLDER/'lower_comparison/period_matrix.json';low=read(lowpath)
    pi=read(FOLDER/'zero_monodromy/transport.json')['pi_certificate']
    pl,pu=F(pi['lower']),F(pi['upper'])
    n=120;l=2*sum(F(1,(2*k+1)*3**(2*k+1)) for k in range(n))
    u=l+2*F(1,(2*n+1)*3**(2*n+1))*F(9,8)
    # q7=-2*(2*pi*i)^2*log(2)=+8*pi^2*log(2).
    q7l,q7u=8*pl**2*l,8*pu**2*u
    q7=roundv([((q7l+q7u)/2,F(0))])[0]
    q7e=radius_up((q7u-q7l)/2+abs(q7[0]-(q7l+q7u)/2),BITS)
    # q8=(2*pi*i)^3=-8*pi^3*i.
    t3l,t3u=-8*pu**3,-8*pl**3
    t3=roundv([(F(0),(t3l+t3u)/2)])[0]
    t3e=radius_up((t3u-t3l)/2+abs(t3[1]-(t3l+t3u)/2),BITS)
    top=[];errs=[]
    for d in ns+[six]:
        v,e=unpack(d['period']);top.append(v[0]);errs.append(e)
    top += [q7,t3];errs += [q7e,t3e]
    P7=vfdecode(low['matrix_center']);P8=[]
    for i in range(7):P8.extend(P7[7*i:7*i+7]+[Z])
    P8+=top;err=max(F(low['matrix_infinity_error_upper']),sum(errs))
    boundary=read(ROOT/'data/boundary_kummer_flags.json')
    Zm=sp.Matrix(boundary['torus_double_residues_Z']);R=sp.eye(30)[boundary['anchor_rows_zero_based'],:]
    edges=boundary['oriented_edges'];idx=edges.index(['E12','F3'])
    tau=sp.eye(30)[idx,:]*(sp.eye(30)-Zm*R)
    c=sp.Matrix(boundary['mass_flag_matrix'])*sp.Matrix([3,2,1,0])
    assert tau*Zm==sp.zeros(1,3) and (tau*c)[0]==1
    flag=(tau*sp.Matrix(boundary['mass_flag_matrix'])*sp.Matrix([0,2,4,6]))[0]
    assert flag==-2
    out={'status':'actual_rank_eight_relative_comparison_enclosed_pending_independent_replay',
         'code_sha256':sha(__file__),'lower_comparison_sha256':sha(lowpath),
         'top_period_certificate_sha256':{p.name:sha(p) for p in [*(OUT/f'node_{a}.json' for a in [1,9,25,49,81]),OUT/'sixth.json']},
         'bits':BITS,'basepoint':['-64','0'],'masses':[1,2,4,8],'matrix_center':encode(P8),
         'matrix_infinity_error_upper':str(err),'top_row_errors':list(map(str,errs)),
         'log2_interval':{'series':'2*atanh(1/3)','terms':n,'lower':str(l),'upper':str(u)},
         'boundary_top_period':record([q7],q7e),'top_Tate_period':record([t3],t3e),
         'corrected_boundary_torus_coefficients':list(map(str,tau)),
         'boundary_top_log2_coefficient':str(flag),'boundary_data_sha256':sha(ROOT/'data/boundary_kummer_flags.json'),
         'rows':['omega','v1','v2','v3','v4','w_published','beta_c','eta=dlog(x)dlog(y)dlog(z)'],
         'columns':['Gamma1','Gamma2','Gamma3','Gamma4','Gamma5','T0_Gamma1','Gamma_boundary_plus','T3'],
         'scope':'untwisted E_act(-3); geometric rational relative basis, not an integral saturated basis',
         'physical_Betti_marking':'intersection with positive simplex: (0,0,0,0,0,0,0,1); see RELATIVE_TOP_COMPARISON.md'}
    save_json(OUT/'period_matrix.json',out)
    print(json.dumps({'P8_error':float(err),'dyadic_bound':str(dyadic_upper(err)),'q7':float(q7[0])}),flush=True)
def main():
    OUT.mkdir(exist_ok=True);tr=trunks();ns=[node(a,tr) for a in [1,9,25,49,81]]
    six=sixth(tr);assemble(ns,six)
if __name__=='__main__':main()

