"""Complete the seventh-row comparison using archived zero transport only."""
from certify_lower_regulator import *
OUT=FOLDER/'lower_comparison'

def read(p):return json.loads(p.read_text())
def unpack(v):return vfdecode(v['center']),F(v['error_upper'])
def fixed(d,key,bits=None):
    return [gs(F(1,1<<(d['bits'] if bits is None else bits)),v) for v in decode(d[key])]
def action(M,me,v,ve):
    raw=mv(M,v);out=roundv(raw)
    e=(mnorm(M)+me)*ve+me*vnorm(v)+vnorm([ga(a,gn(b)) for a,b in zip(raw,out)])
    return out,radius_up(e,BITS)
def main():
    OUT.mkdir(exist_ok=True)
    seedp=FOLDER/'zero_monodromy/frobenius.json';zp=FOLDER/'zero_monodromy/transport.json'
    trp=DEST/'trunk_integrals.json';n1p=DEST/'node_1.json';ap=FOLDER/'marked_lattice/absolute_periods.json'
    seed=read(seedp);zt=read(zp);tr=read(trp);n1=read(n1p);absolute=read(ap)
    s=sp.Symbol('ss');toD=sp.Matrix(seed['to_D'])
    gqsp=(sp.Matrix([[sp.Rational(x) for x,y in G]])*toD).applyfunc(sp.cancel)
    assert all(s not in x.free_symbols for x in gqsp)
    gq=[(rat(x),F(0)) for x in gqsp];N=real_matrix(sp.Matrix(seed['nilpotent']))
    assert mm(N,N)==mz()
    coeff=[vfdecode(A) for A in seed['gauge_coefficients_exact']];h=(F(0),F(1,4))
    exact=[Z]*6
    for n,C in enumerate(coeff):
        r=rowmat(gq,C);rN=rowmat(r,N)
        exact=va(exact,[gm(h,ga(gs(F(1,n+1),x),gn(gs(F(1,(n+1)**2),y)))) for x,y in zip(r,rN)])
    J=roundv(exact)
    raw=ab(h)*rn(gq)*(1+mnorm(N))*F(seed['bounds']['gauge_error_upper'])
    je=radius_up(raw+rn([ga(x,gn(y)) for x,y in zip(exact,J)]),BITS)
    local=record(J,je);logs=[]
    for entry in zt['segments']:
        p=zp.parent/entry['file'];d=read(p);assert sha(p)==entry['sha256']
        R,re=integrated_row(d,gq if d['frame']=='zero_Q' else G)
        C=fixed(d,'input_center_gaussian_integer');ce=F(d['input_error_upper'])
        J,je=append_row(J,je,R,re,C,ce)
        logs.append({'file':entry['file'],'sha256':sha(p),'primitive_row':record(J,je)})
    at_top=record(J,je)
    W=fixed(zt,'based_frame_center_gaussian_integer');we=F(zt['based_frame_error_upper'])
    hor,he=unpack(tr['landmarks']['0'])
    J,je=append_row(J,je,[gn(x) for x in hor],he,W,we)
    at_high=record(J,je)
    bd=read(FOLDER/'node_branches/node_1/transport.json')
    v1=fixed(bd,'based_center_gaussian_integer');e1=F(bd['based_error_upper'])
    v6=fixed(zt,'sixth_vector_center_gaussian_integer');e6=F(zt['sixth_vector_error_upper'])
    dv=va(v6,[gn(v) for v in v1]);de=e6+e1
    Vi=fixed(zt['based_frame_inverse'],'center_gaussian_integer',zt['bits'])
    vierr=F(zt['based_frame_inverse']['error_upper'])
    w,werr=action(Vi,vierr,dv,de)
    bc,bce=unpack(n1['beta_at_high_basepoint'])
    b,be=scalar_update(bc[0],bce,J,je,w,werr)
    atc=record([b],be)
    invp=FOLDER/'high_trunk_inverse.json';inv=read(invp)
    V=fixed(inv,'inverse_center_gaussian_integer');Ve=F(inv['inverse_error_upper'])
    vb,vbe=action(V,Ve,v6,e6)
    jv,jve=unpack(tr['vertical'])
    b,be=scalar_update(b,be,jv,jve,vb,vbe,-1)
    unscaled=record([b],be)
    rho=vfdecode([absolute['rho_center']])[0];rhoe=F(absolute['rho_error_upper'])
    raw=gm(rho,b);bv=roundv([raw])[0]
    bve=radius_up((ab(rho)+rhoe)*be+rhoe*ab(b)+ab(ga(raw,gn(bv))),BITS)
    increment=dot(J,w);increrr=rn(J)*werr+je*(vnorm(w)+werr)
    out={'status':'sixth_actual_beta_period_certified','bits':BITS,'code_sha256':sha(__file__),
         'input_sha256':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in [seedp,zp,trp,n1p,ap,invp]},
         'coupling_zero_Q':encode(gq),'local_primitive':local,'steps':logs,'primitive_at_top':at_top,
         'primitive_at_high_basepoint':at_high,'delta_pure_at_high':record(dv,de),
         'delta_in_zero_constants':record(w,werr),'beta_monodromy_increment_normalized':record([increment],increrr),
         'beta_at_high_basepoint':atc,'pure_at_physical_basepoint':record(vb,vbe),
         'unscaled_beta_at_physical_basepoint':unscaled,'beta_period':record([bv],bve),
         'cycle':'delta6=T0 delta1 in H2(U); includes beta monodromy increment'}
    save_json(OUT/'sixth_period.json',out)
    print(json.dumps({'beta6':[float(x) for x in bv],'error':float(bve),'increment_normalized':[float(x) for x in increment]}),flush=True)
    P6=fixed(absolute,'physical_B_center_gaussian_integer');p6e=F(absolute['physical_B_error_upper'])
    bottom=[];berrs=[];inputs={}
    for a in [1,9,25,49,81]:
        p=DEST/f'node_{a}.json';d=read(p);v,e=unpack(d['beta_period']);bottom.append(v[0]);berrs.append(e);inputs[str(p.relative_to(ROOT)).replace('\\','/')]=sha(p)
    bottom.append(bv);berrs.append(bve)
    pc=F(zt['pi_certificate']['center']);pe=F(zt['pi_certificate']['radius'])
    rawt=(-4*pc*pc,F(0));tate=roundv([rawt])[0]
    te=radius_up(4*(2*abs(pc)*pe+pe**2)+ab(ga(rawt,gn(tate))),BITS)
    P7=[]
    for i in range(6):P7.extend(P6[6*i:6*i+6]+[Z])
    P7+=bottom+[tate];err=max(p6e,sum(berrs)+te)
    flags=read(ROOT/'data/minimality_cycle.json')['flags']
    flag=next(x for x in flags if x['edge']==['E12','F3'])
    assert flag['cycle_coefficient']==1
    inputs.update(out['input_sha256'])
    inputs['data/minimality_cycle.json']=sha(ROOT/'data/minimality_cycle.json')
    inputs['data/banana_pure_ray/lower_comparison/sixth_period.json']=sha(OUT/'sixth_period.json')
    final={'status':'seven_dimensional_lower_comparison_certified','bits':BITS,'code_sha256':sha(__file__),
           'input_sha256':inputs,'basepoint':['-64','0'],'masses':[1,2,4,8],
           'rows':['omega','v1','v2','v3','v4','w_published','beta_c'],
           'columns':['delta1','delta2','delta3','delta4','delta5','T0_delta1','tau_E12_F3'],
           'matrix_center':encode(P7),'matrix_infinity_error_upper':str(err),
           'beta_entry_error_upper':list(map(str,berrs)),'boundary_tate':record([tate],te),
           'boundary_flag':flag,'pure_block_certificate_sha256':sha(ap),
           'basis_scope':'rational Betti basis of V_c dual, represented by actual integral cycles; not a saturated integral basis',
           'Tate_twist':'untwisted V_c; physical lower N=V_c(3)',
           'full_eight_dimensional_comparison_certified':False}
    save_json(OUT/'period_matrix.json',final)
    print(json.dumps({'P7_error':float(err),'dyadic_upper':str(dyadic_upper(err)),'tate':float(tate[0])}),flush=True)
if __name__=='__main__':main()

