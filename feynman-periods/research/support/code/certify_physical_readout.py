"""Certified coordinates of nu=e_F-sigma and physical readout in the balanced frame."""
from certify_lower_comparison import *
from certify_horizontal_landmarks import inverse
OUT=FOLDER/'relative_top'
def main():
    ap=FOLDER/'marked_lattice/absolute_periods.json';a=read(ap)
    pp=OUT/'period_matrix.json';p=read(pp);lo=read(FOLDER/'lower_comparison/period_matrix.json')
    D=[F(-64)]*4+[F(4096),F(-1,64)]
    PD=decode(a['physical_D_center_gaussian_integer'])
    PT=[PD[6*j+i] for i in range(6) for j in range(6)]
    pe=6*F(a['physical_D_error_upper']);iv=inverse(PT,pe,BITS)
    V=fixed(iv,'center_gaussian_integer',BITS);ve=F(iv['error_upper'])
    lc=(F(p['log2_interval']['lower'])+F(p['log2_interval']['upper']))/2
    le=(F(p['log2_interval']['upper'])-F(p['log2_interval']['lower']))/2
    k=-2*lc;ke=2*le
    P8=vfdecode(p['matrix_center']);P7=vfdecode(lo['matrix_center'])
    r=[];ers=[]
    for j in range(6):
        b=P7[42+j];q=P8[56+j];qe=F(p['top_row_errors'][j]);be=F(lo['beta_entry_error_upper'][j])
        raw=ga(q,gn(gs(k,b)));v=roundv([raw])[0]
        r.append(v);ers.append(qe+(abs(k)+ke)*be+ke*ab(b)+ab(ga(raw,gn(v))))
    re=max(ers)
    v,veout=action(V,ve,r,re)
    nuB=[gs(D[i],v[i]) for i in range(6)]+[(k,F(0))]
    errors=[abs(x)*veout for x in D]+[ke]
    # Normalized relative duality gives <Omega,delta w>=+2.
    I=gs(-2,nuB[5]);Ie=2*errors[5]
    assert abs(I[1])<=Ie and I[0]-Ie>0
    interval=[F(73554702266148,10**14),F(73554702266149,10**14)]
    # A coarser display interval can be used if the uniform inverse estimate is wider.
    if not(interval[0]<I[0]-Ie and I[0]+Ie<interval[1]):
        interval=[F(73554702266,10**11),F(73554702267,10**11)]
    assert interval[0]<I[0]-Ie and I[0]+Ie<interval[1]
    out={'status':'nu_coordinates_and_physical_readout_enclosed_pending_independent_replay',
      'code_sha256':sha(__file__),'top_matrix_sha256':sha(pp),'pure_absolute_sha256':sha(ap),
      'lower_matrix_sha256':sha(FOLDER/'lower_comparison/period_matrix.json'),
      'transpose_pure_error_upper':str(pe),'transpose_pure_inverse':iv,
      'rhs':record(r,re),'nu_balanced_D':record(v,veout),'nu_B_center':encode(nuB),
      'nu_B_entry_error_upper':list(map(str,errors)),'physical_integral':record([I],Ie),
      'physical_integral_strict_real_interval':list(map(str,interval)),
      'normalization':'I=integral_positive_simplex Omega/F; nu=e_F-sigma; I=-2*nu_B[w_published]',
      'physical_chain_coordinates_in_relative_dual_basis':[0]*7+[1],
      'geometric_sign_proof':'RELATIVE_TOP_COMPARISON.md: ordered-intersection and annulus calculation',
      'float_proposal_used_only_for_inverse':True}
    save_json(OUT/'physical_readout.json',out)
    print(json.dumps({'I':[float(x) for x in I],'error':float(Ie),'interval':out['physical_integral_strict_real_interval']}),flush=True)
if __name__=='__main__':main()

