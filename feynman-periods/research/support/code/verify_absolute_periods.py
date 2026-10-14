"""Independent exact verification of the Morse calibration and period enclosure."""
from node_frobenius import *
from verify_zero_transport import check_pi,check_product


def main():
    dest=FOLDER/'marked_lattice';p=dest/'absolute_periods.json';d=json.loads(p.read_text())
    bp=dest/'marked_zero.json';b=json.loads(bp.read_text());seedp=FOLDER/'node_frobenius/node_1.json';seed=json.loads(seedp.read_text())
    zp=FOLDER/'zero_monodromy/transport.json';z=json.loads(zp.read_text());ip=FOLDER/'high_trunk_inverse.json';inv=json.loads(ip.read_text())
    assert d['code_sha256']==sha(ROOT/'code/certify_absolute_periods.py')
    for key,path in [('basis_certificate_sha256',bp),('node_seed_sha256',seedp),('zero_transport_sha256',zp),('reverse_trunk_sha256',ip),('system_sha256',FOLDER/'polarized_system.json')]:
        assert d[key]==sha(path)
    assert json.loads(bp.with_name('marked_zero_replay.json').read_text())['certificate_sha256']==sha(bp)
    bits=d['bits'];S=1<<bits
    # Direct analytic Hessian formula, independent of symbolic differentiation.
    r=list(map(sp.Integer,d['morse_signed_masses']));assert r==[-1,-2,-4,8] and sum(r)==1
    H=2*(sum(r)*sp.diag(*r[:3])-sp.Matrix(r[:3])*sp.Matrix(r[:3]).T)
    assert H==sp.Matrix(d['morse_hessian']) and H.det()==-512
    Q=sp.Matrix(seed['Q']);v0=Q*sp.Matrix([sp.Rational(x)+sp.I*sp.Rational(y) for x,y in seed['frobenius_vector_coefficients_exact'][0]])
    v1=Q*sp.Matrix([sp.Rational(x)+sp.I*sp.Rational(y) for x,y in seed['frobenius_vector_coefficients_exact'][1]])
    assert v0[0]==0 and v1[0]==sp.I/128==sp.sympify(d['first_Frobenius_coefficient'])
    sys=json.loads((FOLDER/'polarized_system.json').read_text());ss=sp.Symbol('ss')
    assert sp.Matrix(sys['gauge_P'])[0,:]==sp.eye(6)[0,:]
    Di=sp.diag(*[1/sp.Rational(x) for x in D_VALUES]);K=Di*sp.Matrix(sys['inverse_pairing'])*Di
    assert (v0.T*K.subs(ss,1)*v0)[0]==0
    norm=sp.simplify(sum(v1[i]*K[i,j].subs(ss,1)*v0[j]+v0[i]*K[i,j].subs(ss,1)*v1[j]+sp.I/4*v0[i]*sp.diff(K[i,j],ss).subs(ss,1)*v0[j] for i in range(6) for j in range(6)))
    assert norm==sp.sympify(d['local_vector_pairing'])==sp.I/65536
    # Sphere residue integral is 2*pi*sqrt(t), with squared Morse Jacobian 8/det H.
    rho2=sp.simplify((-64)**2*4*(sp.I/4)*8/(H.det()*v1[0]**2))
    assert rho2==sp.sympify(d['rho_squared_over_pi_squared'])==1048576*sp.I
    assert sp.simplify(-rho2*norm/2)==8
    pc,pe=check_pi(z['pi_certificate']);sl,su=map(F,d['sqrt2_interval'])
    assert 0<sl<su and sl*sl<2<su*su and su-sl==F(1,S)
    sc=(sl+su)/2;se=(su-sl)/2;a=(512*pc*sc,512*pc*sc)
    ae=1024*((abs(pc)+pe)*se+pe*sc)
    assert encode([a])[0]==d['rho_center'] and ae==F(d['rho_error_upper'])
    P=decode(b['basis_center_gaussian_integer']);e=F(b['basis_error_upper'])
    high=decode(d['high_D_center_gaussian_integer']);he=F(d['high_D_error_upper'])
    exact=[gm(a,gs(F(1,S),v)) for v in P]
    rd=max(sum(ab(ga(exact[6*i+j],gn(gs(F(1,S),high[6*i+j])))) for j in range(6)) for i in range(6))
    raw=(ab(a)+ae)*e+ae*F(mnorm(P),S)+rd
    assert raw<=he<raw+F(1,S) and rd==F(d['scaling_rounding_error'])
    # Recheck the archived inverse bound, allowing its coarser power-of-two rounding.
    trunkp=FOLDER/'high_trunk_transport.json';tr=json.loads(trunkp.read_text())
    assert inv['transport_sha256']==sha(trunkp)
    ib=inv['bits'];IS=1<<ib;C=decode(tr['endpoint_center_gaussian_integer']);V0=decode(inv['inverse_center_gaussian_integer']);ce=F(tr['error_upper'])
    R=[ga((IS*IS if i==j else 0,0),gn(sumg(gm(V0[6*i+k],C[6*k+j]) for k in range(6)))) for i in range(6) for j in range(6)]
    rn=F(mnorm(R),IS*IS);vn=F(mnorm(V0),IS);eta=rn+vn*ce
    assert str(rn)==inv['residual_upper'] and str(eta)==inv['eta_upper'] and eta<1
    ie=F(inv['inverse_error_upper']);assert eta*vn/(1-eta)<=ie
    V=[gs(1<<(bits-ib),v) for v in V0]
    baseD=decode(d['physical_D_center_gaussian_integer']);be=F(d['physical_D_error_upper'])
    assert check_product(V,ie,high,he,baseD,be,bits)==F(d['rebasing_rounding_error'])
    base=decode(d['physical_B_center_gaussian_integer']);ee=F(d['physical_B_error_upper'])
    exact=[gs(rat(Di[i,i])/S,baseD[6*i+j]) for i in range(6) for j in range(6)]
    rd=mnorm([ga(v,gn(gs(F(1,S),w))) for v,w in zip(exact,base)])
    raw=mnorm(real_matrix(Di))*be+rd
    assert raw<=ee<raw+F(1,S) and rd==F(d['unbalancing_rounding_error'])
    out={'status':'Morse_calibration_and_absolute_six_cycle_period_enclosure_independently_replayed',
         'certificate_sha256':sha(p),'verifier_sha256':sha(__file__),
         'physical_B_error_dyadic':str(dyadic_upper(ee)),'high_D_error_dyadic':str(dyadic_upper(he)),
         'exact_local_norm':'I/65536','rho_squared_over_pi_squared':'1048576*I',
         'period_pairing_constant_over_pi_squared':8,'floating_point_used':False,
         'geometric_input':'Morse residue sphere integral and identification of first form as the specified bare Poincare residue.',
         'full_integral_Betti_basis_certified':False,'old_scalar_theta_basis_compared':False}
    save_json(p.with_name('absolute_periods_replay.json'),out);print(json.dumps(out))


if __name__=='__main__':main()
