"""Morse residue normalization of the six-cycle period matrix."""
from node_frobenius import *
from refined_transport import compose_fine
from math import isqrt


def scale_matrix(C,e,a,ae,bits):
    S=1<<bits;exact=[gm(a,gs(F(1,S),v)) for v in C]
    new=[(dyadic_round(x,bits),dyadic_round(y,bits)) for x,y in exact]
    rd=mnorm([ga(v,gn(gs(F(1,S),w))) for v,w in zip(exact,new)])
    return new,radius_up((ab(a)+ae)*e+ae*F(mnorm(C),S)+rd,bits),rd


def main():
    dest=FOLDER/'marked_lattice';bp=dest/'marked_zero.json';b=json.loads(bp.read_text())
    assert json.loads(bp.with_name('marked_zero_replay.json').read_text())['certificate_sha256']==sha(bp)
    spath=FOLDER/'node_frobenius/node_1.json';seed=json.loads(spath.read_text())
    Q=sp.Matrix(seed['Q']);p=[Q*sp.Matrix([sp.Rational(x)+sp.I*sp.Rational(y) for x,y in v]) for v in seed['frobenius_vector_coefficients_exact'][:2]]
    assert p[0][0]==0 and p[1][0]==sp.I/128
    x=sp.symbols('x1:4');r=[-1,-2,-4,8];Ssum=sum(r)
    f=(sum(r[i]*x[i] for i in range(3))+r[3])*(sum(r[i]/x[i] for i in range(3))+r[3])
    at=dict.fromkeys(x,1);H=sp.hessian(f,x).subs(at)
    assert f.subs(at)==1 and all(sp.diff(f,v).subs(at)==0 for v in x) and H.det()==-512
    data=json.loads((FOLDER/'polarized_system.json').read_text());ss=sp.Symbol('ss')
    assert data['gauge_P'][0]==['1','0','0','0','0','0']
    D=sp.diag(*map(sp.Rational,D_VALUES));K=D.inv()*sp.Matrix(data['inverse_pairing'])*D.inv();K0=K.subs(ss,1)
    leading_pairing=(p[0].T*K0*p[0])[0];assert leading_pairing==0
    norm=sp.simplify((2*p[1].T*K0*p[0]+sp.I/4*p[0].T*K.diff(ss).subs(ss,1)*p[0])[0])
    assert norm==sp.I/65536
    rho2=sp.simplify((-64)**2*4*(sp.I/4)/(p[1][0]**2*Ssum**2*sp.prod(r)))
    assert rho2==1048576*sp.I and sp.simplify(-rho2*norm/2)==8
    zp=FOLDER/'zero_monodromy/transport.json';z=json.loads(zp.read_text());pi=z['pi_certificate']
    pc=F(pi['center']);pe=F(pi['radius']);bits=b['bits'];scale=1<<bits
    k=isqrt(2*scale*scale);sl=F(k,scale);su=F(k+1,scale);assert sl*sl<2<su*su
    sc=(sl+su)/2;se=(su-sl)/2
    ac=512*pc*sc;ae=1024*((abs(pc)+pe)*se+pe*abs(sc)) # usual complex absolute error <= l1 bound
    rho=(ac,ac)
    P=decode(b['basis_center_gaussian_integer']);e=F(b['basis_error_upper'])
    high,he,hrd=scale_matrix(P,e,rho,ae,bits)
    ip=FOLDER/'high_trunk_inverse.json';inv=json.loads(ip.read_text())
    V=[gs(1<<(bits-inv['bits']),v) for v in decode(inv['inverse_center_gaussian_integer'])]
    ie=F(inv['inverse_error_upper'])
    baseD,be,brd=compose_fine(V,ie,high,he,bits)
    Di=sp.diag(*[1/sp.Rational(v) for v in D_VALUES])
    # Exact left diagonal change to the published-pairing polarized B basis.
    exact=[gs(rat(Di[i,i])/scale,baseD[6*i+j]) for i in range(6) for j in range(6)]
    base=[(dyadic_round(x,bits),dyadic_round(y,bits)) for x,y in exact]
    rd=mnorm([ga(v,gn(gs(F(1,scale),w))) for v,w in zip(exact,base)])
    base_error=radius_up(mnorm(real_matrix(Di))*be+rd,bits)
    out={'status':'six_geometric_cycle_absolute_period_matrix_certified_in_explicit_residue_convention',
         'code_sha256':sha(__file__),'basis_certificate_sha256':sha(bp),'node_seed_sha256':sha(spath),
         'system_sha256':sha(FOLDER/'polarized_system.json'),'zero_transport_sha256':sha(zp),
         'reverse_trunk_sha256':sha(ip),'bits':bits,'morse_signed_masses':r,'morse_hessian':[[str(v) for v in row] for row in H.tolist()],
         'first_Frobenius_coefficient':str(p[1][0]),'local_vector_pairing':str(norm),
         'rho_squared_over_pi_squared':str(rho2),'rho_formula':'512*sqrt(2)*pi*(1+i)',
         'residue_convention':'Res(dlog(x1)^dlog(x2)^dlog(x3)/(f-s)); no factors of 2*pi*i in the 2-form',
         'orientation':'node1 leading period +2*pi*sqrt(s-1)/(8*i), sqrt(i/4)=exp(i*pi/4)/2',
         'sqrt2_interval':[str(sl),str(su)],'rho_center':encode([rho])[0],'rho_error_upper':str(ae),
         'high_basepoint':['-64','64'],'high_D_center_gaussian_integer':encode(high),'high_D_error_upper':str(he),
         'scaling_rounding_error':str(hrd),'physical_basepoint':['-64','0'],
         'physical_D_center_gaussian_integer':encode(baseD),'physical_D_error_upper':str(be),'rebasing_rounding_error':str(brd),
         'physical_B_center_gaussian_integer':encode(base),'physical_B_error_upper':str(base_error),'unbalancing_rounding_error':str(rd),
         'period_pairing_identity':'P_D^t K_D P_D = 8*pi^2 Q6',
         'algebraic_trace_relative_to_published_pairing':-2,
         'full_integral_Betti_basis_certified':False,
         'basis_scope':'Six geometrically marked integral cycles, index two in T, hence a rational Betti basis.',
         'old_scalar_theta_basis_compared':False,'fixed_regulator_nonzero_certified':False}
    save_json(dest/'absolute_periods.json',out)
    print(json.dumps({'status':out['status'],'physical_B_error_dyadic':str(dyadic_upper(base_error)),
                      'high_D_error_dyadic':str(dyadic_upper(he)),'rho':out['rho_formula'],'pairing':out['period_pairing_identity']}))


if __name__=='__main__':main()
