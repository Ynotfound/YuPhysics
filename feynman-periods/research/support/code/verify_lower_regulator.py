"""Independent scalar replay of the lower regulator obstruction.

Does not import the generator or call any ODE/Taylor/Frobenius generator.
It checks the new quadratures against hash-linked archived pure certificates.
"""
from node_frobenius import *
from itertools import permutations
from copy import deepcopy

DEST=FOLDER/'lower_regulator'
GG=[(F(0),F(0)),(F(-3,64),F(0)),(F(-2,64),F(0)),(F(-1,64),F(0)),(F(0),F(0)),(F(0),F(0))]

def read(p):return json.loads(p.read_text())
def rnorm(v):return sum(abs(x)+abs(y) for x,y in v)
def vals(d):return vfdecode(d['center']),F(d['error_upper'])
def check_radius(raw,stored):
    assert raw<=stored<raw+F(1,1<<320)
def scalar_sum(A,B):
    x=y=F(0)
    for (a,b),(c,d) in zip(A,B):x+=a*c-b*d;y+=a*d+b*c
    return x,y
def distance(a,b):return abs(a[0]-b[0])+abs(a[1]-b[1])
def calc_integral(d,g,first=False):
    h=(F(d['end'][0])-F(d['start'][0]),F(d['end'][1])-F(d['start'][1]))
    poly=d['witness_coefficients_ascending'] if first else d['polynomial_coefficients_gaussian_integer']
    out=[]
    for j in range(6):
        value=(F(0),F(0))
        for n,M in enumerate(poly):
            coeff=[]
            for i in range(6):
                v=(F(M[i][j]),F(0)) if first else (F(M[6*i+j][0])/(1<<d['bits']),F(M[6*i+j][1])/(1<<d['bits']))
                if first:
                    for _ in range(n):v=(-v[1],v[0])
                coeff.append(v)
            term=scalar_sum(g,coeff)
            value=(value[0]+term[0]/(n+1),value[1]+term[1]/(n+1))
        out.append((h[0]*value[0]-h[1]*value[1],h[0]*value[1]+h[1]*value[0]))
    eps=(abs(h[0])+abs(h[1]))*rnorm(g)*F(d['bounds']['error_norm_upper'])
    # Generator saves 320-bit intermediate rows. Reconstruct independently.
    rounded=[]
    for a,b in out:
        pair=[]
        for x in [a,b]:
            z=x*(1<<320);pair.append(F((2*z.numerator+z.denominator)//(2*z.denominator),1<<320))
        rounded.append(tuple(pair))
    raw=eps+sum(distance(a,b) for a,b in zip(out,rounded))
    eps=F(-(-raw*(1<<320)//1),1<<320)
    return rounded,eps

def check_scalar(old,olde,row,re,v,ve,out,oute,sign=1):
    t=scalar_sum(row,v);exact=(old[0]+sign*t[0],old[1]+sign*t[1])
    raw=olde+rnorm(row)*ve+re*(max(abs(a)+abs(b) for a,b in v)+ve)+distance(exact,out)
    check_radius(raw,oute)

def check_row(old,olde,row,re,C,ce,out,oute):
    error=F(0)
    for j in range(6):
        t=scalar_sum(row,[C[6*i+j] for i in range(6)])
        error+=distance((old[j][0]+t[0],old[j][1]+t[1]),out[j])
    cn=max(sum(abs(C[6*i+j][0])+abs(C[6*i+j][1]) for j in range(6)) for i in range(6))
    check_radius(olde+rnorm(row)*ce+re*(cn+ce)+error,oute)

def main():
    trp=DEST/'trunk_integrals.json';tr=read(trp);obsp=DEST/'obstruction.json';obs=read(obsp)
    assert tr['code_sha256']==obs['code_sha256']==sha(ROOT/'code/certify_lower_regulator.py')
    assert vfdecode(tr['coupling_D'])==GG
    sys=read(FOLDER/'polarized_system.json');ss=sp.Symbol('ss')
    Pi=sp.Matrix(sys['inverse_gauge']);Di=sp.diag(*[1/sp.Rational(v) for v in D_VALUES])
    original=sp.Matrix(read(FOLDER/'active_lower_system.json')['matrix'])[6,:6]
    assert (original*Pi*Di).applyfunc(sp.cancel)==sp.Matrix([[sp.Rational(a) for a,b in GG]])
    for data,replay,key in [('first_trunk_step.json','first_trunk_step_replay.json','certificate_sha256'),
                            ('high_trunk_transport.json','high_trunk_replay.json','certificate_sha256'),
                            ('horizontal_trunk/transport.json','horizontal_trunk/replay.json','transport_sha256'),
                            ('marked_lattice/absolute_periods.json','marked_lattice/absolute_periods_replay.json','certificate_sha256')]:
        assert read(FOLDER/replay)[key]==sha(FOLDER/data)
    archival={}
    ht=read(FOLDER/'high_trunk_transport.json')
    hr=read(FOLDER/'high_trunk_replay.json')
    for e in ht['segments']:archival[e['file']]=e['sha256']
    for e in read(FOLDER/'horizontal_trunk/transport.json')['segments']:archival['horizontal_trunk/'+e['file']]=e['sha256']
    J=[Z]*6;je=F(0);C=mi();ce=F(0)
    for n,e in enumerate(tr['vertical_steps']):
        p=FOLDER/e['file'];d=read(p);assert e['sha256']==sha(p)
        if n:assert archival[e['file'].replace('\\','/')]==sha(p)
        R,re=calc_integral(d,GG,n==0);out,oe=vals(e['integral']);check_row(J,je,R,re,C,ce,out,oe);J,je=out,oe
        if n==0:
            C=[(F(a),F(b)) for rr,ii in zip(d['endpoint_center_real'],d['endpoint_center_imag']) for a,b in zip(rr,ii)];ce=F(d['bounds']['error_norm_upper'])
        else:C=[gs(F(1,1<<d['bits']),v) for v in decode(d['propagated_center_gaussian_integer'])];ce=F(d['propagated_error_upper'])
    assert (J,je)==vals(tr['vertical'])
    J=[Z]*6;je=F(0);C=mi();ce=F(0)
    for e in tr['horizontal_steps']:
        p=FOLDER/e['file'];d=read(p);assert e['sha256']==archival[e['file'].replace('\\','/')]==sha(p)
        R,re=calc_integral(d,GG);out,oe=vals(e['integral']);check_row(J,je,R,re,C,ce,out,oe);J,je=out,oe
        C=[gs(F(1,1<<d['bits']),v) for v in decode(d['high_center_gaussian_integer'])];ce=F(d['high_error_upper'])
        if d['landmark_pole'] is not None:assert (J,je)==vals(tr['landmarks'][str(d['landmark_pole'])])
    absd=read(FOLDER/'marked_lattice/absolute_periods.json');rho0=vfdecode([absd['rho_center']])[0];rho0e=F(absd['rho_error_upper'])
    inv=read(FOLDER/'high_trunk_inverse.json');V=[gs(F(1,1<<inv['bits']),v) for v in decode(inv['inverse_center_gaussian_integer'])];ve=F(inv['inverse_error_upper'])
    rsdata={1:[-1,-2,-4,8],9:[1,-2,-4,8],25:[-1,2,-4,8],49:[1,2,-4,8],81:[-1,-2,4,8]}
    records={};tail_tamper_rejected=False
    for a in [1,9,25,49,81]:
        dp=DEST/f'node_{a}.json';d=read(dp);spth=FOLDER/f'node_frobenius/node_{a}.json';seed=read(spth)
        bp=FOLDER/f'node_branches/node_{a}/transport.json';bd=read(bp)
        for k,p in [('seed_sha256',spth),('branch_sha256',bp),('trunk_integrals_sha256',trp),('inverse_sha256',FOLDER/'high_trunk_inverse.json'),('absolute_scale_sha256',FOLDER/'marked_lattice/absolute_periods.json')]:assert d[k]==sha(p)
        assert read(bp.with_name('replay.json'))['transport_sha256']==sha(bp)
        assert obs['node_certificate_sha256'][str(a)]==sha(dp)
        Q=real_matrix(sp.Matrix(seed['Q']));gq=[scalar_sum(GG,[Q[6*i+j] for i in range(6)]) for j in range(6)]
        coeff=[vfdecode(v) for v in seed['frobenius_vector_coefficients_exact']]
        rr=ii=F(0)
        for n,v in enumerate(coeff):
            x,y=scalar_sum(gq,v);rr-=y/(4*(F(n)+F(1,2)));ii+=x/(4*(F(n)+F(1,2)))
        bvec,be=vals(d['local_beta']);b=bvec[0]
        raw=F(1,2)*rnorm(gq)*F(seed['bounds']['vector_error_Q_basis'])+distance((rr,ii),b);check_radius(raw,be)
        if a==1:
            try:check_radius(raw,F(0))
            except AssertionError:tail_tamper_rejected=True
        first=scalar_sum(Q[:6],coeff[1]);beta_lead=scalar_sum(gq,coeff[0]);beta_lead=(-beta_lead[1]/2,beta_lead[0]/2)
        om_lead=gs(F(-1,64),first);ratio=gd(beta_lead,om_lead)
        rs=rsdata[a];S=sum(rs);assert ratio==(-S*sum(x*y for x,y in zip(rs,[3,2,1,0])),0)
        assert encode([ratio])[0]==d['local_beta_omega_leading_ratio']
        for ee,src in zip(d['steps'],bd['segments']):
            pp=bp.parent/ee['file'];st=read(pp);assert ee['sha256']==src['sha256']==sha(pp)
            R,re=calc_integral(st,gq if st['frame']=='Q' else GG)
            v=[gs(F(1,1<<st['bits']),z) for z in decode(st['input_center_gaussian_integer'])]
            out,oe=vals(ee['beta']);check_scalar(b,be,R,re,v,F(st['input_error_upper']),out[0],oe);b,be=out[0],oe
        assert len(d['steps'])==64
        cv=[gs(F(1,1<<bd['bits']),v) for v in decode(bd['based_center_gaussian_integer'])];ce=F(bd['based_error_upper'])
        R,re=vals(tr['landmarks'][str(a)]);out,oe=vals(d['beta_at_high_basepoint']);check_scalar(b,be,R,re,cv,ce,out[0],oe,-1);b,be=out[0],oe
        bv,bve=vals(d['pure_at_physical_basepoint']);rd=F(0)
        for i in range(6):rd=max(rd,distance(scalar_sum(V[6*i:6*i+6],cv),bv[i]))
        check_radius((mnorm(V)+ve)*ce+ve*vnorm(cv)+rd,bve)
        R,re=vals(tr['vertical']);out,oe=vals(d['unscaled_beta_at_physical_basepoint']);check_scalar(b,be,R,re,bv,bve,out[0],oe,-1);b,be=out[0],oe
        # Local Morse leading coefficient, divided by the saved first pure coefficient.
        prod=1
        for x in rs:prod*=x
        assert abs(prod)==64
        jacfactor=(F(1),F(0)) if prod>0 else (F(0),F(-1))
        factor=gd(gs(F(-1,128*S),jacfactor),first)
        rhov,rhoe=vals(d['rho']);rho=rhov[0];assert rho==gm(factor,rho0) and rhoe==ab(factor)*rho0e
        for vv,ev,key in [(gs(F(-1,64),bv[0]),bve/64,'omega_period'),(b,be,'beta_period')]:
            out,oe=vals(d[key]);exact=gm(rho,vv)
            check_radius((ab(rho)+rhoe)*ev+rhoe*ab(vv)+distance(exact,out[0]),oe)
        records[a]=d
    assert tail_tamper_rejected
    for minor in obs['minors']:
        matrix=[];radii=[]
        for a in minor['nodes']:
            d=records[a];p,pe=vals(d['omega_period']);b,be=vals(d['beta_period']);matrix.append([p[0][0],p[0][1],b[0][1]]);radii.append([pe,pe,be])
        assert [[str(x) for x in r] for r in matrix]==minor['rows'];assert [[str(x) for x in r] for r in radii]==minor['radii']
        a,b,c=matrix;det=a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
        error=F(0)
        # Expand the seven nonempty error subsets of every determinant monomial.
        for perm in permutations(range(3)):
            for mask in range(1,8):
                term=F(1)
                for i,j in enumerate(perm):term*=radii[i][j] if (mask>>i)&1 else abs(matrix[i][j])
                error+=term
        assert det==F(minor['det_center']) and error==F(minor['det_error_upper'])
        assert minor['excludes_zero']==(abs(det)>error)
    first=obs['minors'][0];det=F(first['det_center']);err=F(first['det_error_upper'])
    assert F(877,100)<det-err<det+err<F(878,100)
    assert err<F(12,10**28)
    out={'status':'fixed_minus64_real_regulator_obstruction_independently_replayed',
         'obstruction_sha256':sha(obsp),'trunks_sha256':sha(trp),'verifier_sha256':sha(__file__),
         'node_certificates':obs['node_certificate_sha256'],'primary_nodes':[1,9,25],
         'primary_determinant_rational_interval':['877/100','439/50'],
         'primary_determinant_error_less_than':'3/2500000000000000000000000000',
         'all_ten_minors_exclude_zero':True,'coupling_checked_against_actual_geometric_lift':True,
         'negative_control_zero_tail_rejected':tail_tamper_rejected,
         'arithmetic':'exact rational scalar quadrature and separate seven-subset determinant bound',
         'pure_transport_restarted':False,'floating_point_used':False,
         'geometric_proof':'LOWER_REGULATOR.md: real splitting forces Im(beta+g omega)=0 on these spheres',
         'lower_real_Hodge_extension_nonzero':True,'fixed_activity_minimal_rank':8,
         'complete_rank_eight_comparison_matrix_certified':False}
    save_json(DEST/'replay.json',out);print(json.dumps(out))

if __name__=='__main__':main()
