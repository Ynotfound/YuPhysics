"""Integrate the actual lower lift on archived certified pure paths.

No ODE is restarted: all matrix polynomials and bounds are existing certificates.
The new scalar integral has its own uniform residual/input error bound.
"""
from node_branch_transport import propagate_vector
from node_frobenius import *
from itertools import combinations
import time

DEST=FOLDER/'lower_regulator'
BITS=320
G=[(F(x),F(0)) for x in [0,F(-3,64),F(-2,64),F(-1,64),0,0]]

def dot(r,v):return sumg(gm(x,y) for x,y in zip(r,v))
def rowmat(r,M):return [sumg(gm(r[i],M[6*i+j]) for i in range(6)) for j in range(6)]
def rn(r):return sum(map(ab,r))
def roundv(v):
    w=[(dyadic_round(x,BITS),dyadic_round(y,BITS)) for x,y in v]
    return [gs(F(1,1<<BITS),z) for z in w]
def record(v,e):return {'center':encode(v),'error_upper':str(e)}
def scalar_update(b,be,r,re,v,ve,sign=1):
    exact=ga(b,gs(sign,dot(r,v)));new=roundv([exact])[0]
    er=be+rn(r)*ve+re*(vnorm(v)+ve)+ab(ga(exact,gn(new)))
    return new,radius_up(er,BITS)

def integrated_row(step,g=G,first=False):
    h=tuple(F(y)-F(x) for x,y in zip(step['start'],step['end']))
    if first:
        P=[[gm(gp((0,1),n),(F(x),F(0))) for row in C for x in row]
           for n,C in enumerate(step['witness_coefficients_ascending'])]
    else:
        P=[[gs(F(1,1<<step['bits']),z) for z in decode(C)]
           for C in step['polynomial_coefficients_gaussian_integer']]
    A=[sumg(gs(F(1,n+1),M[j]) for n,M in enumerate(P)) for j in range(36)]
    exact=[gm(h,z) for z in rowmat(g,A)];r=roundv(exact)
    error=ab(h)*rn(g)*F(step['bounds']['error_norm_upper'])+rn([ga(x,gn(y)) for x,y in zip(exact,r)])
    return r,radius_up(error,BITS)

def append_row(J,je,R,re,C,ce):
    exact=va(J,rowmat(R,C));out=roundv(exact)
    err=je+rn(R)*ce+re*(mnorm(C)+ce)+rn([ga(x,gn(y)) for x,y in zip(exact,out)])
    return out,radius_up(err,BITS)

def trunks():
    fp=FOLDER/'first_trunk_step.json';f=json.loads(fp.read_text())
    J,je=integrated_row(f,first=True)
    C=[(F(x),F(y)) for rr,ii in zip(f['endpoint_center_real'],f['endpoint_center_imag']) for x,y in zip(rr,ii)]
    ce=F(f['bounds']['error_norm_upper']);logs=[{'file':str(fp.relative_to(FOLDER)),'sha256':sha(fp),'integral':record(J,je)}]
    for p in sorted((FOLDER/'high_trunk').glob('step_*.json')):
        d=json.loads(p.read_text());r,re=integrated_row(d);J,je=append_row(J,je,r,re,C,ce)
        C=[gs(F(1,1<<d['bits']),v) for v in decode(d['propagated_center_gaussian_integer'])];ce=F(d['propagated_error_upper'])
        logs.append({'file':str(p.relative_to(FOLDER)),'sha256':sha(p),'integral':record(J,je)})
    vertical=record(J,je)
    J=[Z]*6;je=F(0);C=mi();ce=F(0);landmarks={};hl=[]
    for p in sorted((FOLDER/'horizontal_trunk').glob('step_*.json')):
        d=json.loads(p.read_text());r,re=integrated_row(d);J,je=append_row(J,je,r,re,C,ce)
        C=[gs(F(1,1<<d['bits']),v) for v in decode(d['high_center_gaussian_integer'])];ce=F(d['high_error_upper'])
        hl.append({'file':str(p.relative_to(FOLDER)),'sha256':sha(p),'integral':record(J,je)})
        if d['landmark_pole'] is not None:landmarks[str(d['landmark_pole'])]=record(J,je)
    out={'status':'lower_lift_integrals_along_archived_trunks_certified','code_sha256':sha(__file__),
         'coupling_D':encode(G),'bits':BITS,'vertical':vertical,'landmarks':landmarks,'vertical_steps':logs,'horizontal_steps':hl}
    save_json(DEST/'trunk_integrals.json',out);return out

def node(a,tr):
    seedp=FOLDER/f'node_frobenius/node_{a}.json';d=json.loads(seedp.read_text());Q=real_matrix(sp.Matrix(d['Q']))
    P=[vfdecode(v) for v in d['frobenius_vector_coefficients_exact']];gq=rowmat(G,Q);h=(F(0),F(1,4))
    exact=gm(h,sumg(gs(1/(F(n)+F(1,2)),dot(gq,v)) for n,v in enumerate(P)))
    b=roundv([exact])[0];be=radius_up(2*ab(h)*rn(gq)*F(d['bounds']['vector_error_Q_basis'])+ab(ga(exact,gn(b))),BITS)
    # Independent local geometric consistency: beta/omega = -S * sum r_i c_i.
    rdata={1:[-1,-2,-4,8],9:[1,-2,-4,8],25:[-1,2,-4,8],49:[1,2,-4,8],81:[-1,-2,4,8]}
    rs=rdata[a];leadbeta=gs(2,gm(h,dot(gq,P[0])));leadomega=gs(F(-1,64),mv(Q,P[1])[0])
    ratio=gd(leadbeta,leadomega);assert ratio==(-sum(rs)*sum(x*y for x,y in zip(rs,[3,2,1,0])),0)
    init=record([b],be);logs=[]
    for p in sorted((FOLDER/f'node_branches/node_{a}').glob('step_*.json')):
        d=json.loads(p.read_text());g=gq if d['frame']=='Q' else G
        rr,re=integrated_row(d,g);v=[gs(F(1,1<<d['bits']),z) for z in decode(d['input_center_gaussian_integer'])]
        b,be=scalar_update(b,be,rr,re,v,F(d['input_error_upper']))
        logs.append({'file':p.name,'sha256':sha(p),'beta':record([b],be)})
    bp=FOLDER/f'node_branches/node_{a}/transport.json';bd=json.loads(bp.read_text());bits=bd['bits']
    cv=[gs(F(1,1<<bits),v) for v in decode(bd['based_center_gaussian_integer'])];ce=F(bd['based_error_upper'])
    ht=tr['landmarks'][str(a)];b,be=scalar_update(b,be,vfdecode(ht['center']),F(ht['error_upper']),cv,ce,-1)
    atc=record([b],be)
    invp=FOLDER/'high_trunk_inverse.json';inv=json.loads(invp.read_text());ib=inv['bits'];S=1<<BITS
    vi=[(dyadic_round(x,BITS),dyadic_round(y,BITS)) for x,y in cv]
    V=[gs(1<<(BITS-ib),v) for v in decode(inv['inverse_center_gaussian_integer'])]
    bv,bve,_=propagate_vector(V,F(inv['inverse_error_upper']),vi,ce,BITS)
    bv=[gs(F(1,S),v) for v in bv]
    jt=tr['vertical'];b,be=scalar_update(b,be,vfdecode(jt['center']),F(jt['error_upper']),bv,bve,-1)
    absolute=json.loads((FOLDER/'marked_lattice/absolute_periods.json').read_text())
    rho=vfdecode([absolute['rho_center']])[0];rhoe=F(absolute['rho_error_upper'])
    root=int(a**.5);assert root*root==a
    if a in [9,25,81]:rho=gm((F(0),F(-1)),rho)
    rho=gs(F(1,root),rho);rhoe/=root
    vals=[gs(F(-1,64),bv[0]),b];errs=[bve/64,be];cent=[];ee=[]
    for v,e in zip(vals,errs):
        exact=gm(rho,v);val=roundv([exact])[0]
        er=(ab(rho)+rhoe)*e+rhoe*ab(v)+ab(ga(exact,gn(val)))
        cent.append(val);ee.append(radius_up(er,BITS))
    out={'status':'actual_omega_and_beta_periods_on_transported_vanishing_sphere_certified',
         'node':a,'code_sha256':sha(__file__),'seed_sha256':sha(seedp),'branch_sha256':sha(bp),
         'trunk_integrals_sha256':sha(DEST/'trunk_integrals.json'),'inverse_sha256':sha(invp),
         'absolute_scale_sha256':sha(FOLDER/'marked_lattice/absolute_periods.json'),
         'local_beta_omega_leading_ratio':encode([ratio])[0],'local_beta':init,'steps':logs,
         'beta_at_high_basepoint':atc,'pure_at_physical_basepoint':record(bv,bve),
         'unscaled_beta_at_physical_basepoint':record([b],be),'rho':record([rho],rhoe),
         'omega_period':record([cent[0]],ee[0]),'beta_period':record([cent[1]],ee[1]),
         'basepoint':['-64','0'],'cycle_orientation':'Morse square roots as in ABSOLUTE_COMPARISON'}
    save_json(DEST/f'node_{a}.json',out)
    print(json.dumps({'node':a,'omega':[float(x) for x in cent[0]],'beta':[float(x) for x in cent[1]],
                      'errors':[float(x) for x in ee]}),flush=True);return out

def main():
    DEST.mkdir(exist_ok=True);t=time.time();tr=trunks();print('trunk scalar integrals saved',flush=True)
    ns=[node(a,tr) for a in [1,9,25,49,81]];minors=[]
    for combi in combinations(ns,3):
        rows=[];rs=[]
        for d in combi:
            p=vfdecode(d['omega_period']['center'])[0];b=vfdecode(d['beta_period']['center'])[0]
            rows.append([p[0],p[1],b[1]]);rs.append([F(d['omega_period']['error_upper'])]*2+[F(d['beta_period']['error_upper'])])
        det=F(sp.det(sp.Matrix(rows)));err=F(0)
        from itertools import permutations
        for perm in permutations(range(3)):
            u=F(1);v=F(1)
            for i,j in enumerate(perm):u*=abs(rows[i][j])+rs[i][j];v*=abs(rows[i][j])
            err+=u-v
        minors.append({'nodes':[d['node'] for d in combi],'rows':[[str(x) for x in row] for row in rows],
                       'radii':[[str(x) for x in row] for row in rs],'det_center':str(det),'det_error_upper':str(err),
                       'excludes_zero':abs(det)>err})
        print(json.dumps({'nodes':minors[-1]['nodes'],'det':float(det),'error':float(err),'excludes_zero':abs(det)>err}),flush=True)
    out={'status':'real_Hodge_splitting_obstruction_test','code_sha256':sha(__file__),
         'node_certificate_sha256':{str(d['node']):sha(DEST/f"node_{d['node']}.json") for d in ns},
         'minors':minors,'at_least_one_minor_nonzero':any(x['excludes_zero'] for x in minors),
         'independent_replay_complete':False,'geometric_argument':'see LOWER_REGULATOR.md; necessary splitting condition only',
         'elapsed_seconds':round(time.time()-t,2)}
    save_json(DEST/'obstruction.json',out)

if __name__=='__main__':main()
