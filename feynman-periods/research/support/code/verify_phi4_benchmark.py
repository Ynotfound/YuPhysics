"""Wick-contraction replay and finite two-loop scalar benchmark computations."""
from pathlib import Path
from itertools import permutations,product
from collections import defaultdict
from fractions import Fraction
from math import factorial
import json,hashlib,time
import mpmath as mp
import numpy as np
import sympy as sp
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/phi4_benchmark'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def key(V,ext,edges):
    strings=[]
    for perm in permutations(range(V)):
        a=','.join(str(perm[x]) for x in ext)
        es=sorted(tuple(sorted((perm[u],perm[v]))) for u,v in edges)
        strings.append(a+'|'+','.join(str(u)+str(v) for u,v in es))
    return min(strings)
def pairings(stubs):
    if not stubs:
        yield [];return
    first=stubs[0]
    for j in range(1,len(stubs)):
        for tail in pairings(stubs[1:j]+stubs[j+1:]):
            yield [(first[0],stubs[j][0])]+tail
def connected(V,es):
    reach={0}
    while True:
        new=reach|{v for u,v in es if u in reach}|{u for u,v in es if v in reach}
        if new==reach:return len(reach)==V
        reach=new
def replay_graphs(d):
    reports={}
    for V in (1,2,3):
        weights=defaultdict(Fraction);paircount=0
        for ext in product(range(V),repeat=4):
            degrees=[ext.count(v) for v in range(V)]
            factor=Fraction(1,factorial(V)*24**V)
            for k in degrees:factor*=Fraction(24,factorial(4-k))
            stubs=[(v,j) for v,k in enumerate(degrees) for j in range(k,4)]
            for es in pairings(stubs):
                paircount+=1
                if connected(V,es):weights[key(V,ext,es)]+=factor
        rows=[r for r in d['graphs'] if r['vertices']==V]
        assert len(weights)==len(rows)
        for r in rows:
            assert weights[key(V,r['external_vertices'],r['internal_edges'])]==Fraction(r['symmetry_weight'])
            es=r['internal_edges']
            assert r['one_particle_irreducible']==all(connected(V,es[:j]+es[j+1:]) for j in range(len(es)))
        reports[str(V)]={'connected_classes':len(rows),'unweighted_stub_pairings_examined':paircount}
    return reports

mp.mp.dps=60
x0=mp.mpf(4)/3
def B(x):
    x=mp.mpf(x)
    if not x:return mp.mpf(0)
    if x==4:return mp.mpf(2)
    if x<0:
        r=mp.sqrt(1-4/x);return 2-r*mp.log((r+1)/(r-1))
    if x<4:
        r=mp.sqrt(4/x-1);return 2-2*r*mp.atan(1/r)
    r=mp.sqrt(1-4/x);return 2-r*mp.log((1+r)/(1-r))+mp.j*mp.pi*r
B0=B(x0)
def A(v):
    v=mp.mpf(v)
    if not v:return mp.mpf(0)
    e=mp.atanh(v)
    return 3-(1-v*v)*e*e/(v*v)-2*e/v
def J(x):
    """Once-subtracted wineglass function, spectral representation."""
    x=mp.mpf(x)
    if x==x0:return mp.mpf(0)
    def F(v):
        if not v:return mp.mpf(0)
        if v==1:return mp.mpf(0) # endpoint value immaterial for integrable logarithm
        return 8*(x-x0)*v*v*(A(v)-B0)/(4-x0+x0*v*v)
    if x<4:
        return mp.quad(lambda v:F(v)/(4-x+x*v*v),[0,mp.mpf('.5'),mp.mpf('.9'),1])
    if x==4:
        return mp.quad(lambda v:2*(x-x0)*(A(v)-B0)/(4-x0+x0*v*v),[0,.5,.9,1])
    v0=mp.sqrt(1-4/x);f0=F(v0)
    def reg(v):
        if v==v0:return mp.diff(F,v0)/(2*x*v0)
        return (F(v)-f0)/(x*(v*v-v0*v0))
    re=mp.quad(reg,[0,v0,.9,1])+f0/(2*x*v0)*mp.log((1-v0)/(1+v0))
    return re+mp.j*mp.pi*v0*(A(v0)-B0)
def H(x):return -(B(x)-B0)**2-4*J(x)
def coefficients(s,z):
    t=-(s-4)*(1-z)/2;u=-(s-4)*(1+z)/2
    return sum(B(x)-B0 for x in [s,t,u])/(32*mp.pi**2),sum(H(x) for x in [s,t,u])/(1024*mp.pi**4)
def fmt(x):return [mp.nstr(mp.re(x),48),mp.nstr(mp.im(x),48)]
def direct_J_gauss(x,N):
    """Independent forest-subtracted Schwinger parameters (real x<4)."""
    nodes,weights=np.polynomial.legendre.leggauss(N)
    nodes=(nodes+1)/2;weights=weights/2
    u=nodes[:,None,None];y=nodes[None,:,None];z=nodes[None,None,:]
    yy=y*(1-y);t=z/(1-z)
    D0=1-(4/3)*u*(1-u);dx=(4/3-float(x))*u*(1-u);C0=1-(4/3)*yy
    integrand=(np.log1p(dx/(D0+t*C0))-np.log1p(dx/(D0+t+t*t*yy))/(1+t*yy)**2)/(z*(1-z))
    return float(np.einsum('i,j,k,ijk->',weights,weights,weights,integrand))
def main():
    started=time.monotonic()
    d=json.loads((OUT/'graphs.json').read_text())
    assert d['generator_sha256']==sha(ROOT/'code/certify_phi4_graphs.py')
    report=replay_graphs(d);print('Wick replay passed '+json.dumps(report),flush=True)
    # Symbolic graph denominator and finite scheme identities.
    a,b,c,e,s=sp.symbols('a b c e s')
    U=c*e+(a+b)*(c+e)
    F=sp.expand((a+b+c+e)*U-(a+b)*c*e-s*a*b*(c+e))
    factor=(c+e)*((a+b)**2+(a+b)*(c+e)+c*e-s*a*b)
    assert sp.expand(F-factor)==0
    q,aa,bb,X,Y,Z=sp.symbols('q aa bb X Y Z',real=True)
    g=q+aa*q*q+(2*aa*aa+bb)*q**3
    lhs=sp.series(g*g-2*g**3*X+g**4*(X*X+Y*Y-2*Z),q,0,5).removeO()
    rhs=q*q-2*q**3*(X-aa)+q**4*((X-aa)**2+Y*Y-2*(Z+2*aa*X-2*aa*aa-bb))
    assert sp.expand(lhs-rhs)==0
    checks=[]
    for xx in [-7,-1,mp.mpf(4)/3,3,5,8,12]:
        if xx<4:
            integral=-mp.quad(lambda z:mp.log(1-xx*z*(1-z)),[0,.5,1])
        else:
            beta=mp.sqrt(1-4/mp.mpf(xx));l=(1-beta)/2;h=(1+beta)/2
            integral=-mp.quad(lambda z:mp.log(abs(1-xx*z*(1-z))),[0,l,h,1])+mp.j*mp.pi*beta
        err=abs(integral-B(xx));assert err<mp.mpf('1e-50'),(xx,err)
        checks.append({'x':str(xx),'bubble_direct_residual':mp.nstr(err,4)})
    cuts=[]
    for s in [5,8,12,20]:
        v=mp.sqrt(1-4/mp.mpf(s))
        angular=mp.quad(lambda z:B(-(s-4)*(1-z)/2),[-1,0,1])/2
        assert abs(angular-A(v))<mp.mpf('1e-50')
        K=mp.re(B(s))+2*A(v)-3*B0
        rho=-v*K/(512*mp.pi**3)
        direct=-v*(mp.re(B(s))+2*angular-3*B0)/(512*mp.pi**3)
        assert abs(rho-direct)<mp.mpf('1e-50')
        cuts.append({'s_over_m2':s,'Im_M2_over_g3':mp.nstr(rho,48),'angular_residual':mp.nstr(abs(angular-A(v)),4)})
    direct=[]
    for xx in [-1,0,3]:
        target=float(J(xx));vals=[direct_J_gauss(xx,N) for N in [24,48,96,160]]
        assert abs(vals[-1]-target)<2e-7,(xx,vals,target)
        direct.append({'x':xx,'spectral_J':mp.nstr(J(xx),48),'gauss_orders':[24,48,96,160],
                       'parameter_integrals':vals,'last_residual':abs(vals[-1]-target),
                       'status':'numerical_crosscheck_not_interval_certificate'})
        print('Independent forest integral '+json.dumps(direct[-1]),flush=True)
    samples=[]
    for s,z in [(5,0),(8,0),(8,mp.mpf(1)/2),(12,0)]:
        one,two=coefficients(mp.mpf(s),mp.mpf(z))
        sigma=[mp.mpf(1),-2*mp.re(one),abs(one)**2-2*mp.re(two)]
        assert abs(mp.im(two)-next(mp.mpf(x['Im_M2_over_g3']) for x in cuts if x['s_over_m2']==s))<mp.mpf('1e-45')
        samples.append({'s_over_m2':s,'cos_theta':str(z),'A1':fmt(one),'A2':fmt(two),
                        'squared_amplitude_g2_g3_g4':[mp.nstr(x,48) for x in sigma]})
    ath,bth=coefficients(mp.mpf(4),mp.mpf(0))
    result={'status':'graphs_exact_and_two_loop_dispersive_benchmark_computed',
      'graphs_sha256':sha(OUT/'graphs.json'),'verifier_sha256':sha(__file__),
      'wick_replay':report,'symanzik_factorization_exact':True,'finite_scheme_identity_exact':True,
      'subtraction_point_s_over_m2':'4/3','B_subtraction':mp.nstr(B0,48),
      'bubble_integral_checks':checks,'absorptive_checks':cuts,'independent_forest_checks':direct,
      'physical_samples':samples,'threshold_coupling_conversion':{'a':fmt(ath),'b':fmt(bth),'equation':'g_threshold=g-a*g^2-b*g^3+O(g^4)'},
      'numerical_status':'High-precision quadrature and independent finite-parameter checks; not rigorous interval enclosures.',
      'scope':'Massive equal-mass real phi4; pole mass, unit residue; elastic amplitude to g^3, differential cross section to g^4 below four-particle threshold.',
      'not_done':['integrated_cross_section','2to4_phase_space_above_threshold','interval_certification_of_numeric_amplitudes'],
      'runtime_seconds':round(time.monotonic()-started,3)}
    (OUT/'benchmark.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'runtime_seconds':result['runtime_seconds'],'samples':len(samples)}))
if __name__=='__main__':main()

