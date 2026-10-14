"""Second positive-cell enclosure, using endpoint factor bounds (no generator imports)."""
from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'runtime/python-flint'))
from flint import arb,ctx
import sympy as sp
ctx.prec=220;N=128
E=arb(20).sqrt();e=E-4
def U(x):return 2*x**arb('1.5')/3-4*x**arb('2.5')/5+2*x**arb('3.5')/7
def V(x):
    if x==0:return -arb.pi()/16
    if x==1:return arb.pi()/16
    return ((2*x-1)*(x*(1-x)).sqrt()*2+(2*x-1).asin())/8
lower=arb(0);upper=arb(0)
for i in reversed(range(N)):
    u0=arb(i)/N;u1=arb(i+1)/N
    a0=2+e*u0;a1=2+e*u1
    wu=U(u1)-U(u0)
    for j in reversed(range(N)):
        v0=arb(j)/N;v1=arb(j+1)/N
        b0=2+e*(1-u1)*v0;b1=2+e*(1-u0)*v1
        lo=((E+a0+b0)*(E-a1+b0)*(E+a0-b1)*(a0+2)*(b0+2)).sqrt()
        hi=((E+a1+b1)*(E-a0+b1)*(E+a1-b0)*(a1+2)*(b1+2)).sqrt()
        w=wu*(V(v1)-V(v0))
        lower+=w*lo;upper+=w*hi
fac=e**arb('3.5')/(512*arb.pi()**5*20)
L=fac*lower;H=fac*upper
assert L>arb('7.58e-8') and H<arb('7.63e-8')
flux=48*20*(arb(4)/5).sqrt()
assert L/(100*flux)>arb('8.8e-13')
assert H*(3-4/(20-2*E))**2/flux<arb('6.2e-10')
# Independent symbolic normalization / antiderivative checks.
u,v,ep=sp.symbols('u v ep',positive=True)
AA=2+ep*u;BB=2+ep*(1-u)*v;EE=4+ep
raw=16*AA**2*BB**2*ep**4*(1-u)**2*((EE-AA-BB)*(EE+AA+BB)*(EE-AA+BB)*(EE+AA-BB))/EE**4*(AA**2-4)/AA**2*(BB**2-4)/BB**2
factored=16*ep**7*(1-u)**4*u*v*(1-v)*(EE+AA+BB)*(EE-AA+BB)*(EE+AA-BB)*(AA+2)*(BB+2)/EE**4
assert sp.expand(raw-factored)==0
assert sp.simplify(sp.diff(2*u**sp.Rational(3,2)/3-4*u**sp.Rational(5,2)/5+2*u**sp.Rational(7,2)/7,u)-(1-u)**2*sp.sqrt(u))==0
assert sp.simplify(sp.diff((2*v-1)*sp.sqrt(v*(1-v))/4+sp.asin(2*v-1)/8,v)-sp.sqrt(v*(1-v)))==0
S,a,b=sp.symbols('S a b',positive=True)
massless=sp.integrate(sp.integrate((1-a/S)*(1-b/a),(b,0,a)),(a,0,S))
assert massless==S**2/12
beta_u=sp.Rational(2,3)-sp.Rational(4,5)+sp.Rational(2,7)
assert beta_u==sp.Rational(16,105)
threshold=32*sp.sqrt(2)*beta_u*(sp.pi/8)/(512*sp.pi**5*16)
assert sp.simplify(threshold-sp.sqrt(2)/(13440*sp.pi**4))==0
out=dict(status='independent_positive_endpoint_enclosure_passed',N=N,precision=ctx.prec,
         rigorous_phase_window=['7.58e-8','7.63e-8'],rigorous_sigma_over_g4_window=['8.8e-13','6.2e-10'],
         phase_lower=L.str(25),phase_upper=H.str(25),symbolic_factorization=True,
         primitive_derivatives=True,massless_volume='s^2/(24576*pi^5)',
         threshold_phase='sqrt(2)*epsilon^(7/2)/(13440*pi^4)',
         threshold_sigma='sqrt(2/3)*epsilon^(7/2)/(5160960*pi^4)',
         bounds_sha256=hashlib.sha256((ROOT/'data/phi4_four_body/bounds.json').read_bytes()).hexdigest(),
         script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(ROOT/'data/phi4_four_body/bounds_replay.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))

