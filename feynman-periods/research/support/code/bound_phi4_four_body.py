"""Rigorous weighted Darboux enclosures of Phi4(20), then analytic amplitude bounds."""
from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'runtime/python-flint'))
from flint import arb,ctx
OUT=ROOT/'data/phi4_four_body'
ctx.prec=180
def endpoints(x):
    def exact(a):
        n,e=a.man_exp()
        return [str(n),int(e)]
    return dict(lower=exact(x.lower()),upper=exact(x.upper()),display=x.str(30))
def Fu(x):return 2*x*x.sqrt()/3-4*x*x*x.sqrt()/5+2*x*x*x*x.sqrt()/7
def Fv(x):
    if x==0:return -arb.pi()/16
    if x==1:return arb.pi()/16
    return (2*x-1)*(x*(1-x)).sqrt()/4+(2*x-1).asin()/8
def calc(N):
    E=arb(20).sqrt();eps=E-4;S=arb(0)
    wu=[Fu(arb(i+1)/N)-Fu(arb(i)/N) for i in range(N)]
    wv=[Fv(arb(i+1)/N)-Fv(arb(i)/N) for i in range(N)]
    for i in range(N):
        u=arb(arb(2*i+1)/(2*N),arb(1)/(2*N))
        A=2+eps*u
        for j in range(N):
            v=arb(arb(2*j+1)/(2*N),arb(1)/(2*N))
            B=2+eps*(1-u)*v
            f=((E+A+B)*(E-A+B)*(E+A-B)*(A+2)*(B+2)).sqrt()
            S+=wu[i]*wv[j]*f
    phi=eps**3*eps.sqrt()/(512*arb.pi()**5*20)*S
    # For s=20, each mixed denominator lies in [-10,-2].
    ampmax=3-4/(20-2*E)
    fluxsym=48*20*(arb(1)-arb(4)/20).sqrt()
    lower=phi.lower()/100/fluxsym
    upper=phi.upper()*ampmax**2/fluxsym
    sigma=lower.union(upper)
    return dict(N=N,phase_volume=endpoints(phi),sigma_over_g4=endpoints(sigma),
                explicit_phase_window_verified=bool(phi>arb('7.58e-8') and phi<arb('7.63e-8')),
                explicit_sigma_window_verified=bool(sigma>arb('8.8e-13') and sigma<arb('6.2e-10')))
def main():
    data=dict(status='rigorous_coarse_Darboux_and_amplitude_bounds',precision=ctx.prec,
              script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              rows=[calc(64),calc(128)],not_a_precision_certificate_for_QMC=True,
              integrand='epsilon^(7/2)/(512*pi^5*s)*(1-u)^2*sqrt(u)*sqrt(v*(1-v))*sqrt((E+A+B)*(E-A+B)*(E+A-B)*(A+2)*(B+2))')
    assert data['rows'][-1]['explicit_phase_window_verified']
    assert data['rows'][-1]['explicit_sigma_window_verified']
    (OUT/'bounds.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
    print(json.dumps(data,indent=2))
if __name__=='__main__':main()

