"""Total elastic cross section: angular quadrature versus reordered dispersion."""
from pathlib import Path
import json,hashlib,time
import mpmath as mp
from verify_phi4_benchmark import B,B0,J,H,A,coefficients,x0
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/phi4_benchmark'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def avg_J(s):
    d=s-4
    def integrand(v):
        if v==1:return mp.mpf(0)
        z=1-v*v;k=d*z/4
        if abs(k)<mp.mpf('1e-12'):
            stable=sum((-1)**(j+1)*d**(j-1)*z**(j-2)/(j*4**j) for j in range(2,8))
        else:stable=(mp.log1p(k)-k)/(d*z*z)
        return 8*v*v*(A(v)-B0)*(stable-x0/(4*(4-x0*z)))
    return mp.quad(integrand,[0,.5,.9,1])
def main():
    start=time.monotonic();rows=[]
    for ss in [5,8,12,20]:
        s=mp.mpf(ss);d=s-4;beta=mp.sqrt(1-4/s)
        avg_b2=mp.quad(lambda z:(B(-d*(1-z)/2)-B0)**2,[-1,0,1])/2
        meanJ=avg_J(s)
        avg_a2=(H(s)+2*(-avg_b2-4*meanJ))/(1024*mp.pi**4)
        def one(z):
            return (B(s)+B(-d*(1-z)/2)+B(-d*(1+z)/2)-3*B0)/(32*mp.pi**2)
        avg_a1=mp.quad(lambda z:mp.re(one(z)),[-1,0,1])/2
        avg_square=mp.quad(lambda z:abs(one(z))**2,[-1,0,1])/2
        direct=[]
        for N in [12,24]:
            nodes,weights=mp.gauss_quadrature(N,'legendre')
            value=sum(weights[i]*J(-d*(1-nodes[i])/2) for i in range(N))/2
            direct.append(value)
        assert abs(direct[-1]-meanJ)<mp.mpf('1e-20')
        optical=-beta*avg_a1/(16*mp.pi)
        assert abs(mp.im(avg_a2)-optical)<mp.mpf('1e-50')
        rows.append({'s_over_m2':ss,'mean_J_t':mp.nstr(meanJ,44),
          'nested_quadrature_errors':[mp.nstr(abs(x-meanJ),6) for x in direct],
          'nested_orders':[12,24],
          'sigma_coefficients_inside_1_over_32pi_s':[mp.nstr(x,44) for x in [1,-2*avg_a1,avg_square-2*mp.re(avg_a2)]],
          'optical_g3_residual':mp.nstr(abs(mp.im(avg_a2)-optical),6),
          'inclusive_scope':'elastic equals total through g4 in this energy range' if ss<16 else 'elastic only; add 2to4 at g4'})
        print(json.dumps(rows[-1]),flush=True)
    result={'status':'total_elastic_cross_section_computed','benchmark_sha256':sha(OUT/'benchmark.json'),
      'script_sha256':sha(__file__),'rows':rows,
      'numerical_status':'High-precision quadrature consistency, not rigorous interval certification.',
      'g4_optical_theorem_independently_certified':False,'runtime_seconds':round(time.monotonic()-start,3)}
    (OUT/'elastic_total.json').write_text(json.dumps(result,indent=2)+'\n')
    print('total elastic quadratures finished')
if __name__=='__main__':main()

