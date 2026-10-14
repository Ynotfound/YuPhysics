"""Exact independent models for the weighted-defect lemma and rank-one traces.

This checks proof formulas; it does NOT produce a banana-node certificate.
"""
from pathlib import Path
import json,hashlib
import sympy as sp

def main():
    root=Path(__file__).resolve().parents[1];t=sp.Symbol('t')
    D=sp.diag(0,sp.Rational(1,2),sp.Rational(1,2),sp.Rational(1,2),sp.Rational(1,2),sp.Rational(1,2))
    q=(4-t)**2;N=sp.zeros(6);N[1,0]=8+2*t
    g=sp.Matrix([1,t/(1-t/4),0,0,0,0])
    assert (t*g.diff(t)-D*g-t*N*g/q).applyfunc(sp.cancel)==sp.zeros(6,1)
    P=sp.Matrix([1,t+t*t/4+t**3/16,0,0,0,0]);m=4
    R=(q*t*P.diff(t)-q*D*P-t*N*P).applyfunc(sp.expand)
    assert all(sp.Poly(v,t).nth(j)==0 for v in R for j in range(m))
    S=(R/t**m).applyfunc(sp.cancel)
    norm=lambda V:max(sum(abs(c) for c in sp.Poly(v,t).all_coeffs()) for v in V)
    lower=sp.Rational(9);M=sp.Rational(10,9);delta=norm(S)/lower
    kappa=M/(m+sp.Rational(1,2));error=delta/((m-sp.Rational(1,2))*(1-kappa))
    actual=max(abs(v) for v in (g-P).subs(t,1))
    assert 0<actual<=error and 0<kappa<1
    # The exact model really uses the +1/2 eigenspace and nonzero coupling.
    K=sp.diag(1,1,-1,-1,-1,-1)
    vectors=[sp.Matrix(v) for v in [[1,0,0,0,0,0],[1,1,2,0,0,0],[2,-1,1,0,0,0]]]
    d=[(v.T*K*v)[0] for v in vectors]
    Rs=[-2*v*v.T*K/a for v,a in zip(vectors,d)]
    for R in Rs:
        T=sp.eye(6)+R
        assert T*T==sp.eye(6) and T.T*K*T==K and T.det()==-1
    pair=(vectors[0].T*K*vectors[1])[0]
    assert sp.trace(Rs[0]*Rs[1])==4*pair**2/(d[0]*d[1])
    triple=sp.prod((vectors[i].T*K*vectors[(i+1)%3])[0] for i in range(3))
    assert sp.trace(Rs[0]*Rs[1]*Rs[2])==-8*triple/sp.prod(d)
    # Verify the reflection enclosure on an independent exact perturbed vector.
    V=vectors[0];eps=sp.Rational(1,100)
    v=V+sp.Matrix([sp.I*eps,eps,0,0,0,0])
    dc=(V.T*K*V)[0];dt=(v.T*K*v)[0]
    Delta=6*eps*(2+eps);d0=abs(dc);low=d0-Delta;assert low>0
    center=sp.eye(6)-2*V*V.T*K/dc;true=sp.eye(6)-2*v*v.T*K/dt
    upper=2*(Delta/low+Delta/(d0*low))
    diff=(true-center).applyfunc(sp.expand_complex)
    measured=max(sum(abs(sp.re(diff[i,j]))+abs(sp.im(diff[i,j])) for j in range(6)) for i in range(6))
    assert measured<=upper
    out={'status':'independent_frobenius_and_reflection_models_passed',
         'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'weighted_model':'g=(1,t/(1-t/4),0,0,0,0), D0=diag(0,1/2,...,1/2)',
         'm':m,'M':str(M),'delta':str(delta),'kappa':str(kappa),
         'exact_endpoint_error':str(actual),'certified_error':str(error),
         'reflection_involution_and_polarization':True,'pair_and_triple_trace_signs':True,
         'perturbed_vector_enclosure_checked':True,'banana_node_frobenius_certified':False,
         'global_loop_certified':False}
    (root/'data/frobenius_reflection_lemma_check.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out),flush=True)

if __name__=='__main__': main()
