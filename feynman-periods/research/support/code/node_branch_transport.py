"""Piecewise constant frames and certified vector propagation on node branches."""
from node_frobenius import *
from certify_high_trunk import proposal

SWITCH_HEIGHT={1:F(8),9:F(16),25:F(32),49:F(16),81:F(16),121:F(16),169:F(16),225:F(16)}


def transformed_source(source,Q,Qi):
    return (source[0],[[rat(v) for v in Qi*sp.Matrix(6,6,list(map(sp.Rational,M)))*Q] for M in source[1]])


def initial_Q_vector(seed,bits):
    P=[vfdecode(v) for v in seed['frobenius_vector_coefficients_exact']]
    exact=[sumg(v[i] for v in P) for i in range(6)]
    S=1<<bits;C=[(dyadic_round(x,bits),dyadic_round(y,bits)) for x,y in exact]
    rounding=vnorm([ga(x,gn(gs(F(1,S),y))) for x,y in zip(exact,C)])
    e=radius_up(F(seed['bounds']['vector_error_Q_basis'])+rounding,bits)
    return C,e,rounding


def switch_vector(C,e,Q,bits):
    S=1<<bits;exact=mv(real_matrix(Q),[gs(F(1,S),v) for v in C])
    new=[(dyadic_round(x,bits),dyadic_round(y,bits)) for x,y in exact]
    rd=vnorm([ga(x,gn(gs(F(1,S),y))) for x,y in zip(exact,new)])
    return new,radius_up(mnorm(real_matrix(Q))*e+rd,bits),rd


def propagate_vector(L,le,C,ce,bits):
    S=1<<bits;product=mv(L,C)
    new=[(dyadic_round(F(x,S*S),bits),dyadic_round(F(y,S*S),bits)) for x,y in product]
    rd=max(F(ab(ga(x,gn(gs(S,y)))),S*S) for x,y in zip(product,new))
    raw=(F(mnorm(L),S)+le)*ce+le*F(vnorm(C),S)+rd
    return new,radius_up(raw,bits),rd


def reflection_at_vector(V,e,K,bits):
    # Reuse the already independently checked bound, with no new tail:
    # P is constant, Q is identity, and e is the transported vector enclosure.
    return reflection_from_vector([V],sp.eye(6),e,K,bits)
