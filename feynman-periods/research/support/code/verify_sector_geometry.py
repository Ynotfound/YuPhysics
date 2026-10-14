"""Exact Hepp-sector identities; does not certify motives or period matrices."""
from itertools import permutations, combinations
import json
from pathlib import Path
import sympy as sp

root = Path(__file__).resolve().parents[1]
alpha = sp.symbols("a b c d")
mu = sp.symbols("mu1:5")
t1, t2, t3, tau = sp.symbols("t1 t2 t3 tau")
t = (t1, t2, t3)
q = t1*t2*t3
u = sum(sp.prod(alpha[i] for i in triplet)
        for triplet in combinations(range(4), 3))
f_euclidean = sum(mu[i]*alpha[i] for i in range(4))*u + tau*sp.prod(alpha)
records = []
for p in permutations(range(4)):
    coords = (sp.Integer(1), t1, t1*t2, q)
    substitution = {alpha[p[i]]: coords[i] for i in range(4)}
    pulled_f = sp.expand(f_euclidean.subs(substitution, simultaneous=True))
    jacobian = sp.Matrix(coords[1:]).jacobian(t).det()
    assert jacobian == t1**2*t2
    actual_d = sp.cancel(pulled_f / jacobian)
    expected_d = (sum(mu[p[i]]*coords[i] for i in range(4))
                  * (1+t3+t2*t3+q) + tau*q)
    assert sp.expand(actual_d-expected_d) == 0
    # This certifies D >= mu[p[0]] for all nonnegative variables.
    difference = sp.Poly(actual_d-mu[p[0]], *mu, *t, tau)
    assert all(coefficient >= 0 for coefficient in difference.coeffs())
    assert sp.diff(actual_d, tau) == q
    records.append({"permutation": [i+1 for i in p],
                    "factorization": True,
                    "jacobian": "t1^2*t2",
                    "positive_remainder_coefficients": True,
                    "tau_derivative": "t1*t2*t3"})
result = {"all_passed": True, "sectors": len(records),
          "method": "Exact symbolic polynomial identities over Q",
          "scope": "Hepp sector algebra and denominator lower bound only",
          "checks": records}
(root / "data" / "sector_geometry_checks.json").write_text(
    json.dumps(result, indent=2)+"\n", encoding="utf-8")
print(json.dumps({k: v for k, v in result.items() if k != "checks"}))
