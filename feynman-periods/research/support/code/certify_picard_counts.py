"""Exact banana point counts over F_p, F_p^2, F_p^3 and Picard certificates.

Uses integer finite fields, no floating arithmetic. Geometric justification of
the ten-node correction and Frobenius inference is in PICARD_CERTIFICATE.md.
"""
from itertools import product
from pathlib import Path
import hashlib
import json
import platform
import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path('E:/document/feynman_period/k3_minimality_relative_coaction.pdf')


class Field:
    def __init__(self, p, n):
        self.p, self.n, self.q = p, n, p**n
        q = self.q
        self.modulus = [0, 1] if n == 1 else self.find_modulus()
        digits = np.array([[(a // p**i) % p for i in range(n)]
                           for a in range(q)], dtype=np.int32)
        powers = p**np.arange(n, dtype=np.int32)
        self.add = np.zeros((q, q), dtype=np.int32)
        self.mul = np.zeros((q, q), dtype=np.int32)
        for a in range(q):
            self.add[a] = ((digits[a] + digits) % p) @ powers
            coefficients = np.zeros((q, 2*n-1), dtype=np.int32)
            for i in range(n):
                for j in range(n):
                    coefficients[:, i+j] += digits[a, i] * digits[:, j]
            coefficients %= p
            for k in range(2*n-2, n-1, -1):
                for i in range(n):
                    coefficients[:, k-n+i] -= coefficients[:, k] * self.modulus[i]
                coefficients %= p
            self.mul[a] = coefficients[:, :n] @ powers
        self.neg = ((-digits) % p) @ powers
        self.chi = -np.ones(q, dtype=np.int32)
        nonzero = np.arange(1, q)
        self.chi[self.mul[nonzero, nonzero]] = 1
        self.chi[0] = 0
        assert np.count_nonzero(self.chi == 1) == (q-1)//2
        assert np.array_equal(self.mul, self.mul.T)
        assert np.all(self.mul[1] == np.arange(q))
        assert np.all(self.add[np.arange(q), self.neg] == 0)
        assert np.all(np.sort(self.mul[1:, 1:], axis=1) == nonzero)

    def find_modulus(self):
        # Degrees 2 and 3: a monic polynomial is irreducible iff it has no root.
        for coeffs in product(range(self.p), repeat=self.n):
            if not coeffs[0]:
                continue
            f = list(coeffs) + [1]
            if all(sum(c*pow(a, i, self.p) for i, c in enumerate(f)) % self.p
                   for a in range(self.p)):
                return f
        raise AssertionError('No irreducible modulus')

    def count(self, masses, s):
        p, q = self.p, self.q
        mu = [m*m % p for m in masses]
        add, mul, neg = self.add, self.mul, self.neg
        b = np.arange(q)
        total, linear, identically_zero = 0, 0, 0
        for a in range(q):
            ab = mul[a, b]
            A = add[add[mul[mu[0], a], mul[mu[1], b]], mu[3]]
            B = add[add[ab, a], b]
            c2 = mul[mu[2], B]
            c1 = add[mul[A, B], mul[(mu[2]-s) % p, ab]]
            c0 = mul[A, ab]
            disc = add[mul[c1, c1], neg[mul[4 % p, mul[c2, c0]]]]
            quadratic = c2 != 0
            nonzero_linear = (~quadratic) & (c1 != 0)
            zero_poly = (~quadratic) & (c1 == 0) & (c0 == 0)
            total += int(np.sum(1+self.chi[disc[quadratic]], dtype=np.int64))
            total += int(np.count_nonzero(nonzero_linear))
            total += q*int(np.count_nonzero(zero_poly))
            linear += int(np.count_nonzero(nonzero_linear))
            identically_zero += int(np.count_nonzero(zero_poly))
        return {'n': self.n, 'q': q, 'modulus_low_to_high': self.modulus,
                'affine_quartic_count': total, 'infinity_count': 4*q-2,
                'ten_node_correction': 10*q, 'smooth_count': total+14*q-2,
                'linear_fibres': linear, 'zero_polynomial_fibres': identically_zero}

    def brute_extension(self, masses, s):
        # Independent evaluation of the original quartic for all triples.
        # Use at degree 2, not degree 3 (where quadratic fibres avoid q^3 work).
        add, mul = self.add, self.mul
        mu = [m*m % self.p for m in masses]
        b = np.arange(self.q)[:, None]
        c = np.arange(self.q)[None, :]
        total = 0
        for a in range(self.q):
            ab, ac, bc = mul[a, b], mul[a, c], mul[b, c]
            abc = mul[ab, c]
            U = add[add[abc, ab], add[ac, bc]]
            M = add[add[mul[mu[0], a], mul[mu[1], b]],
                    add[mul[mu[2], c], mu[3]]]
            total += int(np.count_nonzero(mul[M, U] == mul[s % self.p, abc]))
        return total


def brute_prime(masses, s, p):
    mu = [m*m for m in masses]
    return sum(((mu[0]*a+mu[1]*b+mu[2]*c+mu[3])
                * (a*b*c+a*b+a*c+b*c)-s*a*b*c) % p == 0
               for a, b, c in product(range(p), repeat=3))


def good_reduction(masses, s, p):
    signed = sorted({sum(e*m for e, m in zip(signs, masses))**2 % p
                     for signs in product((-1, 1), repeat=4)})
    assert p > 2 and all(m % p for m in masses)
    assert s % p != 0 and s % p not in signed
    return {'nonzero_masses': True, 's_mod_p': s % p,
            'signed_mass_squares_mod_p': signed, 'torus_singular_locus_avoided': True}


def geometry_identities():
    a, b, c, d, s = sp.symbols('a b c d s')
    u, v, w, z = sp.symbols('u v w z')
    mu = sp.symbols('mu1:5')
    F = sum(m*x for m, x in zip(mu, (a, b, c, d)))*(a*b*c+a*b*d+a*c*d+b*c*d)-s*a*b*c*d
    A, B = mu[0]*a+mu[1]*b+mu[3], a*b+a+b
    assert sp.expand(F.subs(d, 1)-(mu[2]*B*c*c+(A*B+(mu[2]-s)*a*b)*c+A*a*b)) == 0
    vertex = mu[0]*(v*w+v*z+w*z)
    assert sp.det(sp.hessian(vertex, (v, w, z))) == 2*mu[0]**3
    edge = u*(mu[0]*v*(w+z)+mu[2]*w*w+mu[3]*z*z+(mu[2]+mu[3]-s)*w*z)
    assert sp.factor(sp.det(sp.hessian(edge, (v, w, z)))) == -2*mu[0]**2*s*u**3
    badF = F.subs(dict(zip(mu, (1, 4, 16, 64)))).subs(s, 1)
    badpoint = {a: 8, b: 4, c: 2, d: -1}
    assert badF.subs(badpoint) == 0
    assert all(sp.diff(badF, x).subs(badpoint) == 0 for x in (a, b, c, d))
    return {'quadratic_fibre_identity': True, 'vertex_hessian_determinant': '2*mu_i^3',
            'edge_hessian_determinant': '-2*mu_i^2*s*(-mu_j/mu_i)^3',
            'wrong_s_plus_one_negative_control_detected': True}


def frobenius(p, records):
    X = sp.Symbol('X')
    traces = [r['smooth_count'] - 1 - r['q']**2 - 16*r['q'] for r in records]
    cs = [sp.Integer(1)]
    for k in range(1, 4):
        cs.append(-sum(cs[k-i]*traces[i-1] for i in range(1, k+1))/k)
    assert all(c.q == 1 for c in cs)
    # Orthogonal reciprocal spectrum: determinant sign - forces middle c3=0.
    # Nonzero c3 therefore determines the + sign, without assuming it.
    assert cs[3] != 0, 'Three counts do not determine reciprocity sign here'
    coeffs = [int(c) for c in cs] + [int(cs[2]*p*p), int(cs[1]*p**4), p**6]
    P = sp.Poly.from_list(coeffs, X)
    f = sp.Poly(P, X, modulus=5)
    x = sp.Poly(X, X, modulus=5)
    def powmod(exponent):
        ans, base = sp.Poly(1, X, modulus=5), x
        while exponent:
            if exponent & 1:
                ans = (ans*base).rem(f)
            base = (base*base).rem(f)
            exponent >>= 1
        return ans
    g2 = sp.gcd(f, powmod(5**2)-x)
    g3 = sp.gcd(f, powmod(5**3)-x)
    r6 = (powmod(5**6)-x).rem(f)
    irreducible = g2.degree() == g3.degree() == 0 and r6.is_zero
    normalized = sp.Poly(P.as_expr().subs(X, p*X)/p**6, X)
    nonintegral = any(c.q != 1 for c in normalized.all_coeffs())
    return {'residual_traces': traces, 'residual_polynomial_coefficients': coeffs,
            'residual_polynomial': str(P.as_expr()), 'middle_coefficient_nonzero': True,
            'mod_5_rabin': {'gcd_x25_minus_x': str(g2.as_expr()),
                            'gcd_x125_minus_x': str(g3.as_expr()),
                            'remainder_x15625_minus_x': str(r6.as_expr()),
                            'irreducible': bool(irreducible)},
            'normalized_polynomial': str(normalized.as_expr()),
            'normalized_has_nonintegral_coefficient': nonintegral,
            'no_residual_root_is_p_times_root_of_unity': bool(irreducible and nonintegral),
            'picard_upper_bound_16_given_geometric_model': bool(irreducible and nonintegral)}


def main():
    samples = [([1, 2, 4, 8], -1, 11, [296, 16472, 1788824]),
               ([1, 3, 7, 15], -2, 13, [376, 31192, 4856104]),
               ([1, 2, 4, 8], -64, 11, None)]
    result = {'status': 'running', 'method': 'exact finite field quadratic fibre enumeration',
              'python': platform.python_version(), 'numpy': np.__version__,
              'sympy': sp.__version__, 'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'geometry_identity_checks': geometry_identities(),
              'source_path': str(SOURCE), 'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              'original_pdf_pages_visually_checked': [2, 4], 'samples': []}
    output = ROOT / 'data/picard_point_counts.json'
    for masses, s, p, expected in samples:
        record = {'masses': masses, 's': s, 'p': p, 'old_reported_counts': expected,
                  'good_reduction_checks': good_reduction(masses, s, p), 'extensions': []}
        for n in (1, 2, 3):
            field = Field(p, n)
            count = field.count(masses, s)
            if n == 1:
                assert count['affine_quartic_count'] == brute_prime(masses, s, p)
                count['independent_original_cubic_loop_checked'] = True
            if n == 2:
                assert count['affine_quartic_count'] == field.brute_extension(masses, s)
                count['independent_original_cubic_loop_checked'] = True
            if expected:
                assert count['smooth_count'] == expected[n-1], (masses, s, n, count)
            record['extensions'].append(count)
            print(json.dumps({'masses': masses, 's': s, 'p': p, **count}), flush=True)
            del field
        record.update(frobenius(p, record['extensions']))
        result['samples'].append(record)
        output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
        print(json.dumps({k: v for k, v in record.items() if k != 'extensions'}), flush=True)
    result['status'] = 'completed'
    output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
