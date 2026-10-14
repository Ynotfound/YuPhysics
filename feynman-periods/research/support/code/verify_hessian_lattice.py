"""Exact banana/Hessian coordinate and marked-configuration comparison.

Primitivity is a theorem input plus a smooth-family argument, not an output
of this finite calculation. No geometric marking of candidate Q6 is asserted.
"""
from pathlib import Path
from itertools import combinations, product
import hashlib
import json
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

ROOT = Path(__file__).resolve().parents[1]
boundary_path = ROOT / 'data/boundary_projector.json'
boundary = json.loads(boundary_path.read_text(encoding='utf-8'))
labels = boundary['labels']
alpha = sp.symbols('a b c d')
mu = sp.symbols('mu1:5', nonzero=True)
s = sp.Symbol('s', nonzero=True)
M = sum(u*x for u, x in zip(mu, alpha))
U = sum(sp.prod(alpha[j] for j in range(4) if j != i) for i in range(4))
F = sp.expand(M*U - s*sp.prod(alpha))
y = [mu[i]*alpha[i] for i in range(4)] + [-M]
lam = list(mu) + [s]
Q = sum(lam[i]*sp.prod(y[j] for j in range(5) if j != i) for i in range(5))
assert sp.expand(Q + sp.prod(mu)*F) == 0

# Independent rank-one determinant identity used for the cubic Hessian.
d = sp.symbols('d1:5')
t = sp.Symbol('t')
H = sp.diag(*d) + t*sp.ones(4)
expected = sp.prod(d) + t*sum(sp.prod(d[j] for j in range(4) if j != i) for i in range(4))
assert sp.expand(H.det(method='domain-ge') - expected) == 0

# Source labels use 0,...,4; our translated source labels use 1,...,5.
def hessian_label(label):
    ids = set(map(int, label[1:]))
    if label[0] == 'V':
        return tuple(sorted(set(range(1, 5)) - ids))
    if label[0] == 'E':
        return tuple(sorted((set(range(1, 5)) - ids) | {5}))
    if label[0] == 'C':
        return tuple(sorted(set(range(1, 5)) - ids))
    if label[0] == 'F':
        return tuple(sorted(ids | {5}))
    raise ValueError(label)

marked = [hessian_label(label) for label in labels]
assert len(set(marked)) == 20
assert set(marked) == set(combinations(range(1, 6), 2)) | set(combinations(range(1, 6), 3))
source_gram = -2*sp.eye(20)
for i, j in combinations(range(20), 2):
    a, b = set(marked[i]), set(marked[j])
    if len(a) != len(b) and (a < b or b < a):
        source_gram[i, j] = source_gram[j, i] = 1
G = sp.Matrix(boundary['gram'])
coordinates = sp.Matrix(boundary['all_divisors_basis_coordinates'])
assert source_gram == coordinates.T*G*coordinates
omitted = set(labels) - set(boundary['basis_labels'])
# N234, N14, N23, N24 in original zero-based source labels.
source_omitted = {(3, 4, 5), (2, 5), (3, 4), (3, 5)}
assert {hessian_label(label) for label in omitted} == source_omitted

def snf(matrix):
    return [abs(int(v)) for v in smith_normal_form(matrix, domain=sp.ZZ).diagonal()]

hyperbolic = sp.Matrix([[0, 1], [1, 0]])
T = sp.diag(hyperbolic, 2*hyperbolic, sp.Matrix([[-4, 2], [2, -4]]))
assert G.det() == -48 and T.det() == 48
assert snf(G) == [1]*12 + [2, 2, 2, 6]
assert snf(T) == [1, 1, 2, 2, 2, 6]

samples = [([1, 2, 4, 8], -1), ([1, 3, 7, 15], -2), ([1, 2, 4, 8], -64)]
sample_checks = []
for masses, parameter in samples:
    thresholds = sorted({sum(e*m for e, m in zip(signs, masses))**2
                         for signs in product((-1, 1), repeat=4)})
    assert parameter != 0 and parameter not in thresholds
    sample_checks.append({'masses': masses, 's': parameter,
                          'threshold_squares': thresholds,
                          'smooth_sylvester_cubic': True})

source_pdf = ROOT / 'tmp/pdfs/dardanelli-van-geemen-0409322.pdf'
output = {
    'status': 'exact_coordinate_and_configuration_comparison_passed',
    'source': {'title': 'Hessians and the moduli space of cubic surfaces',
               'authors': 'Elisa Dardanelli and Bert van Geemen',
               'url': 'https://arxiv.org/pdf/math/0409322',
               'sections': ['1.5', '1.6'], 'pdf_pages_one_based': [2, 3],
               'local_pdf_sha256': hashlib.sha256(source_pdf.read_bytes()).hexdigest()},
    'boundary_input_sha256': hashlib.sha256(boundary_path.read_bytes()).hexdigest(),
    'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'coordinate_identity': 'Q_lambda(mu_i*alpha_i,-M) = -prod(mu_i)*F',
    'cubic_coefficients': ['1/mu1', '1/mu2', '1/mu3', '1/mu4', '1/s'],
    'cubic_hessian_identity': 'det Hess(sum y_i^3/lambda_i) = 6^4 Q_lambda/prod(lambda_i)',
    'curve_dictionary': dict(zip(labels, ['N'+''.join(map(str, x)) for x in marked])),
    'same_omitted_basis_curves': sorted(omitted),
    'boundary_gram_determinant': str(G.det()),
    'boundary_snf': snf(G),
    'transcendental_gram': T.tolist(),
    'transcendental_snf': snf(T),
    'transcendental_determinant': str(T.det()),
    'samples': sample_checks,
    'primitivity_proof': '../LATTICE_SATURATION.md',
    'primitivity_proved_by_this_script': False,
    'candidate_Q6_geometric_marking_certified': False,
    'absolute_period_matrix_certified': False
}
# Convert only exact SymPy integers in matrix entries.
output['transcendental_gram'] = [[int(v) for v in row] for row in T.tolist()]
(ROOT/'data/hessian_lattice.json').write_text(json.dumps(output, indent=2)+'\n', encoding='utf-8')
print(json.dumps({'status': output['status'], 'curves_matched': 20,
                  'boundary_det': -48, 'transcendental_det': 48,
                  'same_source_basis': True, 'geometric_Q6_marking': False}))
