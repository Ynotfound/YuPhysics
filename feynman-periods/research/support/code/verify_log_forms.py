"""Logarithmic form frame, constant residues, and algebraic deformation checks."""
from itertools import combinations
from pathlib import Path
import hashlib
import json
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
paths = [ROOT/'data/boundary_kummer_flags.json', ROOT/'data/regulator_cup_splitting.json']
d, cups = [json.loads(p.read_text(encoding='utf-8')) for p in paths]
labels, edges = d['labels'], d['oriented_edges']
V = sp.Matrix(d['coordinate_unit_valuations']).row_join(sp.zeros(20, 1))
Z = sp.Matrix(d['torus_double_residues_Z'])
W = sp.Matrix(cups['extra_unit_cup_flags'])
pairs = list(combinations(range(4), 2))

# q_ij^- = mu_i alpha_i/alpha_j - mu_j alpha_j/alpha_i.
# Columns: dlog x1^x2, x1^x3, x2^x3, x1^M, x2^M, x3^M.
T = sp.Matrix([[0, 1, 1, 0, 0, -1], [-1, 0, 1, 0, 1, 0],
               [0, 0, 1, 0, 0, 0], [-1, -1, 0, -1, 0, 0],
               [0, -1, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0]])
Gminus = Z.row_join(W)*T.inv()
Gplus = sp.zeros(30, 6)
for e, (a, b) in enumerate(edges):
    for col, (i, j) in enumerate(pairs):
        orders = [V[labels.index(v), i]-V[labels.index(v), j] for v in (a, b)]
        if orders == [-1, -1]:
            Gplus[e, col] = Gminus[e, col]
        elif orders == [1, 1]:
            Gplus[e, col] = -Gminus[e, col]
        else:
            assert Gminus[e, col] == 0
partial = sp.Matrix(d['incidence'])
assert partial*Gplus == sp.zeros(20, 6)
assert Gplus*sp.ones(6, 1) == sp.zeros(30, 1)
assert Gplus.rank() == 5
assert Gminus.row_join(Gplus[:, :5]).rank() == 11
mass_residues = sp.zeros(30, 4)
for col, (i, j) in enumerate(pairs):
    mass_residues[:, i] -= (Gplus[:, col]+Gminus[:, col])/2
    mass_residues[:, j] -= (Gplus[:, col]-Gminus[:, col])/2
assert mass_residues == sp.Matrix(d['mass_flag_matrix'])

x, y, z, s = sp.symbols('x y z s')
mu = sp.symbols('mu1:5')
alpha = [x, y, z, sp.Integer(1)]
M = sum(m*a for m, a in zip(mu, alpha))
N = 1/x+1/y+1/z+1
f = sp.expand(x*y*z*(M*N-s))
fx, fy, fz = [sp.diff(f, a) for a in (x, y, z)]
qm = sp.Matrix([mu[i]*alpha[i]/alpha[j]-mu[j]*alpha[j]/alpha[i] for i, j in pairs])
actual_cups = sp.Matrix([fz/(x*y), -fy/(x*z), fx/(y*z),
    (mu[1]*fz-mu[2]*fy)/(x*M), (-mu[0]*fz+mu[2]*fx)/(y*M),
    (mu[0]*fy-mu[1]*fx)/(z*M)])
field = sp.QQ.frac_field(*mu, s)
for value in actual_cups-T.T*qm:
    numerator = sp.fraction(sp.cancel(value))[0]
    assert sp.reduced(numerator, [f], x, y, z, domain=field)[1] == 0
for i in range(4):
    assert sp.cancel(sp.diff(-mu[i]*alpha[i]*N/f, s)-mu[i]*sp.diff(1/f, mu[i])) == 0
assert sp.cancel(sum(mu[i]*sp.diff(1/f, mu[i]) for i in range(4))
                 + s*sp.diff(1/f, s)+1/f) == 0

# Formal coefficients in 1 and the 12 weighted ordered ratios.
# Adding the defining Laurent relation to the 12 proposed basis columns
# must give an invertible 13x13 matrix.
ordered = [(i, j) for i in range(4) for j in range(4) if i != j]
columns = [sp.eye(13)[:, 0]]
for parity, subset in [(-1, pairs), (1, pairs[:5])]:
    for i, j in subset:
        c = sp.zeros(13, 1)
        c[1+ordered.index((i, j))] = 1
        c[1+ordered.index((j, i))] = parity
        columns.append(c)
relation = sp.Matrix([sum(mu)-s]+[1]*12)
basis_det = sp.factor(sp.Matrix.hstack(*columns, relation).det())
assert basis_det != 0 and not basis_det.free_symbols

# Infinitesimal projective stabilizer of the ten fixed Sylvester nodes.
nodes = [sp.eye(4)[:, i] for i in range(4)]
nodes += [sp.eye(4)[:, i]-sp.eye(4)[:, j] for i, j in combinations(range(4), 2)]
entries = sp.symbols('a0:16')
A = sp.Matrix(4, 4, entries)
equations = [(A*p)[i]*p[j]-(A*p)[j]*p[i]
             for p in nodes for i, j in combinations(range(4), 2)]
matrix, _ = sp.linear_eq_to_matrix(equations, entries)
assert matrix.rank() == 15
null = matrix.nullspace()[0]
assert sp.Matrix(4, 4, null) == sp.eye(4)

def strings(a):
    return [[str(v) for v in row] for row in a.tolist()]

result = {
    'status': 'exact_log_form_and_deformation_checks_passed',
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'input_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
    'pair_order_one_based': [[i+1, j+1] for i, j in pairs],
    'cup_to_antisymmetric_form_matrix': strings(T),
    'antisymmetric_residues': strings(Gminus),
    'symmetric_residues': strings(Gplus),
    'symmetric_residue_rank': 5, 'total_residue_rank': 11,
    'twelve_form_basis_plus_relation_determinant': str(basis_det),
    'mass_lift_residues_equal_mass_flags': True,
    'mass_lift_formula': 'beta_i = -mu_i*alpha_i*(sum_j 1/alpha_j)*omega = -s*mu_i*alpha_i/M*omega on X',
    'connection_identity': 'nabla_s beta_i = mu_i nabla_mu_i omega',
    'active_lift_at_m1248_sminus64': '64*(3*a+8*b+16*c)/(a+4*b+16*c+64*d)*omega',
    'projective_node_stabilizer_linear_system_rank': 15,
    'stabilizer_kernel': 'scalar matrices only',
    'very_general_minimality_proof': '../DE_RHAM_MINIMALITY.md',
    'fixed_basepoint_regulator_nonvanishing_certified': False,
}
(ROOT/'data/log_form_frame.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
print(json.dumps({k: result[k] for k in
    ('status', 'symmetric_residue_rank', 'total_residue_rank',
     'twelve_form_basis_plus_relation_determinant', 'mass_lift_residues_equal_mass_flags',
     'projective_node_stabilizer_linear_system_rank', 'fixed_basepoint_regulator_nonvanishing_certified')}))
