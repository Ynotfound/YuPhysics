"""Exact unit cup residues and rational splitting frames, not regulator values."""
from itertools import combinations
from pathlib import Path
import hashlib
import json
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
source = ROOT/'data/boundary_kummer_flags.json'
d = json.loads(source.read_text(encoding='utf-8'))
labels, edges = d['labels'], d['oriented_edges']
V = sp.Matrix(d['coordinate_unit_valuations'])
P = sp.Matrix(d['Tate_basis_columns'])
S = sp.Matrix(d['Tate_basis_left_inverse'])
Z = sp.Matrix(d['torus_double_residues_Z'])
F = sp.Matrix(d['mass_flag_matrix'])
partial = sp.Matrix(d['incidence'])
R = sp.eye(30)[d['anchor_rows_zero_based'], :]

# Fourth unit is M/d. At E and F, ord M=1; at V and C, ord M=0.
fourth = []
for name in labels:
    ids = list(map(int, name[1:]))
    alpha = [int(i in ids) if name[0] == 'F' else int(i not in ids)
             for i in range(1, 5)]
    fourth.append(int(name[0] in ('E', 'F'))-alpha[3])
W = sp.Matrix([[V[labels.index(a), j]*fourth[labels.index(b)]
                - V[labels.index(b), j]*fourth[labels.index(a)]
                for j in range(3)] for a, b in edges])
assert partial*W == sp.zeros(20, 3)
assert Z.row_join(W).rank() == 6
Wbar = (sp.eye(30)-Z*R)*W
C = S*Wbar
expected = sp.Matrix([[0, -1, -1], [1, 0, 1], [-1, -1, 0],
                      [-1, 1, 0], [0, 1, 1], [2, 0, 0],
                      [1, 0, 1], [1, 1, 0]])
assert C == expected and P*C == Wbar
assert R*Wbar == sp.zeros(3, 3) and C.rank() == 3
mass = sp.eye(8)[:, :3]
assert mass.row_join(C).rank() == 6
active = sp.Matrix([3, 2, 1, 0, 0, 0, 0, 0])
assert C.row_join(active).rank() == 4

# Rational frame: first 3 mass directions, 2 complements, last 3 split cups.
options = []
for i, j in combinations(range(3, 8), 2):
    K = mass.row_join(sp.eye(8)[:, [i, j]]).row_join(C)
    if K.det():
        options.append((abs(K.det()), i, j, K))
_, i, j, K = min(options, key=lambda x: x[:3])
Ki = K.inv()
Jproj = C*Ki[5:8, :]
assert Jproj**2 == Jproj
assert Jproj*C == C and Jproj*mass == sp.zeros(8, 3)
assert Ki*active == sp.Matrix([3, 2, 1, 0, 0, 0, 0, 0])
assert K.det() != 0

def strings(M):
    return [[str(x) for x in row] for row in M.tolist()]

output = {
    'status': 'exact_unit_cup_and_rational_frame_checks_passed',
    'scope': 'Finite residues and frames. Derived regulator comparison is proved separately; no nonvanishing computed.',
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'boundary_input_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
    'units': ['a/d', 'b/d', 'c/d', 'M/d'],
    'fourth_unit_valuations': fourth,
    'extra_unit_cup_flags': strings(W),
    'ambient_torus_correction_R_W': strings(R*W),
    'projected_cup_flags': strings(Wbar),
    'projected_cup_coordinates_in_P': strings(C),
    'all_six_cup_rank': 6,
    'residual_cup_rank': 3,
    'mass_plane_intersection_dimension': 0,
    'complement_standard_basis_indices_one_based': [i+1, j+1],
    'new_basis_columns_in_P': strings(K),
    'new_basis_determinant': str(K.det()),
    'new_basis_inverse': strings(Ki),
    'projection_to_split_cups_in_P': strings(Jproj),
    'projection_kills_mass_plane': True,
    'integral_splitting_claimed': False,
    'regulator_nonzero_proved': False,
}
(ROOT/'data/regulator_cup_splitting.json').write_text(
    json.dumps(output, indent=2)+'\n', encoding='utf-8')
print(json.dumps({k: output[k] for k in
                  ('status', 'all_six_cup_rank', 'residual_cup_rank',
                   'mass_plane_intersection_dimension',
                   'complement_standard_basis_indices_one_based',
                   'new_basis_determinant', 'regulator_nonzero_proved')}))
