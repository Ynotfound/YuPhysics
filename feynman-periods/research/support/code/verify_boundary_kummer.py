"""Exact boundary flags and Tate frames; does not certify an extension comparison."""
from itertools import combinations
from pathlib import Path
import hashlib
import json
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
vertices = tuple(range(1, 5))
pairs = tuple(combinations(vertices, 2))
labels = ([f"V{i}" for i in vertices]
          + [f"E{i}{j}" for i, j in pairs]
          + [f"C{i}{j}" for i, j in pairs]
          + [f"F{i}" for i in vertices])
adjacent = set()
for i, j in pairs:
    links = [(f"C{i}{j}", f"V{i}"), (f"C{i}{j}", f"V{j}"),
             (f"C{i}{j}", f"E{i}{j}")]
    links += [(f"E{i}{j}", f"F{k}") for k in vertices if k not in (i, j)]
    adjacent.update(tuple(sorted((labels.index(a), labels.index(b)))) for a, b in links)
edges = sorted(adjacent)
incidence = sp.zeros(20, 30)
for e, (a, b) in enumerate(edges):
    incidence[a, e], incidence[b, e] = -1, 1

valuations = []
for label in labels:
    ids = tuple(int(i) for i in label[1:])
    v = [int(i in ids) if label[0] == "F" else int(i not in ids)
         for i in vertices]
    valuations.append([v[i]-v[3] for i in range(3)])
Z = sp.Matrix([[valuations[a][i]*valuations[b][j]
                - valuations[a][j]*valuations[b][i]
                for i, j in combinations(range(3), 2)] for a, b in edges])
anchor_labels = [("V1", "C12"), ("V1", "C13"), ("V2", "C23")]
anchor_rows = [edges.index((labels.index(a), labels.index(b))) for a, b in anchor_labels]
R = sp.eye(30)[anchor_rows, :]
assert incidence.rank() == 19 and incidence*Z == sp.zeros(20, 3)
assert R*Z == sp.eye(3)
projection = sp.eye(30)-Z*R
assert projection**2 == projection and projection*Z == sp.zeros(30, 3)

# On Fi, j<k<l are the other labels.  Face-first residues at Ekl, Ejl, Ejk
# are epsilon_i*(lambda_l-lambda_k, lambda_j-lambda_l, lambda_k-lambda_j).
# Our ordered edges are E -> F, so the double residue has the opposite sign.
mass_flags = sp.zeros(30, 4)
face_flags = {}
for i in vertices:
    j, k, l = [x for x in vertices if x != i]
    epsilon = (-1)**(i-1)
    residues = [((k, l), l, k), ((j, l), j, l), ((j, k), k, j)]
    face_flags[f"F{i}"] = []
    for (a, b), positive, negative in residues:
        label = f"E{a}{b}"
        row = edges.index((labels.index(label), labels.index(f"F{i}")))
        face_residue = sp.zeros(1, 4)
        face_residue[0, positive-1] = epsilon
        face_residue[0, negative-1] = -epsilon
        mass_flags[row, :] = -face_residue
        face_flags[f"F{i}"].append({"other_component": label,
                                    "face_first_coefficients": list(map(int, face_residue))})

assert incidence*mass_flags == sp.zeros(20, 4)
assert R*mass_flags == sp.zeros(3, 4)
assert mass_flags*sp.ones(4, 1) == sp.zeros(30, 1)
assert mass_flags.rank() == 3
assert sum(any(mass_flags[r, j] for j in range(4)) for r in range(30)) == 12

# Extend the three mass columns to a rational basis of ker(incidence) cap ker(R).
constraint = incidence.col_join(R)
assert constraint.rank() == 22
basis = mass_flags[:, :3]
for column in constraint.nullspace():
    candidate = basis.row_join(column)
    if candidate.rank() > basis.cols:
        basis = candidate
assert basis.cols == 8 and constraint*basis == sp.zeros(23, 8)
coordinate_rows = list(basis.T.rref()[1])
minor = basis[coordinate_rows, :]
assert all(x.q == 1 for x in basis) and abs(minor.det()) == 1
left_inverse = minor.inv()*sp.eye(30)[coordinate_rows, :]
assert left_inverse*basis == sp.eye(8)
coordinates = left_inverse*mass_flags
expected = sp.zeros(8, 4)
expected[:3, :3] = sp.eye(3)
expected[:3, 3] = -sp.ones(3, 1)
assert coordinates == expected and basis*coordinates == mass_flags
cycle_basis = sp.Matrix.hstack(*incidence.nullspace())
assert basis*left_inverse*projection*cycle_basis == projection*cycle_basis

# Independent corruption check: one reversed nonzero flag breaks reciprocity.
bad = mass_flags.copy()
row = next(r for r in range(30) if any(bad[r, j] for j in range(4)))
bad[row, :] = -bad[row, :]
assert incidence*bad != sp.zeros(20, 4)

def strings(matrix):
    return [[str(x) for x in line] for line in matrix.tolist()]

output = {
    "status": "exact_boundary_flag_and_frame_checks_passed",
    "scope": "Finite residue/cycle identities. No equality of the actual relative extension with the symbol regulator is asserted by this computation.",
    "sympy_version": sp.__version__,
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "labels": labels,
    "oriented_edges": [[labels[a], labels[b]] for a, b in edges],
    "orientation": "lower-index component first; tame_{point} tame_{first component}",
    "incidence": strings(incidence),
    "coordinate_unit_valuations": strings(sp.Matrix(valuations)),
    "torus_double_residues_Z": strings(Z),
    "anchor_edges": anchor_labels,
    "anchor_rows_zero_based": anchor_rows,
    "face_first_residues": face_flags,
    "mass_flag_matrix": strings(mass_flags),
    "Tate_basis_columns": strings(basis),
    "Tate_basis_minor_rows_zero_based": coordinate_rows,
    "Tate_basis_minor_determinant": str(minor.det()),
    "Tate_basis_left_inverse": strings(left_inverse),
    "mass_coordinates_in_Tate_basis": strings(coordinates),
    "checks": {"cycle_dimension": 11, "ambient_torus_dimension": 3,
               "relative_Tate_dimension": 8, "mass_flag_rank": 3,
               "nonzero_flag_rows": 12, "common_mass_rescaling_killed": True,
               "sign_corruption_detected": True,
               "actual_extension_comparison_certified": False}
}
(ROOT/'data/boundary_kummer_flags.json').write_text(
    json.dumps(output, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
print(json.dumps({"status": output["status"], "cycle_dimension": 11,
                  "relative_Tate_dimension": 8, "mass_flag_rank": 3,
                  "basis_minor_determinant": str(minor.det()),
                  "actual_extension_comparison_certified": False}))
