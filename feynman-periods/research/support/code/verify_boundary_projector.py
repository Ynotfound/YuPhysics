"""Rational coefficients for the boundary Chow projector, with relation witnesses."""
from itertools import combinations
from pathlib import Path
import json
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
vertices = tuple(range(1, 5))
pairs = tuple(combinations(vertices, 2))
labels = ([f"V{i}" for i in vertices]
          + [f"E{i}{j}" for i, j in pairs]
          + [f"C{i}{j}" for i, j in pairs]
          + [f"F{i}" for i in vertices])
Q = -2*sp.eye(20)
for i, j in pairs:
    adjacencies = [(f"C{i}{j}", f"V{i}"), (f"C{i}{j}", f"V{j}"),
                   (f"C{i}{j}", f"E{i}{j}")]
    adjacencies += [(f"E{i}{j}", f"F{k}") for k in vertices if k not in (i,j)]
    for a, b in adjacencies:
        Q[labels.index(a), labels.index(b)] = 1
        Q[labels.index(b), labels.index(a)] = 1

omitted = {"E12", "F2", "C12", "F3"}
basis_indices = [i for i, label in enumerate(labels) if label not in omitted]
S = sp.eye(20)[:, basis_indices]
G = S.T*Q*S
assert G.det() == -48
Ginv = G.inv()
coordinates = Ginv*S.T*Q
assert all(x.q == 1 for x in coordinates)
assert coordinates*S == sp.eye(16)
assert coordinates.T*G*coordinates == Q

# Matrix of coefficients of sum C_ij D_i x D_j, on the full list of divisors.
C = S*Ginv*S.T
assert C == C.T and C*Q*C == C and Q*C*Q == Q
presentation_projector = C*Q
assert presentation_projector**2 == presentation_projector
assert presentation_projector.rank() == 16

# Principal divisors of a/d, b/d, c/d, and M/d.
rows = []
for label in labels:
    ids = tuple(int(i) for i in label[1:])
    kind = label[0]
    valuations = [int(i in ids) if kind == "F" else int(i not in ids)
                  for i in vertices]
    m_valuation = int(kind in ("E", "F"))
    rows.append([valuations[i]-valuations[3] for i in range(3)]
                + [m_valuation-valuations[3]])
V = sp.Matrix(rows)
assert V.rank() == 4 and Q*V == sp.zeros(20,4)
relations = sp.eye(20)-S*coordinates
witnesses, parameters = V.gauss_jordan_solve(relations)
assert parameters.rows == 0 and V*witnesses == relations
assert presentation_projector*V == sp.zeros(20,4)

def matrix_strings(matrix):
    return [[str(x) for x in row] for row in matrix.tolist()]

output = {
    "status": "exact_boundary_projector_identities_passed",
    "scope": "Boundary correspondence coefficients and principal-relation witnesses only; no numerical periods or Picard-rank certificate.",
    "labels": labels,
    "basis_labels": [labels[i] for i in basis_indices],
    "gram_determinant": str(G.det()),
    "gram": matrix_strings(G),
    "gram_inverse": matrix_strings(Ginv),
    "inverse_denominator_lcm": int(sp.ilcm(*[int(x.q) for x in Ginv])),
    "all_divisors_basis_coordinates": matrix_strings(coordinates),
    "projector_correspondence_coefficients": matrix_strings(C),
    "principal_unit_order": ["a/d", "b/d", "c/d", "M/d"],
    "principal_divisor_matrix": matrix_strings(V),
    "relation_witnesses": matrix_strings(witnesses),
    "checks": {"C_Q_C_equals_C": True, "Q_C_Q_equals_Q": True,
               "all_relations_are_principal_combinations": True,
               "presentation_projector_rank": 16,
               "presentation_kernel_is_relation_space_dimension": 4},
    "warning": "The 4-dimensional presentation kernel is NOT the 6-dimensional K3 summand."
}
(ROOT/'data/boundary_projector.json').write_text(
    json.dumps(output, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
print(json.dumps({"status": output["status"], "basis_size":16,
                  "gram_determinant":str(G.det()),
                  "inverse_denominator_lcm":output['inverse_denominator_lcm'],
                  "principal_relations_dimension":4}))
