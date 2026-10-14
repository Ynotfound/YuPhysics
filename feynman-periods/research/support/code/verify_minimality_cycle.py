"""Check divisors of the proposed CH^2(X,1) cycle, not regulator nonvanishing."""
from pathlib import Path
import hashlib
import json
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
source = ROOT / 'data/boundary_kummer_flags.json'
data = json.loads(source.read_text(encoding='utf-8'))
P = sp.Matrix(data['Tate_basis_columns'])
incidence = sp.Matrix(data['incidence'])
edges = [tuple(e) for e in data['oriented_edges']]
labels = data['labels']
c = P * sp.Matrix([3, 2, 1, 0, 0, 0, 0, 0])
assert incidence*c == sp.zeros(20, 1)
faces = {1: (-1, 2), 2: (1, -3), 3: (-2, 3), 4: (1, -2)}
edge_powers = {'E12': 1, 'E13': -2, 'E14': 3, 'E23': 1, 'E24': -2, 'E34': 1}
divisors = {label: {} for label in labels}
for i, (zero, one) in faces.items():
    j, k, l = [a for a in range(1, 5) if a != i]
    divisors[f'F{i}'] = {f'E{k}{l}': zero, f'E{j}{l}': one,
                         f'E{j}{k}': -zero-one}
for label, power in edge_powers.items():
    k, l = [a for a in range(1, 5) if str(a) not in label[1:]]
    divisors[label] = {f'F{k}': power, f'F{l}': -power}
flags = []
for e, (a, b) in enumerate(edges):
    da, db = divisors[a].get(b, 0), divisors[b].get(a, 0)
    assert da == int(-incidence[labels.index(a), e]*c[e])
    assert db == int(-incidence[labels.index(b), e]*c[e])
    assert da+db == 0
    flags.append({'edge': [a, b], 'cycle_coefficient': int(c[e]),
                  'divisor_on_first': da, 'divisor_on_second': db})
assert all(sum(ds.values()) == 0 for ds in divisors.values())

# Reversing a single nonzero divisor without its partner must break closure.
bad = [f['divisor_on_first']+f['divisor_on_second'] for f in flags]
idx = next(i for i, f in enumerate(flags) if f['divisor_on_first'])
bad[idx] -= 2*flags[idx]['divisor_on_first']
assert any(bad)

result = {'status': 'divisor_cancellation_verified',
          'scope': 'A cycle representative, not a proof its class or regulator is nonzero.',
          'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'boundary_input_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
          'direction_in_P_basis': [3, 2, 1, 0, 0, 0, 0, 0],
          'face_exponents_u_and_one_minus_u': faces,
          'edge_exponents_v': edge_powers, 'flags': flags,
          'number_of_nonzero_flags': sum(f['cycle_coefficient'] != 0 for f in flags),
          'regulator_comparison_certified_by_this_script': False,
          'regulator_comparison_proof_draft': '../REGULATOR_COMPARISON.md',
          'very_general_nonvanishing_proof_draft': '../DE_RHAM_MINIMALITY.md',
          'regulator_nonzero_proved_at_s_minus_64': False,
          'rank_eight_minimality_proved': False, 'sign_corruption_detected': True}
(ROOT/'data/minimality_cycle.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
print(json.dumps({k: v for k, v in result.items() if k not in
                  ('flags', 'face_exponents_u_and_one_minus_u', 'edge_exponents_v')}))
