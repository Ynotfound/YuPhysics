"""Independent finite algebra checks of transcribed theory inputs; no transport."""
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import re
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

ROOT = Path(__file__).resolve().parents[1]


def snf(matrix):
    return [abs(int(x)) for x in smith_normal_form(matrix, domain=sp.ZZ).diagonal()]


def integral(matrix):
    return all(x.q == 1 for x in matrix)


def verify_lattice(data):
    q = sp.Matrix(data["Q6"])
    source = Path(data["sources"][0]["path"]).read_text(encoding="utf-8")
    block = source.split(r"\label{eq:Qcandidate}", 1)[1].split(r"\end{pmatrix}", 1)[0]
    block = block.split(r"\begin{pmatrix}", 1)[1]
    rows = [list(map(int, re.findall(r"-?\d+", row))) for row in block.split(r"\\") if row.strip()]
    assert q == sp.Matrix(rows), "PDF transcription differs from TeX"
    assert q == q.T and q.det() == 192
    assert snf(q) == [2, 2, 2, 2, 2, 6]
    identity = sp.eye(6)
    matrices = list(map(sp.Matrix, data["mass_matrices"]))
    t0 = sp.Matrix(data["T0"])
    for m in matrices + [t0]:
        assert m.T*q*m == q and (m-identity)**2 == sp.zeros(6)
    for a, b in combinations(matrices, 2):
        assert a*b == b*a
    isotropic, fixed = [], []
    for bits in product(range(2), repeat=6):
        if not any(bits):
            continue
        v = sp.Matrix(bits)
        if (v.T*q*v)[0] % 8 == 0:
            isotropic.append(list(bits))
        if all(all(x % 2 == 0 for x in (m-identity)*v) for m in matrices):
            fixed.append(list(bits))
    assert len(isotropic) == 11 and fixed == [data["half_cycle_numerator"]]
    assert fixed[0] in isotropic
    c = sp.eye(6)
    c[:, 0] = sp.Matrix(fixed[0])/2
    h = c.T*q*c
    p = sp.Matrix(data["standard_isometry"])
    u = sp.Matrix([[0,1],[1,0]])
    standard = sp.diag(u, 2*u, sp.Matrix([[-4,2],[2,-4]]))
    assert c.det() == sp.Rational(1,2) and h.det() == 48
    assert p.det() == 1 and p.T*h*p == standard
    assert snf(h) == [1,1,2,2,2,6]
    assert all(integral(c.inv()*m*c) for m in matrices + [t0])
    reflectors = [identity[:, i]*identity[:, i].T*q for i in range(6)]
    for r in reflectors:
        assert r*r == -2*r and (identity+r)**2 == identity
    pairs = [{"indices": [i+1,j+1], "target": int(sp.trace(reflectors[i]*reflectors[j]))}
             for i,j in combinations(range(6),2)]
    triples = [{"indices": [1,i+1,j+1], "target": int(sp.trace(reflectors[0]*reflectors[i]*reflectors[j]))}
               for i,j in combinations(range(1,6),2) if q[i,j] != 0]
    assert [x["target"] for x in triples] == [-16,-16,-16,-16,-64,-16]
    return {"det_Q6": 192, "snf_Q6": snf(q), "signature_by_exact_isometry": [2,4],
            "isotropic_order_two_classes": len(isotropic), "unique_mass_fixed_nonzero_class": fixed[0],
            "det_completed_lattice": 48, "snf_completed_lattice": snf(h),
            "standard_isometry_verified": True, "mass_and_T0_identities_verified": True,
            "pair_trace_targets": pairs, "triple_trace_targets": triples}


def verify_boundary():
    vertices = range(1,5)
    pairs = list(combinations(vertices,2))
    labels = [("V",(i,)) for i in vertices] + [("E",p) for p in pairs] + [("C",p) for p in pairs] + [("F",(i,)) for i in vertices]
    q = -2*sp.eye(20)
    for p in pairs:
        ci, ei = labels.index(("C",p)), labels.index(("E",p))
        edges = [(ci, labels.index(("V",(i,)))) for i in p] + [(ci,ei)]
        edges += [(ei,labels.index(("F",(i,)))) for i in vertices if i not in p]
        for i,j in edges:
            q[i,j] = q[j,i] = 1
    edges = [(i,j) for i,j in combinations(range(20),2) if q[i,j] == 1]
    incidence = sp.zeros(20,len(edges))
    for k,(i,j) in enumerate(edges):
        incidence[i,k], incidence[j,k] = -1,1
    vals = []
    for kind, ids in labels:
        v = [int((i not in ids) if kind in ("V","C","E") else (i in ids)) for i in vertices]
        v_m = int(kind in ("E","F"))
        vals.append([v[i]-v[3] for i in range(3)] + [v_m-v[3]])
    valuation = sp.Matrix(vals)
    assert q*valuation == sp.zeros(20,4) and valuation.rank() == 4
    assert q.rank() == 16 and len(edges) == 30 and incidence.rank() == 19
    residues = sp.Matrix([[vals[i][p]*vals[j][r]-vals[i][r]*vals[j][p]
                          for p,r in combinations(range(3),2)] for i,j in edges])
    assert incidence*residues == sp.zeros(20,3) and residues.rank() == 3
    anchors = [(("V",(1,)),("C",(1,2))), (("V",(1,)),("C",(1,3))), (("V",(2,)),("C",(2,3)))]
    anchor_rows = [edges.index((labels.index(a),labels.index(b))) for a,b in anchors]
    assert residues[anchor_rows,:] == sp.eye(3)
    omitted = [("E",(1,2)),("F",(2,)),("C",(1,2)),("F",(3,))]
    basis = [i for i,label in enumerate(labels) if label not in omitted]
    gram = q.extract(basis,basis)
    coordinates = gram.inv()*q.extract(basis,range(20))
    assert gram.det() == -48 and integral(coordinates)
    assert coordinates.T*gram*coordinates == q
    return {"boundary_components":20, "crossings":30, "incidence_rank":19,
            "dual_graph_b1":11, "divisor_lattice_rank":16, "divisor_basis_det":-48,
            "valuation_relations_rank":4, "torus_residue_rank":3,
            "residue_anchor_minor":"identity_3", "remaining_weight_four_rank":8,
            "ambient_relative_rank_from_Hodge_argument":15}


def main():
    input_path = ROOT/"data/recovered_lattice.json"
    data = json.loads(input_path.read_text(encoding="utf-8"))
    result = {"status":"exact_algebra_passed", "sympy_version":sp.__version__,
              "scope":"Finite identities only; no operator validation, global trace enclosure, motivic comparison or geometric marking.",
              "input_sha256":hashlib.sha256(input_path.read_bytes()).hexdigest(),
              "lattice":verify_lattice(data), "boundary":verify_boundary()}
    (ROOT/"data/recovered_algebra_checks.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"det_Q6":192,"det_completion":48,"boundary":result["boundary"]}))


if __name__ == "__main__":
    main()
