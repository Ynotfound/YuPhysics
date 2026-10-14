# Feynman periods: research collection

**Yu Zhehao — School of Physics, Peking University**

Current research snapshot: **U41, 14 October 2026**. The documents are predominantly in Chinese, with mathematical notation and English references. They are preserved as exported research drafts, with their existing title pages and version-specific statements.

The **329-page unified research volume** is the comprehensive, survey-style document. The **13-page classification table** is the shortest route into the results. The **191-page structural manuscript** develops the main theoretical framework. The remaining documents provide focused proofs and audits; they overlap with material integrated into the unified volume.

The [English long paper](../english/feynman_periods_long_paper_en.pdf) remains the 130-page editorial version of 13 October. It has **not yet incorporated all U36–U41 additions**. The current research collection and the submission-oriented English paper are separate editions.

## Manuscripts

| Document | Files | Reading guide |
| --- | --- | --- |
| Structures, Visibility and Physical Realizations of Feynman Periods | [PDF (329 pp.)](manuscripts/feynman_periods_research_unified.pdf) · [LaTeX](manuscripts/feynman_periods_research_unified.tex) | Current comprehensive research volume; begin here for the full survey-style treatment. |
| Graph-family classification | [PDF (13 pp.)](manuscripts/GRAPH_FAMILY_CLASSIFICATION.pdf) · [LaTeX](manuscripts/GRAPH_FAMILY_CLASSIFICATION.tex) | 23 classification rows with hypotheses, evidence and open cases. |
| Unified structural framework | [PDF (191 pp.)](manuscripts/UNIFIED_STRUCTURE_THEOREM.pdf) · [LaTeX](manuscripts/UNIFIED_STRUCTURE_THEOREM.tex) | Current theoretical main line, including the U41 additions. |
| Massless graphs and MZV prediction | [PDF (6 pp.)](manuscripts/MASSLESS_MZV_PREDICTION.pdf) · [LaTeX](manuscripts/MASSLESS_MZV_PREDICTION.tex) | Series-parallel sufficient criteria, banana formulae and physical readouts. |
| Geometric comparison and Hypothesis 11.2 | [PDF (9 pp.)](manuscripts/HYPOTHESIS_11_2_AUDIT.pdf) · [LaTeX](manuscripts/HYPOTHESIS_11_2_AUDIT.tex) | Relative objects, physical framings and the scope of the fixed-fiber conclusions. |
| Mass subdivision and confluent classification | [PDF (5 pp.)](manuscripts/SUBDIVISION_COMPARISON_AUDIT.pdf) · [LaTeX](manuscripts/SUBDIVISION_COMPARISON_AUDIT.tex) | Repeated-mass derivative rules, counterexamples and direct literature comparisons. |
| Physical certification and literature audit | [PDF (9 pp.)](manuscripts/PHYSICS_CERTIFICATION_LITERATURE_AUDIT.pdf) · [LaTeX](manuscripts/PHYSICS_CERTIFICATION_LITERATURE_AUDIT.tex) | Angular truncation bounds, the remaining integration problem and source comparison. |
| Differential visibility of the banana | [PDF (9 pp.)](manuscripts/BANANA_DIFFERENTIAL_VISIBILITY.pdf) · [LaTeX](manuscripts/BANANA_DIFFERENTIAL_VISIBILITY.tex) | Parameter derivatives and visible geometric directions. |
| Operational motivic interface | [PDF (10 pp.)](manuscripts/OPERATIONAL_MOTIVIC_INTERFACE.pdf) · [LaTeX](manuscripts/OPERATIONAL_MOTIVIC_INTERFACE.tex) | Actual maps, physical operations and compatible minimalization. |
| Differential descent for forests | [PDF (8 pp.)](manuscripts/FOREST_DIFFERENTIAL_DESCENT.pdf) · [LaTeX](manuscripts/FOREST_DIFFERENTIAL_DESCENT.tex) | Endpoint descent and the scope of mass-insertion rules. |
| Finite MOM comparison | [PDF (8 pp.)](manuscripts/FINITE_MOM_MOTIVIC_COMPARISON.pdf) · [LaTeX](manuscripts/FINITE_MOM_MOTIVIC_COMPARISON.tex) | Compatible subtraction points and formal-period scheme relations. |
| Quadratic subtraction of the sunset | [PDF (8 pp.)](manuscripts/SUNSET_JET_MOTIVE.pdf) · [LaTeX](manuscripts/SUNSET_JET_MOTIVE.tex) | Taylor subtraction, relative objects and endpoint descent. |
| Two-loop MOM renormalization group | [PDF (7 pp.)](manuscripts/MOM_RG_ENDPOINTS.pdf) · [LaTeX](manuscripts/MOM_RG_ENDPOINTS.tex) | Finite-mass endpoint kernels and ultraviolet limits. |
| Earlier research compilation | [PDF (97 pp.)](manuscripts/feynman_periods_research.pdf) · [LaTeX](manuscripts/feynman_periods_research.tex) | Historical 13 October compilation; superseded by the current unified volume. |

## Current results and open questions

The snapshot adds explicit MZV sufficient criteria for accepted massless series-parallel graphs; a graph-family classification table with 23 rows and 57 main-text references; a geometric and physical-framing comparison for the specified banana model; and repeated-mass derivative rules with a counterexample to extending zero-order grid tests across mass collisions.

The following qualifications remain essential:

- Membership of a numerical value in the MZV algebra is distinct from a motivic or differential classification. No general numerical non-MZV test is claimed.
- The banana comparison has a specified geometric domain. Fixed-fiber minimality and group statements retain their separate geometric and arithmetic certificate inputs.
- The confluent-mass rule provides a finite derivative representation and an elliptic-carrier upper bound. General minimality and cancellation in the full sum remain to be established.
- At `s/m^2 = 20`, the bound below `6 × 10^-17` controls the angular truncation contribution only. The remaining integration, a narrow interval for the full physical coefficient, and an independent three-loop forward optical check remain open.
- The literature comparison identifies prior results and candidate differences; it is not an exhaustive novelty or priority certification.

Earlier standalone audits may still describe their integration into the unified volume as pending. That editorial integration is complete in U41; the substantive open problems above remain.

## Supporting material

[Proof notes, research code and certificate archive](support/README.md) accompany the manuscripts. The archive preserves 840 JSON files (about 242 MiB uncompressed), with per-file checksums. The Python scripts and scientific notes are separately browsable. Operational chat ledgers and local application configuration are outside this collection.

[manifest.json](manifest.json) records the published source and PDF hashes, PDF page counts, and supporting-file identities. Source paths in this manifest are relative to the research workspace, not dependencies required to compile the manuscripts.

## Building a manuscript

Every `.tex` file in `manuscripts/` is self-contained, with its bibliography embedded. Use XeLaTeX or Tectonic with Chinese TeX support. For example:

```sh
cd manuscripts
tectonic feynman_periods_research_unified.tex
tectonic GRAPH_FAMILY_CLASSIFICATION.tex
```

If using XeLaTeX, repeat compilation until references and the table of contents stabilize. The structural manuscript, classification table and unified volume were re-exported from U41 sources whose hashes match the completed U41 receipt, because the live workspace had already begun U42. The other PDFs are unchanged copies of their compiled editions. Publication and compilation checks are not mathematical recertification.
