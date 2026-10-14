# Supporting proofs, code and certificates

Research snapshot: **U41, 14 October 2026**. These files accompany the [research manuscripts](../README.md).

## Contents

- Scientific proof notes and source audits are available below in their original Markdown form.
- `code/` contains the original project Python scripts, including certificate generation, verification and document-update utilities.
- [certificates.zip](certificates.zip) contains all 840 JSON files from the research `data/` directory, with the original `data/` paths preserved.
- [certificate-index.json](certificate-index.json) records each uncompressed file's size and SHA-256 hash.

Extract the certificate archive **into this support directory**, so that `code/` and `data/` are siblings. For example:

```sh
python -m zipfile -e certificates.zip .
```

## Scope and reproduction

The code and certificate bytes are preserved, rather than rewritten during publication. This upload checks file identity and archive integrity; it does not rerun the research calculations or independently establish their mathematical claims.

The scripts are archived research code, not a packaged software release. Depending on the entry point, they use SymPy, NumPy, mpmath, SciPy, python-flint and pypdf. Several scripts also refer to the original local runtime, wheel hashes, source PDFs or downloaded public datasets. Those dependencies need to be supplied or their locations configured before replay. Runtime binaries, downloaded literature PDFs, caches and temporary build directories are not included. The scripts and saved JSON record the relevant versions, locations and provenance where available; consult [SOURCE_AUDIT.md](SOURCE_AUDIT.md) and [BANANA_PAIRING.md](BANANA_PAIRING.md) for external inputs.

Some scripts generate data, download sources or update notes. Read the intended entry point before running it. Historical paths and historical status statements remain part of the preserved evidence. For current conclusions and qualifications, use the unified U41 volume and classification table rather than an isolated old status field.

## Scientific notes

- [ABSOLUTE_COMPARISON.md](ABSOLUTE_COMPARISON.md)
- [BANANA_PAIRING.md](BANANA_PAIRING.md)
- [BASED_NODE_TRANSPORT.md](BASED_NODE_TRANSPORT.md)
- [BOUNDED_GRAPH_SUPPORT.md](BOUNDED_GRAPH_SUPPORT.md)
- [CONTRACTION_FOREST_GEOMETRY.md](CONTRACTION_FOREST_GEOMETRY.md)
- [CUT_RESIDUE_VISIBILITY.md](CUT_RESIDUE_VISIBILITY.md)
- [DE_RHAM_MINIMALITY.md](DE_RHAM_MINIMALITY.md)
- [DEGENERATE_EXTENSIONS.md](DEGENERATE_EXTENSIONS.md)
- [FIXED_PURE_GROUP.md](FIXED_PURE_GROUP.md)
- [FOREST_RENORMALIZED_MOTIVE.md](FOREST_RENORMALIZED_MOTIVE.md)
- [FOUR_PARAMETER_CONNECTION.md](FOUR_PARAMETER_CONNECTION.md)
- [FROBENIUS_REFLECTION.md](FROBENIUS_REFLECTION.md)
- [GEOMETRIC_COMPARISON.md](GEOMETRIC_COMPARISON.md)
- [GEOMETRIC_MONODROMY_GROUP.md](GEOMETRIC_MONODROMY_GROUP.md)
- [HIGH_TRUNK_TRANSPORT.md](HIGH_TRUNK_TRANSPORT.md)
- [HODGE.md](HODGE.md)
- [HORIZONTAL_TRANSPORT.md](HORIZONTAL_TRANSPORT.md)
- [K3_KERNEL.md](K3_KERNEL.md)
- [K3_SATURATION_INTAKE.md](K3_SATURATION_INTAKE.md)
- [KUMMER_BOUNDARY.md](KUMMER_BOUNDARY.md)
- [LATTICE_SATURATION.md](LATTICE_SATURATION.md)
- [LOWER_COMPARISON.md](LOWER_COMPARISON.md)
- [LOWER_REGULATOR.md](LOWER_REGULATOR.md)
- [MARKED_LATTICE.md](MARKED_LATTICE.md)
- [MARKED_SATURATION.md](MARKED_SATURATION.md)
- [MINIMALITY.md](MINIMALITY.md)
- [MIXED_PARAMETER_CONNECTION.md](MIXED_PARAMETER_CONNECTION.md)
- [MIXED_UNIPOTENT_GROUP.md](MIXED_UNIPOTENT_GROUP.md)
- [NODE_REFLECTION_CERTIFICATES.md](NODE_REFLECTION_CERTIFICATES.md)
- [PHI4_BENCHMARK.md](PHI4_BENCHMARK.md)
- [PHI4_FOUR_BODY.md](PHI4_FOUR_BODY.md)
- [PICARD_CERTIFICATE.md](PICARD_CERTIFICATE.md)
- [POLARIZED_TRANSPORT.md](POLARIZED_TRANSPORT.md)
- [PRIMITIVE_RENORMALIZED_MOTIVE.md](PRIMITIVE_RENORMALIZED_MOTIVE.md)
- [PURE_CONNECTION.md](PURE_CONNECTION.md)
- [PURE_RAY_SYSTEM.md](PURE_RAY_SYSTEM.md)
- [Q6.md](Q6.md)
- [Q6_GEOMETRIC_CERTIFICATE.md](Q6_GEOMETRIC_CERTIFICATE.md)
- [REGULATOR_COMPARISON.md](REGULATOR_COMPARISON.md)
- [REGULATOR_RANK_AND_GROUP.md](REGULATOR_RANK_AND_GROUP.md)
- [RELATIVE_TOP_AUDIT.md](RELATIVE_TOP_AUDIT.md)
- [RELATIVE_TOP_COMPARISON.md](RELATIVE_TOP_COMPARISON.md)
- [RESEARCH_OVERVIEW_20261013.md](RESEARCH_OVERVIEW_20261013.md)
- [RG_CRITICAL_VISIBILITY.md](RG_CRITICAL_VISIBILITY.md)
- [SOURCE_AUDIT.md](SOURCE_AUDIT.md)
- [THEORY_MANUSCRIPT.md](THEORY_MANUSCRIPT.md)
- [VISIBILITY_AND_CORADICAL.md](VISIBILITY_AND_CORADICAL.md)
- [ZERO_MOMENTUM_FOREST_MOTIVE.md](ZERO_MOMENTUM_FOREST_MOTIVE.md)
- [ZERO_MONODROMY.md](ZERO_MONODROMY.md)
