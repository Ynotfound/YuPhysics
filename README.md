# YuPhysics

Research manuscripts by **Yu Zhehao**, School of Physics, Peking University.

## Feynman periods

**Visible Structures and Renormalization of Feynman Periods**

The paper studies the geometric structures visible in Feynman periods and their behavior under cuts, parameter differentiation and renormalization. The four-propagator banana, forest endpoint descent, the quadratically subtracted sunset and a finite-mass two-loop momentum-subtraction renormalization group provide its principal examples.

The English manuscript is the submission-oriented long-paper edition. The Chinese edition is retained for comparison. The main text and supplement form one paper. These editions date from 13 October and do not yet incorporate all the later research updates linked below.

| Edition | Main paper | Supplement |
| --- | --- | --- |
| English | [PDF, 130 pages](feynman-periods/english/feynman_periods_long_paper_en.pdf) · [LaTeX](feynman-periods/english/feynman_periods_long_paper_en.tex) | [PDF, 45 pages](feynman-periods/english/feynman_periods_supplement_en.pdf) · [LaTeX](feynman-periods/english/feynman_periods_supplement_en.tex) |
| 中文对照版 | [PDF，121 页](feynman-periods/chinese/feynman_periods_long_paper.pdf) · [LaTeX](feynman-periods/chinese/feynman_periods_long_paper.tex) | [PDF，43 页](feynman-periods/chinese/feynman_periods_supplement.pdf) · [LaTeX](feynman-periods/chinese/feynman_periods_supplement.tex) |

**Manuscript version:** 13 October 2026. **Repository upload:** 14 October 2026.

## Current research collection / 最新研究材料

The completed **U41 research snapshot (14 October 2026)** is now available with 14 PDFs, their standalone LaTeX sources, supporting proof notes, research code and a compressed certificate archive.

- **[Comprehensive research volume / 统一研究总稿 — 329 pages](feynman-periods/research/manuscripts/feynman_periods_research_unified.pdf)**: the long, survey-style treatment of the project.
- **[Graph-family classification / 图族分类表 — 13 pages](feynman-periods/research/manuscripts/GRAPH_FAMILY_CLASSIFICATION.pdf)**: 23 rows with assumptions, evidence, physical visibility and outstanding questions.
- **[Unified structural framework / 统一理论主线 — 191 pages](feynman-periods/research/manuscripts/UNIFIED_STRUCTURE_THEOREM.pdf)**: the main theoretical development.
- **[All research manuscripts and sources / 全部专题文稿及源码](feynman-periods/research/README.md)**: MZV prediction, the geometric comparison, mass collisions, forests, MOM, the sunset and physical audits.
- **[Proof notes, code and certificates / 证明记录、代码与证书](feynman-periods/research/support/README.md)**: scientific notes, 107 Python scripts and an archive containing 840 JSON data files.

The research collection is predominantly Chinese and is broader than the English long paper. New research results will require further review and editorial integration before entering that English edition. The source and PDF snapshot is identified in the [publication manifest](feynman-periods/research/manifest.json).

## Build

The long-paper `.tex` files are standalone documents, with their references included in the sources. The research collection also provides standalone sources and [its own build instructions](feynman-periods/research/README.md#building-a-manuscript). Precompiled PDFs are provided above.

With Tectonic, build the English edition from the repository root:

```sh
cd feynman-periods/english
tectonic feynman_periods_long_paper_en.tex
tectonic feynman_periods_supplement_en.tex
```

Alternatively, use XeLaTeX and rerun until cross-references stabilize. The English edition uses `article`, Latin Modern and standard AMS/TikZ packages. The Chinese edition uses `ctexart` and requires a suitable Chinese TeX installation. Keep each main PDF and its supplement in the same directory so that their relative links work.

## Research status

These are research drafts. Geometric assumptions and finite-certificate inputs are stated in the manuscripts. Archived code and certificate data are now included, but this upload does not rerun those calculations or independently recertify the mathematics. Narrow rigorous bounds for the full physical benchmark, an independent same-order three-loop forward optical check, and a general numerical non-MZV test remain open. The angular truncation bound in the latest audit is not a narrow bound on the full cross section.

No arXiv upload or journal submission has been made for this edition as of 14 October 2026.
