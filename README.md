# YuPhysics

Research manuscripts by **Yu Zhehao**, School of Physics, Peking University.

## Feynman periods

**Visible Structures and Renormalization of Feynman Periods**

The paper studies the geometric structures visible in Feynman periods and their behavior under cuts, parameter differentiation and renormalization. The four-propagator banana, forest endpoint descent, the quadratically subtracted sunset and a finite-mass two-loop momentum-subtraction renormalization group provide its principal examples.

The English manuscript is the primary edition. The Chinese edition is retained for comparison. The main text and supplement form one paper.

| Edition | Main paper | Supplement |
| --- | --- | --- |
| English | [PDF, 130 pages](feynman-periods/english/feynman_periods_long_paper_en.pdf) · [LaTeX](feynman-periods/english/feynman_periods_long_paper_en.tex) | [PDF, 45 pages](feynman-periods/english/feynman_periods_supplement_en.pdf) · [LaTeX](feynman-periods/english/feynman_periods_supplement_en.tex) |
| 中文对照版 | [PDF，121 页](feynman-periods/chinese/feynman_periods_long_paper.pdf) · [LaTeX](feynman-periods/chinese/feynman_periods_long_paper.tex) | [PDF，43 页](feynman-periods/chinese/feynman_periods_supplement.pdf) · [LaTeX](feynman-periods/chinese/feynman_periods_supplement.tex) |

**Manuscript version:** 13 October 2026. **Repository upload:** 14 October 2026.

## Build

The four `.tex` files are standalone documents, with their references included in the sources. Precompiled PDFs are provided above.

With Tectonic, build the English edition from the repository root:

```sh
cd feynman-periods/english
tectonic feynman_periods_long_paper_en.tex
tectonic feynman_periods_supplement_en.tex
```

Alternatively, use XeLaTeX and rerun until cross-references stabilize. The English edition uses `article`, Latin Modern and standard AMS/TikZ packages. The Chinese edition uses `ctexart` and requires a suitable Chinese TeX installation. Keep each main PDF and its supplement in the same directory so that their relative links work.

## Research status

This is a research draft. Its geometric assumptions and inherited finite certificates are stated in the manuscript. The fixed-fiber numerical certificates have not been independently rerun for this edition; narrow rigorous error bounds for the physical benchmark remain open. This repository contains the manuscripts and their sources, rather than the computational certificate archive.

No arXiv upload or journal submission has been made for this edition as of 14 October 2026.
