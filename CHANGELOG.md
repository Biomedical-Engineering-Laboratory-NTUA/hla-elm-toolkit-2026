# Changelog

## 1.1.0 — October 2026

Released alongside the revised manuscript (genes-4557658, Genes/MDPI).

### Changed
- `elm/base_elm.py`: output weights are now solved with the unregularized
  Moore-Penrose pseudo-inverse (`np.linalg.pinv(H) @ T`), matching the
  article's Sections 2.4-2.5 and Algorithm S1. The previous ridge-regularized
  solve contradicted both the article and this package's own README.
  `reg_lambda` is retained for API compatibility and is now unused.
- `elm/base_elm.py`: docstrings now state explicitly where this reference
  implementation departs from the production code — single-pass training
  versus the article's multi-epoch accumulation, and the softmax over all
  compatible diplotypes versus the article's Hardy-Weinberg posterior over
  the score-selected candidate set.
- All READMEs: article cross-references updated to the numbering of the
  revised manuscript (Supplementary Sections S2-S4, Algorithms S1-S5,
  Tables S1-S9).
- Root and package READMEs: scope statement corrected. The previous text
  pointed readers to this same repository for "the original implementation",
  which is circular; the production code is not distributed here.
- `requirements.txt`: dependency versions pinned to the environment used for
  the reported analyses, as the article's Data Availability Statement now
  states.

### Removed
- Internal correspondence from the package README (a discussion of Python 3.2
  compatibility, and a section framed around the manuscript's internal review).

### Added
- `CITATION.cff` and this changelog.
