# hla-elm-toolkit

Reference Python implementations of **every method compared** in:

> Kepentzis, S.; Chatzistamatiou, T.; Digalakis, J.; Petropoulou, O.;
> Matsopoulos, G.K.; Koutsouris, D. *Improving the Usability of Donor Data from Registries Using
> an Extreme Learning Machine Approach to Upgrade Low/Mid- to
> High-Resolution HLA Data.* Genes (MDPI),
> manuscript ID genes-4557658, under review.

This package implements the base Extreme Learning Machine (ELM) and its
three extensions (KELM, WELM, Ensemble ELM), the three EM/Bayesian
comparator baselines the article specifies in pseudocode (GRIMM-style,
HaploStats-style, and Hapl-o-Mat-style; Supplementary Section S2,
Algorithms S3-S5),
and the MLP / gradient-boosted-tree comparators used in the article's
framework-level benchmarking (Supplementary Section S3).

**Python version.** This package targets **Python 3.9+** and is
developed and tested on **Python 3.12**, the version used for the
article's own implementation (Section 2.7, "Titanas" server
specification).

---

## 1. Scope and honesty about limitations

This is a **from-scratch reference implementation** of the published
algorithms, written to accompany the article. It is **not** the
registry-scale production code that produced the article's reported
figures, and it is not expected to reproduce them. The production
implementation is bound to the registry data environment and is not
distributed here; requests relating to it should be addressed to the
corresponding author. Three things follow:

1. **No real registry data is included or reachable.** The HTO, ORAM,
   and GRPT donor registries used in the article are not distributable
   (donor privacy; Hellenic Transplant Organization data-sharing
   restrictions). `hla_elm_toolkit.data.make_synthetic_population()`
   generates a small synthetic reference population with a Zipf-like,
   rare-allele-dominated frequency spectrum (loosely mirroring the
   pattern quantified in the article's Tables 4 and 5) purely for
   demonstration and unit testing. **Do not** interpret any
   accuracy/call-rate number produced by this package's demo or tests as
   comparable to the article's reported figures.

2. **The base ELM's domain-structured topology is simplified.** The
   article's network activates input nodes for observed
   alleles/genotype fragments, then propagates through three
   hierarchical hidden layers (HL-1: three-locus partial haplotypes,
   HL-2: five-locus partial haplotypes, HL-3: complete haplotypes) whose
   connectivity is fixed in advance by which partial haplotypes are
   supersets of which (Section 2.4). Reproducing that exact multi-stage
   topology from the article's prose description alone, without the
   original code, is out of scope for this package. `elm/base_elm.py`
   instead implements the same *core ELM mathematics* — randomly
   generated, untrained input-to-hidden weights and a closed-form
   Moore-Penrose-pseudo-inverse solution for the output weights — over a
   **single random hidden layer** applied to a flat multi-hot
   locus/allele input encoding. This preserves the algorithm's defining
   property (no iterative back-propagation) and is sufficient to
   benchmark KELM/WELM/Ensemble/comparators against each other in a
   like-for-like way, but it is not a byte-for-byte reproduction of the
   original network, and its accuracy is not expected to match the
   article's reported figures.

3. **GRIMM-, HaploStats-, and Hapl-o-Mat-style baselines are the
   article's own documented simplifications**, not the published tools.
   The article is explicit about this (Section 2.6 and Supplementary Section S2): "these pseudocode
   summaries omit implementation-specific details (e.g., GRIMM's exact
   graph construction and traversal optimizations, HaploStats'
   race/ethnicity-specific reference tables, and Hapl-o-Mat's
   population-weighting options) ... based on their published method
   descriptions ... rather than on proprietary source code, which we did
   not have access to for any of the three tools." This package
   implements exactly those pseudocode algorithms (Algorithms 3-5) as
   given in the article, faithfully, but they remain approximations of
   the real tools by the article's own account, and running the
   published Hapl-o-Mat tool or the NMDP HaploStats web service itself
   is explicitly listed as outstanding future work in the article
   (Supplementary Sections S4.3-S4.4 and main text Section 5).

The exact production implementation is not available through this
package; see the article's Data Availability Statement.

---

## 2. What corresponds to what

| Article section / algorithm | Module |
|---|---|
| Section 2.2 (resolution levels, ambiguity, missing-locus handling) | `hla_elm_toolkit/data.py` |
| Section 2.3 (accuracy, call rate, posterior probability, top-k) | `hla_elm_toolkit/metrics.py` |
| Section 2.7 (95% CI via Wilson score, chosen over naive bootstrap) | `hla_elm_toolkit/metrics.py::wilson_score_interval` |
| Sections 2.4-2.5, Algorithm S1 (ELM training), Algorithm S2 (inference) | `hla_elm_toolkit/elm/base_elm.py` |
| Section 2.6 and Supplementary Section S3.2, Kernel ELM (KELM) | `hla_elm_toolkit/elm/kelm.py` |
| Section 2.6 and Supplementary Section S3.2, Weighted ELM (WELM) | `hla_elm_toolkit/elm/welm.py` |
| Section 2.6 and Supplementary Section S3.2, Ensemble ELM | `hla_elm_toolkit/elm/ensemble_elm.py` |
| Section 2.6, MLP comparator (Supplementary Table S2) | `hla_elm_toolkit/baselines/mlp_baseline.py` |
| Section 2.6, gradient-boosted-tree comparator (Supplementary Table S2) | `hla_elm_toolkit/baselines/gbt_baseline.py` |
| Supplementary Section S2, Algorithm S3 (GRIMM-style, graph-based) | `hla_elm_toolkit/baselines/grimm_style.py` |
| Supplementary Section S2, Algorithm S4 (HaploStats-style) | `hla_elm_toolkit/baselines/haplostats_style.py` |
| Supplementary Section S2, Algorithm S5 (Hapl-o-Mat-style EM) | `hla_elm_toolkit/baselines/haplo_em.py` |
| Supplementary Section S4.1 (two-locus in-registry EM baseline) | `HaploEM(loci=("A","B"))` — same class, different `loci` argument |
| Supplementary Section S4.4 (five-locus in-house EM baseline) | `HaploEM(loci=("A","B","C","DRB1","DQB1"))` — same class |

---

## 3. Installation

```bash
# From the extracted archive directory:
cd hla_elm_toolkit

# Option A: install as a package (recommended)
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -e .                   # base install (numpy only)
pip install -e ".[all]"            # + scikit-learn (MLP/GBT) + pytest

# Option B: no installation, just add the folder to PYTHONPATH
pip install -r requirements.txt
export PYTHONPATH="$PWD:$PYTHONPATH"
```

---

## 4. Quick start

```python
from hla_elm_toolkit.data import make_synthetic_population
from hla_elm_toolkit.elm import BaseELM
from hla_elm_toolkit.metrics import summarize, format_summary, PredictionResult

# Synthetic reference population (NOT real registry data -- see Section 1 above)
ref = make_synthetic_population(n_donors=1000, n_haplotypes=40, seed=0)
train, test = ref.donors[:800], ref.donors[800:]

model = BaseELM(hidden_size=150, seed=0)
model.fit(train, ref, max_diplotypes=3000)

results = []
for g in test:
    dt, post_p = model.predict_top1(g, ref, theta=0.0)
    called = dt is not None
    correct = None  # compare dt against g.truth yourself, per-locus or jointly
    results.append(PredictionResult(donor_id=g.donor_id, called=called, correct=correct))

print(format_summary("Base ELM", summarize(results)))
```

A complete, runnable, end-to-end comparison of **every** method in the
article (ELM, KELM, WELM, Ensemble ELM, HaploStats-style, Hapl-o-Mat-style
EM, GRIMM-style) is provided in:

```bash
python examples/run_demo_comparison.py --n-donors 1500 --seed 0
```

which prints a per-method summary table in the same
`value [95% CI]` style used throughout the article's tables (e.g. Table 10).

---

## 5. Running the tests

```bash
python tests/test_basic.py     # standalone, no pytest required
# or
pytest tests/                  # if pytest is installed
```

The test suite exercises every method on a small synthetic population and
checks basic invariants (frequencies sum to 1, CIs are well-formed,
predictions are in range, etc.); it does not — and cannot, without the
real registry data — check numerical agreement with the article's
reported accuracy/call-rate figures.

---

## 6. Package layout

```
hla_elm_toolkit/
├── README.md                       (this file)
├── requirements.txt
├── setup.py
├── hla_elm_toolkit/
│   ├── __init__.py
│   ├── data.py                     genotype/haplotype/diplotype types,
│   │                                synthetic reference-population generator
│   ├── metrics.py                  accuracy, call rate, Wilson score 95% CI
│   ├── elm/
│   │   ├── __init__.py
│   │   ├── base_elm.py             Algorithm S1 (training) + Algorithm S2 (inference)
│   │   ├── kelm.py                 Kernel ELM (Section 2.6)
│   │   ├── welm.py                 Weighted ELM (Section 2.6)
│   │   └── ensemble_elm.py         Ensemble ELM (Section 2.6)
│   └── baselines/
│       ├── __init__.py
│       ├── haplo_em.py             Algorithm S5 (Hapl-o-Mat-style EM; used for
│       │                            both the Supplementary Section S4.1 two-locus and
│       │                            Supplementary Section S4.4 five-locus baselines)
│       ├── haplostats_style.py     Algorithm S4 (HaploStats-style)
│       ├── grimm_style.py          Algorithm S3 (GRIMM-style, graph-based)
│       ├── mlp_baseline.py         MLP comparator (Section 2.6; requires scikit-learn)
│       └── gbt_baseline.py         Gradient-boosted-tree comparator (Section 2.6; requires scikit-learn)
├── examples/
│   └── run_demo_comparison.py      end-to-end demo running every method
└── tests/
    └── test_basic.py
```

---

## 7. Citation

If you use this code, please cite the article above. Please also note,
per the article's own acknowledgment, that the published Hapl-o-Mat tool
and the NMDP HaploStats web service were not run directly in the
article's evaluation (no compatible batch interface was available in the
authors' environment); this package's `haplostats_style.py` and
`haplo_em.py` are the article's own open, reproducible approximations of
those tools' core algorithms, not the tools themselves.

## 8. First-field rescoring script

* **`examples/recompute_first_field_accuracy.py`**: the comparators
  (GRIMM, hlaR ImputeHaplo and the in-house EM baseline; main text
  Sections 3.2-3.3 and Supplementary Sections S4.1-S4.4) are scored at
  first-field resolution, while the ELM framework's headline accuracy is
  reported at two-field resolution. This script takes already-computed
  two-field ELM predictions (exported to a simple CSV; see the script's
  docstring for the exact format), truncates both predictions and truth
  to first-field, and rescores them, so that a like-for-like first-field
  figure can be reported alongside the comparators (main text Table 10
  and Table 11). No new model training or inference is required.

## 9. License

This reference toolkit (see Section 1 for its scope) is released under
the MIT License (see the `LICENSE` file in this repository). It is
archived at:

- GitHub: https://github.com/Biomedical-Engineering-Laboratory-NTUA/hla-elm-toolkit-2026
- Zenodo (versioned DOI): 10.5281/zenodo.22061923

consistent with the Data Availability Statement in the article.
