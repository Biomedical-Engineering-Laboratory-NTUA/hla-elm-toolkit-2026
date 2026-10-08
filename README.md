# hla-elm-toolkit-2026

Reference implementation accompanying:

> Kepentzis, S.; Chatzistamatiou, T.; Digalakis, J.; Petropoulou, O.;
> Matsopoulos, G.K.; Koutsouris, D. *Improving the Usability of Donor Data from Registries Using
> an Extreme Learning Machine Approach to Upgrade Low-/Mid- to
> High-Resolution HLA Data.* Genes (MDPI),
> Genes 2026, 17, 1239, https://doi.org/10.3390/genes17101239 

This repository provides a from-scratch **reference implementation** of
the core Extreme Learning Machine (ELM) mathematics described in the
article, of its three extensions (KELM, WELM, Ensemble ELM), and of the
comparator algorithms the article specifies in pseudocode (Algorithms
S1-S5, Supplementary Section S2). Python 3.9+ (developed on 3.12) and
C++17 versions are included.

## What this repository is, and is not

**It is** an independent, readable implementation of the published
algorithms, written so that a reader can follow the method and run it on
their own data.

**It is not** the registry-scale production code that produced the
article's reported figures, and it does not reproduce them. Two
differences matter most:

1. The base ELM here uses a **single random hidden layer** over a flat
   multi-hot locus/allele encoding. The article's network uses a
   domain-structured HL-1/HL-2/HL-3 topology (Section 2.4) fixed in
   advance by partial-haplotype containment. The core ELM property -
   random, untuned input-to-hidden weights and a closed-form
   output-weight solution with no back-propagation - is preserved; the
   topology is not.

2. **No registry data is included or reachable.** The HTO, ORAM and
   GRPT donor data are not distributable (donor privacy; Hellenic
   Transplant Organization data-sharing restrictions). The package
   generates a small synthetic population for demonstration and unit
   testing only. No accuracy or call-rate number produced by this code
   is comparable to the article's results.

The production implementation is bound to the registry data environment
and is not distributed here. Requests relating to it should be addressed
to the corresponding author.

## Contents

| Path | Contents |
|---|---|
| `hla_elm_toolkit_python/` | Python package (base ELM, KELM, WELM, Ensemble ELM, all comparator baselines, metrics) |
| `hla_elm_toolkit_cpp/` | C++17 port of the base ELM and the Hapl-o-Mat-style EM baseline |

Each subdirectory has its own README with an article-section mapping,
installation instructions and a quick-start example.

## Archive and citation

- GitHub: https://github.com/Biomedical-Engineering-Laboratory-NTUA/hla-elm-toolkit-2026
- Zenodo (versioned DOI): 10.5281/zenodo.22061923

Released under the MIT License (see `LICENSE`). If you use this code,
please cite the article above together with the Zenodo DOI.
