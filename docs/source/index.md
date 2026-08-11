# DES Y3 Cluster Cosmology

Technical documentation of the DES Y3 cluster-cosmology prediction pipeline:
the scientific model, the numerical recipes, and their implementation in
`y3_cluster_cpp` and its companion Python modules.

The documentation is organised so that, for any quantity entering the
likelihood, a reader can answer four questions:

1. **What** quantity is being computed?
2. What is its **mathematical definition**?
3. How is the expression **evaluated numerically**?
4. **Where** is the calculation implemented in the repository?

The scientific definition of each observable is kept independent of the
implementation; CPU and GPU code paths are described where relevant.

```{toctree}
:maxdepth: 2
:caption: Getting started

overview
installation
```

```{toctree}
:maxdepth: 2
:caption: Mathematical Framework

science/index
systematics/index
```

```{toctree}
:maxdepth: 2
:caption: Numerics and implementation

numerics/index
modules/index
cosmosis/index
```

```{toctree}
:maxdepth: 2
:caption: Data and validation

data/index
validation/index
```

## Archival documents

The LaTeX documents below are the archival, paper-grade record from which
much of this site is ported. Each ported chapter cites its source; where
this site and a PDF disagree, the site is the living reference. The
sources live under `docs/` (PDFs are built locally with `pdflatex`; they
are deliberately not tracked in git):

- [pipeline_modules.tex](https://github.com/estevesjh/y3_cluster_cpp/blob/master/docs/pipeline_modules.tex)
  — wired-pipeline algorithms, DataBlock contracts, timing audit,
  quadrature knob cheat-sheet.
- [projection_lensing_paper.tex](https://github.com/estevesjh/y3_cluster_cpp/blob/master/docs/projection_lensing_paper.tex)
  — the optical-projection lensing model.
- [emulator_validation.tex](https://github.com/estevesjh/y3_cluster_cpp/blob/master/docs/emulator_validation.tex)
  — validation of the `cp_camb` linear-$P(k)$ emulator against CAMB.
- [shear1h_radial_factorization.tex](https://github.com/estevesjh/y3_cluster_cpp/blob/master/docs/shear1h_radial_factorization.tex)
  — factorisation strategies for the one-halo shear mass integral.

The selection-model derivations live in the `RichnessSelection` repository
(`docs/richness_selection_function.tex`, `docs/richness_selection.tex`,
`docs/delta_sigma_prj_derivation.tex`); they are the source of truth for
the model chapters here.
