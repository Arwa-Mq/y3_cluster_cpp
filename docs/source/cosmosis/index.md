# How CosmoSIS works

A brief primer — the reference pipeline itself is documented in
{doc}`../running`; this page is background for readers new to CosmoSIS.

CosmoSIS runs an ordered list of **modules**, each a shared library
(`.so`) or Python script exposing `setup(options)` /
`execute(block, config)` / `cleanup(config)`. `setup` runs once per
process (reads ini options, builds fixed tables); `execute` runs once
per sample (reads/writes the sample's state); `cleanup` runs once at
exit.

The shared state is the **DataBlock**: a per-sample key-value store,
keyed by `(section, name)`. A module reads its inputs with
`block[section, name]` and writes outputs the same way — there is no
other communication channel between modules, so the DataBlock *is* the
pipeline's data contract. A `put_val` on a `(section, name)` that
another module already wrote in the same sample does **not** overwrite
the value — see the note on this in {doc}`../running` about co-running
DES Y1 and DES Y3 modules for comparison.

A **sampler** drives the pipeline over many samples: `test` runs it
once at the values-file starting point (a smoke test); `apriori` draws
uniformly from the prior box; `emcee` and `polychord` do MCMC/nested
sampling over the varied parameters, calling the pipeline (and reading
the final `likelihoods` section) at every step.

Three categories of module in this analysis:

1. **CosmoSIS Standard Library** — `consistency`, `GrowthFactor`,
   `MfTinker` (external repository).
2. **This analysis's own modules** — the C++/CUDA `.so` files and
   Python modules documented per-module from {doc}`../running`.
3. **External dependencies** — CUBA/cubacpp, PAGANI, GSL,
   `cluster_toolkit`, the CUDA toolkit for GPU variants.

Naming convention: C++ modules use CamelCase ini section names, Python
modules snake_case — except modules whose source hardcodes
`module_label()` (`b_sel_marg`, `bsel`, `shear_prj_frozen_physics`),
which keep their lowercase labels regardless.
