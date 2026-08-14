# Testing inventory

This page indexes the tests by source folder. Each inventory records the CTest
target, the test source, the implementation or script under test, what the
test covers, and its current status.

Run the configured suite with:

```bash
ctest -j 6 --output-on-failure
```

Use `ctest -N` for the exact target list in a particular build. The normal
relative tolerance is `1e-3`; tighter checks are used for exact identities and
looser bounds are called out where they measure a scientific approximation.

Current cross-backend exception: `shear1h_cross_backend_test` has one
deliberately failing comparison. Raw radial-series $\Delta\Sigma$ differs from
the C++ `full_ltmz` oracle by 56–86% because
`nfw_profile_family.py` uses a fixed `CONC = 4.0` instead of the per-sample
Child18 concentration. The C++/Python radial-series identity and the other
backend comparisons remain separate checks. See
[`radial_series_vs_full_ltmz_defect.md`](https://github.com/estevesjh/y3_cluster_cpp/blob/docs/sphinx-site/docs/radial_series_vs_full_ltmz_defect.md).

The current model-level exception is `nfw_dsigma_mis_test`: 29 of 30
independent `cluster_toolkit` comparisons pass at `1e-3`, while one interior
table point misses by 0.497%. This is retained as a strict failing test until
the cause is understood; see
[`nfw_dsigma_mis_defect.md`](https://github.com/estevesjh/y3_cluster_cpp/blob/docs/sphinx-site/docs/nfw_dsigma_mis_defect.md).

## Folder inventories

```{toctree}
:maxdepth: 1

testing/test
testing/src_pipelines_des_y3
testing/src_modules
```

Status labels mean:

- **Passing** — last recorded configured run passed.
- **Known failing** — intentionally exposes an unresolved defect; see the
  linked issue note.
- **Characterization** — asserts or records currently observed behavior,
  including a known numerical defect, rather than declaring that behavior
  scientifically correct.
- **Disabled/diagnostic** — not part of the configured CTest suite.
