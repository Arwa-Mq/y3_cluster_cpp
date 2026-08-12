# `src/pipelines/des_y3` — new maintained DES Y3 implementations

This tree is the namespace approved in
[docs/module_reorganization_plan.md](../../../docs/module_reorganization_plan.md)
for *new* implementations belonging to the maintained DES Y3 observable
family. The validation baseline it is measured against is recorded in
[docs/des_y3_maintenance_manifest.md](../../../docs/des_y3_maintenance_manifest.md).

Layout is `observable -> integration strategy -> language/backend`:

```text
des_y3/
├── shared/                      Python replicas of the shared model layer
│                                (HMF_t, DV_DO_DZ_t, OMEGA_Z_DES, SelGLCore),
│                                convention-exact against src/models/*.hh
└── observables/
    ├── number_counts/
    │   ├── fast_mass/python/    exact z contraction, W(lnM) outside the operator
    │   └── full_ltmz/
    │       ├── python/          fixed-GL reference (per-(M,z) lt brackets)
    │       ├── cpp/             NumCountsFullLtmz.so (adaptive Cuhre)
    │       └── cuda/            NumCountsFullLtmzGpu.so (PAGANI)
    ├── shear_1h2h/
    │   ├── fast_mass/python/    exact z contraction + direct GL mass sum
    │   ├── full_ltmz/python/    explicit (lt, lnM, z) x production profile
    │   └── radial_series/
    │       ├── python/          offline U_ell generator + moment evaluator
    │       └── cpp/             Shear1hRadialSeries.so (same data, GSL interp)
    └── shear_projection/
        └── fast_mass/python/    exact-z port of ShearPrjCore (no freeze)
```

Ground rules (from the approved proposal):

- Nothing here moves, wraps, or replaces a production entry point. The
  production stages (`sel_function`, `NumCountsSel`, `Shear1hMisSel`,
  `b_sel_marg`, `bsel`, `shear_prj_frozen_physics`) stay where they are.
- Existing C++ module and integration templates are immutable dependencies.
- Directories exist only when they contain a runnable implementation or a
  substantive design document; empty language placeholders are not created.
- Each implementation documents its integration variables, DataBlock
  contract, composed models, numerical tolerance against its reference, and
  its status (production / reference / experimental / planned) in its own
  `README.md`.

The offline `radial_series` derived data lives under
[`data/radial_series/`](../../../data/radial_series/) and is generated once
by the generator in `observables/shear_1h2h/radial_series/python/`; it is
never regenerated inside an MCMC sample.

## Measured accuracy and cost (2026-08-12, fiducial point, 12 pinned bins)

**Accuracy policy** (plan owner, 2026-08-12): accuracy is quoted
against the **`full_ltmz` fiducial** — the fully explicit calculation,
whose own precision is certified by internal convergence (doubling
every quadrature node moves it by ≤ 3.8e-4 counts / ≤ 3.1e-4 shear;
widening the λ bracket by ≤ 1.4e-4) and by three independent
quadrature strategies agreeing (Python fixed-GL, Cuhre, PAGANI).
Agreement with a production `.so` is a separate *algorithm-identity*
check — it proves the same computation, not correctness.

Counts (fiducial: Python `full_ltmz`, 83 ms):

| Path | Time | Error vs fiducial | Identity vs production |
|---|---|---|---|
| production `NumCountsSel.so` (`fast_mass`) | 6 ms | **7.6e-4** (S_ij tabulation) | — |
| Python `fast_mass` | 5 ms | 7.6e-4 (same algorithm) | 2.4e-15 |
| C++ `full_ltmz` (Cuhre, eps_rel 1e-4) | 3.1 s | 4.9e-4 — inside the fiducial's own convergence band | — |
| CUDA `full_ltmz` (PAGANI, 1×A100) | 2.0 s | 5.1e-4 (6.0e-5 from the C++ twin) | — |

Shear, 12 bins × 10 radii (fiducial: Python `full_ltmz`, 149 ms;
for `radial_series` the fiducial uses that strategy's own
fixed-convention profile, doubled nodes):

| Path | Time | Error vs fiducial | Identity vs production |
|---|---|---|---|
| production `Shear1hMisSel.so` (`fast_mass`) | 9 ms | **8.4e-4** (S_ij tabulation) | — |
| Python `fast_mass` | 74 ms | 8.4e-4 (same algorithm) | 3.1e-15 |
| Python `radial_series` (ℓ≤2) | 6 ms | **3.7e-3 total** (tabulation + truncation + interpolation) | — |
| C++ `Shear1hRadialSeries.so` | 7 ms | as Python + 1.6e-4 interpolation-scheme difference | — |

Projection shear, 180-point wall (production frozen module = 82 ms;
**no `full_ltmz` fiducial exists yet for this observable** — the fully
z-, mass-, and angle-resolved reference is the open matrix cell, and
until it lands the exact-z evaluator is the best available reference):

| Path | Time | Error vs best available reference | Identity check |
|---|---|---|---|
| production `ShearPrjFrozenPhysics.so` | 82 ms | 5.5e-5 (frozen-physics approx., measured) | — |
| Python `fast_mass` (exact z, no freeze) | 270 ms | reference (pending `full_ltmz`) | 1.6e-11 vs exact `DSigmaPrjEvaluator.so` |

Caution on the Cuhre knobs: `eps_rel = 1e-3` is 9× faster (0.36 s) but
silently ~1% wrong in the lowest-richness bins (the near-delta HOD
ridge; see the cpp README) — the 1e-4 setting is the validated one.
Full validation records live in each implementation's README.
