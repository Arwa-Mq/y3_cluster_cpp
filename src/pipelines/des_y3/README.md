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
    │   └── full_ltmz/
    │       ├── python/          fixed-GL reference (per-(M,z) lt brackets)
    │       ├── cpp/             NumCountsFullLtmz.so (adaptive Cuhre)
    │       └── cuda/            NumCountsFullLtmzGpu.so (PAGANI)
    └── shear_1h2h/
        └── radial_series/
            ├── python/          offline U_ell generator + moment evaluator
            └── cpp/             Shear1hRadialSeries.so (same data, GSL interp)
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

## Measured precision and cost (2026-08-12, fiducial point, 12 pinned bins)

Counts, per sample (production `NumCountsSel.so` fast path = 6 ms
baseline; the `full_ltmz` backends are references, exempt from
production timing by the plan):

| Backend | Time | Accuracy (measured) |
|---|---|---|
| Python `full_ltmz` (fixed GL 96×64×32) | 83 ms | 7.6e-4 vs production |
| C++ `full_ltmz` (Cuhre, `eps_rel` 1e-4) | 3.1 s | 4.9e-4 vs Python; reported err 1.0e-4 |
| CUDA `full_ltmz` (PAGANI, 1×A100) | 2.0 s | 6.0e-5 vs C++; reported err 0.7–1.0e-4 |

Shear, 12 bins × 10 radii per sample (production `Shear1hMisSel.so`
exact-GL = 9 ms baseline):

| Backend | Time | Accuracy (measured) |
|---|---|---|
| Python `radial_series` (ℓ≤2) | 6 ms | truncation ≤0.45% vs exact mass sum |
| C++ `Shear1hRadialSeries.so` | 7 ms | 1.6e-4 vs Python (interpolation scheme) |

Caution on the Cuhre knobs: `eps_rel = 1e-3` is 9× faster (0.36 s) but
silently ~1% wrong in the lowest-richness bins (the near-delta HOD
ridge; see the cpp README) — the 1e-4 setting is the validated one.
Full validation records live in each implementation's README.
