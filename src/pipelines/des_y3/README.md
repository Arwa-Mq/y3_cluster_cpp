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
    │   └── full_ltmz/python/    full (lambda_true, lnM, z) reference counts
    └── shear_1h2h/
        └── radial_series/python/  offline U_ell tables + moment evaluator
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
