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
    │   ├── full_ltmz/
    │   │   ├── python/          explicit (lt, lnM, z) x production profile
    │   │   └── cpp/             Shear1hFullLtmz.so (adaptive Cuhre per (bin,R))
    │   └── radial_series/
    │       ├── python/          offline U_ell generator + moment evaluator
    │       └── cpp/             Shear1hRadialSeries.so (same data, GSL interp)
    └── shear_projection/
        └── fast_mass/
            ├── python/          exact-z port of ShearPrjCore (no freeze)
            └── cpp/             ShearPrjFastMass.so (one core, both observables)
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

## The matrix: accuracy and timing per measurement mode (2026-08-12)

**Accuracy policy** (plan owner): the reference is the **adaptive**
`full_ltmz` calculation (`full_ltmz_core.full_ltmz_mass_integral_adaptive`:
vectorised adaptive mass integral, reported error ≤ 1e-6; mass limits
with the lower bound at richness → 0 — the inverse scaling relation is
invalid where scatter dominates — and the upper bound from the inverted
scaling relation at 4·λ_max). A fixed-GL implementation is never the
reference; it is certified against the adaptive one (3.5e-5 counts,
4.9e-5 shear) and then used as the fast stand-in. Production agreement
is a separate *algorithm-identity* check. All numbers: real pipeline,
fiducial widePlanck point, pinned 12-bin wall; times per MCMC sample.

### Number counts (12 bins; production `NumCountsSel.so` = 6 ms)

| Strategy / backend | Time | Error vs fiducial | Notes |
|---|---|---|---|
| `full_ltmz` / Python adaptive | 25 s | **reference** (reported err ≤ 1e-6) | |
| `full_ltmz` / Python (fixed GL) | 83 ms | 3.5e-5 | certified by the adaptive reference |
| `full_ltmz` / C++ (Cuhre, eps 1e-4) | 3.1 s | 4.9e-4 | inside fiducial band |
| `full_ltmz` / CUDA (PAGANI, A100) | 2.0 s | 5.1e-4 | 6.0e-5 from C++ twin |
| `fast_mass` / Python | 5 ms | 7.6e-4 | identity 2.4e-15 vs production |
| `fast_mass` / C++ | 6 ms | 7.6e-4 | **is** production `NumCountsSel.so` |
| `fast_mass` / CUDA | — | — | not warranted (1-D sum) |
| `radial_series` | n/a | n/a | counts need no radial operator (f = 1) |

### One-halo miscentred shear (12 bins × 10 radii; production `Shear1hMisSel.so` = 9 ms)

| Strategy / backend | Time | Error vs fiducial | Notes |
|---|---|---|---|
| `full_ltmz` / Python adaptive | 35 s | **reference** (reported err ≤ 1e-6) | |
| `full_ltmz` / Python (fixed GL) | 149 ms | 4.9e-5 | certified by the adaptive reference |
| `full_ltmz` / C++ (Cuhre, eps 1e-4) | 51 s | 3.3e-4 | 120 adaptive triples; all converged |
| `full_ltmz` / CUDA (PAGANI) | — | — | planned (device haloModel interp needed) |
| `fast_mass` / Python | 74 ms | 8.4e-4 | identity 3.1e-15 vs production |
| `fast_mass` / C++ | 9 ms | 8.4e-4 | **is** production `Shear1hMisSel.so` (`method=exact`) |
| `radial_series` / Python (ℓ≤2) | 6 ms | 3.7e-3 total | tabulation + truncation + interp, same-profile fiducial |
| `radial_series` / C++ (`Shear1hRadialSeries.so`) | 7 ms | 3.7e-3 + 1.6e-4 | interp-scheme difference vs Python |
| `radial_series` / CUDA | — | — | not warranted (3 table lookups per point) |

### Projection shear (180-point wall; production `ShearPrjFrozenPhysics.so` = 82 ms)

| Strategy / backend | Time | Error vs reference | Notes |
|---|---|---|---|
| `full_ltmz` | — | — | **open cell**: needs a design decision (unfreeze ξ_NL's z argument? resolve the richness selection inside instead of the b_sel(θ) plateaus?) |
| `fast_mass` / Python (exact z, no freeze) | 270 ms | best available reference | identity 1.0e-11 vs C++ |
| `fast_mass` / C++ (`ShearPrjFastMass.so`) | 154 ms | best available reference | 9.9e-12 vs exact `dsigma_prj` evaluator (same core); both observables in one pass |
| `fast_mass` / frozen production | 82 ms | 5.5e-5 from exact | the frozen-physics approximation, measured |
| `radial_series` | — | — | planned: U_ℓ(x, x_θ) tables with the θ coordinate retained (plan §radial_series) |

Caution on the Cuhre knobs: `eps_rel = 1e-3` is 9× faster (0.36 s) but
silently ~1% wrong in the lowest-richness bins (the near-delta HOD
ridge; see the cpp README) — the 1e-4 setting is the validated one.
Full validation records live in each implementation's README.
