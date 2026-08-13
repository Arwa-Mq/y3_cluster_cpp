# Testing

The suite is CTest-driven: every C++ test is a Catch2 executable under
`test/`, registered in `test/CMakeLists.txt`; a small number of tests are
Python scripts run as ordinary CTest commands (they exit non-zero on
failure, same contract as the C++ executables). See {doc}`installation`
for environment setup and `BUILDING.md` in the repo root for the full
build recipe.

```{admonition} Tolerance convention
:class: note
Unless a test documents a tighter or looser bound in a comment, the
default relative tolerance is **1e-3**. Cross-implementation checks
(e.g. a fixed-GL evaluator vs an adaptive-Cuhre reference, or a C++
backend vs its Python counterpart) sometimes carry a looser, explicitly
documented bound where two genuinely different quadrature strategies are
being compared over the same physics; algebraic identities and
structural invariants (endpoints, limits, factorisations) are checked at
or near machine precision instead, since they must hold exactly by
construction.
```

## Running the suite

```bash
# CPU-only tests (no CUDA modules) — builds and runs on a login node
cd release-build && ctest -j10 --output-on-failure

# Everything, including the CUDA backends — requires a GPU and the
# pinned toolchain (BUILDING.md): cudatoolkit/12.2 + gcc-native/12.3.
# Keep -j low on a shared GPU: on a card another job already has most
# of the memory on, higher -j can trip a spurious cudaCheckError()
# out-of-memory in table construction rather than a real regression.
cd gpu-build && ctest -j2 --output-on-failure

# One test by name
ctest --test-dir release-build -R shear1h_radial_series_test --output-on-failure

# Failing-test output only
CTEST_OUTPUT_ON_FAILURE=1 ctest --test-dir release-build
```

GPU test binaries `dlopen` their `.so` module at runtime, which needs
`libnvToolsExt.so.1` on the loader path; the pinned cudatoolkit/12.2
module's default `LD_LIBRARY_PATH` does not include it, so add it
explicitly before running `gpu-build` tests:

```bash
export LD_LIBRARY_PATH="/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/cuda/12.2/targets/x86_64-linux/lib:$LD_LIBRARY_PATH"
```

`Y3_CLUSTER_CPP_DIR` must be set in the environment for any test that
reads fixture data under `data/` or `test/` by relative path.

## Core numerics and utilities

Grid construction, integration volumes, interpolation, and small
free-function helpers used across every model and module.

| Test | What it checks |
|---|---|
| `interp_1d_test`, `interp_2d_test` | 1D/2D interpolation primitives (`Interp1D`/`Interp2D`) against known function values |
| `transform_test` | coordinate/variable-transform helpers used inside integrands |
| `polynomial_test` | the small polynomial-evaluation utility |
| `integer_pow_test` | fast integer-power helper vs `std::pow` |
| `make_grid_points_test`, `make_grid_points_cartesian_product_test`, `make_grid_cartesian_product_test`, `make_grid_splatted_test` | grid-point construction from datablock "wall-of-numbers" configuration, including the Cartesian-product expansion |
| `make_integration_volume_test` | per-bin integration-volume construction from configuration |
| `get_grid_axes_from_datablock_test` | reading grid axis definitions out of a `DataBlock` |
| `param_space_explorer_test` | the parameter-space sampling helper |
| `point_3d_test` | the 3D point type used by volume/grid code |
| `sample_variance_test` | the sample-variance accumulator |
| `bessel_analytic_integral_test`, `bessel_integral_equivalence_test`, `sin_cos_analytic_integral_test` | Bessel/trig integral kernels vs closed-form references |
| `multi_dimensional_integrator_test` | the generic N-dimensional integrator wrapper (Cuhre/PAGANI dispatch) |
| `integration_test` | end-to-end integration-module plumbing |
| `mpi_support_test` | MPI rank/size detection support code |
| `timing_sentry_test` | the `timing = T` per-module wall-clock instrumentation |
| `fpsupport_test` | floating-point environment/exception-flag helpers |

## Physics models (`src/models`)

Header-only model/kernel classes, each tested directly against a known
reference (closed-form, tabulated, or cross-checked against an
independent implementation such as `astropy`).

| Test | Model | Reference |
|---|---|---|
| `hmf_t_test` | `HMF_t` (halo mass function) | tabulated `dn/dlnM` fixture (`test_hmf_z0_z03.dat`) |
| `mor_t2_test`, `mor_probability_unity_test`, `mor2_probability_test` | SDSS mass-observable relation | closed-form + probability-normalization checks |
| `sptxdes_mor_1d`, `mor_y3xspt_t` | SPT$\times$DES MOR variants | closed-form |
| `lc_lt_t_test`, `lc_lt_probability_unity_test`, `lc_lt_probability_test` | $P(\lambda_{\rm true}\mid \lambda_{\rm obs})$ kernel | closed-form + unity-normalization |
| `lo_lc_t_test`, `lo_lc_probability_unity_test` | $P(\lambda_{\rm obs}\mid \lambda_{\rm true})$ kernel | closed-form + unity-normalization |
| `zo_zt_des_test`, `int_zo_zt_des_test` | $P(z_{\rm obs}\mid z_{\rm true})$ photo-z kernel + its integral | closed-form |
| `roffset_t_test`, `roffset_probability_unity_test`, `roffset_probability_test` | miscentering offset kernel | closed-form + unity-normalization |
| `dv_do_dz_t_test` | comoving volume element $dV/d\Omega dz$ | closed-form |
| `omega_z_sdss_test`, `omega_z_y3xspt_test` | survey solid angle $\Omega(z)$ (per-survey variants) | closed-form |
| `sigma_crit_inv_test` | $\langle\Sigma_{\rm crit}^{-1}\rangle(z_l)$ | `astropy.cosmology` cross-check |
| `pxizeta_t_test`, `int_pxizeta_test` | $P(\xi\mid\zeta)$ richness-scatter kernel + its integral | closed-form |
| `mz_power_law`, `ln_mez_power_law` | mass/redshift power-law helper traits | closed-form |

## Harness self-tests

`test_Example*IntegrationModule` and `Example*Integrand_test` exercise the
`CosmoSISScalarIntegrationModule`/`VectorIntegrationModule`/`OneDIntegrationModule`
templates and the `DEFINE_COSMOSIS_*_MODULE` macros themselves (the
machinery every real module is built from), not physics.

## `src/pipelines/des_y3` — radial series and NFW decomposition

`shear1h_radial_series_test` (C++, `test/shear1h_radial_series.test.cc`)
and `des_y3_pipeline_python_test` (Python,
`test/des_y3_pipeline.test.py`) cover the offline `U_ell` table
(`RadialSeriesTable`) and the population-moment decomposition it feeds:

- **NFW mass/scale/amplitude factorisation** — $r_s \propto M^{1/3}$,
  $A_0 \propto r_s$ — checked to `1e-12` (exact by construction).
- **Mixture-decomposition endpoints** — `u_mix(f_mis=0)` and
  `u_mix(f_mis=1)` must reproduce the pure centred/miscentred terms
  exactly, and the mixture must be affine in `f_mis` (`1e-13`).
- **Second- and third-order series reconstruction** — a synthetic
  symmetric ($\mu_3=0$) and a synthetic skewed ($\bar\mu=0$,
  $\mu_3\neq0$) mass population, each evaluated two independent ways
  (direct per-mass sum vs the $\mu_2$/$\mu_3$ series expansion), agree to
  the default **1e-3**.
- **NFW shape stability through $x=1$** (Python,
  `nfw_profile_family.u_cen`) vs an `mpmath` (40-digit) Taylor reference,
  **1e-3**.
- **Offline generator** (`generate_radial_series_tables.py`): finite-
  difference weights reproduce polynomial derivatives exactly (to
  `1e-8`); the miscentering gamma-kernel averaging matrix conserves a
  constant profile to `1e-3`.
- `validate_radial_series.py` (not a CTest — a standalone report) goes
  further: it checks the committed table against an independent
  `mpmath` Taylor route, and reports truncation error against the exact
  fixed-GL mass integral over all 12 real DES Y3 bins (`ell<=2`
  tolerance 0.75%, `ell<=3` tolerance 1.00% — looser than 1e-3 because
  this is a genuine truncation error of a 4-term series, not a
  round-off check), plus the disclosed shape gap against production
  `Shear1hMisSel.so` (different centred-profile convention — see
  `docs/dsigma_hh_debug_flag.md`).

## `src/pipelines/des_y3` — Phase 2 observable backends

Every backend added in the 2026-08-12 Phase 2 rollout (`full_ltmz`,
`fast_mass`, and the traditional 1h+2h max model, across
`number_counts`, `shear_1h2h`, and `shear_projection`) gets a unit test
alongside its Python/C++/CUDA implementation, following the same
pattern as the radial-series tests above: an exact structural/wiring
check plus a numeric check against the reference the module's own
README already establishes (the adaptive `full_ltmz` fiducial per
{doc}`pipeline_organization`'s accuracy policy — never a production
`.so` as ground truth; agreement with production is reported
separately as an algorithm-identity check).

### C++ (`release-build`, CPU-only)

| Test | Module | What it checks |
|---|---|---|
| `numcounts_full_ltmz_test` | `NumCountsFullLtmz` | K_i/K_j richness and photo-z kernels vs an independent `scipy.special.erf/erfc` transcription; `MOR_HOD_t` vs an independent `scipy.special.gammaln` transcription; `operator()` equals the documented HMF·dV/dΩdz·Ω(z)·K_j·K_i·P_HOD product on two independent synthetic fixtures (including an Einstein-de Sitter, exactly-bilinear-HMF fixture with zero interpolation error); bin-selection and out-of-range/support-boundary invariants |
| `shear1h_full_ltmz_test` | `Shear1hFullLtmz` | reproduces the documented O_ij integrand from its own production sub-models; miscentering mixture is affine in `f_mis`; falls off monotonically with radius |
| `shear1h_fast_mass_test` | `Shear1hFastMass` | miscentering defaults/`lob_centers`/closed-form `r_mis`; GAMMA-kernel `NFW_DSIGMA_MIS` vs the shared Python replica; `evaluate()` vs an independently-assembled reference recombining the same lower-level production pieces, and vs the production-equivalent `Shear1hMisSelGL` header class on an identical sample |
| `shear1h2h_max_test` | `Shear1h2hMax` | `f_mis` mixture exact endpoints and boundedness; `dSigma_hh` NaN-sanitization (see `docs/dsigma_hh_debug_flag.md`); the `phi_max(1h, b·2h)` composition, including the zero-two-halo limit recovering the pure one-halo stack |
| `shear_prj_fast_mass_test` | `ShearPrjFastMass` | `R_lambda`/`default_lob_centers`; the theta-grid builder (`build_theta_grid`) and LoS-slab exclusion angle (`theta_excl_at_z`) vs the Python `fast_mass` port; the c=4/single-kernel `NFW_DSIGMA_MIS` convention vs the shared Python replica; `ShearPrjCore`'s own wall-grid construction |

### CUDA (`gpu-build` only)

Built with `-DUSE_CUDA=On` and the pinned cudatoolkit/12.2 +
gcc-native/12.3 toolchain (`BUILDING.md`); not part of the CPU-only
`release-build` suite.

| Test | Module |
|---|---|
| `num_counts_full_ltmz_gpu_test` | `NumCountsFullLtmzGpu` |
| `shear1h_full_ltmz_gpu_test` | `Shear1hFullLtmzGpu` |
| `shear1h2h_max_gpu_test` | `Shear1h2hMaxGpu` |
| `shear_prj_frozen_gpu_test` | `ShearPrjFrozenGpu` |
| `dsigma_prj_full_ltmz_gpu_test` | `DSigmaPrjFullLtmzGpu` |

### Python (`des_y3_pipeline_python_test`)

`test/des_y3_pipeline.test.py` gained one `unittest.TestCase` class per
new Python backend (`TestNumCountsFastMass`, `TestNumCountsFullLtmz`,
`TestShear1hFastMass`, `TestShear1hFullLtmz`, `TestShear1h2hMax`,
`TestShearPrjFastMass`), each replaying a real-pipeline dump and
comparing against the production `.so` output already stored in it.
Two dumps are used, both gitignored and regenerated on demand:

- `docs/figs/real_pipeline_extract_output` (`cosmosis
  docs/figs/real_pipeline_extract.ini`) — the base dump: HMF,
  distances, `sel_function`, `numcountssel`, `shear1hmissel`. Guards
  `@unittest.skipUnless(HAS_DUMP, ...)` — skipped, not failed, if
  absent. `TestShear1h2hMax`'s primary check exercises only the
  `dSigma_hh`-defect-free limit on this dump (two-halo term forced to
  zero must reproduce the validated 1h `fast_mass` backend exactly),
  per the same caveat as `shear1h2h_max_test` above.
- `docs/figs/real_pipeline_extract_prj2h_output` (`cosmosis
  docs/figs/real_pipeline_extract_prj2h.ini`) — the sibling dump with
  `halo_model` run at `compute_lensing_2h = T` (so `haloModel/dSigma_hh`
  is populated, NaN-heavy tail and all) and the projection chain
  (`b_sel_marg` → `bsel` → `dsigma_prj` → `shear_prj_frozen_physics`).
  Guards `@unittest.skipUnless(HAS_DUMP_PRJ2H, ...)`. Two tests use it:
  `TestShear1h2hMax.test_full_model_is_finite_and_two_halo_contributes`
  checks the full (real, NaN-sanitized) 1h+2h composition stays finite
  and that the two-halo branch is genuinely selected somewhere on the
  wall — no golden number is pinned, since that would bake in the
  documented table defects; `TestShearPrjFastMass.
  test_matches_exact_evaluator_and_frozen_physics` runs
  `ShearPrjFastMass.set_sample()`/`wall_outputs()` end to end and
  checks machine-precision agreement with the exact `dsigma_prj`
  evaluator plus the documented frozen-physics bound vs
  `dsigma_prj_frozen_physics`.

`TestShearPrjFastMass` also has two dump-free tests
(`test_theta_grid_is_bounded_monotone_and_matches_golden`,
`test_theta_exclusion_angle_matches_golden`) pinning `build_theta_grid`/
`theta_excl_at_z` against both independently-derived closed-form bounds
and the same golden values `shear_prj_fast_mass_test` pins on the C++
side — a cross-language identity check.

## Offline validation scripts (not CTest)

These are deeper diagnostic reports, not pass/fail unit tests — run them
by hand when investigating accuracy, not as part of the regular suite:

- `src/pipelines/des_y3/validate_against_fiducial.py` — the full
  namespace accuracy matrix, every path vs the `full_ltmz` fiducial.
- `src/pipelines/des_y3/observables/*/*/python/validate_vs_production.py`
  — per-backend Python-vs-production-`.so` ratio reports.
- `src/pipelines/des_y3/observables/shear_1h2h/fast_mass/python/validate_shear1h2h_max.py`
  — the traditional 1h+2h max-model validation, including the disclosed
  `dSigma_hh` two-halo-table caveat.
- `src/pipelines/des_y3/observables/shear_1h2h/radial_series/python/validate_radial_series.py`
  — see above.

Most of these read `docs/figs/real_pipeline_extract_output` (regenerate
with `cosmosis docs/figs/real_pipeline_extract.ini`); the projection-shear
and max-model validators need the sibling dump
`docs/figs/real_pipeline_extract_prj2h_output` (`cosmosis
docs/figs/real_pipeline_extract_prj2h.ini`), which additionally runs the
`b_sel_marg`/`bsel`/`dsigma_prj`/`shear_prj_frozen_physics` chain and
`halo_model` with `compute_lensing_2h = T`. All print a `PASS`/error
message; they intentionally are not wired into CTest because they
need that multi-gigabyte real-pipeline dump rather than a small
in-repo fixture.
