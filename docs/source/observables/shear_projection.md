# Shear Projection

`C++` · `y3_cluster_cpp` · `Cluster observable` · module `shear_prj_frozen_physics` · `82 ms/sample`

Computes $\Sigma^{\rm prj}$ — in the paper's language, **the two-halo
term sourced by correlated line-of-sight structure**, carrying the
selection-affected bias $b_{\rm sel}(\theta)$: the surface density and
tangential shear contributed by foreground/background haloes whose
presence also boosts the optical selection. It replaces the conventional
(unselected-bias) two-halo term in the reference shear composition
$\gamma_t^{\rm theory} = \langle\gamma_t^{1h}\rangle + \gamma_t^{\rm prj}$.

## Script

- Model: [`src/models/sigma_prj_frozen_t.hh`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/models/sigma_prj_frozen_t.hh)
  (`y3_cluster::ShearPrjFrozenPhysics`, "Option C" — frozen-physics
  fixed-grid reduction of the full `ShearPrjEvaluator`, which it shares
  the $\theta$-grid and $z$-grid builders with via `sp_detail`).
- Module driver: [`src/modules/sigma_prj_cpu/ShearPrjFrozenPhysics.cc`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/modules/sigma_prj_cpu/ShearPrjFrozenPhysics.cc).
- Compiled library loaded by CosmoSIS:
  `${Y3_CLUSTER_CPP_DIR}/release-build/src/modules/sigma_prj_cpu/ShearPrjFrozenPhysics.so`.
- Disk tables (loaded once at construction, $c = 4$, **single**/delta
  offset kernel — not the gamma kernel `Shear1hMisSel` uses):
  `data/nfw_off_center/table_1000_1e-03_5e+03_single_{logx, logxmis}.txt`,
  `…_log_deltasigma_single.txt`.

## DES Y3 implementations

The production frozen-physics module above remains path-stable. The new
namespace under `src/pipelines/des_y3/observables/shear_projection` contains:

| Strategy | Backend | Implementation and status |
|---|---|---|
| `full_ltmz` | CUDA | `DSigmaPrjFullLtmzGpu.so`; adaptive $(\ln\theta,z,\ln M)$ PAGANI reference, with the innermost-radius convergence study still open |
| `fast_mass` | Python | Exact-$z$ reference port of `ShearPrjCore` |
| `fast_mass` | C++ | `ShearPrjFastMass.so`; exact-$z$ core emitting $\Delta\Sigma$ and shear in one pass |
| `fast_mass` | CUDA | `ShearPrjFrozenGpu.so`; CUDA implementation of the explicitly frozen production machinery |
| `radial_series` | — | Planned, not implemented |

Here `fast_mass` means that the redshift contraction occurs outside the
radial operator; it does not imply frozen physics. The CUDA cell is labelled
frozen because it deliberately reproduces `ShearPrjFrozenPhysics.so`, whereas
the Python and C++ cells retain the exact redshift dependence. See
{doc}`../pipeline_organization` for the validation policy.

## Numerical framework

The full integral — per $(\lambda^{\rm ob}, z^{\rm ob}, R)$ wall point,
the projected surface density of line-of-sight structure:

$$\Delta\Sigma^{\rm prj}(R) = \int dz\, d\ln M\, d\theta\;
w_z(z, z^{\rm ob})\, \frac{dV}{d\Omega\,dz}\, n(M, z)\,
\big[\underbrace{1}_{\rm rnd} + \underbrace{b(M,z)\, b_{\rm sel}(\theta)\,
\xi_{\rm NL}(|\Delta\chi|, z^{\rm ob})}_{\rm cl}\big]\,
\Delta\Sigma_{\rm mis}\big(R \mid M,\, \theta D_A(z^{\rm ob})\big)\,
\mathbb{1}\big[\theta > \theta_{\rm excl}(z)\big],$$

then $\gamma_t^{\rm prj} = \Delta\Sigma^{\rm prj}\,
\langle\Sigma_{\rm crit}^{-1}\rangle(z^{\rm ob})$. The
population-integral framework — bin-averaged lensing predicted jointly
with the number counts from the same halo population and selection
kernels — is that of
[DES Cluster et al. 2023](https://ui.adsabs.harvard.edu/abs/2023arXiv230906593A/abstract)
(arXiv:[2309.06593](https://arxiv.org/abs/2309.06593)), the reference
paper for this software suite; the projection-selection extension is
the model of
[Costanzi et al. 2026, PhRvD 113, 103508](https://ui.adsabs.harvard.edu/abs/2026PhRvD.113j3508C/abstract)
(arXiv:[2604.05833](https://arxiv.org/abs/2604.05833)) built on top of
it: a scale-dependent parametrization of the optical cluster bias whose
small- and large-scale plateaus are set by the amplitude of projection
effects, with the two-halo profile expressed through off-axis
line-of-sight haloes.

The integrand uses the parabolic photo-$z$ kernel $w_z$, the analytic
$b_{\rm sel}(\theta) = B_{\rm small} + (B_{\rm large} - B_{\rm small})\,
\sigma(\theta)$ from the `bsel` plateaus, and the single-offset miscentred
NFW — **neighbouring-halo miscentering** in the paper's language: the
"offset" is the projected halo–halo separation $\theta D_A(z^{\rm ob})$
itself, a geometric ingredient with no free nuisance parameters (unlike
the target-cluster miscentering of {doc}`shear_halo`). The exclusion
$\theta > \theta_{\rm excl}(z)$ is a line-of-sight **slab** cut at the
redMaPPer aperture $R_\lambda(1+z^{\rm ob})$, not a 3-D ball.

Evaluation is entirely fixed grids and dot products — feature-pinned
log-GL in $\theta$, ring + foreground/background wings in $z$, fixed GL
in $\ln M$, with the random channel's redshift sum hoisted exactly and
the clustered channel frozen at $z^{\rm ob}$ with an $r_s(M)$-anchored
drift amplitude ($\sim 3.2\times$ faster than the full evaluator,
$< 0.2\%$ deviation). **The complete step-by-step recipe — grid
construction, exclusion mask, channel contractions, table lookups, cost
— lives in {doc}`../numerics/index`,
§"The shear-projection recipe, step by step".**

Model derivation: {doc}`../science/index` (projection lensing);
selection-bias inputs: {doc}`../selection/bsel`; full-fidelity and
adaptive validation backends (`ShearPrjEvaluator`, `ShearPrjGsl`,
`ShearPrjCuhre`): {doc}`../variants`.

## CosmoSIS setup

```ini
[shear_prj_frozen_physics]
file = ${Y3_CLUSTER_CPP_DIR}/release-build/src/modules/sigma_prj_cpu/ShearPrjFrozenPhysics.so
zt_low      = 0.10
zt_high     = 0.75
lnm_low     = 29.9336
lnm_high    = 35.6814
R_max_cMpch = 35.0
n_lnm       = 16
n_per_seg   = 10
n_zring     = 20
n_zouter    = 20
include_omega_z = 0
lambda_bin  = <180-entry wall: 15 radii per (richness, zob) bin>
zo_low      = <180 entries>
zo_high     = <180 entries>
radii       = <180 entries: 0.0426 … 24.8771 cMpc/h per bin>
```

(The four wall arrays are 180 entries each — 12 bins × 15 radii, bin
slow / radius fast; full arrays in
[the ini](https://github.com/estevesjh/des-nersc-cluster-scripts/blob/9fd24ddc075d394af4e20241bda716ac4d529fcb/cosmosis-models/mock_mcmc_buzzard.ini).)

- Ordering: **after `bsel`** (needs the bias plateaus), `halo_model`,
  `MfTinker`, `cp_camb`, `average_sigma_crit_inv`; before `likelihoods`.
- The ini section name must be `shear_prj_frozen_physics` — it is the
  hard-coded `module_label()`.

## Configuration options

| Option | Meaning | Units | Reference value |
|---|---|---|---|
| `lambda_bin`, `zo_low`, `zo_high`, `radii` | zipped wall grid (4 equal-length arrays) | $R$: cMpc/$h$ | 180 entries |
| `zt_low`, `zt_high` | true-redshift envelope of the line-of-sight integral | — | 0.10, 0.75 |
| `lnm_low`, `lnm_high` | mass GL limits | $\ln(M_\odot/h)$ | 29.9336, 35.6814 |
| `n_lnm` | GL nodes in $\ln M$ (class default 24) | — | 16 |
| `n_per_seg` | log-GL nodes per $\theta$ segment | — | 10 |
| `n_zring`, `n_zouter` | ring-band / per-wing redshift nodes | — | 20, 20 |
| `R_max_cMpch` | sets $\theta_{\max} = R_{\max}/D_A$ (class default 30) | cMpc/$h$ | 35.0 |
| `include_omega_z` | multiply the $z$ integrand by the survey area $\Omega(z)$ — class default **on**, deliberately set to **0** here | — | 0 |
| `lob_centers` | richness centres for $R_\lambda$, $\theta_\lambda$ | — | 25 37.5 52.5 130 (default) |

```{admonition} Why include_omega_z = 0
:class: important
$\Omega(z)$ belongs in cluster-count operators only; for a surface
density it cancels between numerator and normalisation, and the sibling
`ShearPrjEvaluator` hard-excludes it. The class default is *on* (matching
the `ShearPrjGsl` diagnostic), so the ini overrides it. Verified
empirically: with it on, fiducial self-closure breaks
($\log L = -151.7$ instead of $\approx 0$); off, closure holds
($-0.004$). See {doc}`../modules/survey_area`.
```

## DataBlock inputs

| DataBlock input | Meaning | Units / shape | Produced by |
|---|---|---|---|
| `mass_function/*`, `cluster_abundance/{hmf_s, hmf_q}` | HMF via `HMF_t` | — | `MfTinker`, values file |
| `haloModel/{lnM, z, bias}` | halo bias $b(M,z)$ | `(50, 100)` | `halo_model` |
| `xi_nl/{r, z, xi_nl}` | nonlinear correlation function | `(50, 128)` | `halo_model` |
| `distances/{z, d_a}` | volume element via `DV_DO_DZ_t` | Mpc | `cp_camb` |
| `distances/{z, d_c}` | comoving distance ($\times h_0$ → cMpc/$h$) | Mpc | `cp_camb` |
| `b_sel_marginalised/{lob, zob, b_small, b_large}` | selection-bias plateaus, interpolated linearly to the slice $z^{\rm ob}$ | `(4,)`, `(3,)`, `(3,4)`, `(3,4)` | `bsel` |
| `average_sigma_crit_inv/{zlense, sci_average}` | $\langle\Sigma_{\rm crit}^{-1}\rangle$ — **optional**: if absent, all $\gamma_t$ outputs are identically zero | `(50,)` | `average_sigma_crit_inv` |
| `cosmological_parameters/{omega_m, h0}` | $\bar\rho_m$ multiplier; unit conversion | scalars | `consistency` |

The photo-$z$ width $\sigma_z(z)$ comes from the compiled-in table
`src/models/z_kernel_data.hh`, not the DataBlock.

## DataBlock outputs

Nine arrays, each of length 180 (total = `rnd` + `cl`):

| DataBlock output | Meaning | Units / shape | Consumed by |
|---|---|---|---|
| `dsigma_prj_frozen_physics/{vals, rnd, cl}` | $\Delta\Sigma^{\rm prj}$: total, random-point channel, clustered channel | $M_\odot/(h\,\mathrm{pc}^2)$, `(180,)` | diagnostics |
| `shear_prj_frozen_physics/{vals, rnd, cl}` | $\gamma_t^{\rm prj} = \Delta\Sigma^{\rm prj}\,\langle\Sigma_{\rm crit}^{-1}\rangle(z^{\rm ob})$ | dimensionless, `(180,)` | — |
| `shear_prj/{vals, rnd, cl}` | **alias** of the row above, for drop-in compatibility with `likelihood_cp.py` / `generate_mock_dv.py` | `(180,)` | `likelihoods` |

Do not load this module together with `ShearPrjEvaluator` in one
pipeline: both write `shear_prj/*` and the later module wins.
