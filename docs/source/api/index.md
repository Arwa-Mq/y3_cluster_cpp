# API Reference

The complete inventory of the C++/CUDA source tree: every model header
and every CosmoSIS module directory, with main types, role, and status.

(src-models)=
## `src/models`

Every header under
[`src/models/`](https://github.com/estevesjh/y3_cluster_cpp/tree/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/models)
— the header-only model components composed into the CosmoSIS modules
(the "models" of DES Cluster et al. 2023). A `.cuh` twin means a CUDA
device implementation of the same model exists. Grouped by status.

### Reference path (used by `mock_mcmc_buzzard.ini`)

| Header | Main types | Role | CUDA twin |
|---|---|---|---|
| [`mor_hod_t.hh`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/models/mor_hod_t.hh) | `MOR_HOD_t` | shifted-Poisson richness–mass relation $P(\lambda^{\rm tr}\mid M,z)$ ({doc}`../modules/richness_mass`) | — |
| [`richness_kernel_t.hh`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/models/richness_kernel_t.hh) | `RichnessKernel_t`, `richness_zkernel` | observed-richness kernel $\mathcal S_i$ (Gaussian+EMG CDF) and observed-redshift kernel $\mathcal S_j$ | — |
| [`z_kernel_data.hh`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/models/z_kernel_data.hh) | (data tables) | compiled photo-$z$ scatter $\sigma_z(z)$ ({doc}`../modules/redshift_kernel`) | — |
| [`hmf_t.hh`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/models/hmf_t.hh) | `HMF_t` | HMF reader: $\Omega_m-\Omega_\nu$ mass-axis shift + $(s,q)$ nuisance scaling | ✓ |
| [`dv_do_dz_t.hh`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/models/dv_do_dz_t.hh) | `DV_DO_DZ_t` | comoving volume element from `distances/d_a` | ✓ |
| [`ez.hh`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/models/ez.hh) / `ez_sq.hh` | `EZ`, `EZ_sq` | $E(z)$ and $E^2(z)$ expansion factors | ✓ |
| [`omega_z_des.hh`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/models/omega_z_des.hh) | `OMEGA_Z_DES` | hard-coded survey area $\Omega(z)$ ({doc}`../modules/survey_area`) | ✓ |
| [`sel_function_t.hh`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/models/sel_function_t.hh) | `SelFunction_t` | `sel_function/S_stack` reader + per-bin `Interp2D` slice | — |
| [`n_operator_sel_gl_t.hh`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/models/n_operator_sel_gl_t.hh) | `SelGLCore`, `NumCountsSelGL`, `Shear1hMisSelGL` | fixed-GL $N_i[f]$ engine ({doc}`../numerics/index`) | — |
| [`p_operator_t.hh`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/models/p_operator_t.hh) | `P_operator` | the $\mathcal P[X]$ selection-bias operators $(P_1, I_1, J)$ | — |
| [`sigma_prj_t.hh`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/models/sigma_prj_t.hh) | `ShearPrjCore`, `ThetaGrid`, `Sigma/DSigma/ShearPrjEvaluator` | full projection evaluator + shared $\theta$/$z$ grid builders | — |
| [`sigma_prj_frozen_t.hh`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/models/sigma_prj_frozen_t.hh) | `ShearPrjFrozenPhysics` | frozen-physics projection stage of the reference run ({doc}`../observables/shear_projection`) | — |
| [`nfw_dsigma_mis.hh`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/models/nfw_dsigma_mis.hh) / `nfw_sigma_mis.hh` | `NFW_DSIGMA_MIS`, `NFW_SIGMA_MIS` | offset-NFW table lookups (gamma and single kernels, $c=4$) | ✓ |
| [`models.hh`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/models/models.hh) / `default_models.hh` | `Models<...>` | the model-typelist pattern + per-survey aliases | — |

### Validation and diagnostic backends

| Header | Main types | Role |
|---|---|---|
| `sigma_prj_frozen_interp_t.hh` | `ShearPrjFrozenCuhre` | Option E: frozen physics driven by continuous Cuhre |
| `p_operator_cuhre_t.hh` | `P_operator_cuhre`, `BSelMarg*PaganiWeight` | adaptive/PAGANI $\mathcal P[X]$ reference benchmarks |
| `n_operator_sel_t.hh` | `NOperatorSelScalar/Radial` | retired Cuhre $N_i[f]$ engine — still drives `MassWeightedSel`, `BiasWeightedSel`, `Shear1hSel` |
| `projection_y3_b_i_t.hh` | `RmSelFunction_t` | brute-force $B_i$ selection kernels (`NumCountsFullScalarIntegrand`) |

### 1h+2h composition models (variants)

| Header | Main types | Role | CUDA twin |
|---|---|---|---|
| `kappa_max.hh` | `KAPPA_MAX` | $\kappa_{\max}$: $\max(\Sigma_{\rm NFW}, b\,\Sigma_{\rm hh})\,\Sigma_{\rm crit}^{-1}$ ({doc}`../observables/second_halo_term`) | ✓ |
| `gamma_max.hh` | `GAMMA_MAX` | $\Delta\Sigma$ twin of `KAPPA_MAX` | ✓ |
| `sig_max.hh` / `sig_sum.hh` | `SIG_MAX`, `SIG_SUM` | $\Sigma$-level max / sum compositions | ✓ |
| `xi_max.hh` / `xi_sum.hh` | `XI_MAX`, `XI_SUM` | $\xi$-level max / sum compositions | — |
| `del_sig_sum.hh` | `DEL_SIG_TOM` (sum variant) | $\Delta\Sigma$ 1h+2h sum | — |
| `xi_to_sigma_t.hh` | `SIG_y3` | $\xi \to \Sigma$ Abel projection | — |

### Legacy richness / photo-z / MOR kernels (Y1, SDSS, mock eras)

| Header | Main types | Role | CUDA twin |
|---|---|---|---|
| `lc_lt_t.hh` (+`.cc`) / `lc_lt_y1_t.hh` | `LC_LT_t`, `LC_LT_Y1_t` | $P(\lambda^{\rm cen}\mid\lambda^{\rm tr})$ kernels | ✓ |
| `lo_lc_t.hh` | `LO_LC_t` | $P(\lambda^{\rm ob}\mid\lambda^{\rm cen}, R_{\rm mis})$ kernel | ✓ |
| `int_lc_lt_des_t.hh` / `_t2.hh` (+`_t3.cuh`, `_f.cuh`) | `INT_LC_LT_DES_t*` | pre-integrated DES richness-kernel tables | ✓ |
| `int_zo_zt_des_t.hh` / `int_zo_zt_t.hh` / `zo_zt_des_t.hh` | `INT_ZO_ZT_*`, `ZO_ZT_DES_t` | legacy photo-$z$ kernel integrals | ✓ |
| `plob_ltr_emg_t.hh` | `PlobLtrEMG_t` | EMG $P(\lambda^{\rm ob}\mid\lambda^{\rm tr})$ from DataBlock splines (retired loader path) | — |
| `mor_t.hh` / `mor_des_t.hh` / `mor_des_log_t.hh` / `mor_des_2022.hh` / `mor_des_2022_logm.hh` / `mor_exp_t.hh` / `mor_sdss_t.hh` (+`mor_des_triax.cuh`) | `MOR_*` | log-normal / exponential / survey-specific mass–richness variants | ✓ |
| `sigma_photoz_des.hh` | `SIGMA_PHOTOZ_DES_t` | legacy DES photo-$z$ scatter model | ✓ |

### Legacy lensing / profile / power-spectrum models

| Header | Main types | Role | CUDA twin |
|---|---|---|---|
| `del_sig_t.hh` / `del_sig_tom.hh` / `del_sig_y1.hh` | `DEL_SIG_*` | $\Delta\Sigma$ profile variants (Y1/Tom eras) | — |
| `dsigma_misc.hh` (+`dsigma_full.cuh`, `dsigma_proj.cuh`, `sigma_misc.cuh`, `sigma_mis_joint.cuh`) | `DSIGMA_*`, `SIGMA_*` | miscentred/full/projected $\Sigma$, $\Delta\Sigma$ integrands | ✓ |
| `nfw_sigma_t.hh` / `nfw_xi_t.hh` | `nfw_sigma`, `nfw_xi` | analytic NFW $\Sigma$ and $\xi$ | — |
| `xinl.hh` / `pk_nl_t.hh` (+`pk_cluster_damp.cuh`, `wp_cluster.cuh`) | `xi_nl`, `pk_nl`, … | $\xi_{\rm NL}$ / $P_{\rm NL}$ readers, damped-$P_k$ and $w_p$ experiments | ✓ |
| `roffset_t.hh` / `t_cen_t.hh` / `t_mis_t.hh` / `a_cen_t.hh` / `a_mis_t.hh` | `ROFFSET_t`, `T/A_CEN/MIS_t` | legacy miscentering offset distribution + centred/miscentred kernels and amplitudes | partial |

### Misc / retired

| Header | Main types | Role |
|---|---|---|
| `hmb_t.hh` | `HMB_t` | **retired** halo-bias type — $b(M,z)$ now read from `haloModel/bias` ({doc}`../cosmology/halo_bias`) |
| `sigma_crit_inverse_t.hh` / `weighted_sigma_crit_inv.hh` / `average_sci_t.hh` | `sigma_crit_inv`, … | $\Sigma_{\rm crit}^{-1}$ variants (the pipeline now uses the Python module, {doc}`../cosmology/sigma_crit_inv`) |
| `pzsource_t.hh` / `pzsource_gaussian_t.hh` | `PZSOURCE_*` | source $p(z)$ models |
| `angle_to_dist_t.hh` | `ANGLE_TO_DIST_t` | $\theta \to R$ conversion |
| `sample_variance.hh` (+`.cc`) | `SampleVariance_t` | abundance sample-variance model |
| `y_sz.hh` | `Y_SZ` | Compton-$y$ scaling (DES×SPT era) |
| `op_sel_park.hh` | `OP_SEL_PARK` | Park-style selection operator | 
| `triax_richness.cuh` / `interpolation_tables.cuh` | `TRIAX_RICHNESS`, … | CUDA-only triaxiality experiment / device table helpers |

(src-pipelines-des-y3)=
## `src/pipelines/des_y3`

The additive maintained namespace for new DES Y3 implementations. Its layout
is observable → integration strategy → backend; it does not relocate the
production modules listed in [src/modules](#src-modules). See
{doc}`../pipeline_organization` for the strategy definitions, status matrix,
and compatibility rules.

### Shared Python layer

| File | Role |
|---|---|
| `shared/datablock_models.py` | Convention-matched HMF, volume, survey-area, selection-weight, and DataBlock-dump adapters |
| `shared/full_ltmz_core.py` | Explicit-selection `full_ltmz` contraction shared by counts and shear |
| `shared/lensing_profiles.py` | Production-compatible centred and miscentred profile readers |
| `shared/sel_function.py` | Staged maintained selection module; currently identical to the production entry point |
| `shared/sel_kernels.py` | Cached loader and HOD/richness-kernel adapters for the shared selection module |
| `shared/z_kernel.py` | Projection photo-$z$ kernel conventions |

### Observable products

| Directory | Products |
|---|---|
| `observables/number_counts/full_ltmz` | Python reference, `NumCountsFullLtmz.so`, `NumCountsFullLtmzGpu.so` |
| `observables/number_counts/fast_mass` | Python production-algorithm replica |
| `observables/shear_1h2h/full_ltmz` | Python reference, `Shear1hFullLtmz.so`, `Shear1hFullLtmzGpu.so` |
| `observables/shear_1h2h/fast_mass` | Python references, `Shear1hFastMass.so`, `Shear1h2hMax.so`, `Shear1h2hMaxGpu.so` |
| `observables/shear_1h2h/radial_series` | Python offline generator/evaluator and `Shear1hRadialSeries.so` |
| `observables/shear_projection/full_ltmz` | `DSigmaPrjFullLtmzGpu.so` |
| `observables/shear_projection/fast_mass` | Python exact-$z$ reference, `ShearPrjFastMass.so`, `ShearPrjFrozenGpu.so` |

(src-modules)=
## `src/modules`

Every directory under
[`src/modules/`](https://github.com/estevesjh/y3_cluster_cpp/tree/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/modules)
— the CosmoSIS-facing wrappers around the [models](#src-models).
`src/modules/CMakeLists.txt` is the build registry; Python steps are
loaded by file path from the ini and need no registration.

### Reference path

| Directory | Products | Language | Docs |
|---|---|---|---|
| `num_counts_sel` | `NumCountsSel.so`, `Shear1hMisSel.so` (+ `Shear1hSel.so`, `MassWeightedSel.so`, `BiasWeightedSel.so`) | C++ | {doc}`../observables/number_counts`, {doc}`../observables/shear_halo` |
| `b_sel_marg_cpu` | `BSelMargIntegrand.so` (+ `P1/I1/I2PaganiIntegrand.so` benchmarks) | C++ (+CUDA) | {doc}`../selection/bsel` |
| `sigma_prj_cpu` | `ShearPrjFrozenPhysics.so` (+ `ShearPrjEvaluator.so`, `SigmaPrjEvaluator.so`, `DSigmaPrjEvaluator.so`, `ShearPrjGsl.so`, `ShearPrjCuhre.so`, `ShearPrjFrozenCuhre.so`) | C++ | {doc}`../observables/shear_projection`, {doc}`../variants` |
| `cp_camb` | `cp_camb.py` (unregistered — loaded by path) | Python | {doc}`../cosmology/cp_camb` |
| `sel_function` | `sel_function.py` (unregistered) | Python | {doc}`../selection/sel_function` |
| `average_sigma_crit_inv` | `average_sigma_crit_inv.py` (unregistered) | Python | {doc}`../cosmology/sigma_crit_inv` |

(The remaining reference-path Python steps — `halo_model_cosmosis.py`,
`bsel.py`, `likelihood_cp.py`, `prj_params.py` — live in `y3_buzzard/`,
not `src/modules/`.)

### Registered validation / diagnostics

| Directory | Products | Role |
|---|---|---|
| `num_counts_full` | `NumCountsFullScalarIntegrand.so` | brute-force triple-integral counts cross-check (Cuhre) |
| `cluster_abundance_covariance` | abundance covariance module | sample-variance covariance |
| `ExampleScalar` / `ExampleVector` / `ExampleOneD` | example `.so` | macro-pattern templates for new modules |

### Registered historical (earlier analyses, kept building)

| Directory | Era / purpose |
|---|---|
| `gt_mock_cpu` / `gt_mock_gpu` | mock-catalogue $\gamma_t$ suites (avgCent/avgMisc) |
| `gt_card_cpu` / `gt_card_gpu` / `gt_card_triax_gpu` | Cardinal-simulation $\gamma_t$ (+ triaxial GPU) |
| `gt_park_sel_cpu` / `gt_park_sel_gpu` | Park-style selection $\gamma_t$ |
| `sigma_park_y1` | GPU $\Sigma$ (Y1) — **broken**, awaits `lc_lt.cuh` debugging |
| `sigma_buzzard_y3` | GPU $\Sigma$ on Buzzard |
| `sigma_mort_y1` / `mass_mort_y1` / `mass_y1` | Y1 $\Sigma$/mass integrands |
| `y1_analysis` | Y1 analysis integrands (incl. `_mor_2022` variants) |
| `sdss_analysis` | SDSS analysis integrands |
| `snapshotsim` | snapshot-simulation NC/$\Sigma$ integrands |
| `buzzard_test` | `buzzard_sigma_halos` test module |

### Commented out of the registry

`DESxSPTModule`, `sigma_kappa_y1`, `compton_y_sims`, `model_sigmahm`.

### On disk, unregistered orphans

| Directory | Note |
|---|---|
| `cuda` | CUDA experiments |
| `deltasigma` | superseded $\Delta\Sigma$ module |
| `finish` | pipeline finisher experiment |
| `mass_conversion` | mass-definition conversion |
| `n_operator_ratios` | $\langle M\rangle_i$, $\langle b\rangle_i$ finisher (optional diagnostics) |
| `parabola` | toy example |
| `prj_lens_model` | earlier projection-lensing model |
| `red_shear_prj` | legacy home of the projection branch — physics now in `sigma_prj_cpu` |

Status details and the known-broken legacy paths:
{doc}`../modules/historical`.
