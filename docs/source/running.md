# Running the reference pipeline

The reference configuration of the DES Y3 cluster-cosmology analysis is
[`cosmosis-models/mock_mcmc_buzzard.ini`](https://github.com/estevesjh/des-nersc-cluster-scripts/blob/9fd24ddc075d394af4e20241bda716ac4d529fcb/cosmosis-models/mock_mcmc_buzzard.ini)
in the **des-nersc-cluster-scripts** repository
(branch `polychord-widePlanck-logspace-ab`, @ `9fd24dd`): the Buzzard
convergence-test pipeline, fitting 12 number counts and 180 tangential-shear
points against the Buzzard simulation data vector.

The software suite this pipeline is built from — the CosmoSIS module
pattern, the model/integrand separation, and the number-count and
population-averaged lensing forward model — is described in
[DES Cluster et al. 2023](https://ui.adsabs.harvard.edu/abs/2023arXiv230906593A/abstract)
(arXiv:[2309.06593](https://arxiv.org/abs/2309.06593)), the main
reference for this documentation.

## The pipeline

```ini
[pipeline]
modules = consistency GrowthFactor cp_camb MfTinker halo_model
          average_sigma_crit_inv sel_function
          NumCountsSel Shear1hMisSel
          b_sel_marg bsel
          shear_prj_frozen_physics
          likelihoods
values = ${DES_CLUSTER_NERSC_DIR}/cosmosis-models/mock_mcmc_widePlanck_values.ini
likelihoods = likelihoods
```

| # | Module | What it computes | Language · source |
|---|---|---|---|
| 1 | {doc}`consistency <cosmology/consistency>` | completes the cosmological parameter set | Python · CosmoSIS Standard Library |
| 2 | {doc}`GrowthFactor <cosmology/growth_factor>` | linear growth $D(z)$, $f(z)$ | C · CosmoSIS Standard Library |
| 3 | {doc}`cp_camb <cosmology/cp_camb>` | linear $P(k,z)$ (CosmoPower emulator) + distances | Python · `y3_cluster_cpp` |
| 4 | {doc}`MfTinker <cosmology/mf_tinker>` | Tinker halo mass function | Fortran · CosmoSIS Standard Library |
| 5 | {doc}`halo_model <cosmology/halo_model>` | Tinker bias $b(M,z)$, $\xi_{\rm NL}$, NFW lensing tables | Python · `y3_cluster_cpp` |
| 6 | {doc}`average_sigma_crit_inv <cosmology/sigma_crit_inv>` | $\langle\Sigma_{\rm crit}^{-1}\rangle(z_l)$ | Python · `y3_cluster_cpp` |
| 7 | {doc}`sel_function <selection/sel_function>` | selection tensor $S_{ij}(\ln M, z)$ | Python · `y3_cluster_cpp` |
| 8 | {doc}`NumCountsSel <observables/number_counts>` | cluster counts $N_i[1]$ | C++ · `y3_cluster_cpp` |
| 9 | {doc}`Shear1hMisSel <observables/shear_halo>` | one-halo shear with miscentering | C++ · `y3_cluster_cpp` |
| 10 | {doc}`b_sel_marg <selection/bsel>` | selection-bias operators $(P_1, I_1, J)$ | C++ · `y3_cluster_cpp` |
| 11 | {doc}`bsel <selection/bsel>` | bias plateaus $(B_{\rm small}, B_{\rm large})$ | Python · `y3_cluster_cpp` |
| 12 | {doc}`shear_prj_frozen_physics <observables/shear_projection>` | projection shear $\gamma_t^{\rm prj}(R)$ | C++ · `y3_cluster_cpp` |
| 13 | {doc}`likelihoods <observables/likelihood>` | Gaussian $\log L$ | Python · `y3_cluster_cpp` |

Data flow (edge labels are the DataBlock sections passed between
modules; blue = cosmology quantities, orange = selection effects,
green = cluster observables, grey = likelihood):

```{image} _static/img/pipeline_dataflow.png
:alt: Data flow of the mock_mcmc_buzzard.ini reference pipeline
:width: 100%
```

(Source: `docs/figs/pipeline_dataflow.mmd`; regenerate the PNG with
`npx -y @mermaid-js/mermaid-cli -i docs/figs/pipeline_dataflow.mmd -o
docs/source/_static/img/pipeline_dataflow.png -b white -s 2`.)

## Required repositories and environment

| Repository | Role | Env variable |
|---|---|---|
| [y3_cluster_cpp](https://github.com/estevesjh/y3_cluster_cpp) | C++/CUDA modules (built once, {doc}`installation`) + Python modules | `Y3_CLUSTER_CPP_DIR=/pscratch/sd/j/jesteves/y3_cluster_cpp` |
| [des-nersc-cluster-scripts](https://github.com/estevesjh/des-nersc-cluster-scripts) | this ini, values, data vectors, sbatch scripts (deployed as `des-cluster-nersc`) | `DES_CLUSTER_NERSC_DIR=/pscratch/sd/j/jesteves/github/des-cluster-nersc` |
| [cosmosis-standard-library](https://github.com/joezuntz/cosmosis-standard-library) | `consistency`, `GrowthFactor`, `MfTinker` | `COSMOSIS_STANDARD_LIBRARY` |
| [camb-emulator](https://github.com/estevesjh/camb-emulator) | trained CosmoPower emulators read by `cp_camb` | (paths in the `[cp_camb]` section) |

`fast-cpu/setup_env.sh` (sourced by every job script) performs the
Perlmutter module swaps and exports; the CosmoSIS environment comes from
`setup-cosmosis-nersc` with the `y3cl_je` conda env:

```bash
module swap cudatoolkit/12.9 cudatoolkit/12.2
module swap gcc-native/13.2 gcc-native/12.3
export Y3_CLUSTER_CPP_DIR=/pscratch/sd/j/jesteves/y3_cluster_cpp
export DES_CLUSTER_NERSC_DIR=/pscratch/sd/j/jesteves/github/des-cluster-nersc
export COSMOSIS_STANDARD_LIBRARY=/global/common/software/des/jesteves/cosmosis-standard-library
export PYTHONPATH=${Y3_CLUSTER_CPP_DIR}:${PYTHONPATH:-}
export OMP_NUM_THREADS=1
source ${COSMOSIS_REPO_DIR}/setup-cosmosis-nersc \
       /global/common/software/des/common/Conda_Envs/y3cl_je
```

## Values, priors, data vector, covariance

- **Values**:
  [`mock_mcmc_widePlanck_values.ini`](https://github.com/estevesjh/des-nersc-cluster-scripts/blob/9fd24ddc075d394af4e20241bda716ac4d529fcb/cosmosis-models/mock_mcmc_widePlanck_values.ini)
  — 10 varied parameters, flat priors from the `[min start max]` boxes
  (there is **no separate priors file**):

  | Parameter | min | start | max |
  |---|---|---|---|
  | `cosmological_parameters/h0` | 0.473 | 0.6766 | 0.873 |
  | `cosmological_parameters/omega_m` | 0.11 | 0.311049 | 1.0 |
  | `cosmological_parameters/omega_b` | 0.02 | 0.048975 | 0.10 |
  | `cosmological_parameters/n_s` | 0.8 | 0.9665 | 1.15 |
  | `cosmological_parameters/sigma8` | 0.5 | 0.8238 | 1.5 |
  | `cluster_mor/log10_Mmin` | 10.0 | 11.4 | 13.0 |
  | `cluster_mor/log10_ratio` | 1.0 | 1.3 | 1.5 |
  | `cluster_mor/alpha` | 0.4 | 0.86 | 1.4 |
  | `cluster_mor/epsilon` | −1.0 | 0.0 | 1.0 |
  | `cluster_mor/sigma_lambda` | 0.05 | 0.18 | 0.50 |

  Fixed: `mnu = 0` (the emulator's growth rescaling assumes
  scale-independent $D(z)$), `w = -1`, `wa = 0`,
  `cluster_abundance/{hmf_s = 0, hmf_q = 1}`, `photoz/delta_z = 0`.
  `log10_ratio = \log_{10}(M_1/M_{\rm min})` replaces `log10_M1` to break
  the box-prior degeneracy. No `[miscentering]` section →
  `Shear1hMisSel` uses its in-code defaults
  $(f_{\rm mis}, \tau_{\rm mis}) = (0.22, 0.17)$.

- **Data vector**: `${DES_CLUSTER_NERSC_DIR}/data/mock/mock_dv_buzzard.npz`
  — the Buzzard simulation measurement (committed in the repo; its
  builder `validations/build_buzzard_datavector.py` is referenced by the
  driver script but not currently committed). Keys: `data_NC (12,)`,
  `data_Shear (180,)` (12 bins × **15 radii**,
  $R \in [0.0426, 24.877]$ cMpc/$h$, bin-major), `invcov_NC (12, 12)`,
  `invcov_Shear (180, 180)`, plus the embedded fiducial truth
  (`fiducial_param_names/values`).

- **Covariance**: the shear inverse covariance is the DES Y1 WL layout
  (`y1_rerun/data_files/wl_cov.txt`, 180 × 180); the number-count inverse
  covariance is the full Y1 matrix (`Cov_ij_bestfit_DESY1_105.txt`), not
  diagonal Poisson.

## Commands

Smoke test (single sample, `test` sampler — logL is finite but **not**
$\approx 0$: the Buzzard run is a recovery test against an external
simulation, not a self-closure):

```bash
cd ${DES_CLUSTER_NERSC_DIR}
source fast-cpu/setup_env.sh
srun -n 1 cosmosis cosmosis-models/mock_mcmc_buzzard.ini -p runtime.sampler=test
```

The committed end-to-end driver is
[`fast-cpu/build_buzzard_dv.sh`](https://github.com/estevesjh/des-nersc-cluster-scripts/blob/9fd24ddc075d394af4e20241bda716ac4d529fcb/fast-cpu/build_buzzard_dv.sh)
(debug QOS, 25 min): builds the data vector, asserts its shapes, and runs
the smoke test.

Production (PolyChord, 64 MPI ranks, shared QOS, ~9 h — the pattern of
`fast-cpu/mock_polychord.sh` with this ini; the completed production run
is job 56588673):

```bash
OUTDIR=/pscratch/sd/j/jesteves/cluster_lib/chains/mock_mcmc/buzzard/polychord
mkdir -p ${OUTDIR}/clusters
srun -n 64 cosmosis --mpi cosmosis-models/mock_mcmc_buzzard.ini \
     -p runtime.sampler=polychord runtime.resume=F \
        polychord.live_points=500 polychord.base_dir=${OUTDIR} \
        output.filename=${OUTDIR}/chain.txt
```

```{warning}
Three ini gotchas when launching production runs:

1. `[runtime]` defaults to `sampler = emcee`, `resume = T` — a fresh
   PolyChord run must override both on the command line.
2. The `[polychord]` section of the Buzzard ini uses non-CosmoSIS key
   names (`nlive`, `nrepeat`, …); the CosmoSIS names are `live_points`,
   `num_repeats`, `tolerance`, `base_dir`, `random_seed`. Pass them via
   `-p` as above, and `mkdir -p ${OUTDIR}/clusters` first — without
   `base_dir` PolyChord aborts with "Cannot open file ./clusters/…".
3. Under MPI on Lustre add `-p output.lock=F` (all ranks race for the
   POSIX lock on the chain file otherwise).
```

## Sampler and output configuration

```ini
[runtime]
sampler = emcee            ; overridden per run: test | apriori | emcee | polychord
root = ${COSMOSIS_SRC_DIR}
resume = T

[emcee]
walkers = 64               ; pair with `cosmosis --smp=64`
samples = 10000
nsteps = 50

[output]
filename = /pscratch/sd/j/jesteves/cluster_lib/chains/mock_mcmc/buzzard/chain.txt
format = text
```

Chain columns: the 10 varied parameters in values-file order, then
`prior like post weight`. Post-processing notebooks and convergence plots
live in `des-nersc-cluster-scripts/chains/`.

## What the likelihood compares

$$\log L = -\tfrac12\Big[
\delta_{\rm NC}^{\mathsf T} C_{\rm NC}^{-1} \delta_{\rm NC}
+ \delta_{\gamma}^{\mathsf T} C_{\gamma}^{-1} \delta_{\gamma}\Big],
\qquad
\gamma_t^{\rm theory}(R \mid i) =
\frac{N_i[\gamma_t^{1h,\rm full}](R)}{N_i[1]} + \gamma_t^{\rm prj}(R \mid i),$$

12 number counts + 180 shear points. Details: {doc}`observables/likelihood`
— including one open shape-contract caveat: the committed
`likelihood_cp.py` still hard-codes the 10-radii (120-point) layout of
the widePlanck variant and must be set to `_SHEAR_N_R = 15` for the
Buzzard data vector. Pipeline variants (self-closure widePlanck run,
10-radii grid, conventional $1h{+}2h$): {doc}`variants`.
