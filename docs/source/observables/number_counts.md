# Cluster Number Counts

`C++` · `y3_cluster_cpp` · `Cluster observable` · module `NumCountsSel`

Computes the expected cluster count $N_i[1]$ in each of the 12
$(\lambda^{\rm ob}, z^{\rm ob})$ bins by integrating the halo mass function
against the pre-tabulated selection function. Its output is the first block
of the theory vector and the denominator of the stacked one-halo shear.

## Script

- Model: [`src/models/n_operator_sel_gl_t.hh`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/models/n_operator_sel_gl_t.hh)
  (`y3_cluster::NumCountsSelGL` over `nosel_gl_detail::SelGLCore`).
- Module driver: [`src/modules/num_counts_sel/NumCounts.cc`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/modules/num_counts_sel/NumCounts.cc)
  (`DEFINE_COSMOSIS_SCALAR_EVALUATOR_MODULE`).
- Compiled library loaded by CosmoSIS:
  `${Y3_CLUSTER_CPP_DIR}/release-build/src/modules/num_counts_sel/NumCountsSel.so`.

## Numerical framework

The number-count forward model follows
[DES Cluster et al. 2023](https://ui.adsabs.harvard.edu/abs/2023arXiv230906593A/abstract)
(arXiv:[2309.06593](https://arxiv.org/abs/2309.06593)) — the reference
paper for this software suite: cluster counts predicted per
$(\lambda^{\rm ob}, z^{\rm ob})$ bin from the halo mass function, the
survey area, the volume element, and the richness–mass and photo-$z$
kernels (its Eq. 1). This module is the $f = 1$ instance of the
population operator

$$N_i[f] = \int d\ln M \int dz\;
\Omega(z)\,\frac{dV}{d\Omega\,dz}\,\frac{dn}{d\ln M}(M,z)\,
S_{ij}(\ln M, z)\, f(\ln M, z),$$

with $\Omega(z)$ the survey area ({doc}`../selection/survey_area`) and
$S_{ij}$ the tabulated selection function
({doc}`../selection/sel_function`). Evaluation is fixed Gauss–Legendre:
the $z$ axis is contracted once per sample into per-bin mass weights,
and each count is one 1-D mass sum — deterministic 0.02 s per sample vs
the retired adaptive Cuhre path (mean 0.11 s, tail 1 s),
grid-convergence error $< 0.05\%$. **The complete step-by-step recipe
lives in {doc}`../numerics/index`, §"The number-counts and one-halo
lensing recipe, step by step".** Derivation: {doc}`../science/index`.

## CosmoSIS setup

```ini
[NumCountsSel]
file = ${Y3_CLUSTER_CPP_DIR}/release-build/src/modules/num_counts_sel/NumCountsSel.so
algorithm = cuhre
eps_rel = 1.5e-3
eps_abs = 1.0e-12
max_eval = 1000000
use_cartesian_product = T
bin_index = 0 1 2 3 4 5 6 7 8 9 10 11
zt_low  = 0.05
zt_high = 0.80
lnm_low  = 29.9336
lnm_high = 36.7300
```

- Build once on Perlmutter (see {doc}`../installation`); requires
  `Y3_CLUSTER_CPP_DIR` so the `.so` and its data files resolve.
- Ordering: after `sel_function`, `MfTinker`, `cp_camb`.
- `algorithm`/`eps_*`/`max_eval`/`use_cartesian_product` are legacy Cuhre
  knobs, **ignored** by the current fixed Gauss–Legendre evaluator; they
  remain in the ini for backward compatibility.

## Configuration options

| Option | Meaning | Units | Reference value |
|---|---|---|---|
| `bin_index` | wall of bins to evaluate (richness index fastest) | — | `0 … 11` |
| `zt_low`, `zt_high` | true-redshift integration limits | — | 0.05, 0.80 |
| `lnm_low`, `lnm_high` | mass integration limits | $\ln(M_\odot/h)$ | 29.9336, 36.7300 |
| `n_lnm`, `n_z` | GL nodes in $\ln M$ / $z$ | — | 96, 64 (defaults) |

## DataBlock inputs

| DataBlock input | Meaning | Units / shape | Produced by |
|---|---|---|---|
| `sel_function/{lnM, z, S_stack}` | selection tensor $S_{ij}(\ln M, z)$ | `(12, 64, 192)` | `sel_function` |
| `mass_function/{m_h, z, dndlnmh}` | halo mass function (queried through `HMF_t`, which applies the $\Omega_m - \Omega_\nu$ mass-axis shift) | $h^3\,\mathrm{Mpc}^{-3}$ | `MfTinker` |
| `cluster_abundance/{hmf_s, hmf_q}` | HMF nuisance amplitudes | scalars | sampler (values file) |
| `distances/{z, d_a}` | comoving volume element via `DV_DO_DZ_t` | Mpc | `cp_camb` |
| `cosmological_parameters/{omega_m, omega_nu}` | HMF mass-axis shift | scalars | `consistency` |

## DataBlock outputs

| DataBlock output | Meaning | Units / shape | Consumed by |
|---|---|---|---|
| `numcountssel/vals` | expected counts $N_i[1]$ per bin | `(12,)` | `likelihoods` (data block and shear normalisation) |

The output section name is hard-coded in the module (deliberately not an
ini knob: a CosmoSIS `[DEFAULT]` block would propagate an
`output_section` value into every module and silently redirect writes).

