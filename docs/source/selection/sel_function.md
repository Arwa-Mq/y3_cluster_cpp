# Selection Function

`Python` · `y3_cluster_cpp` · `Selection` · module `sel_function` · `197 ms/sample`

Pre-tabulates, once per sample, the joint richness + photo-$z$ selection
function $S_{ij}(\ln M, z)$ for all 12 observed bins on one shared grid.
`NumCountsSel` and `Shear1hMisSel` slice their bin's plane from the packed
tensor and interpolate it inside their population integrals.

## Script

- Source: [`src/modules/sel_function/sel_function.py`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/modules/sel_function/sel_function.py)
  (`y3_cluster_cpp` @ `d7feb75`).
- Loaded by CosmoSIS as a Python module.

## Numerical framework

$$S_{ij}(\ln M, z) = \Big[\textstyle\sum_k W_k\,
K_i(\lambda_k, z)\, P_{\rm HOD}(\lambda_k \mid M, z)\Big]\, K_j(z),$$

where $K_i$ differences the EMG CDF of
$P(\lambda^{\rm ob}\mid\lambda^{\rm tr}, z)$ at the bin edges, $K_j$ is a
Gaussian photo-$z$ bin kernel of width `sigma_z`, and $P_{\rm HOD}$ is the
shifted-Poisson richness–mass relation — both defined in
{doc}`../observables/richness_mass`. In the paper's language: $K_i$ is
the **observed-richness kernel** $\mathcal{S}_i$, $K_j$ the
**observed-redshift kernel** $\mathcal{S}_j$, and the
$\lambda^{\rm tr}$-integrated product is the **richness selection
function** $S_i(M, z^{\rm tr})$ — together the redMaPPer selection
function $\mathcal{S}_{ij}$. The $\lambda^{\rm tr}$ sum uses $N_q$
GL nodes bracketed at $\mu_{\rm eff} \pm L\sigma_{\rm eff}$; the EMG CDF is
evaluated via `erfcx` at the 5 unique bin edges {20, 30, 45, 60, 200} only.
The Python kernels match the C++ models
(`src/models/mor_hod_t.hh`, `src/models/richness_kernel_t.hh`)
line-for-line. Full derivation: {doc}`../science/index`.

## CosmoSIS setup

```ini
[sel_function]
file = ${Y3_CLUSTER_CPP_DIR}/src/modules/sel_function/sel_function.py
lam_min = 20.0  30.0  45.0  60.0   20.0  30.0  45.0  60.0   20.0  30.0  45.0  60.0
lam_max = 30.0  45.0  60.0  200.0  30.0  45.0  60.0  200.0  30.0  45.0  60.0  200.0
zob_min = 0.20  0.20  0.20  0.20   0.35  0.35  0.35  0.35   0.50  0.50  0.50  0.50
zob_max = 0.35  0.35  0.35  0.35   0.50  0.50  0.50  0.50   0.65  0.65  0.65  0.65
sigma_z = 0.03  0.03  0.03  0.03   0.03  0.03  0.03  0.03   0.03  0.03  0.03  0.03
zt_low   = 0.05
zt_high  = 0.80
lnm_low  = 29.9336
lnm_high = 36.8414
n_lnm = 192
n_z   = 20
n_z_shared = 64
L_z    = 6.0
L_lam  = 6.0
N_q    = 32
```

- Requires `Y3_CLUSTER_CPP_DIR`.
- Ordering: after `consistency` (MOR parameters live in the sampled
  sections); before `NumCountsSel` and `Shear1hMisSel`.
- The 12 array entries define the bin wall: 4 richness bins
  $\{[20,30), [30,45), [45,60), [60,200)\}$ × 3 photo-$z$ bins
  $\{[0.20,0.35), [0.35,0.50), [0.50,0.65)\}$, richness index fastest.

## Configuration options

| Option | Meaning | Units | Reference value |
|---|---|---|---|
| `lam_min`, `lam_max` | per-bin observed-richness edges | — | 12-entry wall |
| `zob_min`, `zob_max` | per-bin observed-redshift edges | — | 12-entry wall |
| `sigma_z` | photo-$z$ scatter in the redshift kernel $K_j$ | — | 0.03 (all bins) |
| `zt_low`, `zt_high` | shared true-$z$ grid envelope | — | 0.05, 0.80 |
| `lnm_low`, `lnm_high` | shared $\ln M$ grid envelope | $\ln(M_\odot/h)$ | 29.9336, 36.8414 |
| `n_lnm` | $\ln M$ nodes — whole-pipeline optimum 192; **64 is pathological** (GL resonance, 4.5% drift on counts) | — | 192 |
| `n_z_shared` | shared $z$ nodes | — | 64 |
| `N_q` | Gauss–Legendre nodes of the $\lambda^{\rm tr}$ quadrature | — | 32 |
| `L_lam`, `L_z` | quadrature bracket half-widths | units of $\sigma$ | 6.0, 6.0 |

## DataBlock inputs

| DataBlock input | Meaning | Units / shape | Produced by |
|---|---|---|---|
| `cluster_mor/{log10_Mmin, log10_M1 \| log10_ratio, alpha, epsilon, sigma_lambda}` | shifted-Poisson HOD richness–mass parameters | — | sampler (values file) |
| `plob_ltr_params/*` | EMG projection-kernel coefficient splines (optional; falls back to the table embedded in `y3_buzzard/prj_params.py`) | — | retired `prj_params` module / in-code fallback |

## DataBlock outputs

| DataBlock output | Meaning | Units / shape | Consumed by |
|---|---|---|---|
| `sel_function/lnM` | shared mass grid | `(192,)` | `NumCountsSel`, `Shear1hMisSel` |
| `sel_function/z` | shared redshift grid | `(64,)` | same |
| `sel_function/S_stack` | packed selection tensor $S_{ij}(\ln M, z)$, layout `(bin, z, lnM)`, C-contiguous | `(12, 64, 192)` | same |

