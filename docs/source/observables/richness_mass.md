# Richness scaling relation

`C++ / Python` (model, not a module) · `y3_cluster_cpp` · `Halo population`

The richness–mass relation $P(\lambda^{\rm ob} \mid M, z)$ connects halo
mass to the observed richness. It is **not a pipeline module**: the model
is evaluated inside `sel_function`, `b_sel_marg`, and `bsel`, from
parameters the sampler writes into `cluster_mor/*`. It factorises as

$$P(\lambda^{\rm ob}\mid M, z) = \int d\lambda^{\rm tr}\,
\underbrace{P(\lambda^{\rm ob}\mid\lambda^{\rm tr}, z)}_{\text{projection kernel (EMG)}}\;
\underbrace{P(\lambda^{\rm tr}\mid M, z)}_{\text{intrinsic relation (shifted Poisson)}} .$$

## Script

- Intrinsic relation (C++): [`src/models/mor_hod_t.hh`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/models/mor_hod_t.hh)
  (`MOR_HOD_t`), used by `b_sel_marg`.
- Projection kernel (C++): [`src/models/richness_kernel_t.hh`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/models/richness_kernel_t.hh)
  (EMG CDF `F_EMG`), used by the selection kernels.
- Python mirrors (line-for-line): `_p_hod_scalar` and `_cdf_lob` in
  [`src/modules/sel_function/sel_function.py`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/src/modules/sel_function/sel_function.py);
  `_p_ltr_given_M` / `_p_lob_given_ltr_emg` in
  [`y3_buzzard/bsel.py`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/y3_buzzard/bsel.py).
- EMG coefficient calibration: [`y3_buzzard/prj_params.py`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/y3_buzzard/prj_params.py)
  (`PrjParams`, the Buzzard-calibrated spline table; the standalone
  `prj_params` DataBlock module is retired).

## Parameters (DataBlock)

Sampled parameters, written by the values file into `cluster_mor`:

| DataBlock key | Meaning | Units | Read by |
|---|---|---|---|
| `cluster_mor/log10_Mmin` | mass where $\mu_{\rm sat}$ turns on | $\log_{10} M_\odot/h$ | `sel_function`, `b_sel_marg`, `bsel` |
| `cluster_mor/log10_M1` or `log10_ratio` | satellite normalisation mass (or its ratio to `Mmin`) | $\log_{10} M_\odot/h$ | same |
| `cluster_mor/alpha` | satellite power-law slope | — | same |
| `cluster_mor/epsilon` | redshift evolution exponent | — | same |
| `cluster_mor/sigma_lambda` | intrinsic scatter parameter | — | same |
| `cluster_mor/z_pivot` | redshift pivot (optional) | — | same |

## Science and numerics

**Intrinsic relation — shifted-Poisson HOD.** With
$\nu = \mu_{\rm sat} + \delta$ and $\delta = (\sigma_\lambda \mu_{\rm sat})^2$:

$$P(\lambda^{\rm tr}\mid M, z) =
\exp\!\big[-\nu + (\lambda^{\rm tr} + \delta - 1)\ln\nu
- \ln\Gamma(\lambda^{\rm tr} + \delta)\big],$$

a continuous shifted-Poisson with mean $\mu_{\rm sat}(M, z)$ set by the
HOD parameters above (narrow-Gaussian fallback where
$\mu_{\rm sat} \le 10^{-8}$; `z_pivot` defaults to 0.45). When both
`log10_ratio` and `log10_M1` are present, `log10_ratio` wins:
$\log_{10} M_1 = \log_{10} M_{\rm min} + \texttt{log10\_ratio}$. Two
known small divergences between evaluators: `b_sel_marg` and `bsel` use
the satellite term without the central-galaxy shift $\lambda_{\rm cen}$,
and `bsel` defaults its pivot to 0.4544.

**Projection kernel — EMG.** $P(\lambda^{\rm ob}\mid\lambda^{\rm tr}, z)$
is a Gaussian core plus an exponentially-modified Gaussian (EMG) tail with
projection fraction $f_{\rm prj}$; its four coefficient functions
$(\mu, \sigma, \tau, f_{\rm prj})(\lambda^{\rm tr}, z)$ come from 8
Buzzard-calibrated splines in `PrjParams`. The selection kernels need only
its CDF, evaluated in closed form via `erfcx` at the bin edges.

Where the reference pipeline evaluates it:

- {doc}`sel_function <../selection/sel_function>` builds
  $S_{ij}(\ln M, z)$ from both pieces — this is how it enters
  {doc}`NumCountsSel <number_counts>` and
  {doc}`Shear1hMisSel <shear_halo>`;
- {doc}`b_sel_marg <../selection/bsel>` weights the $P[X]$ operators by
  the shifted-Poisson pdf;
- {doc}`bsel <../selection/bsel>` uses the EMG pdf in the
  $\lambda^{\rm tr}$ marginalisation.

Full derivation (log-normal alternative, EMG CDF closed form, quadrature
placement): {doc}`../science/index`.
