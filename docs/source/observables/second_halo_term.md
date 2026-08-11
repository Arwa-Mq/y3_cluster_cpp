# Second Halo Term

`Python` (producer) · `y3_cluster_cpp` (`y3_buzzard/`) · `Lensing`

The conventional two-halo term — correlated matter around the cluster,
scaled by its plain halo bias — is, in the paper's language, the
**unselected-bias limit** of the selection-affected two-halo term
$\Sigma^{\rm prj}$ used by the reference pipeline
({doc}`shear_projection`): replace $b(M,z)\,b_{\rm sel}(\theta)$ by the
halo-bias aggregate $b_{\rm halo}$ and drop the selection modulation.
Its ingredients are computed by the
{doc}`halo_model <../cosmology/halo_model>` module, gated behind
`compute_lensing_2h`.

```{admonition} Not active in the reference pipeline
:class: important
`mock_mcmc_buzzard.ini` runs `halo_model` with `compute_lensing_2h = F`:
the reference shear composition is one-halo + projection
({doc}`shear_projection`), which never reads the two-halo tables.
Skipping the branch saves 200–300 ms per sample. The conventional
$1h{+}2h$ composition is retained as a comparison variant — see
{doc}`../variants`.
```

## Script

- Producer: [`y3_buzzard/halo_model_cosmosis.py`](https://github.com/estevesjh/y3_cluster_cpp/blob/d7feb7504ed5dfcad84f99a1791af8a55c858aa0/y3_buzzard/halo_model_cosmosis.py)
  with the Hankel-transform backend in `y3_buzzard/haloModel.py`
  (class `ct_2hTerm`, driving the `cluster_toolkit`
  $P \to \xi \to \Sigma \to \Delta\Sigma$ chain).

## Numerical framework

Around a halo of mass $M$, the correlated-matter surface density is the
matter correlation projected along the line of sight, scaled by the halo
bias:

$$\Delta\Sigma_{2h}(R; M, z) = b(M, z)\,
\bar\rho_m \int dz_\parallel\; \xi_{\rm NL}\!\big(\sqrt{R^2 + z_\parallel^2},\, z\big)
\;\longrightarrow\; \Delta\Sigma \text{ via the Abel/Hankel chain},$$

with $\bar\rho_m = \Omega_m \rho_{\rm crit}$. The stacked variant uses the
population-averaged bias $\langle b\rangle_i = N_i[b]/N_i[1]$ from the
`BiasWeightedSel` diagnostic module ({doc}`../variants`).

At the fiducial point the projection-aware reference composition sits
$1.3$–$2.5\times$ **above** the conventional
$\langle b\rangle\,\Delta\Sigma_{2h}$ term at $R = 5\,h^{-1}$cMpc (the
$b_{\rm sel}$-boosted line-of-sight contribution grows toward low richness
and high redshift), and $\approx 0.81\times$ below it at
$R = 0.2\,h^{-1}$cMpc (miscentering suppression):

```{image} ../_static/img/dsigma_compositions.png
:alt: Stacked lensing compositions at the fiducial point
:width: 100%
```

```{warning}
The wired legacy 1h+2h composition modules (`SigmaTotSel` /
`DSigmaTotSel`) are **currently broken** (NaN `dSigma_hh` below
$R \approx 8.6\,h^{-1}$cMpc; interpolation on the wrong radial axis).
Any 1h+2h comparison must assemble the term from `BiasWeightedSel` +
`xi_nl` directly. Details: {doc}`../modules/historical`.
```

Full derivation and the composition comparison:
{doc}`../science/index`.

## CosmoSIS setup

Enable inside the `[halo_model]` section (everything else as in the
reference — see {doc}`../cosmology/halo_model`):

```ini
[halo_model]
compute_lensing_1h = T
compute_lensing_2h = T   ; reference run sets F
```

## DataBlock outputs (when enabled)

| DataBlock output | Meaning | Units / shape | Consumed by |
|---|---|---|---|
| `haloModel/Rp` | radius grid of the 2h tables | cMpc/$h$, `(128,)` | 1h+2h variant assembly |
| `haloModel/Wp_hh` | matter two-point correlation table of the chain (despite the $W_p$ name) | `(50, 128)` | same |
| `haloModel/Sigma_hh` | two-halo surface density, **bias not applied** (consumer multiplies by $b$ or $\langle b\rangle_i$) | `(50, 128)` | same |
| `haloModel/dSigma_hh` | two-halo excess surface density, bias not applied | `(50, 128)` | same |

