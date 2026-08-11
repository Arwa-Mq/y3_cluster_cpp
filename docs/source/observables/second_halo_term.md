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

The full integral: around a halo of mass $M$, the two-halo surface
density is the line-of-sight projection of the halo–matter correlation,

$$\Sigma_{2h}(R; M, z) = \bar\rho_m \int d\chi_\parallel\;
\xi_{\rm hm}\!\big(\sqrt{R^2 + \chi_\parallel^2},\, z \,\big|\, M\big),
\qquad
\xi_{\rm hm} = b(M, z)\,\xi_{\rm NL},$$

$$\Delta\Sigma_{2h}(R) = \bar\Sigma_{2h}(<R) - \Sigma_{2h}(R),
\qquad \bar\rho_m = \Omega_m\,\rho_{\rm crit},$$

and the **reference combination with the one-halo term is the
$\Sigma_{\max}/\kappa_{\max}$ ("max") model** — the halo–matter
correlation is not the sum but the pointwise maximum,
$\xi_{\rm hm} = \max\!\big(\xi_{\rm 1h},\, b\,\xi_{\rm NL}\big)$
(Hayashi & White 2008), the prescription of the DES Y1 lensing analysis
([McClintock et al. 2019, MNRAS 482, 1352](https://ui.adsabs.harvard.edu/abs/2019MNRAS.482.1352M/abstract),
arXiv:[1805.00039](https://arxiv.org/abs/1805.00039)) — applied at the
$\Sigma$/$\kappa$ level in the profile assembly. The stacked variant
uses the population-averaged bias
$\langle b\rangle_i = N_i[b]/N_i[1]$ from the `BiasWeightedSel`
diagnostic module ({doc}`../variants`).

**The recipe** (`ct_2hTerm` in `y3_buzzard/haloModel.py`, run per slice
of the 50-point redshift grid):

1. $\xi_{\rm mm}$ from the power spectrum on a fixed grid
   $R \in [10^{-3}, 10^{3}]\ {\rm cMpc}/h$ with **50 log nodes** — the
   correlation must be tabulated well past the BAO scale; 50 nodes is
   the speed/accuracy sweet spot (0.1%).
2. $\xi_{\rm 2halo} = b\,\xi_{\rm mm}$ with $b = 1$: **the tables are
   published unbiased** and the consumer applies $b(M,z)$ or
   $\langle b\rangle_i$.
3. $\Sigma_{2h}$ on the `Rp` grid via
   `cluster_toolkit.deltasigma.Sigma_at_R`, which extends the
   $\xi$ table below its inner edge with a **dummy NFW halo**
   ($M_d = 10^{14}\,M_\odot/h$, $c_d = 5$ — placeholder values, not
   physics).
4. $\Sigma \to \Delta\Sigma$ needs the full interior mass
   $\bar\Sigma(<R)$, which the pure two-halo table cannot supply at
   small $R$. The **dummy-halo trick**: *add* the dummy's analytic
   $\Sigma_{\rm NFW}$, Hankel/Abel-transform the sum to $\Delta\Sigma$
   (`DeltaSigma_at_R`), then *subtract* the dummy's analytic
   $\Delta\Sigma_{\rm NFW}$ — regularising the interior integral while
   cancelling the dummy exactly (up to numerics).
5. Residual negative values — where the cancellation of step 4 is
   imperfect at small $R$ — are set to **NaN** rather than silently
   folded in. This is the documented NaN region of `dSigma_hh` below
   $R \approx 8.6\,h^{-1}$cMpc (warning below).

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

