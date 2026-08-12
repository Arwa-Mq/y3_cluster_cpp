# Traditional 1h+2h shear — the max model (`shear1h2h_max`)

**Status: reference implementation** (validated 2026-08-12). Provides
the *traditional-shear* arm of the two-pipeline comparison (traditional
vs projection), alongside the 1h-only `shear1h_fast_mass` in this
directory.

Radial operator — the Y1-era `SIG_MAX`/`KAPPA_MAX`/`GAMMA_MAX`
composition (`src/models/{sig,kappa,gamma}_max.hh`) on the modern
haloModel tables:

```text
Φ_max(R, lnM, z | bin) = max( ΔΣ_cl(R, lnM | bin),  b(lnM, z)·ΔΣ_hh(R, z) )
```

with ΔΣ_cl the production miscentred 1-halo mixture
(`include_miscentering = F` gives the pure centred term). Requires
`halo_model` to run with **`compute_lensing_2h = T`**.

**Why z stays inside the mass integral:** unlike the 1h-only observable
(whose profile is z-free — the radial-factorization study's "free
win"), the two-halo term depends on z through both b(M,z) and
ΔΣ_hh(R,z). The fast path therefore keeps the *z-resolved* tabulated
weight W2d(lnM, z) and does the double contraction; the S_ij
tabulation is still what makes it "fast_mass". A z-resolved
`full_ltmz` variant and an adaptive reference share the same structure
(`full_ltmz_core.full_ltmz_mass_z_weights`).

**Two-halo NaNs are expected** at low radii (the Hankel-transform
producer leaves ΔΣ_hh undefined there — in this dump 3850/6400 table
entries, everything below R ≈ 2.5 cMpc/h). They are zero-filled before
interpolation, which is the faithful treatment for a max model: the
1-halo term always wins where the 2-halo term is undefined, so
`max(1h, 0) = 1h`. Never let them propagate — a NaN silently poisons
every bin through the mass sum.

## Validation (2026-08-12, 2h-enabled pipeline at the fiducial point)

`validate_shear1h2h_max.py <dump>`, 12 bins × 10 radii, against the
**adaptive** z-resolved reference (reported error ≤ 1e-6, per the
namespace accuracy policy):

| Path | Error vs adaptive reference |
|---|---|
| `full_ltmz` GL (direct kernels) | **4.9e-5** |
| fast path (S_ij tabulated) | **8.3e-4** (the usual tabulation class) |
| 2h → 0 limit vs the validated 1h `fast_mass` backend | 4.0e-15 |

Cost: 83 ms/sample (12 bins × 10 radii). At the fiducial point the max
model exceeds the 1-halo-only stack by ≤ 0.8% on the production radial
grid (R ≤ 5 cMpc/h) — the two-halo term only wins at the largest radii,
as expected.

Output: `shear1h2h_max/vals` (hardcoded section, bin slow / R fast).
