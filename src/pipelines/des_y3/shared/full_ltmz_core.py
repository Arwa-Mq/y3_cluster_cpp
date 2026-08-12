"""The full_ltmz selection contraction, shared across observables.

Computes the fully explicit (lambda_true, z) contraction of the DES Y3
selection at fixed-GL mass nodes — the quantity every full_ltmz
observable integrates against its own radial operator:

    W_ij(lnM) = int dz int dlt  n(M,z) dV/dOmega/dz(z) Omega(z)
                [Sigma_crit_inv(z)]  K_j(z) K_i(lt, z) P_HOD(lt | M, z)

with the lt integral on the per-(M,z) HOD bracket exactly as the
maintained sel_function machinery defines it (its kernels are imported,
not copied). Number counts contract W against 1; the 1-halo shear
against its z-free radial profile Phi_i(R, lnM).

This is the *reference* counterpart of the fast path's tabulated
version: production tabulates S_ij once on a fixed grid and interpolates
(sel_function.py -> SelGLCore); here every kernel is evaluated at the
quadrature nodes directly.
"""
from __future__ import annotations

import numpy as np

from . import datablock_models as dm
from . import sel_kernels


def full_ltmz_mass_weights(bins, mor, plob_splines, hmf, dv, *,
                           sci=None, zt_low, zt_high, lnm_low, lnm_high,
                           n_lnm=96, n_z=64, n_q=32, l_lam=6.0):
    """W_ij on GL mass nodes for every configured bin.

    ``bins``: dict of equal-length arrays (lam_min, lam_max, zob_min,
    zob_max, sigma_z). ``sci``: optional Sigma_crit_inv(z) callable
    (shear observables fold it into the z contraction; counts pass
    None). Returns (lnm_x, lnm_w, W) with W of shape (n_bins, n_lnm).
    """
    sf = sel_kernels.load()

    z_x, z_w = dm.gl_nodes(zt_low, zt_high, n_z)
    lnm_x, lnm_w = dm.gl_nodes(lnm_low, lnm_high, n_lnm)
    gl_t, gl_w = np.polynomial.legendre.leggauss(n_q)

    lam_k, w_k, p_mz, degenerate = sf._compute_lam_nodes_and_P_HOD(
        lnm_x, z_x, mor, gl_t, gl_w, L=l_lam)

    mu_p, sig_p, tau_p, fprj_p = sf._plob_params(lam_k, z_x, plob_splines)
    lam_min = np.asarray(bins["lam_min"], dtype=float)
    lam_max = np.asarray(bins["lam_max"], dtype=float)
    edges = np.unique(np.concatenate([lam_min, lam_max]))
    cdfs = sf._cdf_lob_stacked(edges, mu_p, sig_p, tau_p, fprj_p)

    zfac = z_w * dv(z_x) * dm.omega_z_des(z_x)
    if sci is not None:
        zfac = zfac * sci(z_x)
    base_kq = hmf(lnm_x[:, None], z_x[None, :]) * zfac[None, :]

    n_bins = lam_min.size
    weights = np.empty((n_bins, lnm_x.size))
    for b in range(n_bins):
        lo = int(np.searchsorted(edges, lam_min[b]))
        hi = int(np.searchsorted(edges, lam_max[b]))
        k_i = cdfs[hi] - cdfs[lo]
        s_kq = np.sum(w_k * k_i * p_mz, axis=-1)
        s_kq = np.where(degenerate, 0.0, s_kq)
        k_j = sf._K_j(z_x, float(bins["zob_min"][b]),
                      float(bins["zob_max"][b]), float(bins["sigma_z"][b]))
        weights[b] = (base_kq * s_kq) @ k_j
    return lnm_x, lnm_w, weights
