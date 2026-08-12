"""Production-exact replicas of the miscentred 1-halo lensing profile.

The Shear1hMisSel mixture (lensing_weights.hh):

    DSigma_cl(R | M; bin) = (1 - f_mis) DSigma_nfw(R, M)
                          + f_mis      DSigma_mis(R, tau_mis R_lambda, M)

with the centred term interpolated from the per-sample
``haloModel/dSigma_nfw`` table (GSL bilinear, clamped) and the
miscentred term from the fixed gamma-kernel look-up table through
NFW_DSIGMA_MIS (bilinear in ln u over (ln x, ln x_mis), clamped, with
the analytic 2 r_s delta_c rho_crit rho_mult amplitude and the 1e-12
Mpc^2 -> pc^2 conversion).

These are interpolation-exact replicas — linear interpolation with
clamped queries is the same arithmetic GSL performs — so Python
implementations that consume them can match the production .so to
near machine precision. (Contrast with
observables/shear_1h2h/radial_series/python/nfw_profile_family.py,
which deliberately uses a *smooth* spline view of the same gamma table
for offline differentiation.)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

CONC = 4.0
RHOC = 2.77533742639e11
DELTA_C = (200.0 * CONC**3 / 3.0) / (np.log(1.0 + CONC) - CONC / (1.0 + CONC))

GAMMA_TABLE = "table_1000_1e-03_5e+03_{}_{}.txt"


def _repo_root():
    for p in Path(__file__).resolve().parents:
        if (p / "data" / "nfw_off_center").is_dir():
            return p
    raise FileNotFoundError("lensing_profiles: cannot locate data/")


def r_lambda(lob):
    """R_lambda(lambda^ob) = (lambda^ob / 100)^0.2 [Mpc/h]."""
    return (np.asarray(lob, dtype=float) / 100.0) ** 0.2


class HaloModelDSigmaNfw:
    """Bilinear clamp interpolator over haloModel/dSigma_nfw(r_sigma, lnM)."""

    def __init__(self, source):
        self._r = source.array("halomodel", "r_sigma")
        self._lnm = source.array("halomodel", "lnm")
        vals = np.asarray(source.array("halomodel", "dsigma_nfw"), dtype=float)
        if vals.shape != (self._lnm.size, self._r.size):
            vals = vals.reshape(self._lnm.size, self._r.size)
        self._interp = RegularGridInterpolator(
            (self._lnm, self._r), vals, method="linear",
            bounds_error=False, fill_value=None)

    def __call__(self, r_perp, lnM):
        r = np.clip(np.asarray(r_perp, dtype=float), self._r[0], self._r[-1])
        m = np.clip(np.asarray(lnM, dtype=float), self._lnm[0], self._lnm[-1])
        m_b, r_b = np.broadcast_arrays(m, r)
        return self._interp(np.stack([m_b, r_b], axis=-1))


class NfwDsigmaMisProduction:
    """NFW_DSIGMA_MIS replica: gamma-kernel table, bilinear ln u, clamped."""

    def __init__(self, kernel="gamma", data_dir=None):
        d = Path(data_dir) if data_dir else _repo_root() / "data" / "nfw_off_center"
        self._lnx = np.loadtxt(d / GAMMA_TABLE.format(kernel, "logx"))
        self._lnxm = np.loadtxt(d / GAMMA_TABLE.format(kernel, "logxmis"))
        logds = np.loadtxt(d / ("table_1000_1e-03_5e+03_log_deltasigma_"
                                + kernel + ".txt"))
        self._interp = RegularGridInterpolator(
            (self._lnxm, self._lnx), logds, method="linear",
            bounds_error=False, fill_value=None)

    def __call__(self, r_perp, r_mis, lnM, rho_mult=1.0):
        lnM = np.asarray(lnM, dtype=float)
        r_200 = np.cbrt(3.0 * np.exp(lnM) / (800.0 * np.pi * RHOC))
        r_s = r_200 / CONC
        lnx = np.clip(np.log(np.asarray(r_perp, dtype=float) / r_s),
                      self._lnx[0], self._lnx[-1])
        lnxm = np.clip(np.log(np.asarray(r_mis, dtype=float) / r_s),
                       self._lnxm[0], self._lnxm[-1])
        lnxm_b, lnx_b = np.broadcast_arrays(lnxm, lnx)
        log_u = self._interp(np.stack([lnxm_b, lnx_b], axis=-1))
        norm = 2.0 * r_s * DELTA_C * RHOC * rho_mult
        return norm * np.exp(log_u) * 1.0e-12


class MisMixtureProfile:
    """The per-bin production mixture Phi_i(R, lnM) with per-sample knobs."""

    def __init__(self, source, *, lob_centers, f_mis, tau_mis, omega_m):
        self._cen = HaloModelDSigmaNfw(source)
        self._mis = NfwDsigmaMisProduction()
        self._lob = np.asarray(lob_centers, dtype=float)
        self.f_mis = float(f_mis)
        self.tau_mis = float(tau_mis)
        self.omega_m = float(omega_m)

    def r_mis(self, bin_index):
        """tau_mis R_lambda for richness bin = bin_index mod len(lob)."""
        return self.tau_mis * float(r_lambda(self._lob[bin_index % self._lob.size]))

    def __call__(self, bin_index, r_perp, lnM):
        d_cen = self._cen(r_perp, lnM)
        d_mis = self._mis(r_perp, self.r_mis(bin_index), lnM,
                          rho_mult=self.omega_m)
        return (1.0 - self.f_mis) * d_cen + self.f_mis * d_mis
