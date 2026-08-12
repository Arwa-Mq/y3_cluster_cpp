// Unit tests for the radial_series C++ backend table + series core
// (src/pipelines/des_y3/observables/shear_1h2h/radial_series/cpp).
//
// Golden values were computed 2026-08-12 from the committed
// data/radial_series text tables with an independent Python bilinear /
// linear interpolator (see the python/ directory of the backend), so
// they pin both the table contents and the GSL interpolation semantics.
// The amplitude cross-check drives the composite A0(y) * u_mis through
// the production NFW_DSIGMA_MIS reader, so the fixed conventions of the
// two paths cannot drift apart silently.
//
// Requires Y3_CLUSTER_CPP_DIR to point at the source tree (data/).
#include "catch2/catch.hpp"

#include "pipelines/des_y3/observables/shear_1h2h/radial_series/cpp/shear1h_radial_series_t.hh"

#include <cmath>

using y3_cluster::des_y3::A0_of_y;
using y3_cluster::des_y3::RadialSeriesTable;
using y3_cluster::des_y3::y_of_lnM;

TEST_CASE("radial_series table interpolation matches golden values")
{
  RadialSeriesTable const tab;

  struct Golden {
    double lnx, lnxm, u0, u2, u3, c0;
  };
  // Bilinear (mis) / linear (cen) on the committed text tables.
  Golden const golden[] = {
    {0.1234, -0.987, 1.352591378821e-01, 1.588740771068e-01,
     5.981286639335e-02, 2.604112758336e-01},
    {1.5, 0.25, 4.439326085706e-02, 9.551480617341e-02,
     5.754488821577e-02, 7.476007404807e-02},
    {-2.0, -1.5, 2.710753422364e-02, 1.790178350799e-02,
     5.566988339746e-03, 4.771010214995e-01},
  };
  for (auto const& g : golden) {
    CHECK(tab.u_mis(0, g.lnx, g.lnxm) == Approx(g.u0).epsilon(1e-9));
    CHECK(tab.u_mis(2, g.lnx, g.lnxm) == Approx(g.u2).epsilon(1e-9));
    CHECK(tab.u_mis(3, g.lnx, g.lnxm) == Approx(g.u3).epsilon(1e-9));
    CHECK(tab.u_cen(0, g.lnx) == Approx(g.c0).epsilon(1e-9));
  }
}

TEST_CASE("radial_series series assembly matches golden values")
{
  RadialSeriesTable const tab;

  // Synthetic but representative moments; goldens from the Python
  // bilinear reference with identical inputs.
  double const R = 0.83625;
  double const r_mis = 0.17 * std::pow(52.5 / 100.0, 0.2);
  double const norm = 141.4;
  double const ybar = std::log(0.30);
  double const mu2 = 0.021;
  double const mu3 = -0.0035;
  double const f_mis = 0.22;
  double const rho_mult = 0.3096;

  CHECK(tab.series(R, r_mis, norm, ybar, mu2, mu3, f_mis, rho_mult, 2) ==
        Approx(1.346639748480e+04).epsilon(1e-9));
  CHECK(tab.series(R, r_mis, norm, ybar, mu2, mu3, f_mis, rho_mult, 3) ==
        Approx(1.342457714365e+04).epsilon(1e-9));
}

TEST_CASE("radial_series amplitude matches the production NFW_DSIGMA_MIS")
{
  RadialSeriesTable const tab;
  y3_cluster::NFW_DSIGMA_MIS dsigma_mis(4.0, 2.77533742639e+11,
                                        y3_cluster::GAMMA);
  double const omega_m = 0.3096;
  dsigma_mis.set_rho_mult(omega_m);

  // rho_mult * A0(y) * U0(ln x, ln x_mis) must reproduce the production
  // reader (which interpolates the original data/nfw_off_center table)
  // at the level the generator measured for U0 fidelity (<= ~5e-4 over
  // the physical window).
  struct Point {
    double R, r_mis, lnM;
  };
  Point const pts[] = {{0.5, 0.14, 32.0}, {1.5, 0.15, 33.5},
                       {3.0, 0.16, 34.5}};
  for (auto const& p : pts) {
    double const y = y_of_lnM(p.lnM);
    double const mine =
      omega_m * A0_of_y(y) *
      tab.u_mis(0, std::log(p.R) - y, std::log(p.r_mis) - y);
    double const prod = dsigma_mis(p.R, p.r_mis, p.lnM);
    CHECK(mine == Approx(prod).epsilon(5e-4));
  }
}

TEST_CASE("radial_series miscentred table recovers the centred limit")
{
  RadialSeriesTable const tab;
  // For x >> x_mis the gamma-averaged profile reduces to the centred
  // one; at the table's smallest x_mis (~0.01) the suppression at
  // x >= 1 is a few 1e-4.
  for (double lnx : {0.0, 1.0, 2.0}) {
    double const ratio = tab.u_mis(0, lnx, -4.625) / tab.u_cen(0, lnx);
    CHECK(ratio == Approx(1.0).epsilon(2e-3));
  }
}
