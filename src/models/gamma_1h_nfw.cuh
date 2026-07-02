#ifndef Y3_CLUSTER_GAMMA_1H_NFW_CUH
#define Y3_CLUSTER_GAMMA_1H_NFW_CUH
//
// This file defines the GPU-compatible 1-halo shear model, including miscentering.
// IN progress: miscentering 
//
#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "common/cuda/Interp2D.cuh"
#include "models/ez.hh"
#include "utils/make_interp_2d.cuh"
#include "utils/primitives.hh"

namespace y3_cuda {

  namespace gamma_1h_detail {
    inline quad::Interp2D make_dsigma_nfw(cosmosis::DataBlock& sample)
    {
      // quad::Interp2D(xs, ys, zs) stores xs→cols, ys→rows.
      // dSigma_nfw shape (n_M, n_r): M rows, r cols. Call as clamp(r, lnM).
      return make_Interp2D(sample, "haloModel", "r_sigma", "lnM", "dSigma_nfw");
    }

    inline quad::Interp2D make_dsigma_mis_or_cen(cosmosis::DataBlock& sample)
    {
      // JOHNNEY NEED TO CHECK
      // If the pipeline publishes a precomputed miscentered DeltaSigma table,use it. 
      //Otherwise fall back to dSigma_nfw so the code still compiles
      // and gives the old centred result when no miscentering table exists.
      
      if (sample.has_val("haloModel", "dSigma_mis")) {
        return make_Interp2D(sample, "haloModel", "r_sigma", "lnM", "dSigma_mis");
      }
      return make_dsigma_nfw(sample);
    }

    inline double read_or_default(cosmosis::DataBlock& sample,
                                  char const* section,
                                  char const* name,
                                  double fallback)
    {
      return sample.has_val(section, name) ? sample.view<double>(section, name) : fallback;
    }
  }

  class GAMMA_1H_NFW {
  private:
    quad::Interp2D _dsigma_nfw;
    quad::Interp2D _dsigma_mis;
    quad::Interp2D _sigma_crit_inv;
    double f_mis_ = 0.22;
    bool has_mis_table_ = false;

  public:
    size_t
    get_device_mem_footprint()
    {
      size_t size = 0;
      size += _dsigma_nfw.get_device_mem_footprint();
      size += _dsigma_mis.get_device_mem_footprint();
      size += _sigma_crit_inv.get_device_mem_footprint();
      return size;
    }

    GAMMA_1H_NFW(quad::Interp2D const& dsigma_nfw,
                 quad::Interp2D const& dsigma_mis,
                 quad::Interp2D const& sigma_crit_inv,
                 double f_mis,
                 bool has_mis_table)
      : _dsigma_nfw(dsigma_nfw),
        _dsigma_mis(dsigma_mis),
        _sigma_crit_inv(sigma_crit_inv),
        f_mis_(f_mis),
        has_mis_table_(has_mis_table)
    {}

    explicit GAMMA_1H_NFW(cosmosis::DataBlock& sample)
      : _dsigma_nfw(gamma_1h_detail::make_dsigma_nfw(sample))
      , _dsigma_mis(gamma_1h_detail::make_dsigma_mis_or_cen(sample))
      // sigma_crit_inv shape (n_z, n_r): z rows, r cols. Call as clamp(r, z).
      , _sigma_crit_inv(make_Interp2D(sample,
                                      "sigmaCritInv",
                                      "r_sigma",
                                      "sigmaCritInv",
                                      "z",
                                      "sigmaCritInv",
                                      "sigma_crit_inv"))
      , f_mis_(gamma_1h_detail::read_or_default(sample, "miscentering", "f_mis", 0.22))
      , has_mis_table_(sample.has_val("haloModel", "dSigma_mis"))
    {}

    __device__ __host__ double
    operator()(double r, double lnM, double zt) const
    /* r in h^-1 Mpc */ /* M in h^-1 M_solar, represented by lnM */
    {
      double const dsigma_cen = _dsigma_nfw.clamp(r, lnM);
      double const dsigma_mis = has_mis_table_ ? _dsigma_mis.clamp(r, lnM) : dsigma_cen;

      // Full 1-halo DeltaSigma mixture:
      // DeltaSigma_cl = (1 - f_mis) DeltaSigma_NFW + f_mis DeltaSigma_mis.
      double const dsigma_1h = (1.0 - f_mis_) * dsigma_cen + f_mis_ * dsigma_mis;

      double const sigc_inv = _sigma_crit_inv.clamp(r, zt);
      return dsigma_1h * sigc_inv;
    }
  };
}

#endif
