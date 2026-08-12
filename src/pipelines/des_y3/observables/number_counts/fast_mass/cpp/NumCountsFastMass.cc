// Cluster number counts via the fast_mass strategy — C++ backend, inline
// selection function.
//
// The exact-z-contraction fast path stated as a des_y3 module, same
// algorithm as the production NumCountsSel.so (NumCountsSelGL in the
// immutable n_operator_sel_gl_t.hh) EXCEPT the selection tensor S_ij is
// computed inline, per sample, inside SelFunctionInlineCore_t's own
// construction (the Gauss-Legendre quadrature over true richness
// lambda^tr, composed from MOR_HOD_t + RichnessKernel_t + PlobLtrEMG_t),
// instead of being read from the separate sel_function.py tabulation
// (sel_function/S_stack, ~197 ms/sample in the reference pipeline).
//
// Namespace label/section: module `NumCountsFastMass`, output
// `numcounts_fast_mass/vals` (same section as the Python fast_mass
// reference; the two are interchangeable and must never run in one
// pipeline — DataBlock sections do not overwrite). Production's own
// section (`numcountssel`) is different, so this driver CAN co-run with
// production `NumCountsSel.so` in the same pipeline for direct
// comparison.
//
// Options: bin_index wall (default = every configured bin), required
// zt_low/zt_high/lnm_low/lnm_high, n_lnm (96), n_z (64); bin definitions
// lam_min/lam_max/zob_min/zob_max/sigma_z (12-entry walls, same values as
// the `[sel_function]` ini section); optional n_q_lam (32), l_lam (6.0).
//
// Status: reference re-expression, inline-selection variant. Production
// remains NumCountsSel.so (tabulated). No existing template, model, or
// module is modified; registry entry added.
#include "models/sel_function_inline_t.hh"
#include "utils/make_grid_points.hh"
#include "utils/module_macros.hh"

#include <array>
#include <stdexcept>

class NumCountsFastMassCpp {
public:
  using grid_t = y3_cluster::grid_t<1>;
  using grid_point_t = grid_t::value_type;
  static constexpr std::size_t n_outputs = 1;

  explicit NumCountsFastMassCpp(cosmosis::DataBlock& cfg)
    : core_(cfg, module_label())
  {}

  void
  set_sample(cosmosis::DataBlock& s)
  {
    core_.build_weights(s, /*include_sci=*/false);
  }

  std::array<double, n_outputs>
  evaluate(grid_point_t const& pt) const
  {
    int const b = static_cast<int>(pt[0]);
    if (b < 0 || static_cast<std::size_t>(b) >= core_.n_bins())
      throw std::out_of_range("NumCountsFastMassCpp: bin_index out of range");
    return {core_.norm(b)};
  }

  static char const* module_label() { return "NumCountsFastMass"; }

  static std::array<char const*, n_outputs>
  output_sections()
  {
    return {"numcounts_fast_mass"};
  }

  static grid_t
  make_grid_points(cosmosis::DataBlock& cfg)
  {
    return y3_cluster::make_grid_points_wall_of_numbers(
      cfg, module_label(), "bin_index");
  }

private:
  y3_cluster::SelFunctionInlineCore_t core_;
};

DEFINE_COSMOSIS_SCALAR_EVALUATOR_MODULE(NumCountsFastMassCpp)
