// Miscentred one-halo shear via the fast_mass strategy — C++ backend,
// inline selection function.
//
// Identical to Shear1hFastMass.cc (SelGLCore + production miscentred
// mixture) EXCEPT the selection tensor S_ij is computed inline, per
// sample, inside SelFunctionInlineCore_t's own construction (the
// Gauss-Legendre quadrature over true richness lambda^tr, composed from
// MOR_HOD_t + RichnessKernel_t + PlobLtrEMG_t), instead of being read
// from the separate sel_function.py tabulation (sel_function/S_stack,
// ~197 ms/sample in the reference pipeline). Shear1hFastMass.cc itself
// is left untouched (it stays validated bitwise-identical to production
// via the shared tabulated table); this is an additive sibling so both
// remain independently comparable.
//
// Namespace label/section: module `Shear1hFastMassInline`, output
// `shear1h_fast_mass_inline/vals` — distinct from both production
// (`shear1hmissel`) and the existing tabulated fast_mass driver
// (`shear1h_fast_mass`), so all three can co-run in one pipeline for a
// direct 3-way comparison.
//
// Options: bin_index x r_perp cartesian grid (bin slow / R fast),
// lob_centers (default 25 37.5 52.5 130), required zt_low/zt_high/
// lnm_low/lnm_high, n_lnm (96), n_z (64); bin definitions
// lam_min/lam_max/zob_min/zob_max/sigma_z (12-entry walls, same values as
// the `[sel_function]` ini section); optional n_q_lam (32), l_lam (6.0).
// Status: reference re-expression, inline-selection variant. Production
// remains Shear1hMisSel.so. No existing template, model, or module is
// modified; registry entry added.
#include "models/sel_function_inline_t.hh"
#include "utils/module_macros.hh"

#include <array>
#include <optional>
#include <stdexcept>
#include <vector>

class Shear1hFastMassInlineCpp {
public:
  using grid_t = y3_cluster::grid_t<2>;
  using grid_point_t = grid_t::value_type;
  static constexpr std::size_t n_outputs = 1;

  explicit Shear1hFastMassInlineCpp(cosmosis::DataBlock& cfg)
    : core_(cfg, module_label())
    , dsigma_mis_(4.0, 2.77533742639e+11, y3_cluster::GAMMA)
  {
    lob_centers_ =
      y3_cluster_sel_weights::mis_detail::read_lob_centers(cfg,
                                                           module_label());
    if (lob_centers_.empty())
      throw std::runtime_error("Shear1hFastMassInlineCpp: lob_centers is empty");
  }

  void
  set_sample(cosmosis::DataBlock& s)
  {
    namespace w = y3_cluster_sel_weights;
    core_.build_weights(s, /*include_sci=*/true);
    dsigma_nfw_.emplace(y3_cluster::make_Interp2D(
      s, "haloModel", "r_sigma", "lnM", "dSigma_nfw"));
    f_mis_ = w::mis_detail::read_mis_param(s, "f_mis",
                                           w::mis_detail::F_MIS_DEFAULT);
    double const tau_mis = w::mis_detail::read_mis_param(
      s, "tau_mis", w::mis_detail::TAU_MIS_DEFAULT);
    dsigma_mis_.set_rho_mult(
      s.view<double>("cosmological_parameters", "omega_M"));
    r_mis_.assign(core_.n_bins(), 0.0);
    for (std::size_t b = 0; b != core_.n_bins(); ++b)
      r_mis_[b] = tau_mis * w::mis_detail::R_lambda(
                              lob_centers_[b % lob_centers_.size()]);
  }

  std::array<double, n_outputs>
  evaluate(grid_point_t const& pt) const
  {
    int const b = static_cast<int>(pt[0]);
    double const R = pt[1];
    if (b < 0 || static_cast<std::size_t>(b) >= core_.n_bins())
      throw std::out_of_range("Shear1hFastMassInlineCpp: bin_index");
    auto const& wb = core_.weights(b);
    auto const& xs = core_.lnm_x();
    auto const& ws = core_.lnm_w();
    double acc = 0.0;
    for (std::size_t k = 0; k != xs.size(); ++k) {
      double const d_cen = dsigma_nfw_->clamp(R, xs[k]);
      double const d_mis = dsigma_mis_(R, r_mis_[b], xs[k]);
      acc += ws[k] * wb[k] *
             ((1.0 - f_mis_) * d_cen + f_mis_ * d_mis);
    }
    return {acc};
  }

  static char const* module_label() { return "Shear1hFastMassInline"; }

  static std::array<char const*, n_outputs>
  output_sections()
  {
    return {"shear1h_fast_mass_inline"};
  }

  static grid_t
  make_grid_points(cosmosis::DataBlock& cfg)
  {
    return y3_cluster::make_grid_points_cartesian_product(
      cfg, module_label(), "bin_index", "r_perp");
  }

private:
  y3_cluster::SelFunctionInlineCore_t core_;
  std::vector<double> lob_centers_;
  std::optional<y3_cluster::Interp2D> dsigma_nfw_;
  y3_cluster::NFW_DSIGMA_MIS dsigma_mis_;
  double f_mis_{y3_cluster_sel_weights::mis_detail::F_MIS_DEFAULT};
  std::vector<double> r_mis_;
};

DEFINE_COSMOSIS_SCALAR_EVALUATOR_MODULE(Shear1hFastMassInlineCpp)
