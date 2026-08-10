// Option E (frozen-physics benchmark backend): keeps the same frozen
// rnd-hoist / a_b(z)-drift reduction as Option C but replaces the explicit
// N_theta x N_M fixed-grid dot product with a continuous 2-D Cuhre integral
// over (theta, lnM), using Interp1D-tabulated w_rnd(lnM)/w_cl(lnM)/Psi(theta).
// See RichnessSelection/docs/richness_selection_frozen.tex section
// "Extension: frozen physics for <DeltaSigma_prj>" and
// richness_selection.FrozenDeltaSigmaPrj for the Python reference this
// ports, and the plan doc "Port frozen-physics b_sel/DeltaSigma_prj recipe
// into the cosmosis pipeline" (Part III, Option E) for the design
// rationale. Diagnostic/benchmark backend -- not wired into the production
// [pipeline] modules list until benchmarked.
#include "models/sigma_prj_frozen_interp_t.hh"
#include "utils/module_macros.hh"

using ShearPrjFrozenCuhre = y3_cluster::ShearPrjFrozenCuhre;

DEFINE_COSMOSIS_SCALAR_EVALUATOR_MODULE(ShearPrjFrozenCuhre)
