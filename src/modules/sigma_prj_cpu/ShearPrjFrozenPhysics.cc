// Option C (frozen-physics benchmark backend): same frozen rnd-hoist /
// a_b(z)-drift reduction as Option E (ShearPrjFrozenCuhre), but the
// (theta, lnM) assembly is an explicit fixed N_theta x N_M grid + dot
// product (sp_detail::build_theta_grid, the same log-GL per-R-breakpoint
// grid production's own ShearPrjEvaluator uses), not a continuous Cuhre
// integral. See RichnessSelection/docs/richness_selection_frozen.tex
// section "Extension: frozen physics for <DeltaSigma_prj>" and
// richness_selection.FrozenDeltaSigmaPrj for the Python reference this
// ports, and the plan doc "Port frozen-physics b_sel/DeltaSigma_prj recipe
// into the cosmosis pipeline" (Part II/Part III, Option C) for the design
// rationale. Diagnostic/benchmark backend -- not wired into the production
// [pipeline] modules list until benchmarked.
#include "models/sigma_prj_frozen_t.hh"
#include "utils/module_macros.hh"

using ShearPrjFrozenPhysics = y3_cluster::ShearPrjFrozenPhysics;

DEFINE_COSMOSIS_SCALAR_EVALUATOR_MODULE(ShearPrjFrozenPhysics)
