// Fixed-GL evaluator for the projected surface mass density Sigma_prj.
// Thin wrapper around sp_detail::ShearPrjCore; emits sigma_prj/{vals,rnd,cl}.
//
// Module label: "sigma_prj".  See src/models/sigma_prj_t.hh for the
// shared cache + integrand recipe (Costanzi-2026, May-2026 RichnessSelection).
#include "models/sigma_prj_t.hh"
#include "utils/module_macros.hh"

using SigmaPrjEvaluator = y3_cluster::SigmaPrjEvaluator;

DEFINE_COSMOSIS_SCALAR_EVALUATOR_MODULE(SigmaPrjEvaluator)
