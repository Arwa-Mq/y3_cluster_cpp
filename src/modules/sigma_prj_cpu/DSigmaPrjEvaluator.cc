// Fixed-GL evaluator for the projected differential surface density
// DSigma_prj.  Thin wrapper around sp_detail::ShearPrjCore; emits
// dsigma_prj/{vals,rnd,cl}.
//
// Module label: "dsigma_prj".  See src/models/sigma_prj_t.hh.
#include "models/sigma_prj_t.hh"
#include "utils/module_macros.hh"

using DSigmaPrjEvaluator = y3_cluster::DSigmaPrjEvaluator;

DEFINE_COSMOSIS_SCALAR_EVALUATOR_MODULE(DSigmaPrjEvaluator)
