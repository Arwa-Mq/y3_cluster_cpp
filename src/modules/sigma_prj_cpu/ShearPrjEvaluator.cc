// Fixed-GL evaluator for the projected tangential shear
// gamma_t^prj = DSigma_prj * Sigma_crit^-1 (linear; reduced-shear
// denominator retired 2026-05-11).  Thin wrapper around
// sp_detail::ShearPrjCore; emits shear_prj/{vals,rnd,cl}.
//
// Module label: "shear_prj".  See src/models/sigma_prj_t.hh.
#include "models/sigma_prj_t.hh"
#include "utils/module_macros.hh"

using ShearPrjEvaluator = y3_cluster::ShearPrjEvaluator;

DEFINE_COSMOSIS_SCALAR_EVALUATOR_MODULE(ShearPrjEvaluator)
