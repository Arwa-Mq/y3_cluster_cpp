
// projection-shear module wrapper
// It directly evaluates the projection shear following sigma_prj_t.hh / ShearPrjEvaluator.cc.

#include "models/sigma_prj_gpu_t.cuh"
#include "utils/cuda_module_macros.cuh"

using ShearPrjEvaluator = y3_cuda::SigmaPrjGPU<y3_cuda::DSigmaRndWeight>;

DEFINE_COSMOSIS_CUDA_INTEGRATION_MODULE(ShearPrjEvaluator);