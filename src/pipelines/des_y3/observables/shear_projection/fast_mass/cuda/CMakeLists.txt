include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

# fast_mass projection CUDA/PAGANI backend (des_y3 Phase 2), standard
# CUDA integration-module template.
if (USE_CUDA)
  find_library(NVTX_LIBRARY nvToolsExt PATHS ENV LD_LIBRARY_PATH)
  add_library(DSigmaPrjFastMassGpu MODULE DSigmaPrjFastMassGpu.cu)
  target_link_libraries(DSigmaPrjFastMassGpu PRIVATE cosmosis
                        ${NVTX_LIBRARY} utils cuda_models cuda_utils
                        fmt::fmt)
  target_include_directories(DSigmaPrjFastMassGpu PRIVATE
    ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${PAGANI_DIR})
  set_target_properties(DSigmaPrjFastMassGpu PROPERTIES
    POSITION_INDEPENDENT_CODE ON CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH})
  target_compile_options(DSigmaPrjFastMassGpu PRIVATE
    "--expt-relaxed-constexpr" "-lineinfo" "-DCURAND")
endif ()
