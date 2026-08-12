include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

# fast_mass shear C++ backend (des_y3 Phase 2): SelGLCore + production
# mixture, namespace label/section.
add_library(Shear1hFastMass MODULE Shear1hFastMass.cc)
target_link_libraries(Shear1hFastMass PRIVATE cosmosis utils models
                      ${GSL_LIBRARIES} ${CUBA_LIBRARIES} fmt::fmt)
target_include_directories(Shear1hFastMass PRIVATE
  ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})
