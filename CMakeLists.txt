cmake_minimum_required(VERSION 2.8)
project("BeamAdapter")

set(BEAMADAPTER_MAJOR_VERSION 1)
set(BEAMADAPTER_MINOR_VERSION 0)
set(BEAMADAPTER_VERSION ${BEAMADAPTER_MAJOR_VERSION}.${BEAMADAPTER_MINOR_VERSION})

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

set(HEADER_FILES
  constraint/AdaptiveBeamConstraint.h
  constraint/AdaptiveBeamConstraint.inl
  constraint/AdaptiveBeamLengthConstraint.h
  constraint/AdaptiveBeamLengthConstraint.inl
  constraint/ImplicitSurfaceAdaptiveConstraint.h
  constraint/ImplicitSurfaceAdaptiveConstraint.inl

  controller/AdaptiveBeamController.h
  controller/AdaptiveBeamController.inl
  controller/InterventionalRadiologyController.h
  controller/InterventionalRadiologyController.inl
  controller/SutureController.h
  controller/SutureController.inl

  engine/WireRestShape.h
  engine/WireRestShape.inl
  engine/SteerableCatheter.h
  engine/SteerableCatheter.inl

  forcefield/AdaptiveBeamForceFieldAndMass.h
  forcefield/AdaptiveBeamForceFieldAndMass.inl

  mapping/AdaptiveBeamMapping.h
  mapping/AdaptiveBeamMapping.inl
  mapping/MultiAdaptiveBeamMapping.h
  mapping/MultiAdaptiveBeamMapping.inl

  solver/UnbuiltGenericConstraintSolver.h

  #BaseRestShape.h
  #BaseRestShape.inl
  initBeamAdapter.h

  BeamInterpolation.h
  BeamInterpolation.inl
  WireBeamInterpolation.h
  WireBeamInterpolation.inl
  )

set(SOURCE_FILES
  constraint/AdaptiveBeamConstraint.cpp
  constraint/AdaptiveBeamLengthConstraint.cpp
  constraint/ImplicitSurfaceAdaptiveConstraint.cpp

  controller/AdaptiveBeamController.cpp
  controller/InterventionalRadiologyController.cpp
  controller/SutureController.cpp

  engine/WireRestShape.cpp
  engine/SteerableCatheter.cpp

  forcefield/AdaptiveBeamForceFieldAndMass.cpp

  mapping/AdaptiveBeamMapping.cpp
  mapping/MultiAdaptiveBeamMapping.cpp

  solver/UnbuiltGenericConstraintSolver.cpp

  #BaseRestShape.cpp
  initBeamAdapter.cpp
  BeamInterpolation.cpp
  WireBeamInterpolation.cpp
  )

set(README_FILE "BeamAdapter.txt" )

if(PLUGIN_SOFAADVANCEDCONSTRAINT)
set(HEADER_FILES ${HEADER_FILES}
  AdaptiveBeamContactMapper.h
  AdaptiveBeamContactMapper.inl
  AdaptiveBeamFrictionContact.h
  AdaptiveBeamFrictionContact.inl
  MultiAdaptiveBeamContactMapper.h
  MultiAdaptiveBeamContactMapper.inl
)

set(SOURCE_FILES ${SOURCE_FILES}
  AdaptiveBeamContactMapper.cpp
  AdaptiveBeamFrictionContact.cpp
  MultiAdaptiveBeamContactMapper.cpp
)
        #find_package(CUDA REQUIRED)
        #cuda_include_directories(${SOFA_CUDA_DIR})
        #cuda_include_directories(${SOFA_FRAMEWORK_DIR})
        #cuda_add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${CUDA_SOURCES})
endif()

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${README_FILES})


set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-DSOFA_BUILD_BEAMADAPTER")
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${HEADER_FILES}")

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include>")

find_package(SofaFramework REQUIRED)

target_link_libraries(${PROJECT_NAME} SofaCore SofaConstraint SofaBaseMechanics SofaUserInteraction)

## Install rules for the library and headers; CMake package configurations files
sofa_create_package(BeamAdapter ${BEAMADAPTER_VERSION} ${PROJECT_NAME} BeamAdapter)

# if only this was possible...
if(SOFA_BUILD_TESTS)
    add_subdirectory(BeamAdapter_test)
endif()

# this paragraph has not been ported to the new build system yet...
#if(SOFA-PLUGIN_SOFAADVANCEDCONSTRAINT)
    #set(LINKER_DEPENDENCIES SofaAdvancedConstraint)
    #list(APPEND LINKER_DEPENDENCIES SofaCUDA)
    #list(APPEND LINKER_DEPENDENCIES SofaCUDADev)
    #set(LINKER_DEPENDENCIES SofaAdvancedConstraint)
#endif()

#
# Organize projects into folders
#

SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
SET_PROPERTY(TARGET ${PROJECT_NAME}	PROPERTY FOLDER "plugins")
