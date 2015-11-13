cmake_minimum_required(VERSION 2.8)
project("BeamAdapter")

set(BEAMADAPTER_MAJOR_VERSION 1)
set(BEAMADAPTER_MINOR_VERSION 0)
set(BEAMADAPTER_VERSION ${BEAMADAPTER_MAJOR_VERSION}.${BEAMADAPTER_MINOR_VERSION})

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

set(HEADER_FILES
  AdaptiveBeamConstraint.h
  AdaptiveBeamConstraint.inl
  AdaptiveBeamController.h
  AdaptiveBeamController.inl
  AdaptiveBeamForceFieldAndMass.h
  AdaptiveBeamForceFieldAndMass.inl
  AdaptiveBeamLengthConstraint.h
  AdaptiveBeamLengthConstraint.inl
  AdaptiveBeamMapping.h
  AdaptiveBeamMapping.inl
  BaseRestShape.h
  BaseRestShape.inl
  BeamInterpolation.h
  BeamInterpolation.inl
  ImplicitSurfaceAdaptiveConstraint.h
  ImplicitSurfaceAdaptiveConstraint.inl
  initBeamAdapter.h
  InterventionalRadiologyController.h
  InterventionalRadiologyController.inl
  MultiAdaptiveBeamMapping.h
  MultiAdaptiveBeamMapping.inl
  SutureController.h
  SutureController.inl
  UnbuiltGenericConstraintSolver.h
  WireBeamInterpolation.h
  WireBeamInterpolation.inl
  WireRestShape.h
  WireRestShape.inl
  SteerableCatheter.h
  SteerableCatheter.inl
  )

set(SOURCE_FILES
  AdaptiveBeamConstraint.cpp
  AdaptiveBeamController.cpp
  AdaptiveBeamForceFieldAndMass.cpp
  AdaptiveBeamLengthConstraint.cpp
  AdaptiveBeamMapping.cpp
  BaseRestShape.cpp
  BeamInterpolation.cpp
  ImplicitSurfaceAdaptiveConstraint.cpp
  initBeamAdapter.cpp
  InterventionalRadiologyController.cpp
  MultiAdaptiveBeamMapping.cpp
  SutureController.cpp
  UnbuiltGenericConstraintSolver.cpp
  WireBeamInterpolation.cpp
  WireRestShape.cpp
  SteerableCatheter.cpp
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


if(${SOFA_FLOATING_POINT_TYPE} STREQUAL float)
    add_definitions(-DSOFA_WITH_FLOAT)
endif()
if(${SOFA_FLOATING_POINT_TYPE} STREQUAL double)
    add_definitions(-DSOFA_WITH_DOUBLE)
endif()
if(${SOFA_FLOATING_POINT_TYPE} STREQUAL both)
    add_definitions(-DSOFA_WITH_FLOAT)
    add_definitions(-DSOFA_WITH_DOUBLE)
endif()

## Install rules for the library and headers; CMake package configurations files
sofa_create_package(BeamAdapter ${BEAMADAPTER_VERSION} ${PROJECT_NAME} BeamAdapter)

# if only this was possible...
#if(SOFA_BUILD_TESTS)
#    add_subdirectory(BeamAdapter_test)
#endif()

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
