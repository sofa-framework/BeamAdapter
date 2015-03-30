cmake_minimum_required(VERSION 2.8)

project("BeamAdapter")

include(${SOFA_CMAKE_DIR}/preProject.cmake)

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

#message(status "BeamAdapter")


if(SOFA-PLUGIN_SOFAADVANCEDCONSTRAINT)
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

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} )

set(COMPILER_DEFINES "SOFA_BUILD_BEAMADAPTER")
set(LINKER_DEPENDENCIES SofaConstraint SofaComponentAdvanced SofaVolumetricData )


# Define the SOFA_WITH_FLOAT and SOFA_WITH_DOUBLE macro
# these two macro can be used to activate the compilation
# of float or double specific version of the components.
#
# the template to use is:
# ifdef SOFA_WITH_FLOAT
#    .... your "float" version of your code
# endif // SOFA_WITH_FLOAT
# ifdef SOFA_WITH_DOUBLE
#    .... your "double" version of your code
# endif // SOFA_WITH_DOUBLE
#
# This two macros should replace the SOFA_DOUBLE/SOFA_FLOAT
# that are very confusing because they were working
# in an "inverse" way:
#
# ifndef SOFA_FLOAT
#    .... the "double" version of the code
# endif // SOFA_FLOAT
# ifndef SOFA_DOUBLE
#    .... the "float" version of the code
# endif // SOFA_DOUBLE
if(SOFA-MISC_USE_FLOAT)
    AddCompilerDefinitions("SOFA_WITH_FLOAT")
endif(SOFA-MISC_USE_FLOAT)
if(SOFA-MISC_USE_DOUBLE)
    AddCompilerDefinitions("SOFA_WITH_DOUBLE")
endif(SOFA-MISC_USE_DOUBLE)
if(NOT SOFA-MISC_USE_FLOAT AND NOT SOFA-MISC_USE_DOUBLE)
    AddCompilerDefinitions("SOFA_WITH_FLOAT")
    AddCompilerDefinitions("SOFA_WITH_DOUBLE")
endif(NOT SOFA-MISC_USE_FLOAT AND NOT SOFA-MISC_USE_DOUBLE)


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


include(${SOFA_CMAKE_DIR}/postProject.cmake)
