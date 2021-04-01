cmake_minimum_required(VERSION 3.1)
project(BeamAdapter VERSION 1.0)

include(cmake/environment.cmake)

find_package(SofaConstraint REQUIRED)
find_package(SofaBaseMechanics REQUIRED)
find_package(SofaUserInteraction REQUIRED)
find_package(SofaLoader REQUIRED)

set(BEAMADAPTER_SRC "src/BeamAdapter")

set(HEADER_FILES
    ${BEAMADAPTER_SRC}/initBeamAdapter.h

    ${BEAMADAPTER_SRC}/component/BeamInterpolation.h
    ${BEAMADAPTER_SRC}/component/BeamInterpolation.inl
    ${BEAMADAPTER_SRC}/component/WireBeamInterpolation.h
    ${BEAMADAPTER_SRC}/component/WireBeamInterpolation.inl

    ${BEAMADAPTER_SRC}/component/constraint/AdaptiveBeamSlidingConstraint.h
    ${BEAMADAPTER_SRC}/component/constraint/AdaptiveBeamSlidingConstraint.inl
    ${BEAMADAPTER_SRC}/component/constraint/AdaptiveBeamLengthConstraint.h
    ${BEAMADAPTER_SRC}/component/constraint/AdaptiveBeamLengthConstraint.inl

    ${BEAMADAPTER_SRC}/component/controller/AdaptiveBeamController.h
    ${BEAMADAPTER_SRC}/component/controller/AdaptiveBeamController.inl
    ${BEAMADAPTER_SRC}/component/controller/InterventionalRadiologyController.h
    ${BEAMADAPTER_SRC}/component/controller/InterventionalRadiologyController.inl
    ${BEAMADAPTER_SRC}/component/controller/SutureController.h
    ${BEAMADAPTER_SRC}/component/controller/SutureController.inl

    ${BEAMADAPTER_SRC}/component/engine/WireRestShape.h
    ${BEAMADAPTER_SRC}/component/engine/WireRestShape.inl
    ${BEAMADAPTER_SRC}/component/engine/SteerableCatheter.h
    ${BEAMADAPTER_SRC}/component/engine/SteerableCatheter.inl

    ${BEAMADAPTER_SRC}/component/forcefield/AdaptiveBeamForceFieldAndMass.h
    ${BEAMADAPTER_SRC}/component/forcefield/AdaptiveBeamForceFieldAndMass.inl
    ${BEAMADAPTER_SRC}/component/forcefield/AdaptiveInflatableBeamForceField.h
    ${BEAMADAPTER_SRC}/component/forcefield/AdaptiveInflatableBeamForceField.inl

    ${BEAMADAPTER_SRC}/component/mapping/AdaptiveBeamMapping.h
    ${BEAMADAPTER_SRC}/component/mapping/AdaptiveBeamMapping.inl
    ${BEAMADAPTER_SRC}/component/mapping/MultiAdaptiveBeamMapping.h
    ${BEAMADAPTER_SRC}/component/mapping/MultiAdaptiveBeamMapping.inl
    )

set(SOURCE_FILES
    ${BEAMADAPTER_SRC}/initBeamAdapter.cpp

    ${BEAMADAPTER_SRC}/component/BeamInterpolation.cpp
    ${BEAMADAPTER_SRC}/component/WireBeamInterpolation.cpp

    ${BEAMADAPTER_SRC}/component/constraint/AdaptiveBeamSlidingConstraint.cpp
    ${BEAMADAPTER_SRC}/component/constraint/AdaptiveBeamLengthConstraint.cpp

    ${BEAMADAPTER_SRC}/component/controller/AdaptiveBeamController.cpp
    ${BEAMADAPTER_SRC}/component/controller/InterventionalRadiologyController.cpp
    ${BEAMADAPTER_SRC}/component/controller/SutureController.cpp

    ${BEAMADAPTER_SRC}/component/engine/WireRestShape.cpp
    ${BEAMADAPTER_SRC}/component/engine/SteerableCatheter.cpp

    ${BEAMADAPTER_SRC}/component/forcefield/AdaptiveBeamForceFieldAndMass.cpp
    ${BEAMADAPTER_SRC}/component/forcefield/AdaptiveInflatableBeamForceField.cpp

    ${BEAMADAPTER_SRC}/component/mapping/AdaptiveBeamMapping.cpp
    ${BEAMADAPTER_SRC}/component/mapping/MultiAdaptiveBeamMapping.cpp
    )

set(README_FILE "BeamAdapter.txt" )

find_package(SofaAdvancedConstraint QUIET)
if(SofaAdvancedConstraint_FOUND)
    list(APPEND HEADER_FILES
        ${BEAMADAPTER_SRC}/component/AdaptiveBeamContactMapper.h
        ${BEAMADAPTER_SRC}/component/AdaptiveBeamContactMapper.inl
#        ${BEAMADAPTER_SRC}/component/AdaptiveBeamFrictionContact.h
#        ${BEAMADAPTER_SRC}/component/AdaptiveBeamFrictionContact.inl
        ${BEAMADAPTER_SRC}/component/MultiAdaptiveBeamContactMapper.h
        ${BEAMADAPTER_SRC}/component/MultiAdaptiveBeamContactMapper.inl
        )

    list(APPEND SOURCE_FILES
        ${BEAMADAPTER_SRC}/component/AdaptiveBeamContactMapper.cpp
        ${BEAMADAPTER_SRC}/component/AdaptiveBeamFrictionContact.cpp
        ${BEAMADAPTER_SRC}/component/MultiAdaptiveBeamContactMapper.cpp
        )

    #find_package(CUDA REQUIRED)
    #cuda_include_directories(${SOFA_CUDA_DIR})
    #cuda_include_directories(${SOFA_FRAMEWORK_DIR})
    #cuda_add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${CUDA_SOURCES})
endif()

sofa_find_package(SofaImplicitField QUIET) # for ImplicitSurfaceAdaptiveConstraint
if(SofaImplicitField_FOUND)
    list(APPEND HEADER_FILES
        # ${BEAMADAPTER_SRC}/component/constraint/ImplicitSurfaceAdaptiveConstraint.h
        # ${BEAMADAPTER_SRC}/component/constraint/ImplicitSurfaceAdaptiveConstraint.inl
        )
    list(APPEND SOURCE_FILES
        # ${BEAMADAPTER_SRC}/component/constraint/ImplicitSurfaceAdaptiveConstraint.cpp
        )
endif()

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${README_FILES})

target_link_libraries(${PROJECT_NAME} SofaCore SofaConstraint SofaBaseMechanics SofaUserInteraction SofaLoader)
if(SofaImplicitField_FOUND)
    target_link_libraries(${PROJECT_NAME} SofaImplicitField)
endif()

sofa_create_package_with_targets(
    PACKAGE_NAME BeamAdapter
    PACKAGE_VERSION ${PROJECT_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR "BeamAdapter"
    EXAMPLE_INSTALL_DIR "examples"
    RELOCATABLE "plugins"
    )

option(BEAMADAPTER_BUILD_TESTS "Build unit tests" ON)
if(BEAMADAPTER_BUILD_TESTS)
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
SET_PROPERTY(TARGET ${PROJECT_NAME} PROPERTY FOLDER "plugins")

include(cmake/packaging.cmake)
