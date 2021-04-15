cmake_minimum_required(VERSION 3.1)
project(BeamAdapter VERSION 1.0)

include(cmake/environment.cmake)

find_package(SofaFramework REQUIRED)
find_package(SofaImplicitField REQUIRED) # for ImplicitSurfaceAdaptiveConstraint
find_package(SofaAdvancedConstraint QUIET)

set(HEADER_FILES
    component/initBeamAdapter.h

    component/BeamInterpolation.h
    component/BeamInterpolation.inl
    component/WireBeamInterpolation.h
    component/WireBeamInterpolation.inl

    component/constraint/AdaptiveBeamSlidingConstraint.h
    component/constraint/AdaptiveBeamSlidingConstraint.inl
    component/constraint/AdaptiveBeamLengthConstraint.h
    component/constraint/AdaptiveBeamLengthConstraint.inl
#    component/constraint/ImplicitSurfaceAdaptiveConstraint.h
#    component/constraint/ImplicitSurfaceAdaptiveConstraint.inl

    component/controller/AdaptiveBeamController.h
    component/controller/AdaptiveBeamController.inl
    component/controller/InterventionalRadiologyController.h
    component/controller/InterventionalRadiologyController.inl
    component/controller/SutureController.h
    component/controller/SutureController.inl

    component/engine/WireRestShape.h
    component/engine/WireRestShape.inl
    component/engine/SteerableCatheter.h
    component/engine/SteerableCatheter.inl

    component/forcefield/AdaptiveBeamForceFieldAndMass.h
    component/forcefield/AdaptiveBeamForceFieldAndMass.inl
    component/forcefield/AdaptiveInflatableBeamForceField.h
    component/forcefield/AdaptiveInflatableBeamForceField.inl
    component/forcefield/AdaptivePlasticBeamForceField.h
    component/forcefield/AdaptivePlasticBeamForceField.inl

    component/mapping/AdaptiveBeamMapping.h
    component/mapping/AdaptiveBeamMapping.inl
    component/mapping/MultiAdaptiveBeamMapping.h
    component/mapping/MultiAdaptiveBeamMapping.inl
    )

set(SOURCE_FILES
    component/initBeamAdapter.cpp

    component/BeamInterpolation.cpp
    component/WireBeamInterpolation.cpp

    component/constraint/AdaptiveBeamSlidingConstraint.cpp
    component/constraint/AdaptiveBeamLengthConstraint.cpp
#    component/constraint/ImplicitSurfaceAdaptiveConstraint.cpp

    component/controller/AdaptiveBeamController.cpp
    component/controller/InterventionalRadiologyController.cpp
    component/controller/SutureController.cpp

    component/engine/WireRestShape.cpp
    component/engine/SteerableCatheter.cpp

    component/forcefield/AdaptiveBeamForceFieldAndMass.cpp
    component/forcefield/AdaptiveInflatableBeamForceField.cpp
    component/forcefield/AdaptivePlasticBeamForceField.cpp

    component/mapping/AdaptiveBeamMapping.cpp
    component/mapping/MultiAdaptiveBeamMapping.cpp
    )

set(README_FILE "BeamAdapter.txt" )

if(SofaAdvancedConstraint_FOUND)
    set(HEADER_FILES ${HEADER_FILES}
        component/AdaptiveBeamContactMapper.h
        component/AdaptiveBeamContactMapper.inl
#        component/AdaptiveBeamFrictionContact.h
#        component/AdaptiveBeamFrictionContact.inl
        component/MultiAdaptiveBeamContactMapper.h
        component/MultiAdaptiveBeamContactMapper.inl
        )

    set(SOURCE_FILES ${SOURCE_FILES}
        component/AdaptiveBeamContactMapper.cpp
        component/AdaptiveBeamFrictionContact.cpp
        component/MultiAdaptiveBeamContactMapper.cpp
        )

    #find_package(CUDA REQUIRED)
    #cuda_include_directories(${SOFA_CUDA_DIR})
    #cuda_include_directories(${SOFA_FRAMEWORK_DIR})
    #cuda_add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${CUDA_SOURCES})
endif()


add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${README_FILES})

set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-DSOFA_BUILD_BEAMADAPTER")

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include>")

target_link_libraries(${PROJECT_NAME} SofaHelper SofaCore SofaConstraint SofaBaseMechanics SofaUserInteraction SofaImplicitField SofaLoader Sofa.GL)
sofa_create_package_with_targets(
    PACKAGE_NAME BeamAdapter
    PACKAGE_VERSION ${PROJECT_VERSION}
    TARGETS ${PROJECT_NAME}
    INCLUDE_INSTALL_DIR "BeamAdapter"
    EXAMPLE_INSTALL_DIR "examples"
    RELOCATABLE "$plugins"
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
