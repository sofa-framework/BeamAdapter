cmake_minimum_required(VERSION 3.12)
project(BeamAdapter VERSION 1.0)

include(cmake/environment.cmake)

find_package(Sofa.Simulation.Core REQUIRED)

sofa_find_package(Sofa.Component.StateContainer REQUIRED)
sofa_find_package(Sofa.Component.Controller REQUIRED)
sofa_find_package(Sofa.Component.Topology.Container.Dynamic REQUIRED)
sofa_find_package(Sofa.Component.Topology.Mapping REQUIRED)
sofa_find_package(Sofa.Component.Collision.Geometry REQUIRED)
sofa_find_package(Sofa.Component.Constraint.Projective REQUIRED)
sofa_find_package(Sofa.Component.Constraint.Lagrangian REQUIRED)
sofa_find_package(SofaImplicitField QUIET)
sofa_find_package(SofaAdvancedConstraint QUIET)
sofa_find_package(SofaCUDA QUIET)

set(BEAMADAPTER_SRC "src/${PROJECT_NAME}")

set(HEADER_FILES
    ${BEAMADAPTER_SRC}/config.h.in
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
    ${BEAMADAPTER_SRC}/component/controller/BeamAdapterActionController.h
    ${BEAMADAPTER_SRC}/component/controller/BeamAdapterActionController.inl    
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
    ${BEAMADAPTER_SRC}/component/mapping/BeamLengthMapping.h
    ${BEAMADAPTER_SRC}/component/mapping/BeamLengthMapping.inl
    ${BEAMADAPTER_SRC}/component/mapping/MultiAdaptiveBeamMapping.h
    ${BEAMADAPTER_SRC}/component/mapping/MultiAdaptiveBeamMapping.inl

    ${BEAMADAPTER_SRC}/component/model/BaseRodSectionMaterial.h
    ${BEAMADAPTER_SRC}/component/model/BaseRodSectionMaterial.inl
    ${BEAMADAPTER_SRC}/component/model/RodMeshSection.h
    ${BEAMADAPTER_SRC}/component/model/RodMeshSection.inl
    ${BEAMADAPTER_SRC}/component/model/RodSpireSection.h
    ${BEAMADAPTER_SRC}/component/model/RodSpireSection.inl
    ${BEAMADAPTER_SRC}/component/model/RodStraightSection.h
    ${BEAMADAPTER_SRC}/component/model/RodStraightSection.inl

    ${BEAMADAPTER_SRC}/utils/BeamSection.h
    ${BEAMADAPTER_SRC}/utils/BeamActions.h
    ${BEAMADAPTER_SRC}/utils/deprecatedcomponent.h
)

set(SOURCE_FILES
    ${BEAMADAPTER_SRC}/initBeamAdapter.cpp

    ${BEAMADAPTER_SRC}/component/BeamInterpolation.cpp
    ${BEAMADAPTER_SRC}/component/WireBeamInterpolation.cpp

    ${BEAMADAPTER_SRC}/component/constraint/AdaptiveBeamSlidingConstraint.cpp
    ${BEAMADAPTER_SRC}/component/constraint/AdaptiveBeamLengthConstraint.cpp

    ${BEAMADAPTER_SRC}/component/controller/AdaptiveBeamController.cpp
    ${BEAMADAPTER_SRC}/component/controller/BeamAdapterActionController.cpp
    ${BEAMADAPTER_SRC}/component/controller/InterventionalRadiologyController.cpp
    ${BEAMADAPTER_SRC}/component/controller/SutureController.cpp

    ${BEAMADAPTER_SRC}/component/engine/WireRestShape.cpp
    ${BEAMADAPTER_SRC}/component/engine/SteerableCatheter.cpp

    ${BEAMADAPTER_SRC}/component/forcefield/AdaptiveBeamForceFieldAndMass.cpp
    ${BEAMADAPTER_SRC}/component/forcefield/AdaptiveInflatableBeamForceField.cpp

    ${BEAMADAPTER_SRC}/component/mapping/AdaptiveBeamMapping.cpp
    ${BEAMADAPTER_SRC}/component/mapping/BeamLengthMapping.cpp
    ${BEAMADAPTER_SRC}/component/mapping/MultiAdaptiveBeamMapping.cpp
	
    ${BEAMADAPTER_SRC}/component/model/RodMeshSection.cpp
    ${BEAMADAPTER_SRC}/component/model/RodSpireSection.cpp
    ${BEAMADAPTER_SRC}/component/model/RodStraightSection.cpp
)

if(SofaImplicitField_FOUND)
    set(HEADER_FILES ${HEADER_FILES}
        ${BEAMADAPTER_SRC}/component/constraint/ImplicitSurfaceAdaptiveConstraint.h
        ${BEAMADAPTER_SRC}/component/constraint/ImplicitSurfaceAdaptiveConstraint.inl
        )

    set(SOURCE_FILES ${SOURCE_FILES}
        ${BEAMADAPTER_SRC}/component/constraint/ImplicitSurfaceAdaptiveConstraint.cpp
        )
endif()

if(SofaAdvancedConstraint_FOUND)
    set(HEADER_FILES ${HEADER_FILES}
        ${BEAMADAPTER_SRC}/component/AdaptiveBeamContactMapper.h
        ${BEAMADAPTER_SRC}/component/AdaptiveBeamContactMapper.inl
#        ${BEAMADAPTER_SRC}/component/AdaptiveBeamFrictionContact.h
#        ${BEAMADAPTER_SRC}/component/AdaptiveBeamFrictionContact.inl
        ${BEAMADAPTER_SRC}/component/MultiAdaptiveBeamContactMapper.h
        ${BEAMADAPTER_SRC}/component/MultiAdaptiveBeamContactMapper.inl
        )

    set(SOURCE_FILES ${SOURCE_FILES}
        ${BEAMADAPTER_SRC}/component/AdaptiveBeamContactMapper.cpp
        ${BEAMADAPTER_SRC}/component/AdaptiveBeamFrictionContact.cpp
        ${BEAMADAPTER_SRC}/component/MultiAdaptiveBeamContactMapper.cpp
        )

    #find_package(CUDA REQUIRED)
    #cuda_include_directories(${SOFA_CUDA_DIR})
    #cuda_include_directories(${SOFA_FRAMEWORK_DIR})
    #cuda_add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${CUDA_SOURCES})
endif()


if(SofaCUDA_FOUND)
    list (APPEND SOURCE_FILES
        ${BEAMADAPTER_SRC}/gpu/cuda/CudaInstantiations.cpp
    )
    message(STATUS "BeamAdapter: SofaCUDA has been enabled, instantiations for CudaTypes will be enabled.")
endif()

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${README_FILES})

target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Simulation.Core)
target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Component.StateContainer Sofa.Component.Controller Sofa.Component.Topology.Container.Dynamic Sofa.Component.Topology.Mapping)
target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Component.Collision.Geometry Sofa.Component.Constraint.Lagrangian Sofa.Component.Constraint.Projective)

if(SofaImplicitField_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC SofaImplicitField)
endif()
if(SofaAdvancedConstraint_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC SofaAdvancedConstraint)
endif()
if(SofaCUDA_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC SofaCUDA)
endif()

sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${PROJECT_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR ${PROJECT_NAME}
    EXAMPLE_INSTALL_DIR "examples"
    RELOCATABLE "plugins"
    )

# Tests
# If SOFA_BUILD_TESTS exists and is OFF, then these tests will be auto-disabled
cmake_dependent_option(BEAMADAPTER_BUILD_TESTS "Compile the automatic tests" ON "SOFA_BUILD_TESTS OR NOT DEFINED SOFA_BUILD_TESTS" OFF)
if(BEAMADAPTER_BUILD_TESTS)
    enable_testing()
    add_subdirectory(${PROJECT_NAME}_test)
endif()


# this paragraph has not been ported to the new build system yet...
#if(SOFA-PLUGIN_SOFAADVANCEDCONSTRAINT)
    #set(LINKER_DEPENDENCIES SofaAdvancedConstraint)
    #list(APPEND LINKER_DEPENDENCIES SofaCUDA)
    #list(APPEND LINKER_DEPENDENCIES SofaCUDADev)
    #set(LINKER_DEPENDENCIES SofaAdvancedConstraint)
#endif()

# Organize projects into folders
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
SET_PROPERTY(TARGET ${PROJECT_NAME} PROPERTY FOLDER "plugins")

# Config files and install rules for pythons scripts
sofa_install_pythonscripts(PLUGIN_NAME ${PROJECT_NAME} PYTHONSCRIPTS_SOURCE_DIR "python")

include(cmake/packaging.cmake)
