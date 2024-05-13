cmake_minimum_required(VERSION 3.12)

project(BeamAdapter_test VERSION 1.0)

find_package(Sofa.Testing REQUIRED)
find_package(Threads REQUIRED)

set(SOURCE_FILES
    BeamInterpolation_test.cpp
    component/controller/AdaptiveBeamController_test.cpp
    component/controller/InterventionalRadiologyController_test.cpp
    component/constraint/AdaptiveBeamSlidingConstraintTest.cpp
    component/forcefield/AdaptiveBeamForceFieldAndMassTest.cpp
    component/model/WireRestShape_test.cpp
)

# Sofa.Component.Mapping.Testing does not have its own package configuration
if (TARGET Sofa.Component.Mapping.Testing)
    list(APPEND SOURCE_FILES component/mapping/BeamLengthMappingTest.cpp)
endif()

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>")

target_link_libraries(${PROJECT_NAME} PUBLIC BeamAdapter Sofa.Testing Threads::Threads)
if (TARGET Sofa.Component.Mapping.Testing)
    target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Component.Mapping.Testing)
endif()

add_definitions("-DBEAMADAPTER_TEST_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/\"")

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
