cmake_minimum_required(VERSION 3.1)

project(BeamAdapter_test VERSION 1.0)

#find_package(BeamAdapter REQUIRED)
find_package(Sofa.Testing REQUIRED)
find_package(SofaBaseUtils REQUIRED)
find_package(Threads REQUIRED)

set(SOURCE_FILES
    BeamInterpolation_test.cpp
    component/constraint/AdaptiveBeamSlidingConstraintTest.cpp
    component/forcefield/AdaptiveBeamForceFieldAndMassTest.cpp
    )

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>")

target_link_libraries(${PROJECT_NAME} BeamAdapter Sofa.Testing SofaCore SofaConstraint SofaBaseMechanics SofaUserInteraction SofaBaseUtils Threads::Threads)

add_definitions("-DBEAMADAPTER_TEST_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/\"")

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
