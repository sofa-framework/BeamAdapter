cmake_minimum_required(VERSION 2.8.12)

project(BeamAdapter_test)

set(SOURCE_FILES
    BeamInterpolation_test.cpp
    WireBeamInterpolation_test.cpp)

set(HEADER_FILES )

include(component/forcefield/ForceFieldTest.cmake)
include(component/constraint/ConstraintTest.cmake)
include(component/mapping/MappingTest.cmake)

add_executable(${PROJECT_NAME} ${SOURCE_FILES} ${RC_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>")

find_package(SofaFramework REQUIRED)
find_package(SofaCommon REQUIRED)
find_package(SofaTest REQUIRED)
find_package(BeamAdapter REQUIRED)

target_link_libraries(${PROJECT_NAME} SofaTest SofaGTestMain BeamAdapter SofaCore SofaConstraint SofaBaseMechanics SofaUserInteraction SofaComponentCommon)

add_definitions("-DBEAMADAPTER_TEST_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/\"")

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
