cmake_minimum_required(VERSION 2.8)
project("BeamAdapter")

set(BEAMADAPTER_MAJOR_VERSION 1)
set(BEAMADAPTER_MINOR_VERSION 0)
set(BEAMADAPTER_VERSION ${BEAMADAPTER_MAJOR_VERSION}.${BEAMADAPTER_MINOR_VERSION})

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

set(HEADER_FILES
  component/initBeamAdapter.h
  component/BeamInterpolation.h
  component/BeamInterpolation.inl
  component/WireBeamInterpolation.h
  component/WireBeamInterpolation.inl
  )

set(SOURCE_FILES
  component/initBeamAdapter.cpp
  component/BeamInterpolation.cpp
  component/WireBeamInterpolation.cpp
  )

include(component/controller/Controller.cmake)
include(component/engine/Engine.cmake)
include(component/forcefield/Forcefield.cmake)
include(component/mapping/Mapping.cmake)

set(README_FILE "BeamAdapter.txt" )


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

SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
SET_PROPERTY(TARGET ${PROJECT_NAME}	PROPERTY FOLDER "plugins")