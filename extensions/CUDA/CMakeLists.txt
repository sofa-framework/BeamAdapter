cmake_minimum_required(VERSION 3.12)
project(BeamAdapter.CUDA)

set(HEADER_FILES
    src/BeamAdapter.CUDA/init.h
    src/BeamAdapter.CUDA/config.h.in
)

set(SOURCE_FILES
    src/BeamAdapter.CUDA/init.cpp
    src/BeamAdapter.CUDA/CudaInstantiations.cpp
)

sofa_find_package(BeamAdapter REQUIRED)
sofa_find_package(SofaCUDA REQUIRED)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} BeamAdapter)
target_link_libraries(${PROJECT_NAME} SofaCUDA)

sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${Sofa_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR "${PROJECT_NAME}"
    RELOCATABLE "plugins"
)
