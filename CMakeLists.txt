cmake_minimum_required(VERSION 3.15)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True) # Ensure the exact standard is required
set(CMAKE_CXX_EXTENSIONS OFF)         # Avoid GNU-specific extensions

# Default to Release build type if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()


set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_DEBUG "-g --coverage")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_MODULE_PATH} )

# Find LEMON headers dynamically in the Conda environment
set(LIBLEMON_INCLUDE_DIR "$ENV{CONDA_PREFIX}/include")
include_directories(${LIBLEMON_INCLUDE_DIR})

# Add Python3 and pybind11
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
add_subdirectory(extern/pybind11)

# Define the Python extension module
pybind11_add_module(clonelib
    clonelib/clonelib.cpp
    clonelib/basematrix.cpp
    clonelib/readmatrix.cpp
    clonelib/basetree.cpp
    clonelib/cnagraph.cpp
    clonelib/cnatree.cpp
    clonelib/genotypetree.cpp
    clonelib/utils.cpp
)

# Install the Python extension module
install(TARGETS clonelib DESTINATION ${SKBUILD_PROJECT_NAME})



