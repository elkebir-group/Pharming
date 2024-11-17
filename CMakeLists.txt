cmake_minimum_required(VERSION 3.15)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
# set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_DEBUG "-g --coverage")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_MODULE_PATH} )

# set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

# include_directories( "${LIBLEMON_ROOT}/include" )
# link_directories("${LIBLEMON_ROOT}/lib" )

# set(LIBLEMON_ROOT "/scratch/data/leah/anaconda3/envs/pharming_test" CACHE PATH "LEMON library directory")

# include_directories(${LIBLEMON_ROOT}/include)
# link_directories(${LIBLEMON_ROOT}/lib)

# Dynamically find LEMON headers and library
set(LIBLEMON_INCLUDE_DIR "$ENV{CONDA_PREFIX}/include")

# Include and link all LEMON files
include_directories(${LIBLEMON_INCLUDE_DIR})
# link_directories(${LIBLEMON_LIBRARY})

# Add Lemon and Boost as dependencies
# find_package(Boost REQUIRED)
# find_package(lemon REQUIRED PATHS ${LEMON_ROOT})

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
# setup PyBind11
set(PYBIND11_FINDPYTHON ON)

# Add pybind11 (from the extern submodule)
add_subdirectory(extern/pybind11)


# Python extension module using pybind11
pybind11_add_module(_clonelib
    clonelib/clonelib.cpp
    clonelib/basematrix.cpp
    clonelib/readmatrix.cpp
    clonelib/basetree.cpp
    clonelib/cnagraph.cpp
    clonelib/cnatree.cpp
    clonelib/genotypetree.cpp
    clonelib/utils.cpp
    WITH_SOABI
)
install(TARGETS _clonelib DESTINATION ${SKBUILD_PROJECT_NAME})

# install(TARGETS clonelib
#         COMPONENT python
#         LIBRARY DESTINATION "lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages"
#         ARCHIVE DESTINATION "lib"
#         RUNTIME DESTINATION "bin")
# target_include_directories(clonelib PRIVATE pharming/clonelib)

# add_subdirectory(pharming/clonelib)

# # Include paths for pybind11, Lemon, and Boost
# target_include_directories(clonelib PRIVATE
#     pharming/clonelib
#     extern/pybind11/include
#     ${Boost_INCLUDE_DIRS}
#     ${LEMON_INCLUDE_DIRS}
# )

# # Link Boost and Lemon libraries
# target_link_libraries(clonelib PRIVATE
#     ${Boost_LIBRARIES}
#     ${LEMON_LIBRARIES}
# )

# # Additional compiler flags
# if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#     target_compile_options(clonelib PRIVATE -Wall -Wextra -Wpedantic)
# endif()

# # Debug options
# if(CMAKE_BUILD_TYPE STREQUAL "Debug")
#     add_compile_definitions(DEBUG)
#     target_compile_options(clonelib PRIVATE -g)
# endif()