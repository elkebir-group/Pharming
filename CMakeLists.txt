cmake_minimum_required(VERSION 3.5)



project(clonelib)

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_MODULE_PATH} )

set(visualizephylo_src
  src/arg_parser.cpp
  src/visualizephylomain.cpp
  src/basetree.cpp
  src/cnatree.cpp
  src/genotypetree.cpp
  src/utils.cpp
  src/phylogeny.cpp
)

set(visualizephylo_hdr
    src/basetree.h
    src/cnatree.h
    src/genotypetree.h
    src/utils.h
    src/phylogeny.h
)

set( mergegenotypetrees_src
    src/mergegenotypetreesmain.cpp
  src/basetree.cpp
  src/genotypetree.cpp
  src/genotypegraph.cpp
  src/cnatree.cpp
  src/cnagraph.cpp
  src/utils.cpp
  src/basematrix.cpp
)

set( mergegenotypetrees_hdr
  src/genotypetree.h
  src/genotypegraph.h
  src/basetree.h
  src/cnatree.h
  src/cnagraph.h
  src/utils.h
  src/basematrix.h
)

set( generategenotypetrees_src
  src/arg_parser.cpp
  src/generategenotypetreesmain.cpp
  src/cnatree.cpp
  src/cnagraph.cpp
  src/basetree.cpp
  src/genotypetree.cpp
  src/genotypegraph.cpp
  src/utils.cpp
  src/basematrix.cpp
)

set( generategenotypetrees_hdr
  src/arg_parser.h
  src/basetree.h
  src/cnatree.h
  src/cnagraph.h
  src/genotypetree.h
  src/genotypegraph.h
  src/utils.h
  src/basematrix.h
)

set( generatecnatrees_src
    src/arg_parser.cpp
    src/generatecnatreesmain.cpp
    src/cnatree.cpp
    src/cnagraph.cpp
    src/genotypetree.cpp
    src/utils.cpp
    src/basematrix.cpp
    src/basetree.cpp
    )

set( generatecnatrees_hdr
    src/arg_parser.h
    src/cnatree.h
    src/cnagraph.h
    src/genotypetree.h
    src/utils.h
    src/basematrix.h
    src/basetree.h
    )

set( simulate_src
  src/arg_parser.cpp
  src/simulatemain.cpp
  src/basetree.cpp
  src/cnatree.cpp
  src/genotypetree.cpp
  src/genotypegraph.cpp
  src/cnagraph.cpp
  src/utils.cpp
  src/basematrix.cpp
  src/phylogeny.cpp
  src/singlecellgeneration.cpp
)

set( simulate_hdr
  src/arg_parser.h
  src/cnatree.h
  src/basetree.h
  src/genotypetree.h
  src/genotypegraph.h
  src/cnagraph.h
  src/utils.h
  src/basematrix.h
  src/phylogeny.h
  src/singlecellgeneration.h
)

set ( generatesinglecells_src
        src/arg_parser.cpp
        src/singlecellgeneration.cpp
        src/utils.cpp
        src/generatesinglecellsmain.cpp
        src/generatesinglecellsmain.cpp
)

set ( generatesinglecells_hdr
        src/arg_parser.h
        src/utils.h
        src/singlecellgeneration.h
        src/generatesinglecellsmain.cpp
)

set ( clonelib_hdr
        src/cnatree.h
        src/basetree.h
        src/genotypetree.h
        src/cnagraph.h
        src/basematrix.h
        src/utils.h
)


set ( clonelib_src
        src/cnatree.cpp
        src/basetree.cpp
        src/genotypetree.cpp
        src/basematrix.cpp
        src/utils.cpp
        src/cnagraph.cpp

)

set (clonelib_py
  src/pybind/gen_trees_pybind.cpp
)

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs" )
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
#set( Boost_USE_STATIC_LIBS ON )

find_package( Boost COMPONENTS random )

MESSAGE( STATUS "found boost:"
  "\nINCLUDE: ${Boost_INCLUDE_DIRS}"
  "\nLIB: ${Boost_LIBRARIES}"
)

include_directories( "${LIBLEMON_ROOT}/include" ${Boost_INCLUDE_DIRS} )
link_directories( ${Boost_LIBRARY_DIRS} "${LIBLEMON_ROOT}/lib" ${GUROBI_LIBRARIES} ${CPLEX_LIB_DIR} ${CONCERT_LIB_DIR} )

set( CommonLibs
  m
  pthread
  ${Boost_LIBRARIES}
)

set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG" )

add_executable( mergegenotypetrees ${mergegenotypetrees_src} ${mergegenotypetrees_hdr} )
target_link_libraries( mergegenotypetrees ${CommonLibs} )

add_executable( generategenotypetrees ${generategenotypetrees_src} ${generategenotypetrees_hdr} )
target_link_libraries( generategenotypetrees ${CommonLibs} )

add_executable( generatecnatrees ${generatecnatrees_src} ${generatecnatrees_hdr} )
target_link_libraries( generatecnatrees ${CommonLibs} )

add_executable( simulate ${simulate_src} ${simulate_hdr} )
target_link_libraries( simulate ${CommonLibs} )

add_executable( visualizephylo ${visualizephylo_src} ${visualizephylo_hdr} )
target_link_libraries( visualizephylo ${CommonLibs} )

add_executable( generatesinglecells ${generatesinglecells_src} ${generatesinglecells_hdr}
        src/generatesinglecellsmain.cpp)

target_link_libraries( generatesinglecells ${CommonLibs})


if(NOT PYTHON_EXECUTABLE)
    set(Python_FIND_VIRTUALENV ONLY)
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
    # set(Python3_EXECUTABLE "${Python3_EXECUTABLE}" CACHE FILEPATH "Path to Python interpreter")
    set(Python3_EXECUTABLE "${PYTHON_EXECUTABLE}" CACHE FILEPATH "Path to Python interpreter")

else()
  set(Python3_EXECUTABLE "${PYTHON_EXECUTABLE}")
  find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
  execute_process(COMMAND dirname "${PYTHON_EXECUTABLE}" OUTPUT_VARIABLE CONDA_ENV_PATH OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "Clonelib python module installation directory: ${CONDA_ENV_PATH}")
endif()


message(STATUS "Python interpreter: ${PYTHON_EXECUTABLE}")
message(STATUS "Python3_VERSION_MAJOR: ${Python3_VERSION_MAJOR}")
message(STATUS "Python3_VERSION_MINOR: ${Python3_VERSION_MINOR}")





# Install the target to the appropriate locations within the Conda environment

find_package(pybind11 REQUIRED)
pybind11_add_module(clonelib ${clonelib_hdr} ${clonelib_src} ${clonelib_py})

install(TARGETS clonelib
    COMPONENT python
    LIBRARY DESTINATION "${CONDA_ENV_PATH}/lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages"
    ARCHIVE DESTINATION "${CONDA_ENV_PATH}/lib"
    RUNTIME DESTINATION "${CONDA_ENV_PATH}/bin"
)

