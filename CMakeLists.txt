
cmake_minimum_required(VERSION 3.18)
#cmake_policy(SET CMP0074 NEW)
#cmake_policy(SET CMP0148 NEW)

project(IPknot CXX C)
option(CXX "enable C++ compilation" ON)
if (CXX)
  enable_language(CXX)
endif()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
add_definitions(-DPACKAGE_VERSION="1.1.0")

include_directories("spdlog/include")
include_directories("cxxopts/include")

find_package(Python COMPONENTS Interpreter Development)
#find_package(Python REQUIRED COMPONENTS Interpreter Development.Embed)
find_package(pybind11 REQUIRED)

find_package(PkgConfig REQUIRED)

pkg_check_modules(VIENNARNA RNAlib2)
if (VIENNARNA_FOUND)
  include_directories(${VIENNARNA_INCLUDE_DIRS})
  add_definitions(-DHAVE_VIENNA18)
  add_definitions(-DHAVE_VIENNA20)
endif()

set(EXTERNAL_LIBS "")
if(ENABLE_CPLEX)
  find_package(CPLEX)
  if (CPLEX_FOUND)
    add_definitions(-DWITH_CPLEX)
    include_directories(${CPLEX_INCLUDE_DIR})
    list(APPEND EXTERNAL_LIBS ${CPLEX_LIBRARY})
  endif()
elseif(ENABLE_GUROBI)
  find_package(GUROBI)
  if (GUROBI_FOUND)
    add_definitions(-DWITH_GUROBI)
    include_directories(${GUROBI_INCLUDE_DIR})
    list(APPEND EXTERNAL_LIBS ${GUROBI_LIBRARY})
  endif()
elseif(ENABLE_SCIP)
  find_package(SCIP)
  if (SCIP_FOUND)
    message(STATUS "Found SCIP version ${SCIP_VERSION}")
    add_definitions(-DWITH_SCIP)
    include_directories(${SCIP_INCLUDE_DIRS})
    list(APPEND EXTERNAL_LIBS ${SCIP_LIBRARIES})
  endif()
elseif(ENABLE_HIGHS)
  find_package(HiGHS)
  find_package(ZLIB 1.2.3)
  if (HiGHS_FOUND)
    add_definitions(-DWITH_HIGHS)
    include_directories(${HiGHS_INCLUDE_DIR})
    list(APPEND EXTERNAL_LIBS ${HiGHS_LIBRARY} ${ZLIB_LIBRARY})
  endif()
else()
  find_package(GLPK)
  if (GLPK_FOUND)
    add_definitions(-DWITH_GLPK)
    include_directories(${GLPK_INCLUDE_DIR})
    list(APPEND EXTERNAL_LIBS ${GLPK_LIBRARY} ${ZLIB_LIBRARY} ${GMP_LIBRARY} ${DL_LIBRARY} ${LTDL_LIBRARY} ${MPFR_LIBRARY})
  endif()
endif()

if(STATIC_BUILD) # static build for linux system
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

add_executable(
  ipknot
  src/main.cpp
  src/ipknot.h src/ipknot.cpp
  src/ip.h src/ip.cpp
  src/fa.cpp 
  src/aln.cpp 
  src/fold.cpp 
  src/mxfold2.cpp
  src/bpseq.h src/bpseq.cpp
  src/boltzmann_param.c 
  src/nupack/nupack.cpp
  src/contrafold/LogSpace.hpp 
  src/contrafold/Defaults.hpp 
  src/contrafold/InferenceEngine.cpp
  src/contrafold/ParameterManager.cpp 
  src/contrafold/SStruct.cpp 
  src/contrafold/Utilities.cpp 
  src/linearpartition/LinearPartition.cpp 
  #src/linearpartition/bpp.cpp
  )
target_link_libraries(ipknot PRIVATE pybind11::embed
  ${VIENNARNA_LDFLAGS} ${EXTERNAL_LIBS} ${CMAKE_DL_LIBS})

target_compile_features(ipknot PUBLIC cxx_std_17)

install(TARGETS ipknot RUNTIME DESTINATION bin)
