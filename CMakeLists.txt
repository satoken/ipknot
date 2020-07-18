
cmake_minimum_required(VERSION 3.1)

#enable_language(cxx)
#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

project(IPknot)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
add_definitions(-DPACKAGE_VERSION="0.0.5")

find_package(PkgConfig REQUIRED)

pkg_check_modules(VIENNARNA REQUIRED RNAlib2)
include_directories(${VIENNARNA_INCLUDE_DIRS})
add_definitions(-DHAVE_VIENNA18)
add_definitions(-DHAVE_VIENNA20)

if(ENABLE_CPLEX)
  find_package(ILOG)
  if (ILOG_FOUND)
    add_definitions(-DWITH_CPLEX)
    include_directories(${ILOG_CPLEX_INCLUDE_DIR})
  endif()
elseif(ENABLE_GUROBI)
  find_package(GUROBI)
  if (GUROBI_FOUND)
    add_definitions(-DWITH_GUROBI)
    include_directories(${GUROBI_INCLUDE_DIR})
  endif()
else()
  find_package(GLPK 4.33)
  if(GLPK_FOUND)
    find_library(ZLIB_LIBRARY z)
    find_library(GMP_LIBRARY gmp)
    find_library(DL_LIBRARY dl)
    find_library(LTDL_LIBRARY ltdl)
    add_definitions(-DWITH_GLPK)
    include_directories(${GLPK_INCLUDE_DIR})
  endif()
endif()

add_executable(
  ipknot
  src/ipknot.cpp
  src/ip.h src/ip.cpp
  src/fa.cpp src/fa.h
  src/aln.cpp src/aln.h
  src/fold.cpp src/fold.h
  src/bpseq.h src/bpseq.cpp
  src/boltzmann_param.c src/boltzmann_param.h
  src/nupack/nupack.h src/nupack/nupack.cpp
  src/nupack/dptable.h src/nupack/def_param.h
	src/contrafold/Config.hpp src/contrafold/LogSpace.hpp 
  # src/contrafold/Defaults.ipp 
  src/contrafold/rand.h 
	src/contrafold/InferenceEngine.hpp #src/contrafold/InferenceEngine.ipp
	src/contrafold/ParameterManager.hpp #src/contrafold/ParameterManager.ipp 
	src/contrafold/SStruct.cpp src/contrafold/SStruct.hpp 
  src/contrafold/Utilities.cpp src/contrafold/Utilities.hpp #src/contrafold/Utilities.ipp 
  )
target_link_libraries(ipknot ${VIENNARNA_LDFLAGS} ${GLPK_LIBRARY} ${ZLIB_LIBRARY} ${GMP_LIBRARY} ${DL_LIBRARY} ${LTDL_LIBRARY} ${ILOG_CPLEX_LIBRARY} ${GUROBI_LIBRARY} ${GUROBI_CPP_LIBRARY})