
cmake_minimum_required(VERSION 3.8)

option(CXX "enable C++ compilation" ON)
if (CXX)
  enable_language(CXX)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILE_FEATURES cxx_std_17)

project(IPknot CXX C)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
add_definitions(-DPACKAGE_VERSION="1.1.0")

find_package(PkgConfig REQUIRED)

pkg_check_modules(VIENNARNA RNAlib2)
if (VIENNARNA_FOUND)
  include_directories(${VIENNARNA_INCLUDE_DIRS})
  add_definitions(-DHAVE_VIENNA18)
  add_definitions(-DHAVE_VIENNA20)
endif()

if(ENABLE_CPLEX)
  find_package(CPLEX)
  if (CPLEX_FOUND)
    add_definitions(-DWITH_CPLEX)
    include_directories(${CPLEX_INCLUDE_DIR})
  endif()
elseif(ENABLE_GUROBI)
  find_package(GUROBI)
  if (GUROBI_FOUND)
    add_definitions(-DWITH_GUROBI)
    include_directories(${GUROBI_INCLUDE_DIR})
  endif()
else()
  find_package(GLPK)
  if (GLPK_FOUND)
    add_definitions(-DWITH_GLPK)
    include_directories(${GLPK_INCLUDE_DIR})
  endif()
endif()

add_executable(
  ipknot
  src/ipknot.cpp
  src/ip.h src/ip.cpp
  src/fa.cpp 
  src/aln.cpp 
  src/fold.cpp 
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
  src/cxxopts.hpp
  )
target_link_libraries(ipknot 
  ${VIENNARNA_LDFLAGS} ${GLPK_LIBRARY} ${ZLIB_LIBRARY} ${GMP_LIBRARY} ${DL_LIBRARY} ${LTDL_LIBRARY} ${MPFR_LIBRARY}
  ${CPLEX_LIBRARY} ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY} ${CMAKE_DL_LIBS})
target_compile_features(ipknot PUBLIC cxx_std_17)

install(TARGETS ipknot RUNTIME DESTINATION bin)
