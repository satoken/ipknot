
cmake_minimum_required(VERSION 3.18)
#cmake_policy(SET CMP0074 NEW)
#cmake_policy(SET CMP0148 NEW)

project(IPknot CXX C)
option(CXX "enable C++ compilation" ON)
option(WITH_MXFOLD2 "Enable MXfold2 integration" OFF)
option(ENABLE_CPLEX "Enable CPLEX backend" OFF)
option(ENABLE_GUROBI "Enable Gurobi backend" OFF)
option(ENABLE_SCIP "Enable SCIP backend" OFF)
option(ENABLE_HIGHS "Enable HiGHS backend" OFF)
option(STATIC_BUILD "Build statically for Linux" OFF)
if (CXX)
  enable_language(CXX)
endif()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(PkgConfig REQUIRED)

if(STATIC_BUILD) # static build for linux system
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

add_executable(ipknot
  src/main.cpp
  src/ipknot.h src/ipknot.cpp
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
)

target_compile_definitions(ipknot PRIVATE PACKAGE_VERSION="1.1.0")
target_include_directories(ipknot PRIVATE "spdlog/include")
target_include_directories(ipknot PRIVATE "cxxopts/include")
target_compile_features(ipknot PUBLIC cxx_std_17)

pkg_check_modules(VIENNARNA RNAlib2)
if (VIENNARNA_FOUND)
  target_include_directories(ipknot PRIVATE ${VIENNARNA_INCLUDE_DIRS})
  target_compile_definitions(ipknot PRIVATE HAVE_VIENNA18 HAVE_VIENNA20)
  target_link_libraries(ipknot PRIVATE ${VIENNARNA_LDFLAGS} ${CMAKE_DL_LIBS})
endif()

if(ENABLE_CPLEX)
  find_package(CPLEX)
  if (CPLEX_FOUND)
    target_compile_definitions(ipknot PRIVATE WITH_CPLEX)
    target_include_directories(ipknot PRIVATE ${CPLEX_INCLUDE_DIR})
    target_link_libraries(ipknot PRIVATE ${CPLEX_LIBRARY})
  endif()
elseif(ENABLE_GUROBI)
  find_package(GUROBI)
  if (GUROBI_FOUND)
    target_compile_definitions(ipknot PRIVATE WITH_GUROBI)
    target_include_directories(ipknot PRIVATE ${GUROBI_INCLUDE_DIR})
    target_link_libraries(ipknot PRIVATE ${GUROBI_LIBRARY})
  endif()
elseif(ENABLE_SCIP)
  find_package(SCIP)
  if (SCIP_FOUND)
    message(STATUS "Found SCIP version ${SCIP_VERSION}")
    target_compile_definitions(ipknot PRIVATE WITH_SCIP)
    target_include_directories(ipknot PRIVATE ${SCIP_INCLUDE_DIRS})
    target_link_libraries(ipknot PRIVATE ${SCIP_LIBRARIES})
  endif()
elseif(ENABLE_HIGHS)
  find_package(HiGHS)
  find_package(ZLIB 1.2.3)
  if (HiGHS_FOUND)
    target_compile_definitions(ipknot PRIVATE WITH_HIGHS)
    target_include_directories(ipknot PRIVATE ${HiGHS_INCLUDE_DIR})
    target_link_libraries(ipknot PRIVATE ${HiGHS_LIBRARY} ${ZLIB_LIBRARY})
  endif()
else()
  find_package(GLPK)
  if (GLPK_FOUND)
    target_compile_definitions(ipknot PRIVATE WITH_GLPK)
    target_include_directories(ipknot PRIVATE ${GLPK_INCLUDE_DIR})
    target_link_libraries(ipknot PRIVATE ${GLPK_LIBRARY} ${ZLIB_LIBRARY} ${GMP_LIBRARY} ${DL_LIBRARY} ${LTDL_LIBRARY} ${MPFR_LIBRARY})
  endif()
endif()

if (WITH_MXFOLD2)
  find_package(Python REQUIRED COMPONENTS Interpreter Development)
  find_package(pybind11 REQUIRED)
  target_sources(ipknot PRIVATE src/mxfold2.cpp)
  target_link_libraries(ipknot PRIVATE pybind11::embed)
  target_compile_definitions(ipknot PRIVATE WITH_MXFOLD2)
endif()

install(TARGETS ipknot RUNTIME DESTINATION bin)

include(CTest)
if(BUILD_TESTING)
  set(NMR_TEST_FASTA "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nmr_constraints.fa")
  set(NMR_COMPACT_TEST_FASTA "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nmr_compact_stack.fa")
  set(NMR_BULGE_TEST_FASTA "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nmr_single_bulge.fa")
  set(NMR_RIGHT_BULGE_TEST_FASTA "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nmr_single_bulge_right.fa")
  set(NMR_MIXED_BULGE_TEST_FASTA "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nmr_bulge_mixed.fa")
  set(NMR_COAXIAL_TEST_FASTA "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nmr_coaxial_multibranch.fa")
  set(NMR_COAXIAL_CHILDREN_FASTA "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nmr_coaxial_children.fa")
  set(NMR_COAXIAL_LAST_FASTA "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nmr_coaxial_last_child.fa")

  add_test(NAME nmr_base_pair_noncanonical
    COMMAND ipknot -r 0 -t 0.01,0.01 --base-pairs UU=1
            --without-canonical-neighbor ${NMR_TEST_FASTA})

  add_test(NAME nmr_base_pair_impossible
    COMMAND ipknot -r 0 -t 0.01,0.01 --base-pairs GC=99 ${NMR_TEST_FASTA})
  set_tests_properties(nmr_base_pair_impossible PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_base_pair_negative
    COMMAND ipknot -r 0 -t 0.01,0.01 --base-pairs UU=-1 ${NMR_TEST_FASTA})
  set_tests_properties(nmr_base_pair_negative PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_base_pair_duplicate
    COMMAND ipknot -r 0 -t 0.01,0.01 --base-pairs GC=1,CG=2 ${NMR_TEST_FASTA})
  set_tests_properties(nmr_base_pair_duplicate PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_stack_space_delimited
    COMMAND ipknot -r 0 -t 0.01,0.01 --stack-constraint "AU AU" ${NMR_TEST_FASTA})

  add_test(NAME nmr_stack_comma_delimited
    COMMAND ipknot -r 0 -t 0.01,0.01 --stack-constraint "AU,AU" ${NMR_TEST_FASTA})

  add_test(NAME nmr_stack_compact_notation
    COMMAND ipknot -r 0 -t 0.01,0.01 --stack-constraint "(GU)-G-U"
            ${NMR_COMPACT_TEST_FASTA})

  add_test(NAME nmr_stack_compact_invalid_base
    COMMAND ipknot -r 0 -t 0.01,0.01 --stack-constraint "(GU)-A-U" ${NMR_TEST_FASTA})
  set_tests_properties(nmr_stack_compact_invalid_base PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_stack_single_bulge
    COMMAND ipknot -r 0 -t 0.01,0.01 --stack-constraint "GC AU"
            ${NMR_BULGE_TEST_FASTA})
  set_tests_properties(nmr_stack_single_bulge PROPERTIES
    PASS_REGULAR_EXPRESSION "\\(\\.\\(\\.\\.\\.\\.\\.\\.\\)\\)")

  add_test(NAME nmr_stack_single_bulge_right
    COMMAND ipknot -r 0 -t 0.01,0.01 --stack-constraint "GC AU"
            ${NMR_RIGHT_BULGE_TEST_FASTA})
  set_tests_properties(nmr_stack_single_bulge_right PROPERTIES
    PASS_REGULAR_EXPRESSION "\\(\\(\\.\\.\\.\\.\\.\\.\\)\\.\\)")

  add_test(NAME nmr_stack_single_bulge_disabled
    COMMAND ipknot -r 0 -t 0.01,0.01 --stack-constraint "GC AU"
            --without-nmr-bulge ${NMR_BULGE_TEST_FASTA})
  set_tests_properties(nmr_stack_single_bulge_disabled PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_stack_single_bulge_mode_none
    COMMAND ipknot -r 0 -t 0.01,0.01 --stack-constraint "GC AU"
            --nmr-bulge-mode none ${NMR_BULGE_TEST_FASTA})
  set_tests_properties(nmr_stack_single_bulge_mode_none PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_stack_single_bulge_mode_fallback
    COMMAND ipknot -r 0 -t 0.01,0.01 --stack-constraint "GC AU"
            --nmr-bulge-mode fallback --loglevel info ${NMR_BULGE_TEST_FASTA})
  set_tests_properties(nmr_stack_single_bulge_mode_fallback PROPERTIES
    PASS_REGULAR_EXPRESSION "has no direct instance; retained [1-9][0-9]* bulged instance")

  add_test(NAME nmr_stack_bulge_mode_fallback_prefers_direct
    COMMAND ipknot -r 0 -t 0.01,0.01 --stack-constraint "GC AU"
            --nmr-bulge-mode fallback --loglevel info ${NMR_MIXED_BULGE_TEST_FASTA})
  set_tests_properties(nmr_stack_bulge_mode_fallback_prefers_direct PROPERTIES
    PASS_REGULAR_EXPRESSION "uses [1-9][0-9]* direct instance.*discarded [1-9][0-9]* bulged instance")

  add_test(NAME nmr_stack_bulge_mode_invalid
    COMMAND ipknot -r 0 -t 0.01,0.01 --stack-constraint "AU AU"
            --nmr-bulge-mode invalid ${NMR_TEST_FASTA})
  set_tests_properties(nmr_stack_bulge_mode_invalid PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_stack_duplicate_observations_require_distinct_instances
    COMMAND ipknot -r 0 -t 0.01,0.01
            --stack-constraint "GC AU" --stack-constraint "GC AU"
            --nmr-bulge-mode all ${NMR_BULGE_TEST_FASTA})
  set_tests_properties(nmr_stack_duplicate_observations_require_distinct_instances
    PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_stack_coaxial_disabled
    COMMAND ipknot -r 0 -t 0,0 --stack-constraint "UU CC"
            ${NMR_COAXIAL_TEST_FASTA})
  set_tests_properties(nmr_stack_coaxial_disabled PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_stack_coaxial_multibranch
    COMMAND ipknot -r 0 -t 0,0 --stack-constraint "UU CC"
            --coaxial-stacking --loglevel info ${NMR_COAXIAL_TEST_FASTA})
  set_tests_properties(nmr_stack_coaxial_multibranch PROPERTIES
    PASS_REGULAR_EXPRESSION "satisfied by flush coaxial stacking")

  add_test(NAME nmr_stack_coaxial_adjacent_children
    COMMAND ipknot -r 0 -t 0,0 --stack-constraint "UU CC"
            --coaxial-stacking --loglevel info ${NMR_COAXIAL_CHILDREN_FASTA})
  set_tests_properties(nmr_stack_coaxial_adjacent_children PROPERTIES
    PASS_REGULAR_EXPRESSION "satisfied by flush coaxial stacking")

  add_test(NAME nmr_stack_coaxial_last_child
    COMMAND ipknot -r 0 -t 0,0 --stack-constraint "UU CC"
            --coaxial-stacking --loglevel info ${NMR_COAXIAL_LAST_FASTA})
  set_tests_properties(nmr_stack_coaxial_last_child PROPERTIES
    PASS_REGULAR_EXPRESSION "satisfied by flush coaxial stacking")

  add_test(NAME nmr_soft_base_pair_satisfied
    COMMAND ipknot -r 0 -t 0.01,0.01 --base-pairs UU=1
            --without-canonical-neighbor --nmr-soft --loglevel info
            ${NMR_TEST_FASTA})
  set_tests_properties(nmr_soft_base_pair_satisfied PROPERTIES
    PASS_REGULAR_EXPRESSION "Soft NMR base-pair constraint UU=1 satisfied")

  add_test(NAME nmr_soft_base_pair_low_penalty
    COMMAND ipknot -r 0 -t 0.01,0.01 --base-pairs UU=1
            --without-canonical-neighbor --nmr-soft
            --nmr-count-penalty 0.0001 --loglevel info ${NMR_TEST_FASTA})
  set_tests_properties(nmr_soft_base_pair_low_penalty PROPERTIES
    PASS_REGULAR_EXPRESSION "Soft NMR base-pair constraint UU violated")

  add_test(NAME nmr_soft_base_pair_impossible
    COMMAND ipknot -r 0 -t 0.01,0.01 --base-pairs GC=99
            --nmr-soft --loglevel info ${NMR_TEST_FASTA})
  set_tests_properties(nmr_soft_base_pair_impossible PROPERTIES
    PASS_REGULAR_EXPRESSION "Soft NMR base-pair constraint GC violated: expected 99")

  add_test(NAME nmr_soft_stack_impossible
    COMMAND ipknot -r 0 -t 0,0 --stack-constraint "UU CC"
            --nmr-soft --loglevel info ${NMR_COAXIAL_TEST_FASTA})
  set_tests_properties(nmr_soft_stack_impossible PROPERTIES
    PASS_REGULAR_EXPRESSION "Soft NMR stack constraint 1 violated")

  add_test(NAME nmr_soft_invalid_penalty
    COMMAND ipknot -r 0 -t 0.01,0.01 --base-pairs UU=1
            --nmr-soft --nmr-count-penalty -1 ${NMR_TEST_FASTA})
  set_tests_properties(nmr_soft_invalid_penalty PROPERTIES WILL_FAIL TRUE)
endif()
