
cmake_minimum_required(VERSION 3.18)
#cmake_policy(SET CMP0074 NEW)
#cmake_policy(SET CMP0148 NEW)

project(IPknot CXX C)
option(CXX "enable C++ compilation" ON)
option(WITH_MXFOLD2 "Enable MXfold2 integration" OFF)
option(ENABLE_ILP "Discover and link an ILP solver (OFF builds only DD)" ON)
option(ENABLE_CPLEX "Enable CPLEX backend" OFF)
option(ENABLE_GUROBI "Enable Gurobi backend" OFF)
option(ENABLE_SCIP "Enable SCIP backend" OFF)
option(ENABLE_HIGHS "Enable HiGHS backend" OFF)
option(STATIC_BUILD "Build statically for Linux" OFF)
if (CXX)
  enable_language(CXX)
endif()

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(CTest)
add_subdirectory(decoder)
option(IPKNOT_BUILD_CLI "Build the IPknot CLI and probability engines" ON)
if(NOT IPKNOT_BUILD_CLI)
  return()
endif()

find_package(PkgConfig REQUIRED)

if(STATIC_BUILD) # static build for linux system
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

add_executable(ipknot
  src/main.cpp
  src/ipknot.h src/ipknot.cpp src/nmr_pair_index.h src/nmr_sequence_index.h src/nmr_candidate_filter.h
  src/dd_constraint_builder.h src/dd_constraint_builder.cpp
  src/dd_constrained.h src/dd_constrained.cpp
  src/dd_noe.h src/dd_noe.cpp
  src/fa.cpp
  src/aln.cpp
  src/fold.cpp
  src/bpseq.h src/bpseq.cpp
  src/boltzmann_param.c
  src/nupack/nupack.cpp
  src/nupack/linear_nupack.cpp
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
target_include_directories(ipknot PRIVATE src decoder/src)
target_link_libraries(ipknot PRIVATE IPknot::decoder)

pkg_check_modules(VIENNARNA RNAlib2)
if (VIENNARNA_FOUND)
  target_include_directories(ipknot PRIVATE ${VIENNARNA_INCLUDE_DIRS})
  target_compile_definitions(ipknot PRIVATE HAVE_VIENNA18 HAVE_VIENNA20)
  target_link_libraries(ipknot PRIVATE ${VIENNARNA_LDFLAGS} ${CMAKE_DL_LIBS})
endif()

# Solver availability is exposed by the backend target; application code uses IP.
get_target_property(IPKNOT_HAS_ILP ipknot_decoder IPKNOT_HAS_ILP)

if (WITH_MXFOLD2)
  find_package(Python REQUIRED COMPONENTS Interpreter Development)
  find_package(pybind11 REQUIRED)
  target_sources(ipknot PRIVATE src/mxfold2.cpp)
  target_link_libraries(ipknot PRIVATE pybind11::embed)
  target_compile_definitions(ipknot PRIVATE WITH_MXFOLD2)
endif()

install(TARGETS ipknot RUNTIME DESTINATION bin)
if(BUILD_TESTING)
  add_executable(nupack_test tests/nupack_test.cpp src/nupack/nupack.cpp src/nupack/linear_nupack.cpp)
  target_include_directories(nupack_test PRIVATE src)
  target_compile_features(nupack_test PRIVATE cxx_std_17)
  if(MSVC)
    target_compile_options(nupack_test PRIVATE /UNDEBUG)
  else()
    target_compile_options(nupack_test PRIVATE -UNDEBUG)
  endif()
  add_test(NAME nupack_exact_and_linear COMMAND nupack_test)
  set_tests_properties(nupack_exact_and_linear PROPERTIES TIMEOUT 120)
  if(UNIX)
    add_executable(nupack_benchmark experiments/nupack/driver.cpp src/nupack/nupack.cpp src/nupack/linear_nupack.cpp)
    target_include_directories(nupack_benchmark PRIVATE src/nupack)
    target_compile_features(nupack_benchmark PRIVATE cxx_std_17)
  endif()
  add_executable(nmr_pair_index_test tests/nmr_pair_index_test.cpp)
  add_executable(nmr_sequence_index_test tests/nmr_sequence_index_test.cpp)
  add_executable(nmr_candidate_filter_test tests/nmr_candidate_filter_test.cpp)
  target_include_directories(nmr_pair_index_test PRIVATE src)
  target_include_directories(nmr_sequence_index_test PRIVATE src)
  target_include_directories(nmr_candidate_filter_test PRIVATE src)
  target_compile_features(nmr_pair_index_test PRIVATE cxx_std_17)
  target_compile_features(nmr_sequence_index_test PRIVATE cxx_std_17)
  target_compile_features(nmr_candidate_filter_test PRIVATE cxx_std_17)
  # Assertions are the test's checks, including in Release builds.
  if(MSVC)
    target_compile_options(nmr_pair_index_test PRIVATE /UNDEBUG)
    target_compile_options(nmr_sequence_index_test PRIVATE /UNDEBUG)
    target_compile_options(nmr_candidate_filter_test PRIVATE /UNDEBUG)
  else()
    target_compile_options(nmr_pair_index_test PRIVATE -UNDEBUG)
    target_compile_options(nmr_sequence_index_test PRIVATE -UNDEBUG)
    target_compile_options(nmr_candidate_filter_test PRIVATE -UNDEBUG)
  endif()
  add_test(NAME nmr_sparse_pair_index COMMAND nmr_pair_index_test)
  add_test(NAME nmr_sparse_sequence_index COMMAND nmr_sequence_index_test)
  add_test(NAME nmr_exact_candidate_filter COMMAND nmr_candidate_filter_test)
  if(IPKNOT_HAS_ILP)
    add_executable(ip_matrix_transfer_test tests/ip_matrix_transfer_test.cpp)
    target_include_directories(ip_matrix_transfer_test PRIVATE src)
    target_compile_features(ip_matrix_transfer_test PRIVATE cxx_std_17)
    target_link_libraries(ip_matrix_transfer_test PRIVATE IPknot::decoder)
    if(MSVC)
      target_compile_options(ip_matrix_transfer_test PRIVATE /UNDEBUG)
    else()
      target_compile_options(ip_matrix_transfer_test PRIVATE -UNDEBUG)
    endif()
    add_test(NAME ip_exact_matrix_transfer COMMAND ip_matrix_transfer_test)
    add_executable(ip_model_solver_test tests/ip_model_solver_test.cpp)
    target_include_directories(ip_model_solver_test PRIVATE src)
    target_link_libraries(ip_model_solver_test PRIVATE IPknot::decoder)
    add_test(NAME ip_model_solver_reoptimization COMMAND ip_model_solver_test)
    set_tests_properties(ip_model_solver_reoptimization PROPERTIES TIMEOUT 60)
  endif()
endif()
if(BUILD_TESTING AND IPKNOT_HAS_ILP)
  find_package(Python3 QUIET COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_test(NAME bpp_roundtrip
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/bpp_roundtrip_test.py
              $<TARGET_FILE:ipknot>)
    set_tests_properties(bpp_roundtrip PROPERTIES TIMEOUT 60)
  endif()

  set(NMR_TEST_FASTA "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nmr_constraints.fa")
  set(NMR_COMPACT_TEST_FASTA "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nmr_compact_stack.fa")
  set(NMR_BULGE_TEST_FASTA "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nmr_single_bulge.fa")
  set(NMR_RIGHT_BULGE_TEST_FASTA "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nmr_single_bulge_right.fa")
  set(NMR_MIXED_BULGE_TEST_FASTA "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nmr_bulge_mixed.fa")
  set(NMR_COAXIAL_TEST_FASTA "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nmr_coaxial_multibranch.fa")
  set(NMR_COAXIAL_CHILDREN_FASTA "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nmr_coaxial_children.fa")
  set(NMR_COAXIAL_LAST_FASTA "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nmr_coaxial_last_child.fa")

  add_test(NAME nmr_threshold_penalty_scale_zero
    COMMAND ipknot --decoder ilp -r 0 -t auto,auto --nmr-soft --base-pairs GC=99
            --nmr-threshold-penalty-scale 0 --loglevel info ${NMR_TEST_FASTA})
  set_tests_properties(nmr_threshold_penalty_scale_zero PROPERTIES
    PASS_REGULAR_EXPRESSION "NMR automatic-threshold penalty scale: 0")
  add_test(NAME nmr_threshold_penalty_scale_dd_zero
    COMMAND ipknot --decoder dd -r 0 -t auto,auto --nmr-soft --base-pairs GC=99
            --nmr-threshold-penalty-scale 0 --loglevel info ${NMR_TEST_FASTA})
  set_tests_properties(nmr_threshold_penalty_scale_dd_zero PROPERTIES
    PASS_REGULAR_EXPRESSION "NMR automatic-threshold penalty scale: 0" TIMEOUT 60)
  foreach(value -1 nan inf)
    add_test(NAME nmr_threshold_penalty_scale_invalid_${value}
      COMMAND ipknot --decoder ilp --nmr-threshold-penalty-scale ${value} ${NMR_TEST_FASTA})
    set_tests_properties(nmr_threshold_penalty_scale_invalid_${value}
      PROPERTIES WILL_FAIL TRUE)
  endforeach()

  add_test(NAME nmr_base_pair_noncanonical
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --base-pairs UU=1
            --without-canonical-neighbor ${NMR_TEST_FASTA})

  add_test(NAME nmr_base_pair_impossible
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --base-pairs GC=99 ${NMR_TEST_FASTA})
  set_tests_properties(nmr_base_pair_impossible PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_base_pair_negative
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --base-pairs UU=-1 ${NMR_TEST_FASTA})
  set_tests_properties(nmr_base_pair_negative PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_base_pair_duplicate
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --base-pairs GC=1,CG=2 ${NMR_TEST_FASTA})
  set_tests_properties(nmr_base_pair_duplicate PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_stack_space_delimited
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --stack-constraint "AU AU" ${NMR_TEST_FASTA})

  add_test(NAME nmr_stack_comma_delimited
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --stack-constraint "AU,AU" ${NMR_TEST_FASTA})

  add_test(NAME nmr_stack_compact_notation
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --stack-constraint "(GU)-G-U"
            ${NMR_COMPACT_TEST_FASTA})

  add_test(NAME nmr_stack_compact_invalid_base
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --stack-constraint "(GU)-A-U" ${NMR_TEST_FASTA})
  set_tests_properties(nmr_stack_compact_invalid_base PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_stack_single_bulge
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --stack-constraint "GC AU"
            ${NMR_BULGE_TEST_FASTA})
  set_tests_properties(nmr_stack_single_bulge PROPERTIES
    PASS_REGULAR_EXPRESSION "\\(\\.\\(\\.\\.\\.\\.\\.\\.\\)\\)")

  add_test(NAME nmr_stack_single_bulge_right
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --stack-constraint "GC AU"
            ${NMR_RIGHT_BULGE_TEST_FASTA})
  set_tests_properties(nmr_stack_single_bulge_right PROPERTIES
    PASS_REGULAR_EXPRESSION "\\(\\(\\.\\.\\.\\.\\.\\.\\)\\.\\)")

  add_test(NAME nmr_stack_single_bulge_disabled
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --stack-constraint "GC AU"
            --without-nmr-bulge ${NMR_BULGE_TEST_FASTA})
  set_tests_properties(nmr_stack_single_bulge_disabled PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_stack_single_bulge_mode_none
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --stack-constraint "GC AU"
            --nmr-bulge-mode none ${NMR_BULGE_TEST_FASTA})
  set_tests_properties(nmr_stack_single_bulge_mode_none PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_stack_mode_none_keeps_adjacent_support
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --stack-constraint "AU AU"
            --nmr-bulge-mode none --loglevel debug ${NMR_TEST_FASTA})
  set_tests_properties(nmr_stack_mode_none_keeps_adjacent_support PROPERTIES
    PASS_REGULAR_EXPRESSION "Stacking neighbor distance: 1")

  add_test(NAME nmr_stack_mode_all_allows_bulged_support
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --stack-constraint "GC AU"
            --nmr-bulge-mode all --loglevel debug ${NMR_BULGE_TEST_FASTA})
  set_tests_properties(nmr_stack_mode_all_allows_bulged_support PROPERTIES
    PASS_REGULAR_EXPRESSION "Stacking neighbor distance: 2")

  add_test(NAME nmr_stack_single_bulge_mode_fallback
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --stack-constraint "GC AU"
            --nmr-bulge-mode fallback --loglevel info ${NMR_BULGE_TEST_FASTA})
  set_tests_properties(nmr_stack_single_bulge_mode_fallback PROPERTIES
    PASS_REGULAR_EXPRESSION "has no direct instance; retained [1-9][0-9]* bulged instance")

  add_test(NAME nmr_stack_bulge_mode_fallback_prefers_direct
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --stack-constraint "GC AU"
            --nmr-bulge-mode fallback --loglevel info ${NMR_MIXED_BULGE_TEST_FASTA})
  set_tests_properties(nmr_stack_bulge_mode_fallback_prefers_direct PROPERTIES
    PASS_REGULAR_EXPRESSION "uses [1-9][0-9]* direct instance.*skipped bulged enumeration")

  add_test(NAME nmr_stack_bulge_mode_invalid
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --stack-constraint "AU AU"
            --nmr-bulge-mode invalid ${NMR_TEST_FASTA})
  set_tests_properties(nmr_stack_bulge_mode_invalid PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_stack_duplicate_observations_require_distinct_instances
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01
            --stack-constraint "GC AU" --stack-constraint "GC AU"
            --nmr-bulge-mode all ${NMR_BULGE_TEST_FASTA})
  set_tests_properties(nmr_stack_duplicate_observations_require_distinct_instances
    PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_stack_coaxial_disabled
    COMMAND ipknot --decoder ilp -r 0 -t 0,0 --stack-constraint "UU CC"
            ${NMR_COAXIAL_TEST_FASTA})
  set_tests_properties(nmr_stack_coaxial_disabled PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_stack_coaxial_multibranch
    COMMAND ipknot --decoder ilp -r 0 -t 0,0 --stack-constraint "UU CC"
            --coaxial-stacking --loglevel info ${NMR_COAXIAL_TEST_FASTA})
  set_tests_properties(nmr_stack_coaxial_multibranch PROPERTIES
    PASS_REGULAR_EXPRESSION "satisfied by flush coaxial stacking")

  add_test(NAME nmr_stack_coaxial_adjacent_children
    COMMAND ipknot --decoder ilp -r 0 -t 0,0 --stack-constraint "UU CC"
            --coaxial-stacking --loglevel info ${NMR_COAXIAL_CHILDREN_FASTA})
  set_tests_properties(nmr_stack_coaxial_adjacent_children PROPERTIES
    PASS_REGULAR_EXPRESSION "satisfied by flush coaxial stacking")

  add_test(NAME nmr_stack_coaxial_last_child
    COMMAND ipknot --decoder ilp -r 0 -t 0,0 --stack-constraint "UU CC"
            --coaxial-stacking --loglevel info ${NMR_COAXIAL_LAST_FASTA})
  set_tests_properties(nmr_stack_coaxial_last_child PROPERTIES
    PASS_REGULAR_EXPRESSION "satisfied by flush coaxial stacking")

  add_test(NAME nmr_coaxial_no_adjacent_bulge_requires_coaxial
    COMMAND ipknot --decoder ilp --nmr-coaxial-no-adjacent-bulge ${NMR_TEST_FASTA})
  set_tests_properties(nmr_coaxial_no_adjacent_bulge_requires_coaxial PROPERTIES WILL_FAIL TRUE)

  set(NMR_COAXIAL_ADJACENT_BULGE_CONSTRAINT
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nmr_coaxial_adjacent_bulge.constraint")
  set(NMR_COAXIAL_ADJACENT_BULGE_FASTA
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nmr_coaxial_bulge_multibranch.fa")
  add_test(NAME nmr_coaxial_adjacent_bulge_default_allowed
    COMMAND ipknot --decoder ilp -r 0 -t 0,0 --nmr-bulge-mode all --coaxial-stacking
            --nmr-allow-shared-stack-pairs
            --stack-constraint "UU CC" --stack-constraint "AC UU"
            --stack-constraint "GC AU" --loglevel info
            -c ${NMR_COAXIAL_ADJACENT_BULGE_CONSTRAINT} ${NMR_COAXIAL_ADJACENT_BULGE_FASTA})
  set_tests_properties(nmr_coaxial_adjacent_bulge_default_allowed PROPERTIES
    PASS_REGULAR_EXPRESSION "All stack constraints are satisfied" TIMEOUT 30)
  add_test(NAME nmr_coaxial_adjacent_bulge_prior_hard
    COMMAND ipknot --decoder ilp -r 0 -t 0,0 --nmr-bulge-mode all --coaxial-stacking
            --nmr-allow-shared-stack-pairs
            --nmr-coaxial-no-adjacent-bulge
            --stack-constraint "UU CC" --stack-constraint "AC UU"
            --stack-constraint "GC AU"
            -c ${NMR_COAXIAL_ADJACENT_BULGE_CONSTRAINT} ${NMR_COAXIAL_ADJACENT_BULGE_FASTA})
  set_tests_properties(nmr_coaxial_adjacent_bulge_prior_hard PROPERTIES WILL_FAIL TRUE TIMEOUT 30)
  add_test(NAME nmr_coaxial_adjacent_bulge_prior_soft
    COMMAND ipknot --decoder ilp -r 0 -t 0,0 --nmr-bulge-mode all --coaxial-stacking --nmr-soft
            --nmr-allow-shared-stack-pairs
            --nmr-coaxial-no-adjacent-bulge
            --stack-constraint "UU CC" --stack-constraint "AC UU"
            --stack-constraint "GC AU" --loglevel info
            -c ${NMR_COAXIAL_ADJACENT_BULGE_CONSTRAINT} ${NMR_COAXIAL_ADJACENT_BULGE_FASTA})
  set_tests_properties(nmr_coaxial_adjacent_bulge_prior_soft PROPERTIES
    PASS_REGULAR_EXPRESSION "Soft NMR stack constraint 1 violated" TIMEOUT 30)

  add_test(NAME nmr_soft_base_pair_satisfied
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --base-pairs UU=1
            --without-canonical-neighbor --nmr-soft --loglevel info
            ${NMR_TEST_FASTA})
  set_tests_properties(nmr_soft_base_pair_satisfied PROPERTIES
    PASS_REGULAR_EXPRESSION "Soft NMR base-pair constraint UU=1 satisfied")

  add_test(NAME nmr_soft_base_pair_low_penalty
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --base-pairs UU=1
            --without-canonical-neighbor --nmr-soft
            --nmr-count-penalty 0.0001 --loglevel info ${NMR_TEST_FASTA})
  set_tests_properties(nmr_soft_base_pair_low_penalty PROPERTIES
    PASS_REGULAR_EXPRESSION "Soft NMR base-pair constraint UU violated")

  add_test(NAME nmr_soft_base_pair_impossible
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --base-pairs GC=99
            --nmr-soft --loglevel info ${NMR_TEST_FASTA})
  set_tests_properties(nmr_soft_base_pair_impossible PROPERTIES
    PASS_REGULAR_EXPRESSION "Soft NMR base-pair constraint GC violated: expected 99")

  add_test(NAME nmr_soft_stack_impossible
    COMMAND ipknot --decoder ilp -r 0 -t 0,0 --stack-constraint "UU CC"
            --nmr-soft --loglevel info ${NMR_COAXIAL_TEST_FASTA})
  set_tests_properties(nmr_soft_stack_impossible PROPERTIES
    PASS_REGULAR_EXPRESSION "Soft NMR stack constraint 1 violated")

  add_test(NAME nmr_soft_invalid_penalty
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --base-pairs UU=1
            --nmr-soft --nmr-count-penalty -1 ${NMR_TEST_FASTA})
  set_tests_properties(nmr_soft_invalid_penalty PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_soft_base_pair_lower_bound
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --base-pairs AU=1
            --nmr-soft --nmr-count-mode lower-bound --loglevel info
            ${NMR_TEST_FASTA})
  set_tests_properties(nmr_soft_base_pair_lower_bound PROPERTIES
    PASS_REGULAR_EXPRESSION "Soft NMR base-pair constraint AU=1 satisfied")

  add_test(NAME nmr_hard_base_pair_lower_bound_allows_excess
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --base-pairs AU=1
            --nmr-count-mode lower-bound --loglevel info
            ${NMR_TEST_FASTA})
  set_tests_properties(nmr_hard_base_pair_lower_bound_allows_excess PROPERTIES
    PASS_REGULAR_EXPRESSION "Base pair constraints are satisfied"
    FAIL_REGULAR_EXPRESSION "Base pair constraints are NOT satisfied")

  add_test(NAME nmr_soft_base_pair_lower_bound_allows_excess
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --base-pairs AU=0
            --nmr-soft --nmr-count-mode lower-bound --loglevel info
            ${NMR_TEST_FASTA})
  set_tests_properties(nmr_soft_base_pair_lower_bound_allows_excess PROPERTIES
    PASS_REGULAR_EXPRESSION "Soft NMR base-pair constraint AU=0 satisfied"
    FAIL_REGULAR_EXPRESSION "Base pair constraints are NOT satisfied")

  add_test(NAME nmr_count_mode_invalid
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01 --base-pairs AU=1
            --nmr-soft --nmr-count-mode invalid ${NMR_TEST_FASTA})
  set_tests_properties(nmr_count_mode_invalid PROPERTIES WILL_FAIL TRUE)

  add_test(NAME nmr_stack_duplicate_observations_allow_shared_pairs
    COMMAND ipknot --decoder ilp -r 0 -t 0.01,0.01
            --stack-constraint "GC AU" --stack-constraint "GC AU"
            --nmr-bulge-mode all --nmr-allow-shared-stack-pairs --loglevel info
            ${NMR_BULGE_TEST_FASTA})
  set_tests_properties(nmr_stack_duplicate_observations_allow_shared_pairs
    PROPERTIES PASS_REGULAR_EXPRESSION "All stack constraints are satisfied")

  # Exercise all flush-coaxial topologies together with a separate one-base
  # bulge. Fix the structure so these are fast, solver-independent checks of
  # the witness formulation, including soft mode and shared-pair assignment.
  foreach(kind multibranch children last_child)
    set(fixture "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nmr_coaxial_bulge_${kind}")
    foreach(mode hard soft)
      foreach(sharing distinct shared)
        set(options)
        if(mode STREQUAL "soft")
          list(APPEND options --nmr-soft)
        endif()
        if(sharing STREQUAL "shared")
          list(APPEND options --nmr-allow-shared-stack-pairs)
        endif()
        set(test_name "nmr_coaxial_bulge_${kind}_${mode}_${sharing}")
        add_test(NAME ${test_name}
          COMMAND ipknot --decoder ilp -r 0 -t 0,0 --nmr-bulge-mode all --coaxial-stacking
                  --stack-constraint "UU CC" --stack-constraint "GC AU"
                  ${options} --loglevel info -c "${fixture}.constraint"
                  "${fixture}.fa")
        set_tests_properties(${test_name} PROPERTIES
          PASS_REGULAR_EXPRESSION "All stack constraints are satisfied"
          FAIL_REGULAR_EXPRESSION "Some stack constraints are NOT satisfied;violated"
          TIMEOUT 30)
      endforeach()
    endforeach()

    add_test(NAME nmr_coaxial_straight_termini_${kind}
      COMMAND ipknot --decoder ilp -r 0 -t 0,0 --nmr-bulge-mode all --coaxial-stacking
              --nmr-coaxial-no-adjacent-bulge
              --stack-constraint "UU CC" --stack-constraint "GC AU"
              --loglevel info -c "${fixture}.constraint" "${fixture}.fa")
    set_tests_properties(nmr_coaxial_straight_termini_${kind} PROPERTIES
      PASS_REGULAR_EXPRESSION "All stack constraints are satisfied"
      FAIL_REGULAR_EXPRESSION "Some stack constraints are NOT satisfied;violated"
      TIMEOUT 30)

    # Sharing base pairs must not permit two observations to use the same
    # coaxial helix face. In soft mode one of them must incur its penalty.
    add_test(NAME nmr_coaxial_face_capacity_${kind}_hard
      COMMAND ipknot --decoder ilp -r 0 -t 0,0 --coaxial-stacking
              --stack-constraint "UU CC" --stack-constraint "UU CC"
              --stack-constraint "GC AU" --nmr-bulge-mode all
              --nmr-allow-shared-stack-pairs -c "${fixture}.constraint"
              "${fixture}.fa")
    set_tests_properties(nmr_coaxial_face_capacity_${kind}_hard
      PROPERTIES WILL_FAIL TRUE TIMEOUT 30)
    add_test(NAME nmr_coaxial_face_capacity_${kind}_soft
      COMMAND ipknot --decoder ilp -r 0 -t 0,0 --coaxial-stacking --nmr-soft
              --stack-constraint "UU CC" --stack-constraint "UU CC"
              --stack-constraint "GC AU" --nmr-bulge-mode all
              --nmr-allow-shared-stack-pairs --loglevel info
              -c "${fixture}.constraint" "${fixture}.fa")
    set_tests_properties(nmr_coaxial_face_capacity_${kind}_soft PROPERTIES
      PASS_REGULAR_EXPRESSION "Soft NMR stack constraint [12] violated"
      TIMEOUT 30)
  endforeach()

endif()

if(BUILD_TESTING)
  add_executable(dd_constrained_test tests/dd_constrained_test.cpp src/dd_constrained.cpp src/dd_noe.cpp)
  target_include_directories(dd_constrained_test PRIVATE src)
  target_compile_features(dd_constrained_test PRIVATE cxx_std_17)
  add_test(NAME dd_constrained_exhaustive COMMAND dd_constrained_test)
  set_tests_properties(dd_constrained_exhaustive PROPERTIES TIMEOUT 60)
  if(IPKNOT_HAS_ILP)
    add_executable(dd_noe_test tests/dd_noe_test.cpp src/dd_noe.cpp src/dd_constrained.cpp)
    target_compile_features(dd_noe_test PRIVATE cxx_std_17)
    target_include_directories(dd_noe_test PRIVATE src)
    target_link_libraries(dd_noe_test PRIVATE IPknot::decoder)
    add_test(NAME dd_noe_factor_exhaustive COMMAND dd_noe_test)
    set_tests_properties(dd_noe_factor_exhaustive PROPERTIES TIMEOUT 60)
  endif()
  add_executable(dual_decomposition_test tests/dual_decomposition_test.cpp)
  add_executable(dd_diagnostic_driver experiments/dual-decomposition/convergence/diagnostic_driver.cpp)
  add_executable(dd_bounds_driver experiments/dual-decomposition/bounds/diagnostic_driver.cpp)
  add_executable(dd_offline_recovery experiments/dual-decomposition/bounds/offline_recovery.cpp)
  target_include_directories(dd_offline_recovery PRIVATE src)
  target_compile_features(dd_offline_recovery PRIVATE cxx_std_17)
  add_executable(dd_integer_driver experiments/dual-decomposition/integer-gap/diagnostic_driver.cpp)
  target_include_directories(dd_integer_driver PRIVATE src)
  target_compile_features(dd_integer_driver PRIVATE cxx_std_17)
  add_executable(dd_joint_bound_test tests/dd_joint_bound_test.cpp)
  target_include_directories(dd_joint_bound_test PRIVATE src)
  target_compile_features(dd_joint_bound_test PRIVATE cxx_std_17)
  add_test(NAME dd_joint_bound_exhaustive COMMAND dd_joint_bound_test)
  add_executable(dd_exchange_test tests/dd_exchange_test.cpp)
  target_include_directories(dd_exchange_test PRIVATE src)
  target_compile_features(dd_exchange_test PRIVATE cxx_std_17)
  add_test(NAME dd_exchange_exhaustive COMMAND dd_exchange_test)
  add_executable(dd_rational_bound_driver experiments/dual-decomposition/bounds/rational_driver.cpp)
  target_include_directories(dd_rational_bound_driver PRIVATE src)
  target_compile_features(dd_rational_bound_driver PRIVATE cxx_std_17)
  target_include_directories(dd_bounds_driver PRIVATE src)
  target_compile_features(dd_bounds_driver PRIVATE cxx_std_17)
  target_include_directories(dd_diagnostic_driver PRIVATE src)
  target_compile_features(dd_diagnostic_driver PRIVATE cxx_std_17)
  add_executable(dd_bounds_test tests/dd_bounds_test.cpp)
  target_include_directories(dd_bounds_test PRIVATE src)
  target_compile_features(dd_bounds_test PRIVATE cxx_std_17)
  add_test(NAME dd_bounds_exhaustive COMMAND dd_bounds_test)
  add_executable(dd_recovery_test tests/dd_recovery_test.cpp)
  target_include_directories(dd_recovery_test PRIVATE src)
  target_compile_features(dd_recovery_test PRIVATE cxx_std_17)
  add_test(NAME dd_recovery_exhaustive COMMAND dd_recovery_test)
  target_include_directories(dual_decomposition_test PRIVATE src)
  target_compile_features(dual_decomposition_test PRIVATE cxx_std_17)
  add_test(NAME dual_decomposition_exhaustive COMMAND dual_decomposition_test)
  set_tests_properties(dual_decomposition_exhaustive PROPERTIES TIMEOUT 60)
  foreach(decoder_test_target IN ITEMS
      dd_constrained_test dual_decomposition_test dd_diagnostic_driver dd_bounds_driver
      dd_offline_recovery dd_integer_driver dd_joint_bound_test dd_exchange_test
      dd_rational_bound_driver dd_bounds_test dd_recovery_test)
    target_link_libraries(${decoder_test_target} PRIVATE IPknot::decoder)
  endforeach()
  find_package(Python3 QUIET COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_test(NAME linear_nupack_cli COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/linear_nupack_cli_test.py $<TARGET_FILE:ipknot>)
    set_tests_properties(linear_nupack_cli PROPERTIES TIMEOUT 120)
    add_test(NAME removed_options_cli COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/removed_options_cli_test.py $<TARGET_FILE:ipknot>)
    set_tests_properties(removed_options_cli PROPERTIES TIMEOUT 60)
    add_test(NAME dual_decomposition_cli COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/dual_decomposition_cli_test.py $<TARGET_FILE:ipknot>)
    set_tests_properties(dual_decomposition_cli PROPERTIES TIMEOUT 60)
    add_test(NAME dd_constraints_cli COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/dd_constraints_cli_test.py $<TARGET_FILE:ipknot>)
    set_tests_properties(dd_constraints_cli PROPERTIES TIMEOUT 120)
    add_test(NAME dd_noe_cli COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/dd_noe_cli_test.py $<TARGET_FILE:ipknot> ${IPKNOT_HAS_ILP})
    set_tests_properties(dd_noe_cli PROPERTIES TIMEOUT 120)
    add_test(NAME dd_constraint_scaling COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/dd_constraint_scaling_test.py $<TARGET_FILE:ipknot>)
    set_tests_properties(dd_constraint_scaling PROPERTIES TIMEOUT 120)
  endif()
endif()

if(BUILD_TESTING)
  add_executable(sparse_posterior_test tests/sparse_posterior_test.cpp)
  target_include_directories(sparse_posterior_test PRIVATE src)
  target_compile_features(sparse_posterior_test PRIVATE cxx_std_17)
  add_test(NAME sparse_posterior_bitwise COMMAND sparse_posterior_test)
endif()
