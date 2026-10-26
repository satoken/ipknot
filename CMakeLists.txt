
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
  src/ipknot.h src/ipknot.cpp
  src/dd_constraint_builder.h src/dd_constraint_builder.cpp
  src/dd_constrained.h src/dd_constrained.cpp
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


endif()

if(BUILD_TESTING)
  add_executable(dd_constrained_test tests/dd_constrained_test.cpp src/dd_constrained.cpp)
  add_executable(dd_performance_test tests/dd_performance_test.cpp src/dd_constrained.cpp)
  target_include_directories(dd_performance_test PRIVATE src)
  target_link_libraries(dd_performance_test PRIVATE IPknot::decoder)
  add_test(NAME dd_sparse_parallel COMMAND dd_performance_test)
  set_tests_properties(dd_sparse_parallel PROPERTIES TIMEOUT 60)
  target_include_directories(dd_constrained_test PRIVATE src)
  target_compile_features(dd_constrained_test PRIVATE cxx_std_17)
  add_test(NAME dd_constrained_exhaustive COMMAND dd_constrained_test)
  set_tests_properties(dd_constrained_exhaustive PROPERTIES TIMEOUT 60)
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
