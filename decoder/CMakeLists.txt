cmake_minimum_required(VERSION 3.18)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(IPknotDecoder VERSION 1.1.0 LANGUAGES CXX)
  include(CTest)
endif()

# The decoder can be added directly, without configuring the IPknot CLI/engines.
option(ENABLE_ILP "Discover and link an ILP solver" OFF)
option(ENABLE_CPLEX "Enable CPLEX backend" OFF)
option(ENABLE_GUROBI "Enable Gurobi backend" OFF)
option(ENABLE_SCIP "Enable SCIP backend" OFF)
option(ENABLE_HIGHS "Enable HiGHS backend" OFF)
option(IPKNOT_DECODER_BUILD_EXAMPLE "Build the probability decoding example" OFF)
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")
set(ipknot_core_dir "${CMAKE_CURRENT_SOURCE_DIR}/../src")
include(GNUInstallDirs)
set(IPKNOT_DECODER_WITH_SCIP OFF)

add_library(ipknot_decoder
  src/decoder.cpp src/ilp_model.cpp
  ${ipknot_core_dir}/dual_decomposition.cpp
  ${ipknot_core_dir}/dd_bounds.cpp
  ${ipknot_core_dir}/dd_recovery.cpp
  ${ipknot_core_dir}/dd_joint_bound.cpp
  ${ipknot_core_dir}/dd_exchange.cpp
  ${ipknot_core_dir}/ip.cpp)
add_library(IPknot::decoder ALIAS ipknot_decoder)
set_target_properties(ipknot_decoder PROPERTIES EXPORT_NAME decoder)
target_compile_features(ipknot_decoder PUBLIC cxx_std_17)
target_include_directories(ipknot_decoder
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE ${ipknot_core_dir})

if(ENABLE_ILP)
  if(ENABLE_CPLEX)
    find_package(CPLEX)
    if (CPLEX_FOUND)
      target_compile_definitions(ipknot_decoder PRIVATE WITH_CPLEX)
      target_include_directories(ipknot_decoder PRIVATE ${CPLEX_INCLUDE_DIR})
      target_link_libraries(ipknot_decoder PRIVATE ${CPLEX_LIBRARY})
    endif()
  elseif(ENABLE_GUROBI)
    find_package(GUROBI)
    if (GUROBI_FOUND)
      target_compile_definitions(ipknot_decoder PRIVATE WITH_GUROBI)
      target_include_directories(ipknot_decoder PRIVATE ${GUROBI_INCLUDE_DIR})
      target_link_libraries(ipknot_decoder PRIVATE ${GUROBI_LIBRARY})
    endif()
  elseif(ENABLE_SCIP)
    find_package(SCIP)
    if (SCIP_FOUND)
      set(IPKNOT_DECODER_WITH_SCIP ON)
      message(STATUS "Found SCIP version ${SCIP_VERSION}")
      target_compile_definitions(ipknot_decoder PRIVATE WITH_SCIP)
      target_include_directories(ipknot_decoder PRIVATE ${SCIP_INCLUDE_DIRS})
      target_link_libraries(ipknot_decoder PRIVATE ${SCIP_LIBRARIES})
    endif()
  elseif(ENABLE_HIGHS)
    find_package(HiGHS)
    find_package(ZLIB 1.2.3)
    if (HiGHS_FOUND)
      target_compile_definitions(ipknot_decoder PRIVATE WITH_HIGHS)
      target_include_directories(ipknot_decoder PRIVATE ${HiGHS_INCLUDE_DIR})
      target_link_libraries(ipknot_decoder PRIVATE ${HiGHS_LIBRARY} ${ZLIB_LIBRARY})
    endif()
  else()
    find_package(GLPK)
    if (GLPK_FOUND)
      target_compile_definitions(ipknot_decoder PRIVATE WITH_GLPK)
      target_include_directories(ipknot_decoder PRIVATE ${GLPK_INCLUDE_DIR})
      target_link_libraries(ipknot_decoder PRIVATE ${GLPK_LIBRARY} ${ZLIB_LIBRARY} ${GMP_LIBRARY} ${DL_LIBRARY} ${LTDL_LIBRARY} ${MPFR_LIBRARY})
    endif()
  endif()

endif() # ENABLE_ILP

set_property(TARGET ipknot_decoder PROPERTY IPKNOT_HAS_ILP OFF)
get_target_property(ipknot_backend_definitions ipknot_decoder COMPILE_DEFINITIONS)
foreach(definition IN LISTS ipknot_backend_definitions)
  if(definition MATCHES "^WITH_(GLPK|CPLEX|GUROBI|SCIP|HIGHS)$")
    set_property(TARGET ipknot_decoder PROPERTY IPKNOT_HAS_ILP ON)
  endif()
endforeach()

include(CMakePackageConfigHelpers)
install(TARGETS ipknot_decoder EXPORT IPknotDecoderTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ../COPYING DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/IPknotDecoder)
install(EXPORT IPknotDecoderTargets NAMESPACE IPknot::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/IPknotDecoder)
configure_package_config_file(cmake/IPknotDecoderConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/IPknotDecoderConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/IPknotDecoder)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/IPknotDecoderConfigVersion.cmake
  VERSION 1.1.0 COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/IPknotDecoderConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/IPknotDecoderConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/IPknotDecoder)

if(IPKNOT_DECODER_BUILD_EXAMPLE)
  add_executable(decode_probabilities examples/decode_probabilities.cpp)
  target_link_libraries(decode_probabilities PRIVATE IPknot::decoder)
endif()
if(BUILD_TESTING)
  add_executable(decoder_api_test tests/decoder_api_test.cpp)
  target_link_libraries(decoder_api_test PRIVATE IPknot::decoder)
  add_test(NAME decoder_api COMMAND decoder_api_test)
  add_test(NAME decoder_installed_consumer COMMAND ${CMAKE_COMMAND}
    -DDECODER_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DDECODER_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
    -DTEST_CONFIG=$<CONFIG>
    -DTEST_GENERATOR=${CMAKE_GENERATOR}
    -DTEST_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/installed_consumer.cmake)
  set_tests_properties(decoder_installed_consumer PROPERTIES TIMEOUT 60)
endif()
