# Exercise only the installed public headers, exported target, and its dependencies.
function(run_checked)
  execute_process(COMMAND ${ARGV} RESULT_VARIABLE result
    OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Command failed: ${ARGV}\n${output}\n${error}")
  endif()
endfunction()

set(consumer_dir "${DECODER_BINARY_DIR}/consumer-test")
set(install_dir "${consumer_dir}/install")
run_checked("${CMAKE_COMMAND}" --install "${DECODER_BINARY_DIR}"
  --prefix "${install_dir}" --config "${TEST_CONFIG}")
file(MAKE_DIRECTORY "${consumer_dir}/src")
configure_file("${DECODER_SOURCE_DIR}/examples/decode_probabilities.cpp"
  "${consumer_dir}/src/main.cpp" COPYONLY)
file(WRITE "${consumer_dir}/src/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.18)
project(DecoderConsumer LANGUAGES CXX)
find_package(IPknotDecoder CONFIG REQUIRED)
add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE IPknot::decoder)
enable_testing()
add_test(NAME run_consumer COMMAND consumer)
]=])
run_checked("${CMAKE_COMMAND}" -S "${consumer_dir}/src" -B "${consumer_dir}/build"
  -G "${TEST_GENERATOR}" "-DCMAKE_PREFIX_PATH=${install_dir}"
  "-DCMAKE_CXX_COMPILER=${TEST_CXX_COMPILER}" "-DCMAKE_BUILD_TYPE=${TEST_CONFIG}")
run_checked("${CMAKE_COMMAND}" --build "${consumer_dir}/build" --config "${TEST_CONFIG}")
run_checked("${CMAKE_CTEST_COMMAND}" --test-dir "${consumer_dir}/build"
  --build-config "${TEST_CONFIG}" --output-on-failure)
