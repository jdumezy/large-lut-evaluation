cmake_minimum_required(VERSION 3.20)
project(large_lut LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(NATIVEOPT "Enable machine-specific optimizations" ON)

find_package(OpenFHE 1.6.0 CONFIG REQUIRED)
message(STATUS "OpenFHE ${BASE_OPENFHE_VERSION}, native size ${OpenFHE_NATIVE_SIZE}")
if(NOT OpenFHE_NATIVE_SIZE EQUAL 64)
    message(FATAL_ERROR "This implementation requires 64-bit OpenFHE")
endif()
find_package(OpenMP REQUIRED)

add_library(large_lut src/setup.cpp src/digit-decompose.cpp src/trees.cpp src/lut.cpp)
target_include_directories(large_lut PUBLIC include)
target_include_directories(large_lut SYSTEM PUBLIC
    ${OpenFHE_INCLUDE} ${OpenFHE_INCLUDE}/core ${OpenFHE_INCLUDE}/pke
    ${OpenFHE_INCLUDE}/binfhe ${OpenFHE_INCLUDE}/third-party/include)
target_compile_definitions(large_lut PUBLIC MATHBACKEND=${OpenFHE_BACKEND})
target_link_libraries(large_lut PUBLIC OPENFHEpke OPENFHEcore OpenMP::OpenMP_CXX)
target_compile_options(large_lut PRIVATE -Wall -Wextra -Wpedantic)
if(NATIVEOPT)
    target_compile_options(large_lut PUBLIC -march=native)
endif()

add_executable(lut-bench src/benchmark.cpp)
target_link_libraries(lut-bench PRIVATE large_lut)
target_compile_options(lut-bench PRIVATE -Wall -Wextra -Wpedantic)

include(CTest)
if(BUILD_TESTING)
    add_test(NAME cleartext COMMAND lut-bench --input-bits 7 --digit-bits 3 --cleaning 1)
    add_test(NAME plaintext COMMAND lut-bench --input-bits 5 --digit-bits 2 --cleaning 1
             --kind plaintext --slot-varying)
    add_test(NAME ciphertext COMMAND lut-bench --input-bits 5 --digit-bits 2 --cleaning 1
             --kind ciphertext --slot-varying)
    set_tests_properties(cleartext plaintext ciphertext PROPERTIES TIMEOUT 120)
endif()
