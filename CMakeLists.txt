cmake_minimum_required(VERSION 3.14)
project(zombit_test)

# GoogleTest requires at least C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS_RELEASE "-O3")
set(CMAKE -DCMAKE_BUILD_TYPE=Debug)


message("${CMAKE_CXX_FLAGS}")

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


find_path (LIB_INCLUDES ClassLibrary.h /path/to/header/files)
include_directories(${CMAKE_SOURCE_DIR}/external/cpp_libs/include)
include_directories(${CMAKE_SOURCE_DIR}/src/include) # zombit
find_library(SDSL-LIB libsdsl "${CMAKE_SOURCE_DIR}/external/cpp_libs/lib")


enable_testing()

add_executable(
  ${PROJECT_NAME}
  ./test/src/zombit_vector_test.cpp
)

target_link_libraries(
  ${PROJECT_NAME}
  gtest_main
  ${CMAKE_SOURCE_DIR}/external/cpp_libs/lib/libsdsl.a
  ${CMAKE_SOURCE_DIR}/external/cpp_libs/lib/libdivsufsort.a
  ${CMAKE_SOURCE_DIR}/external/cpp_libs/lib/libdivsufsort64.a
)

include(GoogleTest)
gtest_discover_tests(${PROJECT_NAME})
