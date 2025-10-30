cmake_minimum_required(VERSION 3.10)

project(Geometry
    LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# compile options
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-g -O0)
    add_definitions(-DEBUG)
else()
    add_compile_options(-O3 -DNDEBUG)
endif()


# configure eigen 5.0.0
set(EIGEN_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/Eigen")

file(GLOB_RECURSE SOURCES
    "src/*.cpp"
    "src/*.cc"
    "src/*.c")

file(GLOB_RECURSE HEADERS
    "include/*.h"
    "include/*.hpp")

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${EIGEN_INCLUDE_DIR})