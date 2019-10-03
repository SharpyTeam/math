cmake_minimum_required(VERSION 3.13)
project(SandwichMath)

set(CMAKE_CXX_STANDARD 17)

set(SANDWICH_MATH_HEADERS
        src/sw/math/matrix.hpp
        src/sw/math/vector.hpp)

set(SANDWICH_MATH_SOURCES
        src/sw/math/matrix.cpp
        src/sw/math/vector.cpp)

add_library(sandwich_math STATIC ${SANDWICH_MATH_HEADERS} ${SANDWICH_MATH_SOURCES})
add_library(sw::math ALIAS sandwich_math)
target_include_directories(sandwich_math PUBLIC src)