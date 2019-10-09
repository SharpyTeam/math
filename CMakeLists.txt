cmake_minimum_required(VERSION 3.13)
project(SandwichMath)

set(CMAKE_CXX_STANDARD 17)

set(SANDWICH_MATH_HEADERS
        src/sw/math/utils.hpp
        src/sw/math/matrix.hpp
        src/sw/math/vector.hpp
        src/sw/math/axis_angle.hpp src/sw/math/simd.hpp)

set(SANDWICH_MATH_SOURCES
        src/sw/math/vector.cpp
        src/sw/math/utils.cpp
        src/sw/math/matrix.cpp
        src/sw/math/axis_angle.cpp)

option(SANDWICH_MATH_INLINE "Compile it header-only inline" OFF)

add_library(sandwich_math STATIC ${SANDWICH_MATH_HEADERS} ${SANDWICH_MATH_SOURCES})
add_library(sw::math ALIAS sandwich_math)
target_include_directories(sandwich_math PUBLIC src)
if (SANDWICH_MATH_INLINE)
    target_compile_definitions(sandwich_math PUBLIC SANDWICH_MATH_INLINE=1)
elseif ()
    target_compile_definitions(sandwich_math PUBLIC SANDWICH_MATH_INLINE=0)
endif ()