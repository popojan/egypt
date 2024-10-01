cmake_minimum_required(VERSION 3.24)
project(sqrt C)

set(CMAKE_C_STANDARD 11)

set(CMAKE_C_FLAGS "-O3 -fno-fast-math -funroll-loops")

set(CMAKE_LIB)

option(STATIC_LINKING "static linking dependencies" ON)

IF(STATIC_LINKING)
    include_directories(${PROJECT_SOURCE_DIR}/deps/)

    link_directories(
        ${PROJECT_SOURCE_DIR}/deps/gmp-6.2.1/.libs
        ${PROJECT_SOURCE_DIR}/deps/flint-3.1.3/
        ${PROJECT_SOURCE_DIR}/deps/mpfr-4.2.0/src/.libs
    )

    link_libraries(
        m
        libflint.a
        libmpfr.a
        libgmp.a
    )
ELSE()
    find_library(GMP gmp)
    find_library(FLINT flint)
    find_library(MPFR mpfr)
    link_libraries(gmp flint mpfr)
ENDIF()

add_executable(sqrt sqrt.c)
