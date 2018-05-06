cmake_minimum_required(VERSION 2.8)

project(arpack Fortran)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

file(GLOB ARPACK_SRC "*.f")

# if("${CMAKE_Fortran_COMPILER}" MATCHES "ifort")
#     set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /MP /fpp /names:lowercase /assume:underscore")
# endif()

add_library(${PROJECT_NAME} STATIC ${ARPACK_SRC})