cmake_minimum_required(VERSION 3.0)

project(suanPan C CXX Fortran)

set(ROOT .)

include(${ROOT}/.config.cmake)

link_directories(${PROJECT_BINARY_DIR}/Libs)

add_executable(${PROJECT_NAME}
        "suanPan.h"
        "suanPan.cpp"
        "Constraint/CMakeLists.txt"
        "Domain/CMakeLists.txt"
        "Load/CMakeLists.txt"
        "Recorder/CMakeLists.txt"
        "Step/CMakeLists.txt"
        "Toolbox/CMakeLists.txt")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU") # GNU GCC COMPILER
    target_sources(${PROJECT_NAME} PRIVATE "Resource/suanPan_gcc.rc")
endif()

add_subdirectory(Converger)
add_subdirectory(Element)
add_subdirectory(Material)
add_subdirectory(Section)
add_subdirectory(Solver)

add_subdirectory(Constraint)
add_subdirectory(Domain)
add_subdirectory(Load)
add_subdirectory(Recorder)
add_subdirectory(Step)
add_subdirectory(Toolbox)

if(BUILD_DLL_EXAMPLE)
    add_subdirectory(Developer/Element)
    add_subdirectory(Developer/Material)
    add_subdirectory(Developer/Section)
    add_subdirectory(Developer/ModuleBundle)
endif()

add_subdirectory(Toolbox/lapack-ext)
add_subdirectory(Toolbox/arpack-src)
# add_subdirectory(Toolbox/superlu-src)

target_link_libraries(${PROJECT_NAME}
    Converger
    Element
    Material
    Section
    Solver
#   superlu
    arpack
    spmm)

if(USE_EXTERNAL_OPENBLAS)
target_link_libraries(${PROJECT_NAME} openblas)
else()
target_link_libraries(${PROJECT_NAME} openblas)
# add_subdirectory(Toolbox/blas-src)
# add_subdirectory(Toolbox/lapack-src)
# target_link_libraries(${PROJECT_NAME} lapack blas)
endif()
