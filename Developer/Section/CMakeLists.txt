cmake_minimum_required(VERSION 3.0)

project(SectionExample)

add_definitions(-DSUANPAN_DLL)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

if(USE_EXTERNAL_OPENBLAS)
link_libraries(openblas)
else()
link_libraries(openblas)
# link_libraries(lapack blas)
endif()

add_library(${PROJECT_NAME} SHARED
    "SectionExample.cpp"
    "../../Domain/Tag.cpp"
    "../../Material/Material.cpp"
    "../../Section/Section.cpp"
    "../../Toolbox/debug.cpp"
    )
