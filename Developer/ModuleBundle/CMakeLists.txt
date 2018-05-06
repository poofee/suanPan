cmake_minimum_required(VERSION 3.0)

project(ModuleBundle)

add_definitions(-DSUANPAN_DLL)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

include_directories(../Element)
include_directories(../Material)
include_directories(../Section)

if(USE_EXTERNAL_OPENBLAS)
link_libraries(openblas)
else()
link_libraries(openblas)
# link_libraries(lapack blas)
endif()

add_library(${PROJECT_NAME} SHARED
    "../../Domain/DomainBase.cpp"
    "../../Domain/Node.cpp"
    "../../Domain/Tag.cpp"
    "../../Element/Element.cpp"
    "../../Material/Material.cpp"
    "../../Section/Section.cpp"
    "../../Toolbox/debug.cpp"
    "../Element/ElementExample.cpp"
    "../Material/MaterialExample.cpp"
    "../Section/SectionExample.cpp"
    )
