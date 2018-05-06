cmake_minimum_required(VERSION 3.0)

project(ElementExample)

add_definitions(-DSUANPAN_DLL)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

if(USE_EXTERNAL_OPENBLAS)
link_libraries(openblas)
else()
link_libraries(openblas)
# link_libraries(lapack blas)
endif()

add_library(${PROJECT_NAME} SHARED
    "ElementExample.cpp"
    "../../Element/Element.cpp"
    "../../Domain/DomainBase.cpp"
    "../../Domain/Node.cpp"
    "../../Domain/Tag.cpp"
    "../../Toolbox/debug.cpp"
    )
