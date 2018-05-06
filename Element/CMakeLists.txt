cmake_minimum_required(VERSION 3.0)

project(Element)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../Libs)

include(Beam/CMakeLists.txt)
include(Cube/CMakeLists.txt)
include(Membrane/CMakeLists.txt)
include(Plate/CMakeLists.txt)
include(Shell/CMakeLists.txt)
include(Special/CMakeLists.txt)
include(Truss/CMakeLists.txt)
include(Utility/CMakeLists.txt)

add_library(${PROJECT_NAME} STATIC
    Element.cpp
    ElementParser.cpp
    ElementTemplate.cpp
    MaterialElement.cpp
    SectionElement.cpp
    ${Beam}
    ${Cube}
    ${Membrane}
    ${Plate}
    ${Shell}
    ${Special}
    ${Truss}
    ${Utility}
    )