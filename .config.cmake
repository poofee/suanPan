include_directories(${ROOT})
include_directories(${ROOT}/Include)
include_directories(${ROOT}/Include/armadillo)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(USE_HDF5 "ADD HDF5 SUPPORT TO DATA STORAGE" ON)
# option(USE_BUNDLED_BLAS "USE BUNDLED NETLIB BLAS" ON)
# option(USE_BUNDLED_LAPACK "USE BUNDLED NETLIB LAPACK" ON)
option(USE_BUNDLED_ARPACK "USE BUNDLED ARPACK" ON)
option(USE_BUNDLED_SUPERLU "USE BUNDLED SUPERLU" ON)
option(USE_EXTERNAL_MAGMA "USE EXTERNAL MAGMA LIBRARY TO UTILIZE GPU" OFF)
option(USE_EXTERNAL_OPENBLAS "USE EXTERNAL OPENBLAS LIBRARY TO UTILIZE CPU" ON)
option(BUILD_MULTITHREAD "BUILD MULTI THREAD VERSION" OFF)
option(BUILD_DLL_EXAMPLE "BUILD DYNAMIC LIBRARY EXAMPLE" ON)

if(USE_BUNDLED_SUPERLU)
set(USE_BUNDLED_SUPERLU OFF)
add_definitions(-DARMA_DONT_USE_SUPERLU)
endif()

if(CMAKE_SYSTEM_NAME MATCHES "Windows") # WINDOWS PLATFORM
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU") # GNU GCC COMPILER
        set(COMPILER_IDENTIFIER "gcc-win")
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC") # MSVC COMPILER
        set(COMPILER_IDENTIFIER "vs")
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel") # INTEL COMPILER
        set(COMPILER_IDENTIFIER "vs")
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang") # CLANG COMPILER
        set(COMPILER_IDENTIFIER "vs")
    else()
        set(COMPILER_IDENTIFIER "unknown")
    endif()
elseif(CMAKE_SYSTEM_NAME MATCHES "Linux") # LINUX PLATFORM
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU") # GNU GCC COMPILER
        set(COMPILER_IDENTIFIER "gcc-linux")
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel") # INTEL COMPILER
        set(COMPILER_IDENTIFIER "unknown")
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang") # CLANG COMPILER
        set(COMPILER_IDENTIFIER "unknown")
    else()
        set(COMPILER_IDENTIFIER "unknown")
    endif()
endif()

if(COMPILER_IDENTIFIER MATCHES "unknown")
    message(FATAL_ERROR "cannot identify the compiler available.")
endif()

link_directories(${ROOT}/Libs/${COMPILER_IDENTIFIER})

if(USE_HDF5)
    add_definitions(-DSUANPAN_HDF5)
    include_directories(${ROOT}/Include/hdf5-common)
    include_directories(${ROOT}/Include/hdf5-${COMPILER_IDENTIFIER})
    link_libraries(hdf5_hl hdf5)
else()
    add_definitions(-DARMA_DONT_USE_HDF5)
endif()

if(USE_MAGMA)
    find_library(MAGMA_LIBRARY magma.lib)
    if(MAGMA_LIBRARY MATCHES "MAGMA_LIBRARY-NOTFOUND")
        message("MAGMA library is not found, please indicate the location.")
    else()
        set(CUDA_ROOT $ENV{CUDA_PATH} CACHE STRING "CUDA LIBRARY ROOT CONTAINS /INCLUDE")
        if(CUDA_ROOT MATCHES "")
            message("CUDA library is not found, please indicate the location.")
        else()
            add_definitions(-DSUANPAN_MAGMA)
            include_directories(${CUDA_ROOT}/include)
        endif()
    endif()
endif()

if(BUILD_MULTITHREAD)
    add_definitions(-DSUANPAN_MT)
    link_libraries(tbb tbbmalloc tbbmalloc_proxy)
endif()

if(COMPILER_IDENTIFIER MATCHES "vs")
    set(CMAKE_CONFIGURATION_TYPES "Release" CACHE STRING "Debug Release RelWithDebInfo MinSizeRel" FORCE)
    unset(TEST_COVERAGE CACHE)

    set(CMAKE_CXX_FLAGS "/MP /EHsc /arch:AVX")
    set(CMAKE_CXX_FLAGS_DEBUG "/MTd /D \"DEBUG\"")
    set(CMAKE_CXX_FLAGS_RELEASE "/MT")

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /O3 /MP /fpp /names:lowercase /assume:underscore")

elseif(COMPILER_IDENTIFIER MATCHES "gcc-win")

elseif(COMPILER_IDENTIFIER MATCHES "gcc-linux")
    link_libraries(dl)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU") # GNU GCC COMPILER
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Debug Release RelWithDebInfo MinSizeRel" FORCE)
    option(TEST_COVERAGE "TEST CODE COVERAGE USING GCOV" OFF)

    link_libraries(pthread gfortran quadmath)

    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fexceptions -mavx")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG")
    if(TEST_COVERAGE) # COVERAGE ONLY AVAILABLE UNDER GCC
        set(CMAKE_CXX_FLAGS "-fprofile-arcs -ftest-coverage")
        link_libraries(gcov)
    endif()

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -cpp")

endif()
