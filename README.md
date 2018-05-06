***PLEASE NOTE THIS PROJECT IS UNDER DEVELOPMENT SO CURRENTLY YOU MAY NOT BE ABLE TO FIND SUFFICIENT INFORMATION***

suanPan
=======

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/7cb47e58d7dc4c1680c2205c4ba02e72)](https://www.codacy.com/app/TLCFEM/suanPan?utm_source=github.com&utm_medium=referral&utm_content=TLCFEM/suanPan&utm_campaign=Badge_Grade)
[![Build status](https://ci.appveyor.com/api/projects/status/fmdt0amjgd6dauf4?svg=true)](https://ci.appveyor.com/project/TLCFEM/suanpan)
[![Build Status](https://travis-ci.org/TLCFEM/suanPan.svg?branch=master)](https://travis-ci.org/TLCFEM/suanPan)
[![codecov](https://codecov.io/gh/TLCFEM/suanPan/branch/master/graph/badge.svg)](https://codecov.io/gh/TLCFEM/suanPan)
[![GitHub issues](https://img.shields.io/github/issues/TLCFEM/suanPan.svg)](https://github.com/TLCFEM/suanPan/issues)

Intro
-----

[**suanPan**](https://tlcfem.github.io/suanPan/) is an finite element method (FEM) simulation platform for applications in solid mechanics, civil/structural/seismic engineering. The name **suanPan** (in some other places also **suPan**) comes from the term *Suan Pan* (算盤), which literally means [Chinese abacus](https://en.wikipedia.org/wiki/Suanpan). **suanPan** is written in high quality C++ codes and is targeted to provide an efficient, concise and reliable FEM simulation platform.

**suanPan** is partially influenced by popular (non-)commercial FEA packages, such as [ABAQUS UNIFIED FEA](https://www.3ds.com/products-services/simulia/products/abaqus/), [ANSYS](http://www.ansys.com/) and [OpenSees](http://opensees.berkeley.edu/).

Check out the [documentation](https://tlcfem.gitbooks.io/suanpan/content/) (under construction).

Features
--------

The highlights of **suanPan** are

-   **suanPan** is fast.
-   **suanPan** is designed based on the [shared memory](https://en.wikipedia.org/wiki/Shared_memory) model and supports parallelism on heterogeneous architectures. For example multi-threaded CPU + GPU.
-   **suanPan** is open source and easy to be expanded to incorporate user-defined elements, materials, etc.
-   **suanPan** separates the FEA model part from the linear algebra operation part, which significantly reduces the complexity of development.
-   **suanPan** utilizes the new language features shipped with the latest standards (C++11 and C++14), such as new STL containers, smart pointers and many other features.

License
-------

**suanPan** is under `GPL v3`, you can click the badge above to check out the detailed clauses.

How to Compile
--------------

### Overview

As **suanPan** uses new language features, such as the `unique_ptr`, please use compilers that supports C++14 standard. For example,

-   **MSVC 14** (Visual Studio 2015) and/or later version,
-   **GNU GCC 5.4** and/or later version.
-   **Clang 3.4** and/or later version.

On Ubuntu, the external libraries are compiled with **GCC 4.8** using **binutils 2.24**, which is the default configuration for **Ubuntu 14.04 LTS**. On Windows, all libraries are compiled as dynamic libraries with **GCC 5.4** in **MSYS**. To avoid any potential linking error, please use later versions which are available from the **MinGW-w64** project.

1. [**GCC 5.4 x86_64-posix-seh**](https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win64/Personal%20Builds/mingw-builds/5.4.0/threads-posix/seh/x86_64-5.4.0-release-posix-seh-rt_v5-rev0.7z)
2. [**GCC 6.4 x86_64-posix-seh**](https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win64/Personal%20Builds/mingw-builds/6.4.0/threads-posix/seh/x86_64-6.4.0-release-posix-seh-rt_v5-rev0.7z)
3. [**GCC 7.3 x86_64-posix-seh**](https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win64/Personal%20Builds/mingw-builds/7.3.0/threads-posix/seh/x86_64-7.3.0-release-posix-seh-rt_v5-rev0.7z)

The **MSYS** project provides a basic compilation environment which can be downloaded from the [**MinGW**](https://sourceforge.net/projects/mingw/files/Installer/mingw-get-setup.exe/download) project.

The program is *deliberately designed to disable the backward compatibility*. **suanPan** uses [**CMake**](https://cmake.org/) (3.0 and later version) to manage builds.

### Windows

The package is tested under Windows with **MSVC++**, **GNU GCC**, **Intel C++ Compiler** and **Clang**. The libraries (only 64-bit, no 32-bit support anymore):

-   **ARPACK** version 0.96 (need a valid Fortran compiler)
-   **OpenBLAS** version 0.2.20 (without dynamic architecture)
-   **TBB** Threading Building Blocks version 2018U3
-   **HDF5** version 1.10.1

are bundled with the source code package.

#### Visual Studio

The default VS solution uses multi-threaded **OpenBLAS**. Simply `Build (F7)` the solution. You can change the linked library to other equivalent libraries such as **Intel MKL**, if those libraries are available on your machine. Since **OpenBLAS** is quite platform dependent, you are recommended to compile your own version instead of using the bundled one, which may probably not work on your machine. To compile the customized version of **OpenBLAS** on Windows, you may need to download [**MSYS**](http://www.mingw.org/wiki/msys).

The compiled program cannot run directly as it depends on other dynamic libraries. Please copy following files to the path that can be found by the program.

```text
/Libs/gcc-win/arpack.dll
/Libs/gcc-win/libgcc_s_seh-1.dll
/Libs/gcc-win/libgfortran-3.dll
/Libs/gcc-win/libopenblas.dll
/Libs/gcc-win/libquadmath-0.dll
/Libs/gcc-win/libstdc++-6.dll
/Libs/gcc-win/libwinpthread-1.dll
/Libs/gcc-win/spmm.dll
```

In addition, the **TBB** libraries shall be copied as well.

```text
/Libs/vs/tbb.dll
/Libs/vs/tbbmalloc.dll
/Libs/vs/tbbmalloc_proxy.dll
```

The VS solution can also be generated via **CMake**.

#### GNU GCC

Use **CMake** to generate Makefiles, assume current folder is the root of the package and your are using **MinGW**, the commands should look like this.

``` bash
# current folder is /suanPan/source/code/path
mkdir cmake-build && cd cmake-build
cmake -G "MinGW Makefiles" ..
make
```

The GUI may be a good tool for beginners to configure the build. There are several options provided to build the software with different configurations.

After successful compilation, the executable file is under `/cmake-build` folder. The dynamic libraries should also be copied.

```bash
# current folder is cmake-build
cp ../Libs/gcc-win/*.dll .
# run the program
./suanPan.exe
```

### Ubuntu

Again, the shipped **OpenBLAS** may not work on your platform, please compile your own library if any error occurs.

Make sure the compilers are installed.

```bash
sudo apt install gcc g++ gfortran binutils cmake cmake-qt-gui
```

Please do check the versions of those tools. A default configuration is enough for most cases. Simply create a build folder next to the source code folder and configure/make the program, such as

``` bash
# current folder is /suanPan/source/code/path
mkdir cmake-build
cd cmake-build
cmake ../suanPan
make
```

The multi-threaded version uses **TBB**, the corresponding path can be added.

```bash
# current folder cmake-build
export LD_LIBRAY_PATH=../Libs/gcc-linux:${LD_LIBRARY_PATH}
# run the program
./suanPan
```

### macOS

No progress yet.

Dependency
----------

Additional tools are used by **suanPan**, they are

-   [UPX](https://upx.github.io/) --- Executable File Packer

Additional libraries that may be used in **suanPan** are

-   [Armadillo](http://arma.sourceforge.net/) --- Linear Algebra Library
-   [MKL](https://software.intel.com/en-us/mkl) --- High Performance Linear Algebra Driver
-   [SuperLU](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/)

Those libraries may depend on other libraries such as

-   [zlib](https://zlib.net/)
-   [Szip](https://support.hdfgroup.org/doc_resource/SZIP/)

