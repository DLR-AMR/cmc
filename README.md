
### Currently, the repository is under construction and some functionality might not work as expected. Moreover, anything might change at any time! 

### Introduction

CMC is a software package providing data compression techniques based on methods of the field of adaptive mesh refinement (AMR).  
It is especially suited for geo-spatial data originating for example from Earth System Model (ESM) simulations.

CMC can be used either as a post-processing tool in order to read and compress data from binary or netCDF files or it can be directly linked to simulation codes in order to perform an online compression - either in a lossy or lossless fashion.

### Dependencies

CMC relies on functionality featured in C++20, therefore, a C++ compiler supporting this standard is needed, e.g. GCC version 11 or newer.

CMC uses [t8code](https://github.com/DLR-AMR/t8code) as its underlying AMR engine, which allows for broad variety of applications, since t8code provides a highly parallel and scalable AMR implementation of various element types.
Like t8code, cmc utilizes MPI for parallelization and although t8code works in a serial configuration (without MPI) as well, cmc restricts itself to the parallel usage of t8code.
Therefore, if cmc is configured with t8code, it should also be configured with MPI.
Moreover, cmc may be linked with netCDF in order to read input data from netCDF files for compression.
Additionally, a netCDF writer and reader is provided to store the compressed data in a custom netCDF format and read it therefrom for decompression, respectively.

It is encouraged to configure cmc with MPI and t8code linkage.
However, all dependencies are optional.
In order to quickly check, whether the compression approaches are suitable for the data at hand, a quick configuration without any dependencies may work to test the (serial) compression of data provided in multi-dimensional arrays.

### Installation

CMC relies on CMake to configure and build the software.
An autotools build process is currently not supported.

A short CMake setup will be described in the following.

CMC needs to be donwloaded, e.g.:
```
git clone https://github.com/DLR-AMR/cmc.git 
```

It is recommended to create an out-of-source build-directory, e.g.:
```
mkdir cmc_build && cd cmc_build
```

With `cmake -LH ../cmc` all possible configuration options can be displayed with their default values.
A sample of the configuration options reads:
| Option                  | Default Value |
| ----------------------- | ------------- |
| CMC_WITH_T8CODE         | ON            |
| CMC_WITH_NETCDF         | ON            |
| CMC_ENABLE_MPI          | ON            |
| CMC_ENABLE_ALL_WARNINGS | OFF           |

Additional `<lib>_DIR`, `<lib>_ROOT` or `<lib>_ROOT_DIR` might help CMake to find the correct dependencies.

A configuration with t8code and MPI may look similar to
```
CC=mpicc CXX=mpicxx cmake -DCMAKE_BUILD_TYPE=Release -DCMC_ENABLE_MPI=ON -DCMC_WITH_T8CODE=ON -DT8CODE_DIR:PATH=<t8code_install_directory>/cmake -DP4EST_DIR:PATH=<t8code_install_directory>/cmake -DSC_DIR:PATH=<t8code_install_directory>/cmake -DCMC_WITH_NETCDF=OFF -DCMAKE_INSTALL_PREFIX=<cmc_install_directory> ../cmc
```

In order to quickly test the patch-based compression approaches, the configuration may look like:
```
cmake -DCMAKE_BUILD_TYPE=Release -DCMC_WITH_T8CODE=OFF -DCMC_ENABLE_MPI=OFF -DCMC_WITH_NETCDF=OFF ../cmc
```

After the configuration the available options are summarized and displayed.

In order to build cmc with CMake
```
cmake --build .
```
might be called, and optionally it might be installed by
```
cmake --install .
```

### Getting Started

The software provides some example use cases in the `example`-directory.
In order to getting started quickly and test the compression approaches on multi-dimensional arrays, there are some CLI applications in the `workflows`-directory that perform the patch-based compression variations working without any of the above dependencies (i.e. t8code, MPI and netCDF).

The example usage and the required arguments of the predefined workflows are displayed by the help messages of the compression approaches, e.g. `./workflows/patch/cmc_lossless_prefix_extraction_compression -h`.
The usage of the corresponding decompression scheme is given by `./workflows/patch/cmc_lossless_decompression -h` as well.
