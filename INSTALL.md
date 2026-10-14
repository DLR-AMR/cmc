
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
| CMC_ENABLE_MPI          | ON            |
| CMC_ENABLE_ALL_WARNINGS | OFF           |

Additional `<lib>_DIR`, `<lib>_ROOT` or `<lib>_ROOT_DIR` might help CMake to find the correct dependencies.

A configuration with t8code and MPI may look similar to
```
CC=mpicc CXX=mpicxx cmake -DCMAKE_BUILD_TYPE=Release -DCMC_ENABLE_MPI=ON -DCMC_WITH_T8CODE=ON -DT8CODE_DIR:PATH=<t8code_install_directory>/cmake -DP4EST_DIR:PATH=<t8code_install_directory>/cmake -DSC_DIR:PATH=<t8code_install_directory>/cmake -DCMAKE_INSTALL_PREFIX=<cmc_install_directory> ../cmc
```

In order to quickly test the patch-based compression approaches, the configuration may look like:
```
cmake -DCMAKE_BUILD_TYPE=Release -DCMC_WITH_T8CODE=OFF -DCMC_ENABLE_MPI=OFF ../cmc
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

Afterwards, the supplied tests may be executed with
```
ctest
```