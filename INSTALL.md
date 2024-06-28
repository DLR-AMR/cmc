
## Installation  

CMC is a software package that uses Autotools as a build system and relies on t8code (as a submodule).  
For parallelization MPI is used and optionally, it may be linked against netCDF.  

Many features are based on a working netCDF installation. Therefore, it is highly recommended to link against netCDF.

### Get the Code  

In order to get started, clone the repository, e.g.  
```
	git clone https://github.com/DLR-AMR/cmc.git
```

Suppose the project has been cloned into the directory `~/cmc`.
For the next steps, we need to change into the projects's directory, e.g. `cd cmc`.

### Initialize the submodule  

CMC relies on [t8code](https://github.com/DLR-AMR/t8code) as a submodule.  
Therefore, it has to be initialized and updated.  
```
	git submodule update --init --recursive
```
If a t8code installation is already available, this submodule-initialization step may be skipped.  
In this case, CMC has to be linked manually against t8code during the configuration (see the option `--with-t8code`).  
However, to ensure compatibility, it is recommended to use t8code as a submodule.  

### Create the configure script  

Afterwards, we need to call the `bootstrap` script.
CMC does not ship with a `configure` script. It needs to be generated on the local machine by calling the `bootstrap` script in the project's directory.

```
	./bootstrap
```

### Configuration of CMC  

Now, we have created a configure script. Move to the directory in which the CMC is ought to be built.  
Suppose we want to build it in the directory `~/cmc_build`.  

Therefore, we need to call the freshly created `configure` script from the directory `~/cmc_build`.  
```
../cmc/configure [OPTIONS]
```
Common options for the configuration are `--enable-debug`, `--enable-mpi`, `--with-netcdf`.  
A quick start may look like 
```
../cmc/configure --enable-mpi --with-netcdf
```
Please take a look at the configuration help for a complete list of available options,  
```
../cmc/configure --help
```

### Building of CMC  

In order to build the software with the previously configured options, we will need to call  
```
make -j
```

### Testing

After the compilation, we are able to verify the 'installed' software by building and running it's 'test'-target
```
make check -j
```

If all checks have passed, the software is good to use, otherwise please verify the availabilty of all dependencies or reach out to us via GitHub [CMC](https://github.com/DLR-AMR/cmc).
