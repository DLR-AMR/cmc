
### !Currently, the repository is under construction and some functionality might not work as expected! 

### Introduction

CMC is a software package providing data compression techniques based on methods of the field of adaptive mesh refinement (AMR).  
It is especially suited for geo-spatial data originating for example from Earth System Model (ESM) simulations.  

CMC can be used either as a post-processing tool in order to read and compress data from netCDF files or it can be directly linked to simulation codes in order to perform an online compression.  

Interfaces for C and Fortran codes as well as several further compression approaches are underway.  

The capabilities of CMC encompass compression based on point-wise absolute and relative error critera. Moreover, splitting of higher dimensional data into several lower dimensional data slices is provided alongside the opportunity to formulate region-wise varying error thresholds for the compression to comply to - that includes in particular nested error domains.  

Besides the opportunity to perform lossy compression, a lossless compression mode is available as well.  

CMC uses [t8code](https://github.com/DLR-AMR/t8code) as its underlying AMR engine, which allows for broad variety of applications, since t8code provides a highly parallel and scalable AMR implementation of various element types.  
  