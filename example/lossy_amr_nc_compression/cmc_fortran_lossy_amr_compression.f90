PROGRAM cmc_fortran_lossy_amr_compression
USE cmc_fortran
USE cmc_netcdf_interface
USE cmc_amr_compression_interface
USE, INTRINSIC :: ISO_C_BINDING

IMPLICIT NONE

include 'mpif.h'


INTEGER(C_INT) :: ncid
TYPE(C_PTR) :: ncdata, amr_data
CHARACTER (30), POINTER :: filename
CHARACTER (3), POINTER :: var1
CHARACTER (4), POINTER :: var2
INTEGER (KIND=8), POINTER, DIMENSION(:) :: start_ptr
INTEGER (KIND=8), POINTER, DIMENSION(:) :: count_ptr
INTEGER, PARAMETER :: comm = MPI_COMM_WORLD
INTEGER, PARAMETER :: compression_mode = 0
TYPE(C_PTR) :: c_filename, c_var1, c_var2, c_start_ptr, c_count_ptr

ALLOCATE(filename)
filename = '../data/ECMWF_ERA-40_subset.nc'
c_filename = C_LOC(filename)

ALLOCATE(var1)
var1 = 'p2t'
c_var1 = C_LOC(var1)
ALLOCATE(var2)
var2 = 'tco3'
c_var2 = C_LOC(var2)

ALLOCATE(start_ptr(1:3))
start_ptr = (/ 0, 0, 0 /)
ALLOCATE(count_ptr(1:3))
count_ptr = (/ 1, 73, 144/)
c_start_ptr = C_LOC(start_ptr)
c_count_ptr = C_LOC(count_ptr)

PRINT *, "Hello, this is a Fortran-Example of CMC's Lossy AMR Compressor"

! Initialize CMC
CALL cmcf_initialize
PRINT *, "[cmc] CMC has been initialized."

! Open a netCDF-File
ncid = cmcf_nc_open(c_filename, LEN(filename), comm)

! Create a netCDF struct
ncdata = cmcf_nc_create(ncid)

! Add several variables from the netCDF File which should be compressed
CALL cmcf_nc_add_variable(ncdata, c_var1, LEN(var1))
CALL cmcf_nc_add_variable(ncdata, c_var2, LEN(var2))

! Inquire the data of these variables
CALL cmcf_nc_inquire_vars(ncdata, c_start_ptr, c_count_ptr)

! Close the netCDF File (all data has been inquired)
CALL cmcf_nc_close(ncid)

! Create an Lossy AMR Compression context
amr_data = cmcf_create_amr_compression_data(ncdata, comm)

! Perform the setup of the compression
CALL cmcf_amr_setup_compression(amr_data, compression_mode)

! Execute the compression setup
CALL cmcf_amr_compress(amr_data, C_NULL_PTR, C_NULL_PTR)

! Destroy the netCDF data struct
CALL cmcf_nc_destroy(ncdata)

! Destroy the Lossy AMR data struct
CALL cmcf_amr_destroy(amr_data)

! Finalize CMC
CALL cmcf_finalize
PRINT *, "[cmc] CMC has been finalized."

END PROGRAM cmc_fortran_lossy_amr_compression