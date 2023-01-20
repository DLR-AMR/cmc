#ifndef CMC_FALLBACKS_H
#ifdef CMC_CONFIG_H
#define CMC_FALLBACKS_H

/* Definition of a standard error type */
#define CMC_ERR 1
typedef int cmc_err_t;

#if CMC_ENABLE_DEBUG
#define cmc_assert(condition) assert(condition)
#else
#define cmc_assert(condition) ((void)0)
#endif

#define cmc_static_assert(condition) static_assert(condition)

/* Definitions to ensure compilation without MPI */
#ifndef CMC_ENABLE_MPI
#define MPI_SUCCESS 0
#define MPI_COMM_WORLD -1
#define MPI_ERR_OTHER -1
#define MPI_INFO_NULL 0
typedef MPI_Comm int;
typedef MPI_Info int;
typedef MPI_File int;
typedef MPI_Offset int;
#endif

/* Definitions to ensure compilation without netCDF */
#ifndef CMC_WITH_NETCDF
#define NC_NOWRITE 0
#define NC_MAX_VAR_DIMS 1024 /* equivalent to NC_MAX_VAR_DIMS (netCDF 4.9.0) */
#define cmc_nc_err(err) (cmc_err_msg("A netCDF error occured.\n"))
#endif

/* Definitions to ensure compilation without parallel netCDF */
#ifndef CMC_WITH_NETCDF_PAR

#endif

/* Definitions to ensure compilation without t8code */
#ifndef CMC_WITH_T8CODE
typedef cmc_err_t t8_forest_t;
typedef cmc_err_t t8_cmesh_t;
#endif

/* Definitions to ensure compilation without t8code */
#ifndef CMC_WITH_MESSY

#endif

#endif /* Check existence of CMC_CONFIG_H */
#endif /* CMC_FALLBACKS_H */