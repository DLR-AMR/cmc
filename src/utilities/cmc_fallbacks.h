#ifndef CMC_FALLBACKS_H
#ifdef CMC_CONFIG_H
#define CMC_FALLBACKS_H

/* Definition of a standard error type */
#define CMC_ERR 1
typedef int cmc_err_t;
typedef void* cmc_err_ptr_t;
typedef void *cmc_err_func_ptr_t(void);

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
typedef MPI_Comm cmc_err_t;
typedef MPI_Info cmc_err_t;
typedef MPI_File cmc_err_t;
typedef MPI_Offset cmc_err_t;
#endif

/* Definitions to ensure compilation without netCDF */
#ifndef CMC_WITH_NETCDF
#define NC_NOWRITE 0
#define NC_MAX_VAR_DIMS 1024 /* equivalent to NC_MAX_VAR_DIMS (netCDF 4.9.0) */
#define cmc_nc_check_err(err) (cmc_err_msg("A netCDF error occured.\n"))
#define nc_strerror(err) "NetCDF error"
#endif

/* Definitions to ensure compilation without parallel netCDF */
#ifndef CMC_WITH_NETCDF_PAR

#endif

/* Definitions to ensure compilation without t8code */
#ifndef CMC_WITH_T8CODE
typedef cmc_err_ptr_t t8_forest_t;
typedef cmc_err_ptr_t t8_cmesh_t;
typedef cmc_err_ptr_t t8_element_t;
typedef cmc_err_ptr_t t8_eclass_scheme_c;
typedef cmc_err_func_ptr_t t8_forest_adapt_t;
typedef cmc_err_func_ptr_t t8_forest_replace_t;
typedef int t8_locidx_t;
typedef long long t8_gloidx_t;
#define t8_forest_unref(err) void(0);
#define t8_forest_ref(err) void(0);
#endif

/* Definitions to ensure compilation without t8code */
#ifndef CMC_WITH_MESSY

#endif

#endif /* Check existence of CMC_CONFIG_H */
#endif /* CMC_FALLBACKS_H */