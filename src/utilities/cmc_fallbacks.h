#ifndef CMC_FALLBACKS_H
#ifdef CMC_CONFIG_H
#define CMC_FALLBACKS_H

namespace cmc
{

/* Definition of a standard error type */
#define CMC_ERR 1
typedef int cmc_err_t;
typedef void* cmc_err_ptr_t;
typedef void *cmc_err_func_ptr_t(void);

#ifdef CMC_ENABLE_DEBUG
#define cmc_assert(condition) assert(condition)
#else
#define cmc_assert(condition) ((void)0)
#endif

#define cmc_static_assert(condition) static_assert(condition)

/* Definitions to ensure compilation without MPI */
#ifndef CMC_ENABLE_MPI
#if 0
#define MPI_SUCCESS 0
#define MPI_COMM_WORLD -1
#define MPI_COMM_NULL -1
#define MPI_ERR_OTHER -1
#define MPI_INFO_NULL 0
#define MPI_COMM_SELF -1
#define MPI_DATATYPE_NULL -1;
typedef cmc_err_t MPI_Comm;
typedef cmc_err_t MPI_Datatype;
typedef cmc_err_t MPI_Request;

//typedef cmc_err_t MPI_Info;
//typedef cmc_err_t MPI_File;
//typedef cmc_err_t MPI_Offset;
//template <typename T> class VariableMessage {};
#endif
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

}

#endif /* Check existence of CMC_CONFIG_H */
#endif /* CMC_FALLBACKS_H */
