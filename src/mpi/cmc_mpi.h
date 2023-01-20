#ifndef CMC_MPI_H
#define CMC_MPI_H

#include "cmc_config.h"
#include "utilities/cmc_constants_definitions.h"

#ifdef CMC_ENABLE_MPI
#include "mpi.h"
#endif


void
cmc_mpi_initialize();

void
cmc_mpi_finalize();

void cmc_mpi_abort(const int _err_code, const char* _location);


#define cmc_mpi_check_err(err) ((err) == MPI_SUCCESS ? (void) 0 : cmc_mpi_abort(err, CMC_FILE_LOCATION))


#endif /* CMC_MPI_H */