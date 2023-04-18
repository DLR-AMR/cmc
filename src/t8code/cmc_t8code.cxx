#include "cmc.h"
#include "cmc_t8code.h"
/* Include t8code dependent header files */
#ifdef CMC_WITH_T8CODE
#include <t8.h>
#endif
#include "utilities/cmc_log_functions.h"


/* Initialize t8code and it's submodules */
void
cmc_t8code_initialize(MPI_Comm comm)
{
    #ifdef CMC_WITH_T8CODE
    int mpi_already_initialized{0};
    int mpiret;
    #if CMC_ENABLE_MPI
    /* Check if MPI has been initialized already */ 
    mpiret = MPI_Initialized(&mpi_already_initialized);
    cmc_mpi_check_err(mpiret);
    #endif
    if(mpi_already_initialized == 0)
    {
        std::cout << "t8code initializes mpi" << std::endl;
        mpiret = sc_MPI_Init(NULL, NULL);
        SC_CHECK_MPI(mpiret);
    }

    sc_init(comm, 1, 1, NULL, SC_LP_ESSENTIAL);

    t8_init(SC_LP_PRODUCTION);

    #else
    cmc_err_msg("CMC is not linked against t8code.");
    #endif
}

/* Finalize t8code and it's submodules */
void
cmc_t8code_finalize(const int flag_finalize_mpi)
{
    #ifdef CMC_WITH_T8CODE
    sc_finalize();
    if (flag_finalize_mpi)
    {
        std::cout << "t8code finalizes mpi"<<std::endl;
        int mpiret{sc_MPI_Finalize()};
        SC_CHECK_MPI(mpiret);
    }
    #endif
}
