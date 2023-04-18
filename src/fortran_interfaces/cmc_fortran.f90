MODULE cmc_fortran

    USE, INTRINSIC :: ISO_C_BINDING

    IMPLICIT NONE

    INTERFACE
    ! Initialize cmc, MPI and t8code using MPI_COMM_WORLD as a default communicator (if MPI is enabled)
    SUBROUTINE cmcf_initialize() BIND(C, NAME='cmc_initialize')
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE
    END SUBROUTINE cmcf_initialize
    END INTERFACE

    INTERFACE
    ! Initialize cmc, MPI and t8code given a specific MPI communicator
    SUBROUTINE cmcf_initialize_mpi_comm(comm_f) BIND(C, NAME='cmcc_initialize_mpi_comm')
    USE, INTRINSIC :: ISO_C_BINDING, only: C_INT
    IMPLICIT NONE
    INTEGER(C_INT), VALUE, INTENT(IN) :: comm_f
    END SUBROUTINE cmcf_initialize_mpi_comm
    END INTERFACE

    INTERFACE
    ! Finalize cmc, MPI and t8code
    SUBROUTINE cmcf_finalize() BIND(C, NAME='cmc_finalize')
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE
    END SUBROUTINE cmcf_finalize
    END INTERFACE

    INTERFACE
    ! Finalize cmc and t8code
    SUBROUTINE cmcf_finalize_without_mpi() BIND(C, NAME='cmc_finalize_without_mpi')
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE
    END SUBROUTINE cmcf_finalize_without_mpi
    END INTERFACE

END MODULE cmc_fortran
