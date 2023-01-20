MODULE cmc_fortran

    USE, INTRINSIC :: ISO_C_BINDING

    IMPLICIT NONE

    INTERFACE
    ! Initialize cmc, MPI and t8code 
    SUBROUTINE cmcf_initialize() BIND(C, NAME='cmc_initialize')
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE
    END SUBROUTINE cmcf_initialize
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