MODULE cmc_amr_compression_interface

    USE, INTRINSIC :: ISO_C_BINDING

    IMPLICIT NONE

    INTERFACE
    ! Create a Lossy AMR Compression context from netCDF data
    TYPE(C_PTR) FUNCTION cmcf_create_amr_compression_data(ncdata, f_comm) &
    BIND(C, NAME='cmcc_create_amr_compression_data')
    USE, INTRINSIC :: ISO_C_BINDING, only: C_PTR, C_INT
    IMPLICIT NONE
    TYPE(C_PTR), VALUE, INTENT(IN)    :: ncdata
    INTEGER(C_INT), VALUE, INTENT(IN) :: f_comm
    END FUNCTION cmcf_create_amr_compression_data
    END INTERFACE

    INTERFACE
    ! Create a Lossy AMR Compression context from MESSy data
    TYPE(C_PTR) FUNCTION cmcf_create_amr_compression_data_messy(messy_data, f_comm) &
    BIND(C, NAME='cmcc_create_amr_compression_data_messy')
    USE, INTRINSIC :: ISO_C_BINDING, only: C_PTR, C_INT
    IMPLICIT NONE
    TYPE(C_PTR), VALUE, INTENT(IN)    :: messy_data
    INTEGER(C_INT), VALUE, INTENT(IN) :: f_comm
    END FUNCTION cmcf_create_amr_compression_data_messy
    END INTERFACE

    INTERFACE
    ! Setup the compression based on a given compression mode
    SUBROUTINE cmcf_amr_setup_compression(amr_data, compression_mode) &
    BIND(C, NAME='cmc_amr_setup_compression')
    USE, INTRINSIC :: ISO_C_BINDING, only: C_PTR, C_INT
    IMPLICIT NONE
    ENUM, BIND(C)
    ENUMERATOR :: CMC_T8_COMPRESSION_UNDEFINED = 0
    ENUMERATOR :: ONE_FOR_ALL_2D
    ENUMERATOR :: ONE_FOR_ONE_2D
    ENUMERATOR :: GROUPED_2D
    ENUMERATOR :: ONE_FOR_ALL_3D
    ENUMERATOR :: ONE_FOR_ONE_3D
    ENUMERATOR :: GROUPED_3D
    END ENUM
    TYPE(C_PTR), VALUE, INTENT(IN)    :: amr_data
    INTEGER(C_INT), VALUE, INTENT(IN) :: compression_mode
    END SUBROUTINE cmcf_amr_setup_compression
    END INTERFACE

    INTERFACE
    ! Perform the Lossy AMR Compression of the data
    ! based on the previously given compression mode
    SUBROUTINE cmcf_amr_compress(amr_data, adapt_func, interpolation_func) &
    BIND(C, NAME='cmc_amr_compress')
    USE, INTRINSIC :: ISO_C_BINDING, only: C_PTR, C_FUNPTR
    IMPLICIT NONE
    TYPE(C_PTR), VALUE, INTENT(IN)    :: amr_data
    TYPE(C_FUNPTR), VALUE, INTENT(IN) :: adapt_func
    TYPE(C_FUNPTR), VALUE, INTENT(IN) :: interpolation_func
    END SUBROUTINE cmcf_amr_compress
    END INTERFACE

    INTERFACE
    ! Destroy the allocated data
    SUBROUTINE cmcf_amr_destroy(amr_data) BIND(C, NAME='cmc_amr_destroy')
    USE, INTRINSIC :: ISO_C_BINDING, only: C_PTR
    IMPLICIT NONE
    TYPE(C_PTR), VALUE, INTENT(IN) :: amr_data
    END SUBROUTINE cmcf_amr_destroy
    END INTERFACE

    INTERFACE
    ! Decompress the lossy AMR compressed data
    SUBROUTINE cmcf_amr_decompress(amr_data) BIND(C, NAME='cmc_amr_decompress')
    USE, INTRINSIC :: ISO_C_BINDING, only: C_PTR
    IMPLICIT NONE
    TYPE(C_PTR), VALUE, INTENT(IN) :: amr_data
    END SUBROUTINE cmcf_amr_decompress
    END INTERFACE

    INTERFACE
    ! Destroy the allocated data
    SUBROUTINE cmcf_amr_write_vtk_file(amr_data, file_prefix, file_prefix_length) &
    BIND(C, NAME='cmcc_amr_write_vtk_file')
    USE, INTRINSIC :: ISO_C_BINDING, only: C_PTR, C_CHAR, C_INT
    IMPLICIT NONE
    TYPE(C_PTR), VALUE, INTENT(IN)       :: amr_data
    CHARACTER(C_CHAR), VALUE, INTENT(IN) :: file_prefix
    INTEGER(C_INT), VALUE, INTENT(IN)    :: file_prefix_length
    END SUBROUTINE cmcf_amr_write_vtk_file
    END INTERFACE

END MODULE cmc_amr_compression_interface