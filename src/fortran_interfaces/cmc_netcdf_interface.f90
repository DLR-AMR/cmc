MODULE cmc_netcdf_interface

    USE, INTRINSIC :: ISO_C_BINDING

    IMPLICIT NONE

    INTERFACE
    ! Create a netCDF data struct
    TYPE(C_PTR) FUNCTION cmcf_nc_create(ncid) BIND(C, NAME='cmc_nc_create')
    USE, INTRINSIC :: ISO_C_BINDING, only: C_PTR, C_INT
    IMPLICIT NONE
    INTEGER(C_INT), VALUE, INTENT(IN) :: ncid
    END FUNCTION cmcf_nc_create
    END INTERFACE

    INTERFACE
    ! Open a netCDF-File and return the id corresponding to this file
    INTEGER(C_INT) FUNCTION cmcf_nc_open(path_to_file, path_length, f_comm) &
    BIND(C, NAME='cmcc_nc_open')
    USE, INTRINSIC :: ISO_C_BINDING, only: C_PTR, C_INT
    IMPLICIT NONE
    TYPE(C_PTR), VALUE, INTENT(IN)    :: path_to_file
    INTEGER(C_INT), VALUE, INTENT(IN) :: path_length
    INTEGER(C_INT), VALUE, INTENT(IN) :: f_comm
    END FUNCTION cmcf_nc_open
    END INTERFACE

    INTERFACE
    ! Close the open netCDF File corresponding to the given id
    SUBROUTINE cmcf_nc_close(ncid) &
    BIND(C, NAME='cmc_nc_close')
    USE, INTRINSIC :: ISO_C_BINDING, only: C_INT
    IMPLICIT NONE
    INTEGER(C_INT), VALUE, INTENT(IN) :: ncid
    END SUBROUTINE cmcf_nc_close
    END INTERFACE

    INTERFACE
    ! Inquire information and data about the supplied varible names
    SUBROUTINE cmcf_nc_add_variable(ncdata, var_name, var_length) &
    BIND(C, NAME='cmcc_nc_add_variable')
    USE, INTRINSIC :: ISO_C_BINDING, only: C_PTR, C_INT
    IMPLICIT NONE
    TYPE(C_PTR), VALUE, INTENT(IN) :: ncdata
    TYPE(C_PTR), VALUE, INTENT(IN) :: var_name
    INTEGER(C_INT), VALUE, INTENT(IN) :: var_length
    END SUBROUTINE cmcf_nc_add_variable
    END INTERFACE

    INTERFACE
    ! Inquire the data corresponding to the previously added variables
    SUBROUTINE cmcf_nc_inquire_vars(ncdata, start_ptr, count_ptr) &
    BIND(C, NAME='cmcc_nc_inquire_added_vars')
    USE, INTRINSIC :: ISO_C_BINDING, only: C_PTR
    IMPLICIT NONE
    TYPE(C_PTR), VALUE, INTENT(IN) :: ncdata
    TYPE(C_PTR), VALUE, INTENT(IN) :: start_ptr
    TYPE(C_PTR), VALUE, INTENT(IN) :: count_ptr
    END SUBROUTINE cmcf_nc_inquire_vars
    END INTERFACE

    INTERFACE
    ! Destroy the allocated data
    SUBROUTINE cmcf_nc_destroy(ncdata) BIND(C, NAME='cmc_nc_destroy')
    USE, INTRINSIC :: ISO_C_BINDING, only: C_PTR
    IMPLICIT NONE
    TYPE(C_PTR), VALUE, INTENT(IN) :: ncdata
    END SUBROUTINE cmcf_nc_destroy
    END INTERFACE

END MODULE cmc_netcdf_interface