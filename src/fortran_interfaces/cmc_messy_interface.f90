MODULE cmc_messy_interface

    USE, INTRINSIC :: ISO_C_BINDING

    IMPLICIT NONE
    
    INTERFACE
    ! Allocate a MESSy_data struct and assign some initial values
    TYPE(C_PTR) FUNCTION cmcf_setup_messy(dimension_sizes, axis_representation, num_variables, &
    missing_value, comm_f) BIND(C, NAME='cmc_setup_messy')
    USE, INTRINSIC :: ISO_C_BINDING, only: C_PTR, C_INT, C_CHAR, C_DOUBLE
    IMPLICIT NONE
    TYPE(C_PTR), VALUE, INTENT(IN)       :: dimension_sizes
    CHARACTER(C_CHAR), VALUE, INTENT(IN) :: axis_representation
    REAL(C_DOUBLE), VALUE, INTENT(IN)    :: missing_value
    INTEGER(C_INT), VALUE, INTENT(IN)    :: num_variables
    INTEGER(C_INT), VALUE, INTENT(IN)    :: comm_f
    END FUNCTION cmcf_setup_messy
    END INTERFACE

    INTERFACE
    ! Set and save a variable and its data in the messy_data struct
    SUBROUTINE cmcf_messy_set_tracer_variable(messy_data, tracer_name, tracer_name_length, data_array) &
    BIND(C, NAME='cmc_messy_set_tracer_variable')
    USE, INTRINSIC :: ISO_C_BINDING, only: C_PTR, C_CHAR, C_INT
    IMPLICIT NONE
    TYPE(C_PTR), VALUE, INTENT(IN)       :: messy_data
    CHARACTER(C_CHAR), VALUE, INTENT(IN) :: tracer_name
    INTEGER(C_INT), VALUE, INTENT(IN)    :: tracer_name_length
    TYPE(C_PTR), VALUE, INTENT(IN)       :: data_array
    END SUBROUTINE cmcf_messy_set_tracer_variable
    END INTERFACE

    INTERFACE
    ! Set and save a refernece variable and its data in the messy_data struct
    ! (this may be the reference used for describing the tracer concentration per cell)
    SUBROUTINE cmcf_messy_set_reference_tracer(messy_data, tracer_name, tracer_name_length, data_array) &
    BIND(C, NAME='cmc_messy_set_reference_tracer')
    USE, INTRINSIC :: ISO_C_BINDING, only: C_PTR, C_CHAR, C_INT
    IMPLICIT NONE
    TYPE(C_PTR), VALUE, INTENT(IN)       :: messy_data
    CHARACTER(C_CHAR), VALUE, INTENT(IN) :: tracer_name
    INTEGER(C_INT), VALUE, INTENT(IN)    :: tracer_name_length
    TYPE(C_PTR), VALUE, INTENT(IN)       :: data_array
    END SUBROUTINE cmcf_messy_set_reference_tracer
    END INTERFACE

    INTERFACE
    ! Free the allocated tracer data
    SUBROUTINE cmcf_messy_free_tracer_variable(tracer_to_free) BIND(C, NAME='cmc_messy_free_tracer_variable')
    USE, INTRINSIC :: ISO_C_BINDING, only: C_PTR
    IMPLICIT NONE
    TYPE(C_PTR), VALUE, INTENT(IN) :: tracer_to_free
    END SUBROUTINE cmcf_messy_free_tracer_variable
    END INTERFACE

    INTERFACE
    ! Free the allocated messy_data struct
    SUBROUTINE cmcf_messy_free_data(messy_data) BIND(C, NAME='cmc_messy_free_data')
    USE, INTRINSIC :: ISO_C_BINDING, only: C_PTR
    IMPLICIT NONE
    TYPE(C_PTR), VALUE, INTENT(IN) :: messy_data
    END SUBROUTINE cmcf_messy_free_data
    END INTERFACE

END MODULE cmc_messy_interface