#include "cmc_messy_data.h"
#if 0
#include "fortran_interfaces/cmcc.h"
#include "utilities/cmc_constants_definitions.h"

#ifdef CMC_ENABLE_MPI
#include "mpi.h"
#endif

struct cmc_tracer_data {
    /* Name of the tracer-variable */
    char* var_name;
    /* Data array of the variable */
    double* var_data;

};


struct CMC_MESSY_DATA {
    /* Dimension of the tracer data (2D, 3D) */
    int data_dimension;
    /* Length of latitiude, longitude, elevation coordinate */
    int* dimension_sizes;

    /* Axis ordering of the data in C row-major layout */
    int* axis_ordering;
    /* Axis representation of the Fortran data */
    char* axis_representation_f;
    /* The ordering scheme to which the data complies */
    enum CMC_DATA_ORDERING_SCHEME data_scheme;

    /* Number of tracers */
    int num_variables;
    /* Pointer to an array of data variables */
    cmc_tracer_data_t** vars;

    /* Reference tracer data for computation of tracer concentration on grid cells */
    cmc_tracer_data_t* reference_var;
    /* Flag indicating whether a refence variable is present or not */
    int reference_var_present;
    
    /* Number of elements */
    size_t num_elements;

    /* A Missing value */
    double missing_value;

    /* The C-communicator to use */
    MPI_Comm comm;

    /* A helper variable */
    int _intern_id_counter;
};

#endif
