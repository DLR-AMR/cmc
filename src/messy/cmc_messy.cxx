#include "cmc_messy.h"
#include "utilities/cmc_constants_definitions.h"
//#include "cmc.h"
#include "utilities/cmc_log_functions.h"
#include "utilities/cmc_util.h"
#include "utilities/cmc_geo_util.h"


struct cmc_messy_data {
    cmc_messy_data(){};
    cmc_messy_data(int _num_variables)
    : num_variables{_num_variables}{
        std::cout << "Created messy data struct with num vars = " << _num_variables << std::endl;
        vars.reserve(_num_variables);
    };
    cmc_messy_data(MPI_Comm _comm)
    : comm{_comm}{};

    /* Number of tracers (excluding the reference variable) */
    int num_variables{0};

    /* Dimension of the tracer data (2D, 3D) */
    int data_dimension{0};

    /* A flag whether a refernce variable is given */
    bool reference_values_present{false};

    /* Length of latitiude, longitude, elevation coordinate */
    std::vector<int> dimension_sizes;

    /* A vector containing all variables (excluding a reference variable) */
    std::vector<cmc_var_t> vars;

    /* A potential reference variable */
    cmc_var_t reference_variable;

    /* Axis ordering of the data in C row-major layout */
    int* axis_ordering;

    /* Axis representation of the Fortran data */
    std::array<char, 4> axis_representation_f;

    /* The C-communicator to use */
    MPI_Comm comm;
};


cmc_messy_data_t
cmc_setup_messy_data(const MPI_Fint comm)
{
    #ifdef CMC_ENABLE_FORTRAN
    std::cout << "Hallo" << std::endl;
    MPI_Comm c_comm = 0;
    #ifdef CMC_ENABLE_MPI
    /** Convert the Fortran MPI Communicator to a 'C-MPI Communicator' */
    c_comm = MPI_Comm_f2c(comm);
    #endif
    cmc_debug_msg("The cmc_messy_data struct will be initialized and passed back to MESSy");
    return new cmc_messy_data(c_comm);
    #endif
}


void
cmc_destroy_messy_data(cmc_messy_data_t messy_data)
{
    #ifdef CMC_ENABLE_FORTRAN
    if (messy_data != nullptr)
    {
        delete messy_data;
    }
    cmc_debug_msg("The messy data struct has been deallocated.");
    #endif
}




#if 0
//TODO: update

/* Opaque pointer to cmc_tracer_data (currently, it resides here, because no public exposure is needed at the moment) */
typedef struct cmc_tracer_data* cmc_tracer_data_t;


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
    cmc_tracer_data_t* vars;

    /* Reference tracer data for computation of tracer concentration on grid cells */
    cmc_tracer_data_t reference_var;
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

static void
cmc_messy_get_c_data_layout_from_fortran_data(CMC_MESSY_DATA_T messy_data)
{
    int counter = 0;

    /* Get the linear ordering of the Fortran data */
    for (int axis_id = 0; axis_id < CMC_NUM_COORD_IDS; axis_id++)
    {
        messy_data->axis_ordering[axis_id] = CMC_COORDINATE_NOT_CONSIDERED;
        if ((messy_data->axis_representation_f)[axis_id] == 'X' )
        {
            messy_data->axis_ordering[counter] = CMC_LON;
            counter++;
        } else if ((messy_data->axis_representation_f)[axis_id] == 'Y' )
        {
            messy_data->axis_ordering[counter] = CMC_LAT;
            counter++;
        } else if ((messy_data->axis_representation_f)[axis_id] == 'Z' && messy_data->data_dimension == 3)
        {
            messy_data->axis_ordering[counter] = CMC_LEV;
            counter++;
        }
    }

    /* Reverse the obtained Fortran ordering in order to obtain the C ordering of this data */
    if (messy_data->data_dimension == 3)
    {
        int tmp = messy_data->axis_ordering[0];
        messy_data->axis_ordering[0] = messy_data->axis_ordering[2];
        messy_data->axis_ordering[2] = tmp;
    } else if (messy_data->data_dimension == 2)
    {
        int tmp = messy_data->axis_ordering[0];
        messy_data->axis_ordering[0] = messy_data->axis_ordering[1];
        messy_data->axis_ordering[1] = tmp;
    } else {
        cmc_err_msg("Only 2D and 3D data (defined on the coordinates latitude, longitude and elevation) is supported.");
    }
}


CMC_MESSY_DATA_T cmc_setup_messy(int* _dimension_sizes, const char* axis_representation, const int _num_variables, const double _missing_value, MPI_Fint comm_f)
{
    #ifdef CMC_ENABLE_FORTRAN
    cmc_debug_msg("Setting up a messy_data struct...");
    /* Allocate memory for a MESSy_Data class */
    CMC_MESSY_DATA_T messy_data = (CMC_MESSY_DATA_T) malloc(sizeof(struct CMC_MESSY_DATA));
    /* Allocate an array holding the dimension sizes */
    messy_data->dimension_sizes = (int*) malloc(sizeof(int) * CMC_NUM_COORD_IDS);
    /* Copy the lengths */
    memcpy(messy_data->dimension_sizes, _dimension_sizes, sizeof(int) * CMC_NUM_COORD_IDS);
    /* Allocate an array holding the Fortran axis_representation */
    messy_data->axis_representation_f = (char*) malloc(sizeof(char) * (CMC_NUM_COORD_IDS +1));
    /* Copy the representation and add a null terminating character */
    memcpy(messy_data->axis_representation_f, _dimension_sizes, sizeof(char) * CMC_NUM_COORD_IDS);
    (messy_data->axis_representation_f)[CMC_NUM_COORD_IDS] = '\0';
    /* Save the amount of variables to consider */
    messy_data->num_variables = _num_variables;
    /* Allocate memory for an array holding the variables */
    messy_data->vars = (cmc_tracer_data_t*) malloc(sizeof(cmc_tracer_data_t) * _num_variables);
    /* Flag indicating whether a reference varibale has been set or not */
    messy_data->reference_var_present = 0;

    /* Convert the Fortran MPI-Communicator to a C-communicator */
    messy_data->comm = MPI_Comm_f2c(comm_f);

    /* Calculate the number of elements */
    messy_data->num_elements = 1;
    /* (Currently,) only latitude, longitude and eleveation data can be considered */
    for (int axis_id = 0; axis_id < CMC_NUM_COORD_IDS; axis_id++)
    {
        if ((messy_data->axis_representation_f)[axis_id] == 'X' ||
            (messy_data->axis_representation_f)[axis_id] == 'Y' )
        {
            messy_data->num_elements *= (size_t) (messy_data->dimension_sizes)[axis_id];
        } else if ((messy_data->axis_representation_f)[axis_id] == 'Z')
        {
            if ((messy_data->dimension_sizes)[axis_id] > 1)
            {
                messy_data->num_elements *= (size_t) (messy_data->dimension_sizes)[axis_id];
                messy_data->data_dimension = 3;
            } else
            {
                messy_data->data_dimension = 2;
            }
        }
    }
    /* Set the counter to zero since no variables are yet created/saved */
    messy_data->_intern_id_counter = 0;

    /* The data is linearily ordered compliant to the coordinate varibales (lat, lon, lev) */
    messy_data->data_scheme = CMC_GEO_DATA_LINEAR;

    /* Obtain the equivalent C row-major data layout from the ordering of the Fortran data layout */
    /* The correspondance between MESSy and cmc coordinates is:
     *  'X' equals CMC_LON
     *  'Y' equals CMC_LAT
     *  'Z' equals CMC_LEV
    */
    messy_data->axis_ordering = (int*) malloc(sizeof(int) * CMC_NUM_COORD_IDS);
    cmc_messy_get_c_data_layout_from_fortran_data(messy_data);

    cmc_debug_msg("The messy_data struct has been set up.");
    return messy_data;

    #endif
}


void
cmc_messy_set_tracer_variable(CMC_MESSY_DATA_T messy_data, const char*  tracer_name, const int tracer_name_length, const double* data_array)
{
    #ifdef CMC_ENABLE_FORTRAN
    /* We receive a double pointer to the tracer data of the given variable data which lays in a contiguous array.
     * We are not able to to just own/move the data pointer, since the allocation and deallocation of this data should remain on the MESSy side.
     * Therefore, we have to copy the data (but since no underlying mesh is defined at this point, we cannot directly reorder compliant to the Morton curve order)
     * We can convert the linear order of the tracer in Fortran colum-major layout to C row-major layout without explicitly reordering the data
     * Since all tracer variables have the same data layout, we are converting the ordering directly during the initialization
    */

    /* Allocate memory for the variable and it's data */
    messy_data->vars[messy_data->_intern_id_counter] = (cmc_tracer_data_t) malloc(sizeof(struct cmc_tracer_data));
    messy_data->vars[messy_data->_intern_id_counter]->var_data = (double*) malloc(sizeof(double) * messy_data->num_elements);

    /* Save the name of the variable */
    messy_data->vars[messy_data->_intern_id_counter]->var_name = (char*) malloc(sizeof(char) * (tracer_name_length + 1));
    /* Copy the 'Fortran'-name */
    memcpy(messy_data->vars[messy_data->_intern_id_counter]->var_name, tracer_name, sizeof(char) * tracer_name_length);
    /* Add a trailing null-terminating character */
    messy_data->vars[messy_data->_intern_id_counter]->var_name[tracer_name_length] = '\0';
    
    /* Copy the data over */
    memcpy(messy_data->vars[messy_data->_intern_id_counter]->var_data, data_array, sizeof(double) * messy_data->num_elements);
    
    /* Update the intern counter */
    (messy_data->_intern_id_counter)++;

    cmc_debug_msg("Tracer varibale ", tracer_name, " has beeen added to messy_data.");
    #endif
}

void
cmc_messy_set_reference_tracer(CMC_MESSY_DATA_T messy_data, const char*  tracer_name, const int tracer_name_length, const double* data_array)
{
    #ifdef CMC_ENABLE_FORTRAN
    if (messy_data->reference_var_present != 0)
    {
        cmc_err_msg("A reference variable for the MESSy data struct has already been added.");
    }

    /* Allocate memory for the reference variable and it's data */
    messy_data->reference_var = (cmc_tracer_data_t) malloc(sizeof(struct cmc_tracer_data));
    messy_data->reference_var->var_data = (double*) malloc(sizeof(double) * messy_data->num_elements);
    
    /* Set the flag inidicating the presence of a reference tracer */
    messy_data->reference_var_present = 1;

    /* Save the name of the variable */
    messy_data->reference_var->var_name = (char*) malloc(sizeof(char) * (tracer_name_length + 1));
    /* Copy the 'Fortran'-name */
    memcpy(messy_data->reference_var->var_name, tracer_name, sizeof(char) * tracer_name_length);
    /* Add a trailing null-terminating character */
    messy_data->reference_var->var_name[tracer_name_length] = '\0';
    
    /* Copy the data over */
    memcpy(messy_data->reference_var->var_data, data_array, sizeof(double) * messy_data->num_elements);

    cmc_debug_msg("Tracer varibale ", tracer_name, " has beeen added to as a reference varibale to messy_data.");
    
    #endif
}

static void
cmc_messy_free_tracer_variable(cmc_tracer_data_t _tracer_to_free)
{
    #ifdef CMC_ENABLE_FORTRAN
    /* Free the name of the tracer */
    free(_tracer_to_free->var_name);
    /* Free the data of the tracer variable */
    free(_tracer_to_free->var_data);
    #endif
}

void
cmc_messy_free_data(CMC_MESSY_DATA_T messy_data)
{
    #ifdef CMC_ENABLE_FORTRAN
    /* Free all allocations of the messy_data struct */
    free(messy_data->dimension_sizes);
    free(messy_data->axis_ordering);
    free(messy_data->axis_representation_f);
    if (messy_data->reference_var_present != 0)
    {
        free(messy_data->reference_var);
    }
    /* Free all single variables */
    for (int id = 0; id < messy_data->_intern_id_counter; ++id)
    {
        cmc_messy_free_tracer_variable(messy_data->vars[id]);
    }
    free(messy_data->vars);
    #endif
}

#ifdef CMC_WITH_T8CODE
void
_cmc_transform_messy_data_to_t8code_data(CMC_MESSY_DATA_T messy_data,  CMC_T8_DATA* t8_data)
{
    #ifdef CMC_ENABLE_FORTRAN
    cmc_debug_msg("Setting up the compression with MESSy data and t8code.");
    /* Store the number of variables considered */
    t8_data->num_variables = (messy_data->reference_var_present != 0) ? messy_data->num_variables +1 : messy_data->num_variables;
    /* Allocate an array holding pointers to the variable data */
    t8_data->var_data = new CMC_T8_VAR_DATA*[t8_data->num_variables];
    /* create a 'geo_data' which holds information needed for the creation of a t8code mesh and calculate the mesh dimension which is induced by the (geo-spatial) netCDF coordinates */
    t8_data->geo_data = new CMC_T8_GEO_MESH_DATA{messy_data->data_dimension};

    /* Assign the coordinate array int8_data.geo_data */
    t8_data->geo_data->data_size[CMC_COORD_IDS::CMC_LAT] = messy_data->dimension_sizes[CMC_COORD_IDS::CMC_LAT];
    t8_data->geo_data->data_size[CMC_COORD_IDS::CMC_LON] = messy_data->dimension_sizes[CMC_COORD_IDS::CMC_LON];
    t8_data->geo_data->data_size[CMC_COORD_IDS::CMC_LEV] = messy_data->dimension_sizes[CMC_COORD_IDS::CMC_LEV];
    t8_data->geo_data->data_size[CMC_COORD_IDS::CMC_TIME] = messy_data->dimension_sizes[CMC_COORD_IDS::CMC_TIME];
    
    /* Loop over all variables */
    for (int var_id{0}; var_id < messy_data->num_variables; ++var_id)
    {
        /* Allocate new t8_data_variable */
        t8_data->var_data[var_id] = new CMC_T8_VAR_DATA{};

        /* Save the name of the variable */
        t8_data->var_data[var_id]->name = std::string(messy_data->vars[var_id]->var_name);
        /* Initialize each data variable in t8_data */
        t8_data->var_data[var_id]->missing_value_present = false;
        t8_data->var_data[var_id]->missing_value = messy_data->missing_value;
        
        /* Set the data ordering scheme */
        t8_data->var_data[var_id]->data_scheme = CMC_DATA_ORDERING_SCHEME::CMC_GEO_DATA_LINEAR;
    
        /* Copy the axis ordering */
        std::memcpy(static_cast<void*>(&(t8_data->var_data[var_id]->data_axis_ordering)[0]), static_cast<void*>(&(messy_data->axis_ordering[0])), CMC_NUM_COORD_IDS * sizeof(int));
        /* Copy the dimension sizes */
        std::memcpy(static_cast<void*>(&(t8_data->var_data[var_id]->variable_dimension_sizes)[0]), static_cast<void*>(&(messy_data->dimension_sizes[0])), CMC_NUM_COORD_IDS * sizeof(int));

        /* Save the pointer to the data array */
        (t8_data->var_data[var_id]->data).push_back(new cmc_array_t<double>{messy_data->num_elements, static_cast<void*>(&(messy_data->vars[var_id]->var_data)[0]), CMC_DOUBLE});
        messy_data->vars[var_id]->var_data = nullptr;
    }
    
    /* Add the reference variable at the end */
    if (messy_data->reference_var_present != 0)
    {
        /* Allocate new t8_data_variable */
        t8_data->var_data[messy_data->num_variables] = new CMC_T8_VAR_DATA{};
        /* Save the name of the variable */
        t8_data->var_data[messy_data->num_variables]->name = std::string(messy_data->reference_var->var_name);
        /* Initialize each data variable in t8_data */
        t8_data->var_data[messy_data->num_variables]->missing_value_present = false;
        t8_data->var_data[messy_data->num_variables]->missing_value = messy_data->missing_value;
        /* Save the data dimension */
        t8_data->var_data[messy_data->num_variables]->data_dim = messy_data->data_dimension;
        /* Set the data ordering scheme */
        t8_data->var_data[messy_data->num_variables]->data_scheme = CMC_DATA_ORDERING_SCHEME::CMC_GEO_DATA_LINEAR;
    
        /* Copy the axis ordering */
        memcpy(static_cast<void*>(&(t8_data->var_data[messy_data->num_variables]->data_axis_ordering)[0]), static_cast<void*>(&(messy_data->axis_ordering[0])), CMC_NUM_COORD_IDS * sizeof(int));
        /* Copy the dimension sizes */
        memcpy(static_cast<void*>(&(t8_data->var_data[messy_data->num_variables]->variable_dimension_sizes)[0]), static_cast<void*>(&(messy_data->dimension_sizes[0])), CMC_NUM_COORD_IDS * sizeof(int));

        /* Save the pointer to the data array */
        (t8_data->var_data[messy_data->num_variables]->data).push_back(new cmc_array_t<double>{messy_data->num_elements, static_cast<void*>(&(messy_data->reference_var->var_data)[0]), CMC_DOUBLE});
        messy_data->reference_var->var_data = nullptr;
    }
    cmc_debug_msg("The data has been transferred from MESSy data classes to t8code classes");

    #endif
}

#endif

#endif