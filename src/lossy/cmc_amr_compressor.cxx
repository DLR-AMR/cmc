#include "cmc_amr_compressor.h"
#include "t8code/cmc_t8code_data.hxx"
#include "utilities/cmc_log_functions.h"
#include "utilities/cmc_container.h"
#include "component_interfaces/cmc_t8_nc.h"
#include "component_interfaces/cmc_t8_messy.h"
#include "netcdf/cmc_netcdf.h"
#include "messy/cmc_messy.h"
#include "t8code/cmc_t8code.h"
#include "t8code/cmc_t8code_geo_data.h"
#include "t8code/cmc_t8code_geo_mesh.h"
#include "t8code/cmc_t8_adapt_callbacks.h"
#include "t8code/cmc_t8_replace_callbacks.h"

#include <iostream>
#include <fstream>

/* Class definitions */
struct cmc_amr_data
{
    cmc_t8_data_t t8_data{nullptr};
    bool compression_applied{false};

    cmc_amr_data(){};
    cmc_amr_data(cmc_t8_data_t initial_t8_data)
    : t8_data{initial_t8_data} {};

    ~cmc_amr_data(){};
};

/* In order to write a netCDF-file of the decompressed data, the data has to be defined on the same geo domain */
#ifdef CMC_WITH_NETCDF
static void
cmc_amr_write_nc_file_decompressed_data(cmc_amr_data_t amr_data, const char* path, const std::vector<int>& var_ids)
{
    #ifdef CMC_WITH_T8CODE
    cmc_assert(amr_data->compression_applied == false);
    cmc_assert(var_ids.size() > 0);

    int ncid;
    int err;
    std::array<int, CMC_NUM_COORD_IDS +1> dim_ptrs;
    std::array<int, CMC_NUM_COORD_IDS> dimids_ptr;
    std::array<int, CMC_NUM_COORD_IDS> coordinate_vars;
    std::vector<int> nc_coord_ids(CMC_NUM_COORD_IDS +1, CMC_VAR_NOT_CONSIDERED);
    std::vector<int> nc_var_ids;
    nc_var_ids.reserve(var_ids.size());
    char ordering[20] = "";

    /* Create a new netCDF File */
    err = nc__create(path, NC_CLOBBER|NC_NETCDF4 , NC_SIZEHINT_DEFAULT, NULL, &ncid);
    cmc_nc_check_err(err);

    /* Count which dimensions are really present in the data */
    int counter{0};

    /* Define netCDF dimensions (if they are considered) and their coordinate variables */
    /** \note The 'cmc_type' is equal to the 'nc_type' */
    if (amr_data->t8_data->geo_data->get_global_coord_length(CMC_COORD_IDS::CMC_LAT) > 1)
    {
        /* Define the coordinate dimension */
        err = nc_def_dim(ncid, "lat", static_cast<size_t>(amr_data->t8_data->geo_data->get_global_coord_length(CMC_COORD_IDS::CMC_LAT)), &(dim_ptrs[CMC_COORD_IDS::CMC_LAT]));
        cmc_nc_check_err(err);
        /* Define the coordinate variable */
        err = nc_def_var(ncid, "lat", (*(amr_data->t8_data->geo_data->coords))[CMC_COORD_IDS::CMC_LAT].get_data_type(), 1, &(dim_ptrs[CMC_COORD_IDS::CMC_LAT]), &(coordinate_vars[CMC_COORD_IDS::CMC_LAT]));
        cmc_nc_check_err(err);
    }
    if (amr_data->t8_data->geo_data->get_global_coord_length(CMC_COORD_IDS::CMC_LON) > 1)
    {
        /* Define the coordinate dimension */
        err = nc_def_dim(ncid, "lon", static_cast<size_t>(amr_data->t8_data->geo_data->get_global_coord_length(CMC_COORD_IDS::CMC_LON)), &(dim_ptrs[CMC_COORD_IDS::CMC_LON]));
        cmc_nc_check_err(err);
        /* Define the coordinate variable */
        err = nc_def_var(ncid, "lon", (*(amr_data->t8_data->geo_data->coords))[CMC_COORD_IDS::CMC_LON].get_data_type(), 1, &(dim_ptrs[CMC_COORD_IDS::CMC_LON]), &(coordinate_vars[CMC_COORD_IDS::CMC_LON]));
        cmc_nc_check_err(err);
    }
    if (amr_data->t8_data->geo_data->get_global_coord_length(CMC_COORD_IDS::CMC_LEV) > 1)
    {
        /* Define the coordinate dimension */
        err = nc_def_dim(ncid, "lev", static_cast<size_t>(amr_data->t8_data->geo_data->get_global_coord_length(CMC_COORD_IDS::CMC_LEV)), &(dim_ptrs[CMC_COORD_IDS::CMC_LEV]));
        cmc_nc_check_err(err);
        /* Define the coordinate variable */
        err = nc_def_var(ncid, "lev", (*(amr_data->t8_data->geo_data->coords))[CMC_COORD_IDS::CMC_LEV].get_data_type(), 1, &(dim_ptrs[CMC_COORD_IDS::CMC_LEV]), &(coordinate_vars[CMC_COORD_IDS::CMC_LEV]));
        cmc_nc_check_err(err);
    }
    if (amr_data->t8_data->geo_data->get_global_coord_length(CMC_COORD_IDS::CMC_TIME) > 1)
    {
        /* Define the coordinate dimension */
        err = nc_def_dim(ncid, "time", static_cast<size_t>(amr_data->t8_data->geo_data->get_global_coord_length(CMC_COORD_IDS::CMC_TIME)), &(dim_ptrs[CMC_COORD_IDS::CMC_TIME]));
        cmc_nc_check_err(err);
        /* Define the coordinate variable */
        err = nc_def_var(ncid, "time", (*(amr_data->t8_data->geo_data->coords))[CMC_COORD_IDS::CMC_TIME].get_data_type(), 1, &(dim_ptrs[CMC_COORD_IDS::CMC_TIME]), &(coordinate_vars[CMC_COORD_IDS::CMC_TIME]));
        cmc_nc_check_err(err);
    }
    
    /* A flag indicating whther a new dimension for zcruve-ordering needs to be added */ 
    bool dim_num_cells_is_present{false};
    size_t num_cells{0};
    for (size_t id{0}; id < var_ids.size(); ++id)
    {
        /* Check the data scheme of the variable */
        if (amr_data->t8_data->vars[var_ids[id]]->var->data_scheme == CMC_DATA_ORDERING_SCHEME::CMC_GEO_DATA_Z_CURVE)
        {
            /* If the variable's data is in z-curve order */
            if (!dim_num_cells_is_present)
            {
                /* Save the actual diemnsion id of 'ncells' in the last position of the vector */
                nc_coord_ids[CMC_NUM_COORD_IDS] = counter;
                /* Since it is assumed (and checked), that only data on the same geo domain are written, we can use the forest of the current variable */
                num_cells = static_cast<size_t>(t8_forest_get_global_num_elements(amr_data->t8_data->vars[var_ids[id]]->assets->forest));
                err = nc_def_dim(ncid, "ncells", num_cells, &(dim_ptrs[CMC_NUM_COORD_IDS]));
                cmc_nc_check_err(err);
                /* Set the flag, that this dimension is already defined */
                dim_num_cells_is_present = true;
            }
            /* Define the variable */
            /* If the variable is in z-curve order, it is dependent on the 'ncells' dimension (since the variable is linearized) */
                err = nc_def_var(ncid, (amr_data->t8_data->vars[var_ids[id]]->var->name).c_str(), amr_data->t8_data->vars[var_ids[id]]->get_type(), 1, &dim_ptrs[CMC_NUM_COORD_IDS], &(nc_var_ids[id]));
                cmc_nc_check_err(err);
                /* Define attributes describing the ordering of the data */
                err = nc_put_att_text(ncid, nc_var_ids[id], "order", 7, "zcurve\0");
                cmc_nc_check_err(err);
                for (int i{0}; i < amr_data->t8_data->vars[var_ids[id]]->var->num_dimensions; ++i)
                {
                    if (i > 0)
                    {
                        /* Add a space between dimensions */
                        strcat(ordering, " ");
                    }
                    /* Add the abbreviation for the dimension */
                    switch ((amr_data->t8_data->vars[var_ids[id]]->var->axis_ordering)[i])
                    {
                        case CMC_COORD_IDS::CMC_LAT :
                            strcat(ordering, "lat");
                            break;
                        case CMC_COORD_IDS::CMC_LON :
                            strcat(ordering, "lon");
                            break;
                        case CMC_COORD_IDS::CMC_LEV :
                            strcat(ordering, "lev");
                            break;
                        case CMC_COORD_IDS::CMC_TIME :
                            strcat(ordering, "time");
                            break;
                    }
                }
                /* Save the precedence of the dimensions as an attribute */
                err = nc_put_att_text(ncid, nc_var_ids[id], "precedence", strlen(ordering) +1, ordering);
                cmc_nc_check_err(err);
                /* Reset the 'ordering' string */
                ordering[0] = '\0';
        }
        else if (amr_data->t8_data->vars[var_ids[id]]->var->data_scheme == CMC_DATA_ORDERING_SCHEME::CMC_GEO_DATA_LINEAR)
        {
            /* Define the variable */
            /* If the data is linearily ordered, it depends on the coordinate variables */
            /* Assign the right axis ordering */
            for (int i{0}; i < amr_data->t8_data->vars[var_ids[id]]->var->num_dimensions; ++i)
            {
                dimids_ptr[i] = dim_ptrs[(amr_data->t8_data->vars[var_ids[id]]->var->axis_ordering)[i]];
            }
            /* Define the variable */
            err = nc_def_var(ncid, (amr_data->t8_data->vars[var_ids[id]]->var->name).c_str(), amr_data->t8_data->vars[var_ids[id]]->get_type(), amr_data->t8_data->vars[var_ids[id]]->var->num_dimensions, &dimids_ptr[0], &(nc_var_ids[id]));
            cmc_nc_check_err(err);
        }
        else
        {
            cmc_err_msg("The data scheme of the variable is either undefined or unsupported to be written as a netCDF file.");
        }

        /* Check if missing values are present */
        if (amr_data->t8_data->vars[var_ids[id]]->var->missing_value_present)
        {
            switch (amr_data->t8_data->vars[var_ids[id]]->get_type())
            {
                case cmc_type::CMC_INT8_T :
                    {
                        int8_t missing_val = cmc_get_universal_data<int8_t>(amr_data->t8_data->vars[var_ids[id]]->var->missing_value);
                        err = nc_put_att(ncid, nc_var_ids[id], "missing_value", NC_BYTE, 1, &missing_val);
                        cmc_nc_check_err(err);
                    }
                    break;
                case cmc_type::CMC_CHAR :
                    {
                        char missing_val = cmc_get_universal_data<char>(amr_data->t8_data->vars[var_ids[id]]->var->missing_value);
                        err = nc_put_att(ncid, nc_var_ids[id], "missing_value", NC_CHAR, 1, &missing_val);
                        cmc_nc_check_err(err);
                    }
                    break;
                case cmc_type::CMC_INT16_T :
                    {
                        int16_t missing_val = cmc_get_universal_data<int16_t>(amr_data->t8_data->vars[var_ids[id]]->var->missing_value);
                        err = nc_put_att(ncid, nc_var_ids[id], "missing_value", NC_SHORT, 1, &missing_val);
                        cmc_nc_check_err(err);
                    }
                    break;
                case cmc_type::CMC_INT32_T :
                    {
                        int32_t missing_val = cmc_get_universal_data<int32_t>(amr_data->t8_data->vars[var_ids[id]]->var->missing_value);
                        err = nc_put_att(ncid, nc_var_ids[id], "missing_value", NC_INT, 1, &missing_val);
                        cmc_nc_check_err(err);
                    }
                    break;
                case cmc_type::CMC_FLOAT :
                    {
                        float missing_val = cmc_get_universal_data<float>(amr_data->t8_data->vars[var_ids[id]]->var->missing_value);
                        err = nc_put_att(ncid, nc_var_ids[id], "missing_value", NC_FLOAT, 1, &missing_val);
                        cmc_nc_check_err(err);
                    }
                    break;
                case cmc_type::CMC_DOUBLE :
                    {
                        double missing_val = cmc_get_universal_data<double>(amr_data->t8_data->vars[var_ids[id]]->var->missing_value);
                        err = nc_put_att(ncid, nc_var_ids[id], "missing_value", NC_DOUBLE, 1, &missing_val);
                        cmc_nc_check_err(err);
                    }
                    break;
                case cmc_type::CMC_UINT8_T :
                    {
                        uint8_t missing_val = cmc_get_universal_data<uint8_t>(amr_data->t8_data->vars[var_ids[id]]->var->missing_value);
                        err = nc_put_att(ncid, nc_var_ids[id], "missing_value", NC_UBYTE, 1, &missing_val);
                        cmc_nc_check_err(err);
                    }
                    break;
                case cmc_type::CMC_UINT16_T :
                    {
                        uint16_t missing_val = cmc_get_universal_data<uint16_t>(amr_data->t8_data->vars[var_ids[id]]->var->missing_value);
                        err = nc_put_att(ncid, nc_var_ids[id], "missing_value", NC_USHORT, 1, &missing_val);
                        cmc_nc_check_err(err);
                    }
                    break;
                case cmc_type::CMC_UINT32_T :
                    {
                        uint32_t missing_val = cmc_get_universal_data<uint32_t>(amr_data->t8_data->vars[var_ids[id]]->var->missing_value);
                        err = nc_put_att(ncid, nc_var_ids[id], "missing_value", NC_UINT, 1, &missing_val);
                        cmc_nc_check_err(err);
                    }
                    break;
                case cmc_type::CMC_INT64_T :
                    {
                        int64_t missing_val = cmc_get_universal_data<int64_t>(amr_data->t8_data->vars[var_ids[id]]->var->missing_value);
                        err = nc_put_att(ncid, nc_var_ids[id], "missing_value", NC_INT64, 1, &missing_val);
                        cmc_nc_check_err(err);
                    }
                    break;
                case cmc_type::CMC_UINT64_T :
                    {
                        uint64_t missing_val = cmc_get_universal_data<uint64_t>(amr_data->t8_data->vars[var_ids[id]]->var->missing_value);
                        err = nc_put_att(ncid, nc_var_ids[id], "missing_value", NC_UINT64, 1, &missing_val);
                        cmc_nc_check_err(err);
                    }
                    break;
                default:
                    cmc_err_msg("The missing value of the variable (name: ", amr_data->t8_data->vars[var_ids[id]]->var->name, "contains an unknown attribute type.");
            }
        }
    }

    /* All dimensions and variables have been defined */
    /* Therefore, we are leaving the 'define-mode' and switch to the data mode */
    err = nc_enddef(ncid);
    cmc_nc_check_err(err);

    /* Write the data of all coordinate variables to the file */
    if (amr_data->t8_data->geo_data->get_global_coord_length(CMC_COORD_IDS::CMC_LAT) > 1)
    {
        err = nc_put_var(ncid, coordinate_vars[CMC_COORD_IDS::CMC_LAT], (*(amr_data->t8_data->geo_data->coords))[CMC_COORD_IDS::CMC_LAT].get_initial_data_ptr());
        cmc_nc_check_err(err);
    }
    if (amr_data->t8_data->geo_data->get_global_coord_length(CMC_COORD_IDS::CMC_LON) > 1)
    {
        err = nc_put_var(ncid, coordinate_vars[CMC_COORD_IDS::CMC_LON], (*(amr_data->t8_data->geo_data->coords))[CMC_COORD_IDS::CMC_LON].get_initial_data_ptr());
        cmc_nc_check_err(err);
    }
    if (amr_data->t8_data->geo_data->get_global_coord_length(CMC_COORD_IDS::CMC_LEV) > 1)
    {
        err = nc_put_var(ncid, coordinate_vars[CMC_COORD_IDS::CMC_LEV], (*(amr_data->t8_data->geo_data->coords))[CMC_COORD_IDS::CMC_LEV].get_initial_data_ptr());
        cmc_nc_check_err(err);
    }
    if (amr_data->t8_data->geo_data->get_global_coord_length(CMC_COORD_IDS::CMC_TIME) > 1)
    {
        err = nc_put_var(ncid, coordinate_vars[CMC_COORD_IDS::CMC_TIME], (*(amr_data->t8_data->geo_data->coords))[CMC_COORD_IDS::CMC_TIME].get_initial_data_ptr());
        cmc_nc_check_err(err);
    }

    /* Write the data of all other variables to the file */
    for (size_t id{0}; id < var_ids.size(); ++id)
    {
        err = nc_put_var(ncid, nc_var_ids[id], amr_data->t8_data->vars[var_ids[id]]->get_initial_data_ptr());
        cmc_nc_check_err(err);
    }
    /* All data has been written. Therefore, the file may be closed */
    err = nc_close(ncid);
    cmc_nc_check_err(err);

    #endif
}
#endif
/** End Static Functions **/

cmc_amr_data_t
cmc_create_amr_compression_data(cmc_nc_data_t nc_data, const MPI_Comm comm)
{
    cmc_amr_data_t amr_data{new cmc_amr_data()};
    amr_data->t8_data = new cmc_t8_data(comm);
    cmc_t8_nc_setup_compression(nc_data, amr_data->t8_data);
    return amr_data;
}

#if 0
//TODO: update messy function
cmc_amr_data_t
cmc_create_amr_compression_data_messy(cmc_messy_data_t messy_data, const MPI_Comm comm)
{
    #ifdef CMC_WITH_T8CODE
    cmc_amr_data_t amr_data{new cmc_amr_data()};
    amr_data->t8_data = new cmc_t8_data{comm};
    cmc_t8_messy_setup_compression(messy_data, amr_data->t8_data);
    return amr_data;
    #endif
}
#endif

/**
 * @brief This function sets a flag that the initial dta will be kept until the end of the compression. If it is set and the decompression
 * is called, the introduced error (from the compression) will be calculated.
 * 
 * @param amr_data The amr data struct pointer holding the context of the comrpession
 * @param flag_keep_initial_data The flag indicating whther tot not to keep the initial data
 */
void
cmc_amr_pre_setup_set_flag_in_order_to_keep_the_initial_data(cmc_amr_data_t amr_data, const int flag_keep_initial_data)
{
    if (flag_keep_initial_data != 0)
    {
        /* If true, we will keep the data */
        amr_data->t8_data->is_initial_data_kept = true;
    } else
    {
        /* if the flag is false, we will abandon the initial data during the compression */
        amr_data->t8_data->is_initial_data_kept = false;
    }
}

/**
 * @brief This functions sets the information that a certain variable (with id @var var_id) will be split before the compression
 * accordingly to the given @var preferred_layout
 * 
 * @param t8_data The pointer to the @struct cmc_t8_data gholding the variables and the information about the compression
 * @param var_id The variable id corersponding to the variable which is ought to be split
 * @param preferred_layout The preferred layout into which the variable will be transformed
 *
 * @note Passing 'CMC_APPLY_TRANSFORMATION_3D_TO_2D_TO_ALL_VARS' as @var var_id will transform all variables compliant to the given layout
 */
void
cmc_amr_pre_setup_split_3D_variable(cmc_amr_data_t amr_data, const int var_id, const DATA_LAYOUT preferred_data_layout)
{
    #ifdef CMC_WITH_T8CODE
    cmc_t8_geo_data_set_split_variable(amr_data->t8_data, var_id, preferred_data_layout);
    #endif
}

void
cmc_amr_pre_setup_set_compression_criterium_relative_error_threshold(cmc_amr_data_t amr_data, const double maximum_error_tolerance)
{
    #ifdef CMC_WITH_T8CODE
    /* Set the error threshold in t8_data */
    cmc_t8_geo_data_set_relative_error_criterium(amr_data->t8_data, maximum_error_tolerance);
    #endif
}

void
cmc_amr_pre_setup_set_compression_criterium_absolute_error_threshold(cmc_amr_data_t amr_data, const double maximum_abs_error_tolerance)
{
    #ifdef CMC_WITH_T8CODE
    /* Set the absolute error criterion */
    cmc_t8_geo_data_set_absolute_error_criterium(amr_data->t8_data, maximum_abs_error_tolerance);
    #endif
}

void
cmc_amr_pre_setup_reset_compression_settings(cmc_amr_data_t amr_data)
{
    #ifdef CMC_WITH_T8CODE
    cmc_t8_geo_data_reset_compression_settigs(amr_data->t8_data);
    #endif
}

#if __cplusplus
void
cmc_amr_pre_setup_set_compression_criterium_exclude_area(cmc_amr_data_t amr_data, const CMC_COORD_IDS coord_id, const cmc_universal_type_t& start_value, const cmc_universal_type_t& end_value)
{
    #ifdef CMC_WITH_T8CODE
    /* Set the exclusion area, in which no coarsening will be applied */
    cmc_t8_geo_data_set_exclude_area(amr_data->t8_data, coord_id, start_value, end_value);
    #endif
}
#else
void
cmc_amr_pre_setup_set_compression_criterium_exclude_area(cmc_amr_data_t amr_data, const enum CMC_COORD_IDS coord_id, const double start_value, const double end_value)
{
    #ifdef CMC_WITH_T8CODE
    /* Transform the start and end values to their corresponding universal_type */
    cmc_universal_type_t start_val = convert_to_universal_type(amr_data->t8_data->geo_data->coords->operator[](static_cast<size_t>(coord_id)).get_data_type(), start_value);
    cmc_universal_type_t end_val = convert_to_universal_type(amr_data->t8_data->geo_data->coords->operator[](static_cast<size_t>(coord_id)).get_data_type(), end_value);
    /* Set the exclusion area, in which no coarsening will be applied */
    cmc_t8_geo_data_set_exclude_area(amr_data->t8_data, coord_id, start_value, end_value);
    #endif
}
#endif

/** Create an initial mesh which will be able to contain the (simulation) data */          
void
cmc_amr_setup_compression(cmc_amr_data_t amr_data, CMC_AMR_COMPRESSION_MODE compression_mode)
{
    #ifdef CMC_WITH_T8CODE
    #ifdef CMC_ENABLE_DEBUG
    if (amr_data->t8_data->compression_mode != CMC_AMR_COMPRESSION_MODE::CMC_T8_COMPRESSION_UNDEFINED
        && amr_data->t8_data->compression_mode != compression_mode)
    {
        cmc_warn_msg("Eventually, the (t8_)variables' data already have been manipulated assuming a different compression approach as now supplied via cmc_amr_compress().");
    }
    #endif

    /* Apply offset and scaling to the variables (curently, this only applies in a netCDF case) */
    cmc_t8_apply_offset_and_scaling(amr_data->t8_data, CMC_APPLY_OFFSET_AND_SCALING_TO_ALL_VARS);

    /* Save the supplied compression approach */
    amr_data->t8_data->compression_mode = compression_mode;

    /* Perform different actions based on the compression approach */
    switch (amr_data->t8_data->compression_mode)
    {
        case CMC_AMR_COMPRESSION_MODE::ONE_FOR_ALL:
            cmc_debug_msg("Lossy AMR 'One for All' compression will be applied...");
            break;
        case CMC_AMR_COMPRESSION_MODE::ONE_FOR_ONE:
            cmc_debug_msg("Lossy AMR 'One for One' compression will be applied...");
            break;
        default:
            cmc_err_msg("An unknown compression_mode was supplied.");
    }

    /* Save the diemensionality of the data */
    amr_data->t8_data->dimension_of_compression_mode = amr_data->t8_data->geo_data->dim;

    /* Check if variables are ought to be split at a certain dimension before the mesh is built and the data is redistributed */
    cmc_t8_geo_data_split_variables(amr_data->t8_data);

    /* Build an enclosing mesh */
    cmc_t8_create_enclosing_geo_mesh(*(amr_data->t8_data));

    /* Apply the Morton order to the data and (in a parallel case,) distribute the data accordingly beforehand */
    cmc_t8_geo_data_distribute_and_apply_ordering((amr_data->t8_data));

    #ifdef CMC_ENABLE_DEBUG
    /* In a DEBUG-configuration we are writing the initial data to '.vtu'-files */
    char file_prefix[55];

    /* Create a file-name dependent on the compression mode */
    switch (compression_mode)
    {
        case CMC_AMR_COMPRESSION_MODE::ONE_FOR_ALL:
        {
            int ret_val{snprintf(file_prefix, 49, "%s%d%s", "cmc_amr_lossy_comp-Initial_Forest_One_For_All_", amr_data->t8_data->dimension_of_compression_mode, "D")};
            if (ret_val < 0)
            {
                cmc_err_msg("An error occured during the creation of a file name for the vtk-output file.");
            }
        }
        break;
        case CMC_AMR_COMPRESSION_MODE::ONE_FOR_ONE:
        {
            int ret_val{snprintf(file_prefix, 49, "%s%d%s", "cmc_amr_lossy_comp-Initial_Forest_One_For_One_", amr_data->t8_data->dimension_of_compression_mode, "D")};
            if (ret_val < 0)
            {
                cmc_err_msg("An error occured during the creation of a file name for the vtk-output file.");
            }
        }
        break;
        default:
            cmc_err_msg("An unknown compression mode ist present.");
    }
    
    /* Write out the forest with the variable's data at the beginning of the compression process */
    cmc_t8_write_forest_all_vars(amr_data->t8_data, file_prefix); 
    #endif
    #endif
}

void
cmc_amr_compress(cmc_amr_data_t amr_data, const t8_forest_adapt_t adapt_function, cmc_t8_forest_interpolation_t interpolation_function)
{
    #ifdef CMC_WITH_T8CODE
    cmc_debug_msg("Compression starts...");

    /* Perform different actions based on the compression approach */
    if (adapt_function != nullptr && interpolation_function != nullptr)
    {
        /* Perform the adaptation/compression with the supplied adapt and replace functions */
        cmc_t8_coarsen_data(amr_data->t8_data, adapt_function, interpolation_function);
    } else {
        switch (amr_data->t8_data->settings.compression_criterium)
        {
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_CRITERIUM_UNDEFINED:
                cmc_err_msg("No compression criterion has been specified. Please call one of the 'cmc_amr_pre_setup_set_...()' functions earlier.");
            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_REL_ERROR_THRESHOLD:
                /* In case a relative error criterion has been chosen */
                cmc_t8_coarsen_data(amr_data->t8_data, cmc_t8_adapt_callback_coarsen_based_on_relative_error_threshold, cmc_t8_compression_interpolation_standard_mean);
            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_ABS_ERROR_THRESHOLD:
                /* In case a relative error criterion has been chosen */
                cmc_t8_coarsen_data(amr_data->t8_data, cmc_t8_adapt_callback_coarsen_based_on_absolute_error_threshold, cmc_t8_compression_interpolation_standard_mean);
            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_EXCLUDE_AREA:
                /* In case an exclude area crierion has been chosen */
                //cmc_t8_coarsen_data(amr_data->t8_data, cmc_t8_adapt_callback_coarsen_exclude_area, cmc_t8_geo_data_interpolate_std_mean);
            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_COMBINED_CRITERION:
                /* In case several compression criteria has been chosen */
                //cmc_t8_coarsen_data(amr_data->t8_data, cmc_t8_adapt_callback_coarsen_combined_criteria, cmc_t8_geo_data_interpolate_error_threshold_adaption);
            break;
            default :
                cmc_err_msg("An unknown compression criterion has bee supplied.");
        }
    }
    cmc_debug_msg("The Lossy AMR compressor has been successfully applied.");
    
    /* Set the compression flag */
    amr_data->compression_applied = true;

    #ifdef CMC_ENABLE_DEBUG
    /* In a DEBUG-configuration we are writing out the comrpessed data */
    char file_prefix[55];

    switch (amr_data->t8_data->compression_mode)
    {
        case CMC_AMR_COMPRESSION_MODE::ONE_FOR_ALL:
        {
            int ret_val{snprintf(file_prefix, 52, "%s%d%s", "cmc_amr_lossy_comp-Compressed_Forest_One_For_All_", amr_data->t8_data->dimension_of_compression_mode, "D")};
            if (ret_val < 0)
            {
                cmc_err_msg("An error occured during the creation of a file name for the vtk-output file.");
            }
        }
        break;
        case CMC_AMR_COMPRESSION_MODE::ONE_FOR_ONE:
        {
            int ret_val{snprintf(file_prefix, 52, "%s%d%s", "cmc_amr_lossy_comp-Compressed_Forest_One_For_One_", amr_data->t8_data->dimension_of_compression_mode, "D")};
            if (ret_val < 0)
            {
                cmc_err_msg("An error occured during the creation of a file name for the vtk-output file.");
            }
        }
        break;
        default:
            cmc_err_msg("An unknown compression mode ist present.");
    }

    #if 0
    /* Save the amount of compressed data points in a file */
    std::ofstream compresseion_results("compression_result_element_count_per_variable.txt");
    compresseion_results << "Elements of the compressed forests\nMax error = " << amr_data->t8_data->settings.max_err << std::endl;
    for (size_t i{0}; i < amr_data->t8_data->vars.size(); ++i)
    {
        compresseion_results << t8_forest_get_global_num_elements(amr_data->t8_data->vars[i]->assets->forest) << " Elements has the forest of variable (name: " << amr_data->t8_data->vars[i]->var->name << ")" << std::endl; 
    }
    compresseion_results.close();

    /* Write out the forest(s) with the compressed data */
    cmc_t8_write_forest_all_vars(amr_data->t8_data, file_prefix);
    #endif

    #endif
    #endif
}

/**
 * @brief Calculate the introduced data inaccuracy after the decompression if the initial data has been kept
 * 
 * @param amr_data The amr_data holding the decompressed data as well as the intial data
 */
static void
cmc_amr_check_inaccuracy_after_decompression(cmc_amr_data_t amr_data)
{
    cmc_assert(!amr_data->compression_applied);
    cmc_assert(amr_data->t8_data->is_initial_data_kept);

    /* Declare some variables needed for the computation */
    double max_err{0.0}, current_err{0.0};
    double approx_value{0.0};
    double exact_value{1.0};

    /* If an absolute error threshold has been used, we will check the absolute error, otherwise we will check the introduced relative error */
    if (amr_data->t8_data->settings.compression_criterium == CMC_T8_COMPRESSION_CRITERIUM::CMC_ABS_ERROR_THRESHOLD)
    {
        cmc_debug_msg("The Lossy AMR compressor introduced a maximum point-wise absolute error of:");
    } else
    {
        cmc_debug_msg("The Lossy AMR compressor introduced a maximum point-wise relative error of:");
    }

    /* Check which mode is used */
    if (amr_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL)
    {
        std::vector<double> max_error(amr_data->t8_data->vars.size(), 0.0);
        
        for (size_t i{0}; i < static_cast<size_t>(t8_forest_get_local_num_elements(amr_data->t8_data->assets->forest)); ++i)
        {
            for (size_t var_id{0}; var_id < amr_data->t8_data->vars.size(); ++var_id)
            {
                exact_value = std::abs(cmc_get_universal_data<double>((amr_data->t8_data->initial_data[var_id]).operator[](i)));
                approx_value = std::abs(cmc_get_universal_data<double>(amr_data->t8_data->vars[var_id]->var->data->operator[](i)));
                /* Calculate the new relative error */
                if (amr_data->t8_data->settings.compression_criterium == CMC_T8_COMPRESSION_CRITERIUM::CMC_ABS_ERROR_THRESHOLD)
                {
                    /* Calculate the new absolute error */
                    current_err = std::abs(approx_value - exact_value);
                    /* Eventually update the maximim error for this variable */
                    if (max_error[var_id] < current_err)
                    {
                        max_error[var_id] = current_err;
                    }
                } else
                {
                    if (exact_value > 0.0 && approx_value > 0.0)
                    {
                        /* Calculate the new relative error */
                        current_err = std::abs((approx_value - exact_value) / exact_value);

                        /* Eventually update the maximim error for this variable */
                        if (max_error[var_id] < current_err)
                        {
                            max_error[var_id] = current_err;
                        }
                    }
                }
            }
        }
        for (size_t var_id{0}; var_id < amr_data->t8_data->vars.size(); ++var_id)
        {
            if (amr_data->t8_data->settings.compression_criterium == CMC_T8_COMPRESSION_CRITERIUM::CMC_ABS_ERROR_THRESHOLD)
            {
                cmc_debug_msg("\tVariable (name: ", amr_data->t8_data->vars[var_id]->var->name, ") has a maximum absolute data inaccuracy of ", max_error[var_id], ".");
            } else
            {
                cmc_debug_msg("\tVariable (name: ", amr_data->t8_data->vars[var_id]->var->name, ") has a maximum relative data inaccuracy of ", max_error[var_id] * 100, "%.");
            }
        }
    }
    else
    {
        /* In a case of a previous 'One for One' compression */
        /* Iterate over all variables and check their introduced relative data inaccuracy */
        for (size_t var_id{0}; var_id < amr_data->t8_data->vars.size(); ++var_id)
        {
            max_err = 0.0;
            cmc_assert((amr_data->t8_data->initial_data[var_id]).size() == amr_data->t8_data->vars[var_id]->var->data->size());

            for (size_t i{0}; i < static_cast<size_t>(t8_forest_get_local_num_elements(amr_data->t8_data->vars[var_id]->assets->forest)); ++i)
            {
                /* Check if a missing value would be at the position of the exact value */
                if (!(amr_data->t8_data->vars[var_id]->var->is_equal_to_missing_value(i)))
                {
                    exact_value = std::abs(cmc_get_universal_data<double>((amr_data->t8_data->initial_data[var_id]).operator[](i)));
                    approx_value = std::abs(cmc_get_universal_data<double>(amr_data->t8_data->vars[var_id]->var->data->operator[](i)));
                    if (amr_data->t8_data->settings.compression_criterium == CMC_T8_COMPRESSION_CRITERIUM::CMC_ABS_ERROR_THRESHOLD)
                    {
                        /* Calculate the new absolute error */
                        current_err = std::abs(approx_value - exact_value);

                        /* Eventually update the maximim error for this variable */
                        if (max_err < current_err)
                        {
                            max_err = current_err;
                        }
                    } else
                    {
                        if (exact_value > 0.0 && approx_value > 0.0)
                        {
                            /* Calculate the new relative error */
                            current_err = std::abs((approx_value - exact_value) / exact_value);

                            /* Eventually update the maximim error for this variable */
                            if (max_err < current_err)
                            {
                                max_err = current_err;
                            }
                        }
                    }
                }
            }
            if (amr_data->t8_data->settings.compression_criterium == CMC_T8_COMPRESSION_CRITERIUM::CMC_ABS_ERROR_THRESHOLD)
            {
                cmc_debug_msg("\tVariable (name: ", amr_data->t8_data->vars[var_id]->var->name, ") has a maximum absolute data inaccuracy of ", max_err, ".");
            } else
            {
                cmc_debug_msg("\tVariable (name: ", amr_data->t8_data->vars[var_id]->var->name, ") has a maximum relative data inaccuracy of ", max_err * 100, "%.");
            }
        }
    }

}

void
cmc_amr_decompress(cmc_amr_data_t amr_data)
{
    #ifdef CMC_WITH_T8CODE
    cmc_assert(amr_data->compression_applied);

    cmc_debug_msg("Decompression of data starts...");

    /* Refine the data of all variables to corresponding intial mesh */
    cmc_t8_refine_to_initial_level(amr_data->t8_data);

    cmc_debug_msg("Decompression has been finished.");
    
    /* Set the (de-)comrpession flag */
    amr_data->compression_applied = false;

    /* If the initial data is kept, we can calculate the maximum introduced data inaccuracy */
    if (amr_data->t8_data->is_initial_data_kept)
    {
        cmc_amr_check_inaccuracy_after_decompression(amr_data);
    }

    #endif
}

void
cmc_amr_destroy(cmc_amr_data_t amr_data)
{
    if (amr_data != nullptr)
    {
        delete amr_data->t8_data;
        delete amr_data;
    }
}


void
cmc_amr_write_netcdf_file(cmc_amr_data_t amr_data, const char* path, const int var_id)
{
    #ifdef CMC_WITH_NETCDF
    std::vector<int> output_variable_ids;
    if (var_id == CMC_AMR_WRITE_ALL_VARS_TO_NETCDF)
    {
        /* Allocate space for the variables */
        output_variable_ids.reserve(amr_data->t8_data->vars.size());
        /* Iterate over all variables */
        for (size_t id{0}; id < amr_data->t8_data->vars.size(); ++id)
        /* Variable's data needs to be in linear or zcurve order */
        if(amr_data->t8_data->vars[id]->var->data_scheme == CMC_DATA_ORDERING_SCHEME::CMC_GEO_DATA_LINEAR ||
           amr_data->t8_data->vars[id]->var->data_scheme == CMC_DATA_ORDERING_SCHEME::CMC_GEO_DATA_Z_CURVE)
        {
            output_variable_ids.push_back(id);
        }
    } else
    {
        /* If only one variable should be written in a netCDF file */
        cmc_assert(var_id >= 0 && var_id < static_cast<int>(amr_data->t8_data->vars.size()));
        /* Onlx save the single id */
        output_variable_ids.push_back(var_id);
    }

    /* Check if compressed or not compressed data is present */
    if (amr_data->compression_applied)
    {
        /* If we have to handle compressed data */
        cmc_err_msg("Writing out compressed data is currently not supported, since it cannot be used without information about the underlying forest...");
    }
    else
    {
        /* If the data is not compressed */
        cmc_amr_write_nc_file_decompressed_data(amr_data, path, output_variable_ids);
    }
    #endif
}


void
cmc_amr_write_vtk_file(cmc_amr_data_t amr_data, const char* file_prefix)
{
    #ifdef CMC_WITH_T8CODE
    cmc_t8_write_forest_all_vars(amr_data->t8_data, file_prefix);
    #endif
}
