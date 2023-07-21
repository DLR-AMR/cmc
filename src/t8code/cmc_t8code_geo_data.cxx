#include "cmc_t8code_geo_data.h"
#include "cmc_t8code_data.hxx"
#include "cmc_t8code_geo_mesh.h"
#include "cmc_t8_adapt_callbacks.h"
#include "cmc_t8_replace_callbacks.h"
#include "utilities/cmc_log_functions.h"
#include "utilities/cmc_container.h"

/** Begin STATIC Functions **/
/****************************/

#ifdef CMC_WITH_T8CODE
/** \note: These functions are considered to be used when linear ordered data is transformed to z-curve compliant data.
  * \note: It is assumed that longitude coordinates are ordered like: lon=0.0, 2.5, 5.0, ..., 357.5
  *        It is assumed that latitude coordinates are ordered like:  lat=90.0, 87.5, 85, ..., -90.0
  *        It is assumed that elevation coordinates are ordered like: lev=0,1,2,3,.... (increasingly ordered)
  *
  * \note: Different geo-spatial coordinate combinations are corresponding to certain x- and y- coordinates. This labeling is consistent througout the 'cmc_t8_...'-functions. 
  *        Data only concerning latitude and longitude will be reordered, s.t. longitude equals x, latitude equals y
  *        Data only concerning latitude and elevation will be reordered, s.t. latitude equals x, elevation equals y
  *        Data only concerning longitude and elevation will be reordered, s.t. longitude equals x, elevation equals y

  * \note: After the data of the variable has been reordered compliant to a z-curve ordering scheme
**/
/* 2D offset functions */
/* Data concerning latitude and longitude will be reordered, s.t. longitude equals x, latitude equals y */
static size_t
cmc_t8_calc_lin_order_offset_for_zcurve_lon_lat(const std::vector<size_t>& dim_lengths, const int x_coord, const int y_coord, const int z_coord)
{
    #ifdef CMC_WITH_T8CODE
    /* In this case the data layout equals 'CMC_2D_LON_LAT' */
    return x_coord * dim_lengths[CMC_COORD_IDS::CMC_LAT] + (dim_lengths[CMC_COORD_IDS::CMC_LAT] -1 -y_coord);
    #else
    return CMC_ERR;
    #endif
}
/* Data concerning latitude and longitude will be reordered, s.t. longitude equals x, latitude equals y */
static size_t
cmc_t8_calc_lin_order_offset_for_zcurve_lat_lon(const std::vector<size_t>& dim_lengths, const int x_coord, const int y_coord, const int z_coord)
{
    #ifdef CMC_WITH_T8CODE
    /* In this case the data layout equals 'CMC_2D_LAT_LON' */
    return x_coord + (dim_lengths[CMC_COORD_IDS::CMC_LAT] -1 - y_coord) * dim_lengths[CMC_COORD_IDS::CMC_LON];
    #else
    return CMC_ERR;
    #endif
}
/* Data concerning latitude and longitude will be reordered, s.t. latitude equals x, elevation equals y */
static size_t
cmc_t8_calc_lin_order_offset_for_zcurve_lat_lev(const std::vector<size_t>& dim_lengths, const int x_coord, const int y_coord, const int z_coord)
{
    #ifdef CMC_WITH_T8CODE
    /* In this case the data layout equals 'CMC_2D_LAT_LEV' */
    return y_coord + (dim_lengths[CMC_COORD_IDS::CMC_LAT] -1 - x_coord) * dim_lengths[CMC_COORD_IDS::CMC_LEV];
    #else
    return CMC_ERR;
    #endif
}
/* Data concerning latitude and longitude will be reordered, s.t. latitude equals x, elevation equals y */
static size_t
cmc_t8_calc_lin_order_offset_for_zcurve_lev_lat(const std::vector<size_t>& dim_lengths, const int x_coord, const int y_coord, const int z_coord)
{
    #ifdef CMC_WITH_T8CODE
    /* In this case the data layout equals 'CMC_2D_LEV_LAT' */
    return (dim_lengths[CMC_COORD_IDS::CMC_LAT] -1 -x_coord) + y_coord * dim_lengths[CMC_COORD_IDS::CMC_LAT];
    #else
    return CMC_ERR;
    #endif
}
/* Data concerning latitude and longitude will be reordered, s.t. longitude equals x, elevation equals y */
static size_t
cmc_t8_calc_lin_order_offset_for_zcurve_lon_lev(const std::vector<size_t>& dim_lengths, const int x_coord, const int y_coord, const int z_coord)
{
    #ifdef CMC_WITH_T8CODE
    /* In this case the data layout equals 'CMC_2D_LON_LEV' */
    return y_coord + x_coord * dim_lengths[CMC_COORD_IDS::CMC_LEV];
    #else
    return CMC_ERR;
    #endif
}
/* Data concerning latitude and longitude will be reordered, s.t. longitude equals x, elevation equals y */
static size_t
cmc_t8_calc_lin_order_offset_for_zcurve_lev_lon(const std::vector<size_t>& dim_lengths, const int x_coord, const int y_coord, const int z_coord)
{
    #ifdef CMC_WITH_T8CODE
    /* In this case the data layout equals 'CMC_2D_LEV_LON' */
    return x_coord + y_coord * dim_lengths[CMC_COORD_IDS::CMC_LON];
    #else
    return CMC_ERR;
    #endif
}
/* 3D offset functions */
/* Data concerning latitude, longitude and elevation will be reordered, s.t. longitude equals x, latitude equals y and elevation equals z */
static size_t
cmc_t8_calc_lin_order_offset_for_zcurve_lon_lat_lev(const std::vector<size_t>& dim_lengths, const int x_coord, const int y_coord, const int z_coord)
{
    #ifdef CMC_WITH_T8CODE
    /* In this case the data layout equals 'CMC_3D_LON_LAT_LEV' */
    return z_coord + dim_lengths[CMC_COORD_IDS::CMC_LEV] * (dim_lengths[CMC_COORD_IDS::CMC_LAT] -1 -y_coord + x_coord * dim_lengths[CMC_COORD_IDS::CMC_LAT]);
    #else
    return CMC_ERR;
    #endif
}
/* Data concerning latitude, longitude and elevation will be reordered, s.t. longitude equals x, latitude equals y and elevation equals z */
static size_t
cmc_t8_calc_lin_order_offset_for_zcurve_lon_lev_lat(const std::vector<size_t>& dim_lengths, const int x_coord, const int y_coord, const int z_coord)
{
    #ifdef CMC_WITH_T8CODE
    /* In this case the data layout equals 'CMC_3D_LON_LEV_LAT' */
    return (dim_lengths[CMC_COORD_IDS::CMC_LAT] -1 -y_coord) + dim_lengths[CMC_COORD_IDS::CMC_LAT] * (z_coord + x_coord * dim_lengths[CMC_COORD_IDS::CMC_LEV]);
    #else
    return CMC_ERR;
    #endif
}
/* Data concerning latitude, longitude and elevation will be reordered, s.t. longitude equals x, latitude equals y and elevation equals z */
static size_t
cmc_t8_calc_lin_order_offset_for_zcurve_lev_lon_lat(const std::vector<size_t>& dim_lengths, const int x_coord, const int y_coord, const int z_coord)
{
    #ifdef CMC_WITH_T8CODE
    /* In this case the data layout equals 'CMC_3D_LEV_LON_LAT' */
    return (dim_lengths[CMC_COORD_IDS::CMC_LAT] -1 -y_coord) + dim_lengths[CMC_COORD_IDS::CMC_LAT] * (x_coord + z_coord * dim_lengths[CMC_COORD_IDS::CMC_LON]);
    #else
    return CMC_ERR;
    #endif
}
/* Data concerning latitude, longitude and elevation will be reordered, s.t. longitude equals x, latitude equals y and elevation equals z */
static size_t
cmc_t8_calc_lin_order_offset_for_zcurve_lev_lat_lon(const std::vector<size_t>& dim_lengths, const int x_coord, const int y_coord, const int z_coord)
{
    #ifdef CMC_WITH_T8CODE
    /* In this case the data layout equals 'CMC_3D_LEV_LAT_LON' */
    return x_coord + dim_lengths[CMC_COORD_IDS::CMC_LON] * (dim_lengths[CMC_COORD_IDS::CMC_LAT] -1 -y_coord + z_coord * dim_lengths[CMC_COORD_IDS::CMC_LAT]);
    #else
    return CMC_ERR;
    #endif
}
/* Data concerning latitude, longitude and elevation will be reordered, s.t. longitude equals x, latitude equals y and elevation equals z */
static size_t
cmc_t8_calc_lin_order_offset_for_zcurve_lat_lev_lon(const std::vector<size_t>& dim_lengths, const int x_coord, const int y_coord, const int z_coord)
{
    #ifdef CMC_WITH_T8CODE
    /* In this case the data layout equals 'CMC_3D_LAT_LEV_LON' */
    return x_coord + dim_lengths[CMC_COORD_IDS::CMC_LON] * (z_coord + (dim_lengths[CMC_COORD_IDS::CMC_LAT] -1 -y_coord) * dim_lengths[CMC_COORD_IDS::CMC_LEV]);
    #else
    return CMC_ERR;
    #endif
}
/* Data concerning latitude, longitude and elevation will be reordered, s.t. longitude equals x, latitude equals y and elevation equals z */
static size_t
cmc_t8_calc_lin_order_offset_for_zcurve_lat_lon_lev(const std::vector<size_t>& dim_lengths, const int x_coord, const int y_coord, const int z_coord)
{
    #ifdef CMC_WITH_T8CODE
    /* In this case the data layout equals 'CMC_3D_LAT_LON_LEV' */
    return z_coord + dim_lengths[CMC_COORD_IDS::CMC_LEV] * (x_coord + (dim_lengths[CMC_COORD_IDS::CMC_LAT] -1 -y_coord) * dim_lengths[CMC_COORD_IDS::CMC_LON]);
    #else
    return CMC_ERR;
    #endif
}
/* A generic error offset function if the deduction of the correct offset function has failed */
static size_t
cmc_t8_calc_lin_order_offset_for_zcurve_offset_err(const std::vector<size_t>& dim_lengths, const int x_coord, const int y_coord, const int z_coord)
{
    cmc_err_msg("The deduction of the correct offset function has failed.");
    return CMC_ERR;
}

/* Array holding all offset functions */
const std::array<std::function<size_t(const std::vector<size_t>&, const int, const int, const int)>, 13> _offset_functions
{
/* 2D Offset funcions */
cmc_t8_calc_lin_order_offset_for_zcurve_lon_lat,
cmc_t8_calc_lin_order_offset_for_zcurve_lat_lon,
cmc_t8_calc_lin_order_offset_for_zcurve_lat_lev,
cmc_t8_calc_lin_order_offset_for_zcurve_lev_lat,
cmc_t8_calc_lin_order_offset_for_zcurve_lon_lev,
cmc_t8_calc_lin_order_offset_for_zcurve_lev_lon,
/* 3D offset functions */
cmc_t8_calc_lin_order_offset_for_zcurve_lon_lat_lev,
cmc_t8_calc_lin_order_offset_for_zcurve_lon_lev_lat,
cmc_t8_calc_lin_order_offset_for_zcurve_lev_lon_lat,
cmc_t8_calc_lin_order_offset_for_zcurve_lev_lat_lon,
cmc_t8_calc_lin_order_offset_for_zcurve_lat_lev_lon,
cmc_t8_calc_lin_order_offset_for_zcurve_lat_lon_lev,
/* Error function */
cmc_t8_calc_lin_order_offset_for_zcurve_offset_err};


/* Return an offset functions based on the data layout of the variable corresponding to the given id */
static
std::function<size_t(const std::vector<size_t>&, const int, const int, const int)>
cmc_t8_get_offset_function_based_on_data_layout(const cmc_t8_data& t8_data, const int var_id)
{
    #ifdef CMC_WITH_T8CODE
    switch (t8_data.vars[var_id]->var->data_layout)
    {
        case DATA_LAYOUT::CMC_2D_LON_LAT:
            return _offset_functions[0];
        break;
        case DATA_LAYOUT::CMC_2D_LAT_LON:
            return _offset_functions[1];
        break;
        case DATA_LAYOUT::CMC_2D_LAT_LEV:
            return _offset_functions[2];
        break;
        case DATA_LAYOUT::CMC_2D_LEV_LAT:
            return _offset_functions[3];
        break;
        case DATA_LAYOUT::CMC_2D_LON_LEV:
            return _offset_functions[4];
        break;
        case DATA_LAYOUT::CMC_2D_LEV_LON:
            return _offset_functions[5];
        break;
        case DATA_LAYOUT::CMC_3D_LON_LAT_LEV:
            return _offset_functions[6];
        break;
        case DATA_LAYOUT::CMC_3D_LON_LEV_LAT:
            return _offset_functions[7];
        break;
        case DATA_LAYOUT::CMC_3D_LEV_LON_LAT:
            return _offset_functions[8];
        break;
        case DATA_LAYOUT::CMC_3D_LEV_LAT_LON:
            return _offset_functions[9];
        break;
        case DATA_LAYOUT::CMC_3D_LAT_LEV_LON:
            return _offset_functions[10];
        break;
        case DATA_LAYOUT::CMC_3D_LAT_LON_LEV:
            return _offset_functions[11];
        break;
        default :
            cmc_err_msg("The variable contains an undefined data layout.");
            return _offset_functions[12];
    }
    #endif
}

//TODO: check if the two functions may be combined
/** Apply the reordering for the given variables compliant to one forest mesh which suites as compression base for all variables */
static void
cmc_t8_apply_zcurve_ordering_one_for_all(cmc_t8_data& t8_data, const std::vector<int>& var_ids, const int dim_compression)
{
    #ifdef CMC_WITH_T8CODE
    cmc_assert(t8_data.vars.size() > 0);
    cmc_assert(dim_compression == 2 || dim_compression == 3);
    cmc_debug_msg("Z-curve reordering of data variables starts...");

    /* Variables for a forest element and it's vertex coordinates */
    t8_element_t *elem;
    int x_coord{0}, y_coord{0}, z_coord{0};
    int vertex_coords[3];

    /* Schemes for quadrilaterals and hexahedrons */
    t8_default_scheme_quad_c scheme_quad;
    t8_default_scheme_hex_c scheme_hex;
    t8_eclass_scheme_c* scheme_eclass;

    /* These pointers will be used in order to copy plain bytes to the reordered new data_array; without explicitly checking the data type, therefore, the use is not encouraged */
    std::byte* data_ptr_src{nullptr};
    std::byte* data_ptr_dest{nullptr};
    size_t offset{0};

    /* Number of elements within the forest */
    const t8_locidx_t num_elems{t8_forest_get_local_num_elements(t8_data.assets->forest)};
    cmc_debug_msg("The forest has ", num_elems, " elements.");

    /* Save the size of the underlying data type of each variable (since it is not required to have only one data type for all variables) */
    std::vector<size_t> size_of_data;
    size_of_data.reserve(var_ids.size());
    
    /* Maximum possible refinement level of the forest, concerning the dimensionality of the forest */ 
    const int max_ref_lvl{t8_data.geo_data->dim == 2 ? P4EST_MAXLEVEL : P8EST_MAXLEVEL};

    /* Get the correct offset functions for each variable (depending on their axis-orderings) */
    std::vector<std::function<size_t(const std::vector<size_t>&, const int, const int, const int)>> offset_functions;
    offset_functions.reserve(var_ids.size());

    /* Save initial data pointers of all variables */
    std::vector<std::byte*> initial_data_ptr;
    initial_data_ptr.reserve(var_ids.size());
    std::vector<std::byte*> initial_data_new_ptr;
    initial_data_new_ptr.reserve(var_ids.size());

    /* Save the data type of the variable */
    std::vector<cmc_type> var_data_type;
    var_data_type.reserve(var_ids.size());

    /* Iterate over all variables and apply the z-curve reordering */
    for (int id{0}; id < static_cast<int>(var_ids.size()); ++id)
    {
        /* Save the data type of each variable */
        var_data_type.push_back(t8_data.vars[var_ids[id]]->get_type());

        /* Allocate a new data array which will hold the reordered element data */
        /** If this functions is called the first time, just after the (e.g. simulation/netCDF) data has been transferred to 't8code-variables', the array size of 'data' correlates to the amount of data points,
          * since the built forest mesh encloses the data points, it holds '#data_new >= #data', therefore 'data_new needs allocate data points equal to the number of mesh elements */
        //t8_data.vars[var_ids[id]]->var = new cmc_var(num_elems, var_data_type[id]);
        t8_data.vars[var_ids[id]]->var->data_new = new var_array_t(static_cast<size_t>(num_elems), var_data_type[id]);
        if (t8_data.vars[var_ids[id]]->var->data_new->size() > t8_data.vars[var_ids[id]]->var->data->size())
        {
            /* In this case dummy elements will be added */
            t8_data.vars[var_ids[id]]->var->missing_value_present = true;
        }

        /* Save the size of the underlying data size of this variable */
        size_of_data.push_back(t8_data.vars[var_ids[id]]->get_data_size());

        /* Save the initial data pointer of the variables */
        initial_data_ptr.push_back(static_cast<std::byte*>(t8_data.vars[var_ids[id]]->get_initial_data_ptr()));
        initial_data_new_ptr.push_back(static_cast<std::byte*>(t8_data.vars[var_ids[id]]->get_initial_data_new_ptr()));

        /* Receive the correct function to calculate the offset */
        if (t8_data.vars[var_ids[id]]->var->data_scheme == CMC_DATA_ORDERING_SCHEME::CMC_GEO_DATA_LINEAR)
        {
            /* Save the offset function corresponding to the variable's data */
            offset_functions.push_back(cmc_t8_get_offset_function_based_on_data_layout(t8_data, var_ids[id]));
        } else
        {
            cmc_err_msg("Other data-schemes are not implemented yet\n");
        }
    }

    cmc_assert(t8_data.vars[0]->var->dim_lengths.size() >= CMC_NUM_COORD_IDS);
    cmc_debug_msg("Dimensions of data to be reordered are: #Lat: ", t8_data.vars[0]->var->dim_lengths[CMC_COORD_IDS::CMC_LAT],
                  ", #LON: ", t8_data.vars[0]->var->dim_lengths[CMC_COORD_IDS::CMC_LON],
                  ", #LEV: ", t8_data.vars[0]->var->dim_lengths[CMC_COORD_IDS::CMC_LEV],
                  ", #TIME: ", t8_data.vars[0]->var->dim_lengths[CMC_COORD_IDS::CMC_TIME]);
    
    /* Iterate over the forest */
    for (t8_locidx_t elem_id{0}; elem_id < num_elems; ++elem_id)
    {
        /* Since the mesh is designed to hold just one tree, the tree_id of every element is 0 */
        elem = t8_forest_get_element_in_tree(t8_data.assets->forest, 0, elem_id);

        /* Get the vertex coordinates of the quad element */
        if (dim_compression == 2)
        {
            scheme_quad.t8_element_vertex_coords(elem, 0, vertex_coords);
            scheme_eclass = static_cast<t8_eclass_scheme_c*>(&scheme_quad);
            /* Transform coordinates into the range of the maximum refinement */
            x_coord = vertex_coords[0] >> (max_ref_lvl - t8_data.assets->initial_refinement_lvl);
            y_coord = vertex_coords[1] >> (max_ref_lvl - t8_data.assets->initial_refinement_lvl);
        } else
        {
            /* Get the vertex coordinates of the hex element */
            scheme_hex.t8_element_vertex_coords(elem, 0, vertex_coords);
            scheme_eclass = static_cast<t8_eclass_scheme_c*>(&scheme_hex);
            /* Transform coordinates into the range of the maximum refinement */
            x_coord = vertex_coords[0] >> (max_ref_lvl - t8_data.assets->initial_refinement_lvl);
            y_coord = vertex_coords[1] >> (max_ref_lvl - t8_data.assets->initial_refinement_lvl);
            z_coord = vertex_coords[2] >> (max_ref_lvl - t8_data.assets->initial_refinement_lvl);
        }
       
        /* Check if the element is in the "lat x lon x lev" mesh (Since all variables are defined in the same mesh, we just retrieve the information from the first variable */
        if (cmc_t8_elem_inside_geo_mesh(elem, scheme_eclass, t8_data, 0) != 0)
        {
            /* Loop over all variables which have to be reordered */
            for (int var_id{0}; var_id < static_cast<int>(var_ids.size()); ++var_id)
            {
                /* Calculate the offset in the linear array */
                offset = offset_functions[var_id](t8_data.vars[var_ids[var_id]]->var->dim_lengths, x_coord, y_coord, z_coord);

                /* Move the pointers to the correct positions, in order to copy the value */
                data_ptr_dest = initial_data_new_ptr[var_id] + elem_id * size_of_data[var_id];
                data_ptr_src = initial_data_ptr[var_id] + offset * size_of_data[var_id];

                /* Copy the value of the unordered data array to the new data array compliant to the z-curve ordering */
                memcpy(static_cast<void*>(data_ptr_dest), static_cast<void*>(data_ptr_src), size_of_data[var_id]);
            }    
        } else
        {
            for (int var_id{0}; var_id < static_cast<int>(var_ids.size()); ++var_id)
            {
                /* If the element is not inside the mesh, a 'missing_value' will be assigned */
                t8_data.vars[var_ids[var_id]]->var->data_new->assign_value(static_cast<size_t>(elem_id), t8_data.vars[var_ids[var_id]]->var->missing_value);
            }
        }
    }

    /* Free the 'unordered' data and assign the newly ordered for each variable */
    for (int var_id{0}; var_id < static_cast<int>(var_ids.size()); ++var_id)
    {
        /* Free the 'unordered' data and assign the newly 'ordered' data for each variable */
        t8_data.vars[var_ids[var_id]]->var->switch_data();
        /* Update the data scheme info */
        t8_data.vars[var_ids[var_id]]->var->data_scheme = CMC_DATA_ORDERING_SCHEME::CMC_GEO_DATA_Z_CURVE;

        cmc_debug_msg("The data of variable ", t8_data.vars[var_ids[var_id]]->var->name, " has been reordered, compliant to the Morton/z-curve order.");
    }
    #endif
}

/** Apply the reordering for the given variables */
static void
cmc_t8_apply_zcurve_ordering_one_for_one(cmc_t8_data& t8_data, const std::vector<int>& var_ids, const int dim_compression)
{
    #ifdef CMC_WITH_T8CODE
    cmc_assert(dim_compression == 2 || dim_compression == 3);
    cmc_debug_msg("Z-curve reordering of data variables starts...");

    /* Variables for a forest element and it's vertex coordinates */
    t8_element_t *elem;
    int x_coord{0}, y_coord{0}, z_coord{0};
    int vertex_coords[3];

    t8_locidx_t num_elems{0};

    /* Schemes for quadrilaterals and hexahedrons */
    t8_default_scheme_quad_c scheme_quad{};
    t8_default_scheme_hex_c scheme_hex{};
    t8_eclass_scheme_c* scheme_eclass{nullptr};

    /* These pointers will be used in order to copy plain bytes to the reordered new data_array; without explicitly checking the data type, therefore, the use is not encouraged */
    std::byte* data_ptr_src{nullptr};
    std::byte* data_ptr_dest{nullptr};
    size_t offset{0};

    /* Used for saving the size of the underlying data type */
    size_t size_of_data{0};

    /* Offset function */
    std::function<size_t(const std::vector<size_t>&, const int, const int, const int)> offset_function;

    /* Maximum possible refinement level of the forest, concerning the dimensionality of the forest */ 
    int max_ref_lvl;

    /* Variables to save initial data pointers */
    std::byte* initial_data_ptr{nullptr};
    std::byte* initial_data_new_ptr{nullptr};

    /* Variable to save the data type */
    cmc_type var_data_type;

    /* Iterate over all variables and apply the z-curve reordering */
    for (int id{0}; id < static_cast<int>(var_ids.size()); ++id)
    {
        /* Get the number of elements of the variable's corresponding forest */
        num_elems = t8_forest_get_local_num_elements(t8_data.vars[var_ids[id]]->assets->forest);
        /* Get the data size */
        size_of_data = t8_data.vars[var_ids[id]]->get_data_size();
        /* Get the correct offset function */
        offset_function = cmc_t8_get_offset_function_based_on_data_layout(t8_data, var_ids[id]);
        /* Save the data type */
        var_data_type = t8_data.vars[var_ids[id]]->get_type();
        /* Allocate a new data array which will hold the reordered element data */
        /** If this functions is called the first time, just after the (e.g. simulation/netCDF) data has been transferred to 't8code-variables', the array size of 'data' correlates to the amount of data points.
          * Since the built forest mesh "encloses" the data points, it holds '#data_new >= #data', therefore 'data_new needs allocate data points equal to the number of mesh elements */
        t8_data.vars[var_ids[id]]->var->data_new = new var_array_t(static_cast<size_t>(num_elems), var_data_type);
        if (t8_data.vars[var_ids[id]]->var->data_new->size() > t8_data.vars[var_ids[id]]->var->data->size())
        {
            /* In this case dummy elements will be added */
            t8_data.vars[var_ids[id]]->var->missing_value_present = true;
        }
        /* Update initial data pointers for each variable */
        initial_data_ptr = static_cast<std::byte*>(t8_data.vars[var_ids[id]]->get_initial_data_ptr());
        initial_data_new_ptr = static_cast<std::byte*>(t8_data.vars[var_ids[id]]->get_initial_data_new_ptr());

        cmc_assert(t8_data.vars[var_ids[id]]->var->dim_lengths.size() >= CMC_NUM_COORD_IDS);
        cmc_debug_msg("Dimensions of data to be reordered are: #Lat: ", t8_data.vars[var_ids[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LAT],
                      ", #LON: ", t8_data.vars[var_ids[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LON],
                      ", #LEV: ", t8_data.vars[var_ids[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LEV],
                      ", #TIME: ", t8_data.vars[var_ids[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_TIME]);
        
        /* Set the maximal refinement level for this variable */
        max_ref_lvl = t8_data.vars[var_ids[id]]->var->num_dimensions == 2 ? P4EST_MAXLEVEL : P8EST_MAXLEVEL;

        /* Iterate over the forest */
        for (t8_locidx_t elem_id{0}; elem_id < num_elems; ++elem_id)
        {
            /* Since the mesh is designed to hold just one tree, the tree_id of every element is 0 */
            elem = t8_forest_get_element_in_tree(t8_data.vars[var_ids[id]]->assets->forest, 0, elem_id);

            /* Get the vertex coordinates of the quad element */
            if (dim_compression == 2)
            {
                scheme_quad.t8_element_vertex_coords(elem, 0, vertex_coords);
                scheme_eclass = static_cast<t8_eclass_scheme_c*>(&scheme_quad);
                /* Transform coordinates into the range of the maximum refinement */
                x_coord = vertex_coords[0] >> (max_ref_lvl - t8_data.vars[var_ids[id]]->assets->initial_refinement_lvl);
                y_coord = vertex_coords[1] >> (max_ref_lvl - t8_data.vars[var_ids[id]]->assets->initial_refinement_lvl);
            } else
            {
                /* Get the vertex coordinates of the hex element */
                scheme_hex.t8_element_vertex_coords(elem, 0, vertex_coords);
                scheme_eclass = static_cast<t8_eclass_scheme_c*>(&scheme_hex);
                /* Transform coordinates into the range of the maximum refinement */
                x_coord = vertex_coords[0] >> (max_ref_lvl - t8_data.vars[var_ids[id]]->assets->initial_refinement_lvl);
                y_coord = vertex_coords[1] >> (max_ref_lvl - t8_data.vars[var_ids[id]]->assets->initial_refinement_lvl);
                z_coord = vertex_coords[2] >> (max_ref_lvl - t8_data.vars[var_ids[id]]->assets->initial_refinement_lvl);
            }

            /* Check if the element is in the "lat x lon x lev" mesh */
            if (cmc_t8_elem_inside_geo_mesh(elem, scheme_eclass, t8_data, var_ids[id]) != 0)
            {
                /* Calculate the offset in the linear array */
                offset = offset_function(t8_data.vars[var_ids[id]]->var->dim_lengths, x_coord, y_coord, z_coord);
                /* Move the pointers to the correct positions, in order to copy the value */
                data_ptr_dest = initial_data_new_ptr + elem_id * size_of_data;
                data_ptr_src = initial_data_ptr + offset * size_of_data;
                /* Copy the value of the unordered data array to the new data array compliant to the z-curve ordering */
                memcpy(static_cast<void*>(data_ptr_dest), static_cast<void*>(data_ptr_src), size_of_data);
            } else
            {
                /* If the element is not inside the mesh, a 'missing_value' will be assigned */
                t8_data.vars[var_ids[id]]->var->data_new->assign(static_cast<size_t>(elem_id), t8_data.vars[var_ids[id]]->var->missing_value);
            }
        }

        /* Free the 'unordered' data and assign the newly 'ordered' data for each variable */
        t8_data.vars[var_ids[id]]->var->switch_data();
        /* Update the data scheme info */
        t8_data.vars[var_ids[id]]->var->data_scheme = CMC_DATA_ORDERING_SCHEME::CMC_GEO_DATA_Z_CURVE;

        cmc_debug_msg("The data of variable ", t8_data.vars[var_ids[id]]->var->name, " has been reordered, compliant to the Morton/z-curve order.");
    }
    #endif
}

/** TODO: Check each case for correctness */
static void
cmc_geo_data_fetch_2d_data(const cmc_t8_var& var3d, cmc_t8_var& var2d, const size_t vanishing_coord_id)
{
    #ifdef CMC_WITH_T8CODE
    /* A linear data scheme is assumed currently */
    cmc_assert(var3d.var->data_scheme == CMC_DATA_ORDERING_SCHEME::CMC_GEO_DATA_LINEAR);
 
    /* We are assuming 3d variables (lat, lon, lev) */
    cmc_assert(var3d.var->num_dimensions == 3);

    /* These pointers will be used in order to copy plain bytes to the reordered new data_array */
    std::byte* data_ptr_src{nullptr};
    std::byte* data_ptr_dest{nullptr};

    /* Get the size in bytes of one array element */
    const size_t data_size{var3d.get_data_size()};

    switch (var3d.var->data_layout)
    {
        case CMC_3D_LAT_LON_LEV:
        {
            switch (var2d.var->data_layout)
            {
                case CMC_2D_LAT_LON:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lat{0}; lat < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                    {
                        for (size_t lon{0}; lon < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (vanishing_coord_id + var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] * (lon + lat * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]));
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lon + lat * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LON_LAT:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lat{0}; lat < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                    {
                        for (size_t lon{0}; lon < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (vanishing_coord_id + var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] * (lon + lat * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]));
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lat + lon * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LAT_LEV:
                {
                    /* This layout constellation enables us to copy a whole elevation chunk at once */
                    for (size_t lat{0}; lat < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                    {
                        data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] * (vanishing_coord_id + lat * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]));
                        data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lat * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]);
                        memcpy(data_ptr_dest, data_ptr_src, data_size * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]); 
                    }
                }
                break;
                case CMC_2D_LEV_LAT:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lat{0}; lat < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                    {
                        for (size_t lev{0}; lev < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] * (vanishing_coord_id + lat * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]));
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lev * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] + lat);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LON_LEV:
                {
                    /* This layout constellation enables us to copy the whole 2D variable at once */
                    data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (vanishing_coord_id * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]);
                    memcpy(var2d.get_initial_data_ptr(), data_ptr_src, data_size * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]);
                }
                break;
                case CMC_2D_LEV_LON:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lon{0}; lon < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                    {
                        for (size_t lev{0}; lev < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (vanishing_coord_id * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] + lon * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] + lev);
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lon + lev * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                default:
                    cmc_err_msg("An unknown preferred (2D) data layout was supplied.");
            }
        }
        break;
        case CMC_3D_LAT_LEV_LON:
        {
            switch (var2d.var->data_layout)
            {
                case CMC_2D_LAT_LON:
                {
                    /* This layout constellation enables us to copy a whole longitude chunk at once */
                    for (size_t lat{0}; lat < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                    {
                        data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * (vanishing_coord_id + lat * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]));
                        data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lat * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]);
                        memcpy(data_ptr_dest, data_ptr_src, data_size * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]);
                    }
                }
                break;
                case CMC_2D_LON_LAT:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lat{0}; lat < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                    {
                        for (size_t lon{0}; lon < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * (vanishing_coord_id + lat * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] + lon));
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lat + lon * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LAT_LEV:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lat{0}; lat < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                    {
                        for (size_t lev{0}; lev < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (vanishing_coord_id + var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * (lev + lat * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]));
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lev + lat * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LEV_LAT:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lat{0}; lat < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                    {
                        for (size_t lev{0}; lev < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (vanishing_coord_id + var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * (lev + lat * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]));
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lat + lev * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LON_LEV:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lev{0}; lev < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                    {
                        for (size_t lon{0}; lon < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * (vanishing_coord_id * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] + lev) + lon);
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lev + lon * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LEV_LON:
                {
                    /* This layout constellation enables us to copy a whole longitude chunk at once */
                    for (size_t lev{0}; lev < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                    {
                        data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * (vanishing_coord_id * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] + lev));
                        data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lev * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]);
                        memcpy(data_ptr_dest, data_ptr_src, data_size * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]);
                    }
                }
                break;
                default:
                    cmc_err_msg("An unknown preferred (2D) data layout was supplied.");
            }
        }
        break;
        case CMC_3D_LEV_LAT_LON:
        {
            switch (var2d.var->data_layout)
            {
                case CMC_2D_LAT_LON:
                {
                    /* This layout constellation enables us to copy the whole 2D variable at once */
                    data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (vanishing_coord_id * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]);
                    memcpy(var2d.get_initial_data_ptr(), data_ptr_src, data_size * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]);
                }
                break;
                case CMC_2D_LON_LAT:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lat{0}; lat < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                    {
                        for (size_t lon{0}; lon < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * (vanishing_coord_id * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] + lat) + lon);
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lat + lon * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LAT_LEV:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lev{0}; lev < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                    {
                        for (size_t lat{0}; lat < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * lev + lat) + vanishing_coord_id);
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lev + lat * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LEV_LAT:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lev{0}; lev < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                    {
                        for (size_t lat{0}; lat < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * lev + lat) + vanishing_coord_id);
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lat + lev * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LON_LEV:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lev{0}; lev < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                    {
                        for (size_t lon{0}; lon < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * lev + vanishing_coord_id) + lon);
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lev + lon * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LEV_LON:
                {
                    /* This layout constellation enables us to copy a whole longitude chunk at once */
                    for (size_t lev{0}; lev < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                    {
                        data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * (vanishing_coord_id + lev * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]));
                        data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lev * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]);
                        memcpy(data_ptr_dest, data_ptr_src, data_size * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]);
                    }
                }
                break;
                default:
                    cmc_err_msg("An unknown preferred (2D) data layout was supplied.");
            }
        }
        break;
        case CMC_3D_LEV_LON_LAT:
        {
            switch (var2d.var->data_layout)
            {
                case CMC_2D_LAT_LON:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lon{0}; lon < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                    {
                        for (size_t lat{0}; lat < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * vanishing_coord_id + lon) + lat);
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lon + lat * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LON_LAT:
                {
                    /* This layout constellation enables us to copy the whole 2D variable at once */
                    data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (vanishing_coord_id * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]);
                    memcpy(var2d.get_initial_data_ptr(), data_ptr_src, data_size * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]);
                }
                break;
                case CMC_2D_LAT_LEV:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lev{0}; lev < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                    {
                        for (size_t lat{0}; lat < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * lev + vanishing_coord_id) + lat);
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lev + lat * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LEV_LAT:
                {
                    /* This layout constellation enables us to copy a latitude chunk at once */
                    for (size_t lev{0}; lev < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                    {
                        data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * lev + vanishing_coord_id));
                        data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lev * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]);
                        memcpy(data_ptr_dest, data_ptr_src, data_size * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]);
                    }
                }
                break;
                case CMC_2D_LON_LEV:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lev{0}; lev < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                    {
                        for (size_t lon{0}; lon < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * lev + lon) + vanishing_coord_id);
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lev + lon * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LEV_LON:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lev{0}; lev < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                    {
                        for (size_t lon{0}; lon < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON] * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * lev + lon) + vanishing_coord_id);
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lon + lev * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                default:
                    cmc_err_msg("An unknown preferred (2D) data layout was supplied.");
            }
        }
        break;
        case CMC_3D_LON_LEV_LAT:
        {
            switch (var2d.var->data_layout)
            {
                case CMC_2D_LAT_LON:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lon{0}; lon < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                    {
                        for (size_t lat{0}; lat < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] * lon + vanishing_coord_id) + lat);
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lon + lat * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LON_LAT:
                {
                    /* This layout constellation enables us to copy a latitude chunk at once */
                    for (size_t lon{0}; lon < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                    {
                        data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] * lon + vanishing_coord_id));
                        data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lon * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]);
                        memcpy(data_ptr_dest, data_ptr_src, data_size * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]);
                    }
                }
                break;
                case CMC_2D_LAT_LEV:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lev{0}; lev < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                    {
                        for (size_t lat{0}; lat < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] * vanishing_coord_id + lev) + lat);
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lev + lat * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LEV_LAT:
                {
                    /* This layout constellation enables us to copy the whole 2D variable at once */
                    data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (vanishing_coord_id * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]);
                    memcpy(var2d.get_initial_data_ptr(), data_ptr_src, data_size * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]);
                }
                break;
                case CMC_2D_LON_LEV:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lon{0}; lon < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                    {
                        for (size_t lev{0}; lev < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] * lon + lev) +vanishing_coord_id);
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lev + lon * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LEV_LON:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lon{0}; lon < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                    {
                        for (size_t lev{0}; lev < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] * lon + lev) +vanishing_coord_id);
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lon + lev * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                default:
                    cmc_err_msg("An unknown preferred (2D) data layout was supplied.");
            }
        }
        break;
        case CMC_3D_LON_LAT_LEV:
        {
            switch (var2d.var->data_layout)
            {
                case CMC_2D_LAT_LON:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lon{0}; lon < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                    {
                        for (size_t lat{0}; lat < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * lon + lat) + vanishing_coord_id);
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lon + lat * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LON_LAT:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lon{0}; lon < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                    {
                        for (size_t lat{0}; lat < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * lon + lat) + vanishing_coord_id);
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lat + lon * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LAT_LEV:
                {
                    /* This layout constellation enables us to copy the whole 2D variable at once */
                    data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (vanishing_coord_id * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]);
                    memcpy(var2d.get_initial_data_ptr(), data_ptr_src, data_size * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]);
                }
                break;
                case CMC_2D_LEV_LAT:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lat{0}; lat < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                    {
                        for (size_t lev{0}; lev < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * vanishing_coord_id + lat) + lev);
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lat + lev * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                case CMC_2D_LON_LEV:
                {
                    /* This layout constellation enables us to copy an elevation chunk at once */
                    for (size_t lon{0}; lon < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                    {
                        data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * lon + vanishing_coord_id));
                        data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lon * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]);
                        memcpy(data_ptr_dest, data_ptr_src, data_size * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]);
                    }
                }
                break;
                case CMC_2D_LEV_LON:
                {
                    /* This layout constellation forces us to copy each value independently */
                    for (size_t lon{0}; lon < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                    {
                        for (size_t lev{0}; lev < var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                        {
                            data_ptr_src = static_cast<std::byte*>(var3d.get_initial_data_ptr()) + data_size * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LEV] * (var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * lon + vanishing_coord_id) + lev);
                            data_ptr_dest = static_cast<std::byte*>(var2d.get_initial_data_ptr()) + data_size * (lon + lev * var3d.var->dim_lengths[CMC_COORD_IDS::CMC_LON]);
                            memcpy(data_ptr_dest, data_ptr_src, data_size);
                        }
                    }
                }
                break;
                default:
                    cmc_err_msg("An unknown preferred (2D) data layout was supplied.");
            }
        }
        break;
        default:
            cmc_err_msg("There was an unknown 3D data layout supplied.");
    }
    #endif
}

/**
 * @brief This functions sets up the adapatation and interpolation structs based on a given compression criterion (Allcoations etc. may be done here)
 * 
 * @param adapt_data The @struct cmc_t8_adapt_data to set up based on the compression criterion
 * @param interpolation_data The @struct cmc_t8_interpolation_data to set up based on the compression criterion
 * @param t8_data A pointer to the @struct cmc_t8_data holding all variables and data
 * @param var_id Eventually a var_id (only used within a 'One For One' compression mode)
 */
static void
cmc_t8_set_up_adapt_data_and_interpolate_data_based_on_compression_settings(cmc_t8_adapt_data& adapt_data, cmc_t8_interpolation_data& interpolation_data, const cmc_t8_data_t t8_data, const int var_id = -1)
{
    #ifdef CMC_WITH_T8CODE
    /* Set a pointer to t8_data if it has not happend before */
    if (adapt_data.t8_data == nullptr)
    {
        adapt_data.t8_data = t8_data;
    }
    /* Set a pointer to t8_data if it has not happend before */
    if (interpolation_data.t8_data == nullptr)
    {
        interpolation_data.t8_data = t8_data;
    }

    /* Set up the data structures based on the compression mode */
    if (t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
        t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D)
    {
        /* If a 'One for All' compression mode is used */
        switch (t8_data->settings.compression_criterium)
        {
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_CRITERIUM_UNDEFINED:
            /* If the criterium is undefined, we use by default an error threshold criterium */
            cmc_warn_msg("No compression criterion has been specified.");
            [[fallthrough]];
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_REL_ERROR_THRESHOLD:
                adapt_data.initial_ref_lvl_ids.emplace_back(std::unordered_map<t8_locidx_t, t8_locidx_t>());
                adapt_data.initial_ref_lvl_ids.back().reserve(static_cast<size_t>(t8_forest_get_local_num_elements(t8_data->assets->forest) / pow(t8_data->geo_data->dim, 2)));
                adapt_data._counter = 0;
                adapt_data._counter_nxt_lvl = 0;
                /* Create a var_vector holding data which may be used by the interpolation */
                adapt_data.adapted_data = new var_vector_t();
                adapt_data.adapted_data->reserve(t8_data->vars.size());

                /* Initialize the deviations vectors */
                adapt_data.associated_max_deviations.insert(adapt_data.associated_max_deviations.begin(), t8_data->vars.size(), std::vector<double>());
                adapt_data.associated_max_deviations_new.insert(adapt_data.associated_max_deviations_new.begin(), t8_data->vars.size(), std::vector<double>());
            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_EXCLUDE_AREA:
                // In this case nothing has to be set up
            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_COMBINED_CRITERION:
                adapt_data.initial_ref_lvl_ids.emplace_back(std::unordered_map<t8_locidx_t, t8_locidx_t>());
                adapt_data.initial_ref_lvl_ids.back().reserve(static_cast<size_t>(t8_forest_get_local_num_elements(t8_data->assets->forest) / pow(t8_data->geo_data->dim, 2)));
                adapt_data._counter = 0;
                adapt_data._counter_nxt_lvl = 0;
                /* Create a var_vector holding data which may be used by the interpolation */
                adapt_data.adapted_data = new var_vector_t();
                adapt_data.adapted_data->reserve(t8_data->vars.size());
            break;
            default:
                cmc_err_msg("The supplied lossy compression criterium is not yet implemented.");
        }
    }
    else
    {
        /* If a 'One for One' compression mode is used */
        switch (t8_data->settings.compression_criterium)
        {
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_CRITERIUM_UNDEFINED:
            /* If the mode is undefined, we use by default an error threshold criterium */
            [[fallthrough]];
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_REL_ERROR_THRESHOLD:
                /* Reset the adapt_data class */
                adapt_data.adapt_step = 0;
                /* Save the current variable ID */
                adapt_data.current_var_id = var_id;
                interpolation_data.current_var_id = var_id;
                /* Save the initial data */
                adapt_data.initial_ref_lvl_ids.emplace_back(std::unordered_map<t8_locidx_t, t8_locidx_t>());
                adapt_data.initial_ref_lvl_ids.back().reserve(static_cast<size_t>(t8_forest_get_local_num_elements(t8_data->vars[var_id]->assets->forest) / pow(t8_data->vars[var_id]->var->num_dimensions, 2)));

                /* Clear all previous vectors in order to start new for the next variable */
                adapt_data.associated_max_deviations.clear();
                adapt_data.associated_max_deviations_new.clear();
                adapt_data.associated_deviations_gelement_id.clear();
                adapt_data.associated_deviations_gelement_id_new.clear();
                /* Initialize the deviations vectors */
                adapt_data.associated_max_deviations.push_back(std::vector<double>());
                adapt_data.associated_max_deviations_new.push_back(std::vector<double>());

            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_EXCLUDE_AREA:
                /* Reset the adapt_data class */
                adapt_data.adapt_step = 0;
                /* Save the current variable ID */
                adapt_data.current_var_id = var_id;
                interpolation_data.current_var_id = var_id;
            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_COMBINED_CRITERION:
                /* Reset the adapt_data class */
                adapt_data.adapt_step = 0;
                /* Save the current variable ID */
                adapt_data.current_var_id = var_id;
                interpolation_data.current_var_id = var_id;
                /* Save the initial data */
                adapt_data.initial_ref_lvl_ids.emplace_back(std::unordered_map<t8_locidx_t, t8_locidx_t>());
                adapt_data.initial_ref_lvl_ids.back().reserve(static_cast<size_t>(t8_forest_get_local_num_elements(t8_data->vars[var_id]->assets->forest) / pow(t8_data->vars[var_id]->var->num_dimensions, 2)));
            
                adapt_data.associated_max_deviations.push_back(std::vector<double>());
                adapt_data.associated_max_deviations_new.push_back(std::vector<double>());
            break;
            default:
                cmc_err_msg("The supplied lossy compression criterium is not yet implemented.");
        }
    }
    #endif
}

static void
cmc_t8_update_adapt_and_interpolation_data_beginning_of_iteration(cmc_t8_adapt_data& adapt_data, cmc_t8_interpolation_data& interpolation_data, const int var_id = -1)
{
    #ifdef CMC_WITH_T8CODE
    /* Set up the data structures based on the compression mode */
    if (adapt_data.t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
        adapt_data.t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D)
    {
        /* If a 'One for All' compression mode is used */
        switch (adapt_data.t8_data->settings.compression_criterium)
        {
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_CRITERIUM_UNDEFINED:
            /* If the criterium is undefined, we use by default an error threshold criterium */
            [[fallthrough]];
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_REL_ERROR_THRESHOLD:
                /* Reset the counters */
                adapt_data._counter = 0;
                adapt_data._counter_nxt_lvl = 0;
                adapt_data.coarsening_counter = 0;
                interpolation_data.coarsening_counter = 0;
                /* Allocate space for each variable in order to save data which may be used by the interpolation */
                for (size_t id{0}; id < adapt_data.t8_data->vars.size(); ++id)
                {
                    adapt_data.adapted_data->push_back(new var_array_t(static_cast<size_t>(t8_forest_get_local_num_elements(adapt_data.t8_data->assets->forest) / (2 * adapt_data.t8_data->vars[id]->var->num_dimensions) +1), adapt_data.t8_data->vars[id]->get_type()));
                }
                /* Save a pointer to the 'adapted_data' in the interpolation struct */
                interpolation_data.adapted_data = adapt_data.adapted_data;
            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_EXCLUDE_AREA:
                //Here has nothing to be done
            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_COMBINED_CRITERION:
                /* Reset the counters */
                adapt_data._counter = 0;
                adapt_data._counter_nxt_lvl = 0;
                adapt_data.coarsening_counter = 0;
                interpolation_data.coarsening_counter = 0;
                /* Allocate space for each variable in order to save data which may be used by the interpolation */
                for (size_t id{0}; id < adapt_data.t8_data->vars.size(); ++id)
                {
                    adapt_data.adapted_data->push_back(new var_array_t(static_cast<size_t>(t8_forest_get_local_num_elements(adapt_data.t8_data->assets->forest) / (2 * adapt_data.t8_data->vars[id]->var->num_dimensions) +1), adapt_data.t8_data->vars[id]->get_type()));
                }
                /* Save a pointer to the 'adapted_data' in the interpolation struct */
                interpolation_data.adapted_data = adapt_data.adapted_data;    
            break;
            default:
                cmc_err_msg("The supplied lossy compression criterium is not yet implemented.");
        }
    }
    else
    {
        /* If a 'One for One' compression mode is used */
        switch (adapt_data.t8_data->settings.compression_criterium)
        {
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_CRITERIUM_UNDEFINED:
            /* If the mode is undefined, we use by default an error threshold criterium */
            [[fallthrough]];
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_REL_ERROR_THRESHOLD:
                /* Reset the adapt counters */
                adapt_data._counter = 0;
                adapt_data._counter_nxt_lvl = 0;
                adapt_data.coarsening_counter = 0;
                adapt_data.adapted_data->push_back(new var_array_t(static_cast<size_t>(t8_forest_get_local_num_elements(adapt_data.t8_data->vars[var_id]->assets->forest) / (2 * adapt_data.t8_data->vars[var_id]->var->num_dimensions) +1), adapt_data.t8_data->vars[var_id]->get_type()));
                /* Save a pointer to the 'adapted_data' for the interpolation */
                interpolation_data.adapted_data = adapt_data.adapted_data;
                interpolation_data.coarsening_counter = 0;
            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_EXCLUDE_AREA:
                //Here has nothing to be done
            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_COMBINED_CRITERION:
                /* Reset the adapt counters */
                adapt_data._counter = 0;
                adapt_data._counter_nxt_lvl = 0;
                adapt_data.coarsening_counter = 0;
                adapt_data.adapted_data->push_back(new var_array_t(static_cast<size_t>(t8_forest_get_local_num_elements(adapt_data.t8_data->vars[var_id]->assets->forest) / (2 * adapt_data.t8_data->vars[var_id]->var->num_dimensions) +1), adapt_data.t8_data->vars[var_id]->get_type()));
                /* Save a pointer to the 'adapted_data' for the interpolation */
                interpolation_data.adapted_data = adapt_data.adapted_data;
                interpolation_data.coarsening_counter = 0;
            break;
            default:
                cmc_err_msg("The supplied lossy compression criterium is not yet implemented.");
        }
    }
    #endif
}
static void
cmc_t8_adapt_struct_allocate_one_for_one(cmc_t8_adapt_data& adapt_data)
{
    #ifdef CMC_WITH_T8CODE
    cmc_assert(adapt_data.t8_data != nullptr);
    /* Allocate space for all unordered maps (for each varibale) */
    adapt_data.initial_ref_lvl_ids.reserve(adapt_data.t8_data->vars.size());
    /* Create a var_vector holding data which may be used by the interpolation */
    adapt_data.adapted_data = new var_vector_t();
    adapt_data.adapted_data->reserve(1);
    #endif
}

struct
cmc_t8_mpi_recv_max_deviations
{
    uint64_t* elem_ids; //!> A vector collecting the received element ids to which the deviaitons corresponds 
    double** deviations; //!< A vector of double vecctor receiving the deviaitons for each variable
    int num_receiving_elems{0}; //!< Number of elements which has been received
};

/**
 * @brief This functions updates the deviations calcultaed for the realtive error criterion. The deviaitons will be associated with the adapted forest
 *        and updated acordintgly. 
 * 
 * @param adapt_data The @struct cmc_t8_adapt_data which holds all information about the data, the deviations and the reference forest
 */
static void
cmc_rel_error_threshold_update_deviations(cmc_t8_adapt_data& adapt_data, int num_preceeding_coarsenings = 0)
{
    #ifdef CMC_WITH_T8CODE

    cmc_assert(adapt_data.forest_reference != nullptr);

    /* Get the scheme of the forest's only tree */
    t8_eclass_scheme_c* ts =  t8_forest_get_eclass_scheme (adapt_data.forest_reference, t8_forest_get_eclass(adapt_data.forest_reference, 0));

    /* Get the childrent the elements refine to */
    const uint64_t num_children = static_cast<uint64_t>(ts->t8_element_num_children(t8_forest_get_element_in_tree(adapt_data.forest_reference, 0, 0)));

    /* Multiply this number by the amount of element reduction per coarsening */
    num_preceeding_coarsenings *= (num_children - 1);

    /* Update the global ids (corresponding to the forest which was not adapted) by subtracting the num_preceeding_coarsenings. This is possible because there is no refinement ever introduced during the compression step */
    for (auto iter = adapt_data.associated_deviations_gelement_id_new.begin(); iter != adapt_data.associated_deviations_gelement_id_new.end(); ++iter)
    {
        /* Subtract the coarsenings */
        *iter = *iter - num_preceeding_coarsenings;
        /* Update the num_preceeding_coarsenings by this coarsening */
        num_preceeding_coarsenings += (num_children - 1);
    }

    #endif
}

/**
 * @brief This function is the serial equivalent to the function 'cmc_rel_error_threshold_update_and_partition_deviations(...)'. The vectors are prepared and switched for the next iteration
 *        and the element ids are updated such that they comply with the adapted forest.
 * 
 * @param adapt_data The @struct cmc_t8_adapt_data which holds all information about the data, the deviations and the reference forest
 */
static void
cmc_rel_error_threshold_update_deviations_serial(cmc_t8_adapt_data& adapt_data)
{
    #ifdef CMC_WITH_T8CODE

    /* Update the deviations accordingly to the reference forest */
    cmc_rel_error_threshold_update_deviations(adapt_data);

    /* Swap the vectors holding the global element ids */
    std::swap(adapt_data.associated_deviations_gelement_id, adapt_data.associated_deviations_gelement_id_new);

    /* Clear the associated_deviations_gelement_id_new in order to fill it within the next adaptation step */
    adapt_data.associated_deviations_gelement_id_new.clear();

    const size_t num_deviation_vars = adapt_data.associated_max_deviations_new.size();

    for(size_t var_iter = 0; var_iter < num_deviation_vars; ++var_iter)
    {
        /* Swap the deviations */
        std::swap(adapt_data.associated_max_deviations[var_iter], adapt_data.associated_max_deviations_new[var_iter]);
        /* Clear the previous deviations */
        adapt_data.associated_max_deviations_new[var_iter].clear();
    }
    #endif
}

/**
 * @brief This functions updates and partitions the deviations calculated for the relative error criterion. The deviaitons will be associated with the (partitined) forest
 *        and communicated acordintgly.
 * 
 * @param adapt_data The @struct cmc_t8_adapt_data which holds all information about the data, the deviations and the reference forest (indicating the partition)
 */
static void
cmc_rel_error_threshold_update_and_partition_deviations(cmc_t8_adapt_data& adapt_data)
{
    #ifdef CMC_WITH_T8CODE
    #ifdef CMC_ENABLE_MPI

    cmc_debug_msg("The maximum relative deviations from the previous adaptation step will be updated and partitioned."); 

    /* Adjust the previous maximum deviations element ids, such that they comply with the new adapted forest */
    int err, mpisize, mpirank;

    /* Get the size of the communicator */
    err = MPI_Comm_size(adapt_data.t8_data->comm, &mpisize);
    cmc_mpi_check_err(err);

    /* Get the (local) rank id */
    err = MPI_Comm_rank(adapt_data.t8_data->comm, &mpirank);
    cmc_mpi_check_err(err);

    /* Allocate memory for the Allgather function call */
    uint64_t* amount_of_coarsenings = new uint64_t[mpisize];

    /* Allocate an offset array holding the SFC index (at the intitial refinement level) of each starting process-local element from each process */
    uint64_t* offsets = new uint64_t[mpisize];

    /* Get the offest of the first local element as a SFC index (at the initial refienment level) */
    uint64_t elem_offset = static_cast<uint64_t>(t8_forest_get_first_local_element_id(adapt_data.forest_reference));

    /* Distribute the offsets between all processes */
    err = MPI_Allgather(&elem_offset, 1, MPI_UINT64_T, offsets, 1, MPI_UINT64_T, adapt_data.t8_data->comm);
    cmc_mpi_check_err(err);

    /* Inquire the amount of coarsenings per process */
    err = MPI_Allgather(&(adapt_data.coarsening_counter), 1, MPI_UINT64_T, amount_of_coarsenings, 1, MPI_UINT64_T, adapt_data.t8_data->comm);
    cmc_mpi_check_err(err);

    uint64_t num_preceeding_coarsenings = 0;

    /* Calculate the coarsenings which have happend up to my local process */
    for (int i = 0; i < mpirank; ++i)
    {
        /* Add up all preceeding coarsenings */
        num_preceeding_coarsenings += amount_of_coarsenings[i];
    }

    /* Update the deviations according the amount of coarsenings that had happend */
    cmc_rel_error_threshold_update_deviations(adapt_data, num_preceeding_coarsenings);

    /* Save the lcoal offset andf length of the data which will stay process-local within this pair (first: local_offset within array; second: length of prospective process-local range) */
    std::pair<int, int> already_process_local_range;

    /* Create a vector of MPI requests */
    std::vector<MPI_Request> send_requests;
    send_requests.reserve((mpisize - 1) * adapt_data.associated_max_deviations_new.size());

    /* A counter for the send requests */
    size_t send_req_idx = 0;

    /* Define an arbitrary tag for the element ids */
    const int mpi_tag_gelem_ids = adapt_data.associated_max_deviations_new.size();

    /* A flag indicating whether some objects are already correctly distributed or if all elements have been sent away to different ranks */
    bool flag_some_data_elements_are_kept_local = false;

    /* Itearte over all element ids (of the deviations) */
    for (auto elem_id_iter = adapt_data.associated_deviations_gelement_id_new.begin(); elem_id_iter != adapt_data.associated_deviations_gelement_id_new.end();)
    {
        int recv_rank = mpisize -1;
        for (int ir{0}; ir < mpisize -1; ++ir)
        {
            if (*elem_id_iter >= offsets[ir] && *elem_id_iter < offsets[ir+1])
            {
                recv_rank = ir;
                break;
            }
        }

        /* Check the length of the contiguous elements which will be sent to this process */
        const int msg_length = std::min(static_cast<int>(std::distance(elem_id_iter, adapt_data.associated_deviations_gelement_id_new.end())), static_cast<int>((recv_rank != mpisize - 1 ?  offsets[recv_rank +1] - *elem_id_iter : t8_forest_get_global_num_elements(adapt_data.forest_reference) - *elem_id_iter)));

        /* Since at least one element has been found before, this assertion should be fullfilled */
        cmc_assert(msg_length > 0);

        if (recv_rank != mpirank)
        {
            /* If the receiver is any other rank than itself, we can send the corresponding data */
            /* We will send the data concerning the element ids first */
            err = MPI_Isend(adapt_data.associated_deviations_gelement_id_new.data() + std::distance(adapt_data.associated_deviations_gelement_id_new.begin(), elem_id_iter), msg_length, MPI_UINT64_T, recv_rank, mpi_tag_gelem_ids, adapt_data.t8_data->comm, &send_requests[send_req_idx]);
            cmc_mpi_check_err(err);

            /* Increment the send_req_idx counter */
            ++send_req_idx;

            /* Afterwards we will send the deviations per variable */
            for (size_t var_iter = 0; var_iter < adapt_data.associated_max_deviations_new.size(); ++var_iter)
            {
                err = MPI_Isend(adapt_data.associated_max_deviations_new.data() + std::distance(adapt_data.associated_deviations_gelement_id_new.begin(), elem_id_iter), msg_length, MPI_DOUBLE, recv_rank, var_iter, adapt_data.t8_data->comm, &send_requests[send_req_idx]);
                cmc_mpi_check_err(err);

                /* Increment the send_req_idx counter */
                ++send_req_idx;
            }
        } else
        {
            /* Save which (already process-local) range belongs to this process */
            already_process_local_range = std::make_pair(std::distance(adapt_data.associated_deviations_gelement_id_new.begin(), elem_id_iter), msg_length);

            /* Set the flag that some elements reside locallly */
            flag_some_data_elements_are_kept_local = true;
        }


        /* Now we can skip all elements which we have sent */
        std::advance(elem_id_iter, msg_length);
    }

    /* Wait until all messages are staged */
    err = MPI_Barrier(adapt_data.t8_data->comm);
    cmc_mpi_check_err(err);


    /* Counter for the amount of data which will be process-local */
    uint64_t counter_max_deviations = already_process_local_range.second;

    /* Create a map saving the received data */
    std::map<int, cmc_t8_mpi_recv_max_deviations> recv_list;

    /* Variable which will hold the amount of elements within the message */
    int received_objs = 0;

    /* MPI Status objects for the send commands */
    MPI_Status recv_stat, actual_recv_stat;

    /* Flags for working through the received messages */
    bool flag_messages_present = true;
    int mpi_msg_flag = 0;

    /* Receive the messages concerning the max deviations */
    /* Since all messages have been staged before, the MPI_Iprobe call finds a matching message in each iteration */
    while (flag_messages_present)
    {
        /* Check for a message containing the global element indices of the deviations or the deviations themself from another rank */
        err = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, adapt_data.t8_data->comm, &mpi_msg_flag, &recv_stat);
        cmc_mpi_check_err(err);

        /* Check if a message from any source with any tag is available */
        if (mpi_msg_flag)
        {
            /* Get the size of the message, respectively the amount of max deviations which will be received */
            /* This number equals the amount of data which will be receied from the partner rank for each variable */
            err = MPI_Get_count(&recv_stat, (recv_stat.MPI_TAG != mpi_tag_gelem_ids ? MPI_DOUBLE : MPI_UINT64_T), &received_objs);
            cmc_mpi_check_err(err);

            /* Check if the rank which has sent the data is already in the 'recv_list' */
            /* Check for the sending rank in the 'recv_list' */
            auto sender = recv_list.find(recv_stat.MPI_SOURCE);

            /* Check if the rank is not yet a sending rank in the 'recv_list' */
            if (sender == recv_list.end())
            {
                /* Create a new cmc_mpi_t8_recv_data in the recv_list */
                recv_list[recv_stat.MPI_SOURCE] = cmc_t8_mpi_recv_max_deviations();
                /* Reserve memory for the incoming bytes resembling the global eleemnt id and the maximum deviation */
                recv_list[recv_stat.MPI_SOURCE].elem_ids = new uint64_t[received_objs];
                /* Reserve memory for the deviations to come */
                recv_list[recv_stat.MPI_SOURCE].deviations = new double*[adapt_data.associated_max_deviations_new.size()];
                for (size_t dev_iter = 0; dev_iter < adapt_data.associated_max_deviations_new.size(); ++dev_iter)
                {
                    recv_list[recv_stat.MPI_SOURCE].deviations[dev_iter] = new double[received_objs];
                }
                /* Save the amount of received objects from this rank */
                recv_list[recv_stat.MPI_SOURCE].num_receiving_elems = received_objs;
            }

            /* Distinguish between global element ids and deviation data */
            if (recv_stat.MPI_TAG == mpi_tag_gelem_ids)
            {
                /* In this case we have received element ids corresponding to the deviation data */
                err = MPI_Recv(recv_list[recv_stat.MPI_SOURCE].elem_ids, received_objs, MPI_UINT64_T, recv_stat.MPI_SOURCE, mpi_tag_gelem_ids, adapt_data.t8_data->comm, &actual_recv_stat);
                cmc_mpi_check_err(err);

                /* Accumulate the incoming elements */
                counter_max_deviations += received_objs;
            } else
            {
                /* In this case we have received deviation data */
                err = MPI_Recv(recv_list[recv_stat.MPI_SOURCE].deviations[recv_stat.MPI_TAG], received_objs, MPI_DOUBLE, recv_stat.MPI_SOURCE, recv_stat.MPI_TAG, adapt_data.t8_data->comm, &actual_recv_stat);
                cmc_mpi_check_err(err);
            }

            cmc_debug_msg("Received ", received_objs, " data elements from rank ", actual_recv_stat.MPI_SOURCE, " with tag ", actual_recv_stat.MPI_TAG);
        } else
        {
            /* If none messages are left */
            flag_messages_present = false;
        }
    }


    /* If all messages are received, we can order the data */
    /* Therefore, we resize the previous associated_deviations_gelement_id array to the correct size and copy the already ordered element ids from the other ranks */
    adapt_data.associated_deviations_gelement_id.resize(counter_max_deviations);

    /* Resize the arrays holding the deviations of all variables */
    for(size_t var_iter = 0; var_iter < adapt_data.associated_max_deviations_new.size(); ++var_iter)
    {
        adapt_data.associated_max_deviations[var_iter].resize(counter_max_deviations);
    }

    /* An offset for copying the data */
    int offset = 0;
    int offset_for_previous_local_data = 0;

    /* Since the map is already ordered by the ranks which have sent the data, we can iteratively draw the data from the map */
    for (auto iter = recv_list.begin(); iter != recv_list.end(); ++iter)  
    {
        /* Copy all element ids over */
        std::copy(iter->second.elem_ids, iter->second.elem_ids + iter->second.num_receiving_elems, adapt_data.associated_deviations_gelement_id.begin() + offset + (iter->first < mpirank ? 0 : already_process_local_range.second));
        /* Copy the data for each variable */
        for(size_t var_iter = 0; var_iter < adapt_data.associated_max_deviations_new.size(); ++var_iter)
        {
            std::copy(iter->second.deviations[var_iter], iter->second.deviations[var_iter] + iter->second.num_receiving_elems, adapt_data.associated_max_deviations[var_iter].begin() + offset + (iter->first < mpirank ? 0 : already_process_local_range.second));
        }
        if (iter->first < mpirank)
        {
            /* Save the offset until the data from this rank will come */
            offset_for_previous_local_data += iter->second.num_receiving_elems;
        }
        /* Save the general offest currently present in the array */
        offset += iter->second.num_receiving_elems;
    }

    /* If some data was already local, we need to copy it as well */
    if (flag_some_data_elements_are_kept_local)
    {
        /* Copy all element ids over */
        std::copy(adapt_data.associated_deviations_gelement_id_new.data() + already_process_local_range.first, adapt_data.associated_deviations_gelement_id_new.data() + already_process_local_range.first + already_process_local_range.second, adapt_data.associated_deviations_gelement_id.begin() + offset_for_previous_local_data);
        /* Copy the data for each variable */
        for(size_t var_iter = 0; var_iter < adapt_data.associated_max_deviations_new.size(); ++var_iter)
        {
            std::copy(adapt_data.associated_max_deviations_new[var_iter].data() + already_process_local_range.first, adapt_data.associated_max_deviations_new[var_iter].data() + already_process_local_range.first + already_process_local_range.second, adapt_data.associated_max_deviations[var_iter].begin() + offset_for_previous_local_data);
        }
    }

    /* If all messages have been sent, we can clear the array associated_max_deviations_new for the next iteration */
    err = MPI_Waitall(send_req_idx, send_requests.data(), MPI_STATUS_IGNORE);
    cmc_mpi_check_err(err);

    /* Clear the arrays and the element indices */
    for (size_t var_iter = 0; var_iter < adapt_data.associated_max_deviations_new.size(); ++var_iter)
    {
        adapt_data.associated_max_deviations_new[var_iter].clear();
    }
    adapt_data.associated_deviations_gelement_id_new.clear();

    /* For convenience, we store the global element offset of the partitined forest in here */
    adapt_data.first_global_elem_id = elem_offset;

    /* Deallocate the array holding the coarsenings per process and the array holding the offsets */
    delete[] amount_of_coarsenings;
    delete[] offsets;

    /* Delete the allocated array within the recv_list */
    for (auto iter = recv_list.begin(); iter != recv_list.end(); ++iter)  
    {
        /* Delete the array containing the received element ids */
        delete[] iter->second.elem_ids;
        /* Delete all arrays holding the deviations from the varibales */
        for (size_t var_iter = 0; var_iter < adapt_data.associated_max_deviations.size(); ++var_iter)
        {
            delete[] iter->second.deviations[var_iter];
        }
        /* Delete the array deviations array itself */
        delete[] iter->second.deviations;
    }

    cmc_debug_msg("The maximum relative deviations from the previous adaptation step were updated and partitioned.");
    #endif
    #endif
}

/**
 * @brief This function updates the members of the adapt_data which are responsible for the relative error criterion. Based on the whether
 *        or not the data is distributed, the deviations will be partitioned after updating.
 * 
 * @param adapt_data The @struct cmc_t8_adapt_data whose member concerning the relative error criterion will be updated
 */
static void
cmc_rel_error_threshold_communicate_and_update_deviations(cmc_t8_adapt_data& adapt_data)
{
    #ifdef CMC_WITH_T8CODE
    if (adapt_data.t8_data->use_distributed_data)
    {
        /* In case of a parallel environment */
        cmc_rel_error_threshold_update_and_partition_deviations(adapt_data);
    } else
    {
        /* In case of a serial environemnt */
        cmc_rel_error_threshold_update_deviations_serial(adapt_data);
    }
    #endif
}

/**
 * @brief This functions updates the adapt and interpolation structs at the beginning of new adaptation step depending on the used compression mode and compression criterion
 * 
 * @param adapt_data The @struct cmc_t8_adapt_data to update/reset
 * @param interpolation_data The @struct cmc_t8_interpolation_data to update/reset
 */
static void
cmc_t8_update_adapt_and_interpolation_data_end_of_iteration(cmc_t8_adapt_data& adapt_data, cmc_t8_interpolation_data& interpolation_data)
{
    #ifdef CMC_WITH_T8CODE
    /* Set up the data structures based on the compression mode */
    if (adapt_data.t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
        adapt_data.t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D)
    {
        /* If a 'One for All' compression mode is used */
        switch (adapt_data.t8_data->settings.compression_criterium)
        {
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_CRITERIUM_UNDEFINED:
            /* If the criterium is undefined, we use by default an error threshold criterium */
            [[fallthrough]];
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_REL_ERROR_THRESHOLD:
                /* Free the data which has been saved during the adaptation */
                adapt_data.adapted_data->clear();
                interpolation_data.adapted_data = nullptr;

                cmc_rel_error_threshold_communicate_and_update_deviations(adapt_data);

                #if 0
                if (adapt_data.t8_data->use_distributed_data)
                {
                    /* In parallel environment, we need to communicate the previous maximum deviations */
                    /* Within this function call, we communicate the maximum deviations, order them and return them within the adapt_data for the next iteration */
                    /* The global elementid of the first local element (of the adapt_data->forest_reference) is inquired there as well and stored withint the adapt_data for the next iteratzion */
                    cmc_rel_error_threshold_communicate_and_update_deviations(adapt_data);
                } else
                {
                    /* In a serial environment, it suffices to just update the deviations */
                    cmc_rel_error_threshold_update_deviations(adapt_data);
                }
                #endif
            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_EXCLUDE_AREA:
                //Here has nothing to be done
            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_COMBINED_CRITERION:
                /* Free the data which has been saved during the adaptation */
                adapt_data.adapted_data->clear();
                interpolation_data.adapted_data = nullptr;
            break;
            default:
                cmc_err_msg("The supplied lossy compression criterium is not yet implemented.");
        }
    }
    else
    {
        /* If a 'One for One' compression mode is used */
        switch (adapt_data.t8_data->settings.compression_criterium)
        {
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_CRITERIUM_UNDEFINED:
            /* If the mode is undefined, we use by default an error threshold criterium */
            [[fallthrough]];
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_REL_ERROR_THRESHOLD:
                /* Free the data which has been saved during the adaptation */
                adapt_data.adapted_data->clear();
                interpolation_data.adapted_data = nullptr;

                cmc_rel_error_threshold_communicate_and_update_deviations(adapt_data);

                #if 0
                if (adapt_data.t8_data->use_distributed_data)
                {
                    /* In parallel environment, we need to communicate the previous maximum deviations */
                    /* Within this function call, we communicate the maximum deviations, order them and return them within the adapt_data for the next iteration */
                    /* The global elementid of the first local element (of the adapt_data->forest_reference) is inquired there as well and stored withint the adapt_data for the next iteratzion */
                    cmc_rel_error_threshold_communicate_and_update_deviations(adapt_data);
                } else
                {
                    /* In a serial environment, it suffices to just update the deviations */
                    cmc_rel_error_threshold_update_deviations(adapt_data);
                }
                #endif

            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_EXCLUDE_AREA:
                //Here has nothing to be done
            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_COMBINED_CRITERION:
                /* Free the data which has been saved during the adaptation */
                adapt_data.adapted_data->clear();
                interpolation_data.adapted_data = nullptr;
            break;
            default:
                cmc_err_msg("The supplied lossy compression criterium is not yet implemented.");
        }
    }
    /* Generally, independet of the mode and criterium, the counter will be incremented */
    ++(adapt_data.adapt_step);
    #endif
}

/**
 * @brief This functions deallocates previous members of the adapt and interpolation struct when all adaptation steps are finished 
 * 
 * @param adapt_data The @struct cmc_t8_adapt_data to update/reset
 * @param interpolation_data The @struct cmc_t8_interpolation_data to update/reset
 */
static void
cmc_t8_deconstruct_adapt_and_interpolate_data(cmc_t8_adapt_data& adapt_data, cmc_t8_interpolation_data& interpolation_data)
{
    #ifdef CMC_WITH_T8CODE
    /* Set up the data structures based on the compression mode */
    if (adapt_data.t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
        adapt_data.t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D)
    {
        /* If a 'One for All' compression mode is used */
        switch (adapt_data.t8_data->settings.compression_criterium)
        {
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_CRITERIUM_UNDEFINED:
            /* If the criterium is undefined, we use by default an error threshold criterium */
            [[fallthrough]];
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_REL_ERROR_THRESHOLD:
                /* Delete the allocation of the var_vector */
                delete adapt_data.adapted_data;
            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_EXCLUDE_AREA:
                //Here has nothing to be done
            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_COMBINED_CRITERION:
                /* Delete the allocation of the var_vector */
                delete adapt_data.adapted_data;
            break;
            default:
                cmc_err_msg("The supplied lossy compression criterium is not yet implemented.");
        }
    }
    else
    {
        /* If a 'One for One' compression mode is used */
        switch (adapt_data.t8_data->settings.compression_criterium)
        {
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_CRITERIUM_UNDEFINED:
            /* If the mode is undefined, we use by default an error threshold criterium */
            [[fallthrough]];
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_REL_ERROR_THRESHOLD:
                /* Delete the allocation of the var_vector */
                delete adapt_data.adapted_data;
            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_EXCLUDE_AREA:
                //Here has nothing to be done
            break;
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_COMBINED_CRITERION:
                /* Delete the allocation of the var_vector */
                delete adapt_data.adapted_data;
            break;
            default:
                cmc_err_msg("The supplied lossy compression criterium is not yet implemented.");
        }
    }
    #endif
}

/**
 * @brief This function performs the repartitioning of the forest and the data during a compression/coarsening step
 * 
 * @param t8_data A pointer to the @struct cmc_t8_data holding all variables and data
 * @param forest The forest which will be partitioned
 * @param var_id Eventually an id indicating to which variable the forest belongs (only in a 'One For One'-mode, otherwise this parameter has no effect)
 * @note The onwership of the forest is taken, if it is not references directly before calling this funciton (@see t8_forest_ref(...))
 * @return t8_forest_t The partitioned forest (The @var t8_data holds now the partitioned data corresponding to this retunred forest)
 */
static
t8_forest_t
cmc_t8_geo_data_repartition_during_compression(cmc_t8_data_t t8_data, t8_forest_t forest, const int var_id = -1)
{
    #ifdef CMC_WITH_T8CODE
    cmc_debug_msg("Repartition forest and data during the compression.");

    /* Declare a new forest variable */
    t8_forest_t forest_partitioned;
    /* Initialize the forest */
    t8_forest_init(&forest_partitioned);

    const int partition_for_coarsening = 0; //Currently, this is not yet available

    /* We need to distinguish between the compression modes */
    if (t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
        t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D)
    {
        /* If a 'One For All' compression is used */
        /* Reference the forest */
        t8_forest_ref(forest);
        /* Set the forest from which the partition will be derived */
        t8_forest_set_partition(forest_partitioned, forest, partition_for_coarsening);
        /* Commit the forest */
        t8_forest_commit(forest_partitioned);

        /* Partition all variables */
        for (size_t var_iter = 0; var_iter < t8_data->vars.size(); ++var_iter)
        {
            /* Create an sc input array from the variable's data */
            sc_array_t in_data;

            /* Set up the sc_array of the input */
            in_data.elem_size = t8_data->vars[var_iter]->get_data_size();
            in_data.elem_count = t8_data->vars[var_iter]->var->data->size();
            in_data.byte_alloc = static_cast<ssize_t>(in_data.elem_size * in_data.elem_count);
            in_data.array = static_cast<char*>(t8_data->vars[var_iter]->var->data->get_initial_data_ptr());

            /* Allcate a new array in data_new which will hold the partitioned data */
            t8_data->vars[var_iter]->var->data_new = new var_array_t(static_cast<size_t>(t8_forest_get_local_num_elements(forest_partitioned)), t8_data->vars[var_iter]->get_type());

            /* Create an sc output array for the partitioned variable's data */
            sc_array_t out_data;

            /* Setup the output sc array */
            out_data.elem_size = t8_data->vars[var_iter]->get_data_size();
            out_data.elem_count = t8_forest_get_local_num_elements(forest_partitioned);
            out_data.byte_alloc = static_cast<ssize_t>(in_data.elem_size * t8_forest_get_local_num_elements(forest_partitioned));
            out_data.array = static_cast<char*>(t8_data->vars[var_iter]->get_initial_data_new_ptr());

            /* The data has to be partitioned as well according to new partitioning scheme of the forest */
            t8_forest_partition_data(forest, forest_partitioned, &in_data, &out_data);

            /* Switch the old with the partitioned data */
            t8_data->vars[var_iter]->var->switch_data();
        }

        /* Dereference the old partitioned forest */
        t8_forest_unref(&forest);

        cmc_debug_msg("The forest and the data has been re-partitioned.");

        /* Return the newly partitioned forest */
        return forest_partitioned;
    }
    else
    {
        /* If a 'One for One' compression mode is used */
        /* Reference the forest */
        t8_forest_ref(forest);
        /* Set the forest from which the partition will be derived */
        t8_forest_set_partition(forest_partitioned, forest, partition_for_coarsening);
        /* Commit the forest */
        t8_forest_commit(forest_partitioned);

        /* Create an sc input array from the variable's data */
        sc_array_t in_data;

        /* Set up the sc_array of the input */
        in_data.elem_size = t8_data->vars[var_id]->get_data_size();
        in_data.elem_count = t8_data->vars[var_id]->var->data->size();
        in_data.byte_alloc = static_cast<ssize_t>(in_data.elem_size * in_data.elem_count);
        in_data.array = static_cast<char*>(t8_data->vars[var_id]->var->data->get_initial_data_ptr());

        /* Allcate a new array in data_new which will hold the partitioned data */
        t8_data->vars[var_id]->var->data_new = new var_array_t(static_cast<size_t>(t8_forest_get_local_num_elements(forest_partitioned)), t8_data->vars[var_id]->get_type());

        /* Create an sc output array for the partitioned variable's data */
        sc_array_t out_data;

        /* Setup the output sc array */
        out_data.elem_size = t8_data->vars[var_id]->get_data_size();
        out_data.elem_count = t8_forest_get_local_num_elements(forest_partitioned);
        out_data.byte_alloc = static_cast<ssize_t>(in_data.elem_size * t8_forest_get_local_num_elements(forest_partitioned));
        out_data.array = static_cast<char*>(t8_data->vars[var_id]->get_initial_data_new_ptr());

        /* The data has to be partitioned as well according to new partitioning scheme of the forest */
        t8_forest_partition_data(forest, forest_partitioned, &in_data, &out_data);

        /* Switch the old with the partitioned data */
        t8_data->vars[var_id]->var->switch_data();

        /* Dereference the old partitioned forest */
        t8_forest_unref(&forest);

        cmc_debug_msg("The forest and the data has been re-partitioned.");

        /* Return the newly partitioned forest */
        return forest_partitioned;

    }
    #endif
}

/**
 * @brief This functions perform the adaptation and interpolation of the forest and the variables in an 'One For All' compression mode.
 *        This means all variables are defined on the same forest and the coarsening/compression happes analogously for all variables.
 * 
 * @param t8_data A pointer to the @struct cmc_t8_data holding the forest and the variables
 * @param adapt_function The adapt_function to use
 * @param interpolation_function The interpolation funcion to use
 */
static void
cmc_t8_adapt_interpolate_data_func_one_for_all(cmc_t8_data_t t8_data, t8_forest_adapt_t adapt_function, t8_forest_replace_t interpolation_function)
{
    #ifdef CMC_WITH_T8CODE
    cmc_debug_msg("Compression Mode: 'One for All'");
    int ref_lvl{t8_data->assets->initial_refinement_lvl};
    t8_gloidx_t num_elems_former_forest{0};

    t8_forest_t forest{t8_data->assets->forest};
    t8_forest_t coarsened_forest;

    /* Create an adapt data struct, and pass eventually supplied compression settings */
    cmc_t8_adapt_data adapt_data{t8_data};
    
    /* Create an interpolation data struct */
    cmc_t8_interpolation_data interpolation_data{t8_data};

    /* Set up the adapt data struct based on the given compression criterium and the compression mode */
    cmc_t8_set_up_adapt_data_and_interpolate_data_based_on_compression_settings(adapt_data, interpolation_data, t8_data);

    /* Reserve enough space for the initial data of all variables */ 
    t8_data->initial_data.reserve(t8_data->vars.size());
    /* Save the initial data of each variable */
    for (size_t id{0}; id < t8_data->vars.size(); ++id)
    {
        t8_data->initial_data.push_back(t8_data->vars[id]->var->data);
    }

    cmc_debug_msg("Adaptation of the forest starts.");

    /* Apply the adaptation/coarsening as often as possible */
    while (ref_lvl > 0 && num_elems_former_forest != t8_forest_get_global_num_elements(forest))
    {
        /* Update the data structurs at the beginning of a iteration */
        cmc_t8_update_adapt_and_interpolation_data_beginning_of_iteration(adapt_data, interpolation_data);

        /* Check if the data is distributed */
        if(t8_data->use_distributed_data)
        {
            /* In case we are using a parallel environment, we store the global element id of the first local element */
            adapt_data.first_global_elem_id = t8_forest_get_first_local_element_id(forest);
        }

        /** Adaptation process starts **/
        /* Keep the 'forest' after the adaptation step */
        t8_forest_ref(forest);

        /* Update the number of elements of the former forest */
        num_elems_former_forest = t8_forest_get_global_num_elements(forest);

        /* Create the adapted forest */
        coarsened_forest = t8_forest_new_adapt(forest, adapt_function, 0, 0, static_cast<void*>(&adapt_data));
        /** Adaptation process ends **/

        /** Interpolation process starts **/
        /* Iterate over all variables */
        for (size_t var_id{0}; var_id < t8_data->vars.size(); ++var_id)
        {
            /* Allocate memory equal to the new elements */
            t8_data->vars[var_id]->var->data_new = new var_array_t(static_cast<size_t>(t8_forest_get_local_num_elements(coarsened_forest)), t8_data->vars[var_id]->get_type());
        }

        /* Set the interpolation data accordingly */
        t8_forest_set_user_data(coarsened_forest, static_cast<void*>(&interpolation_data));
        
        /* Interpolate the element data onto the new coarsened forest */
        t8_forest_iterate_replace(coarsened_forest, forest, interpolation_function);
        
        /* Iterate over all variables */
        for (size_t var_id{0}; var_id < t8_data->vars.size(); ++var_id)
        {
            /* Delete the previous (fine/uncrompressed) data (except in the first iteration, because the intitial data is kept) */
            if (t8_data->assets->initial_refinement_lvl != ref_lvl)
            {
                delete t8_data->vars[var_id]->var->data;
            }
            /* Assign the (coarsened/compressed) data */
            t8_data->vars[var_id]->var->data = t8_data->vars[var_id]->var->data_new;
            t8_data->vars[var_id]->var->data_new = nullptr;
        }
        /** Interpolation process ends **/

        /* Free the former forest */
        t8_forest_unref(&forest);

        /* Switch to coarsened forest */
        forest = coarsened_forest;

        /* Decrement the refinement level */
        --ref_lvl;

        /* In case of a parallel execution with distributed data, we need to repartition the forest, as well as the data at this stage */
        if (t8_data->use_distributed_data)
        {
            /* Repartition the variable's data */
            forest = cmc_t8_geo_data_repartition_during_compression(t8_data, coarsened_forest);
        }

        /* Save the adapated and/or partitioned forest temporarily in case the adapt-interpolate function needs it */
        adapt_data.forest_reference = forest;

        /* Update the data structs at the end of a iteration */
        cmc_t8_update_adapt_and_interpolation_data_end_of_iteration(adapt_data, interpolation_data);
    }

    /* Delete allocations or perform any finalizing steps */
    cmc_t8_deconstruct_adapt_and_interpolate_data(adapt_data, interpolation_data);

    /* Free the former forest and Save the adapted forest */
    t8_data->assets->forest = forest;
    
    cmc_debug_msg("Adaptation/Compression is finished.");

    #endif
}

/**
 * @brief This functions perform the adaptation and interpolation of the forest and the variables in an 'One For One' compression mode.
 *        This means each variable is defined on the seperate forest and the coarsening/compression happes independently of the otehr variables.
 * 
 * @param t8_data A pointer to the @struct cmc_t8_data holding the variables and their forests
 * @param adapt_function The adapt_function to use
 * @param interpolation_function The interpolation funcion to use
 */
static void
cmc_t8_adapt_interpolate_data_func_one_for_one(cmc_t8_data_t t8_data, t8_forest_adapt_t adapt_function, t8_forest_replace_t interpolation_function)
{
    #ifdef CMC_WITH_T8CODE
    cmc_debug_msg("Compression Mode: 'One for One'");

    int ref_lvl;
    t8_gloidx_t num_elems_former_forest{0};

    t8_forest_t forest;
    t8_forest_t coarsened_forest;

    /* Reserve enough space for the initial data of all variables */ 
    t8_data->initial_data.reserve(t8_data->vars.size());

    /* Create an adapt data struct, and pass eventually supplied compression settings */
    cmc_t8_adapt_data adapt_data(t8_data);
    /* Allocate members in the adapt struct */
    cmc_t8_adapt_struct_allocate_one_for_one(adapt_data);

    /* Create an interpolation data struct */
    cmc_t8_interpolation_data interpolation_data(t8_data);

    /* Iterate over all variables */
    for (size_t var_id{0}; var_id < t8_data->vars.size(); ++var_id)
    {
        cmc_debug_msg("Adaptation of the forest of variable ", t8_data->vars[var_id]->var->name, " starts.");

        /* Save the initial data */
        t8_data->initial_data.push_back(t8_data->vars[var_id]->var->data);
        /* Reset to the initial refinement level (this equals the maximum number of possible coarsening steps) */
        ref_lvl = t8_data->vars[var_id]->assets->initial_refinement_lvl;
        /* Get a pointer to the forest of the variable */
        forest = t8_data->vars[var_id]->assets->forest;
        /* Reset the number of former elements in the forest */
        num_elems_former_forest = 0;

        /* Set up the adaptation and interpoaltion data based on the compression mode and on the compression settings */
        cmc_t8_set_up_adapt_data_and_interpolate_data_based_on_compression_settings(adapt_data, interpolation_data, t8_data, var_id);
    
        /* Check if the data is distributed */
        if(t8_data->use_distributed_data)
        {
            /* In case we are using a parallel environment, we store the global element id of the first local element */
            adapt_data.first_global_elem_id = t8_forest_get_first_local_element_id(forest);
        }
        
        /* Reference the initial forest */
        t8_forest_ref(forest);

        /* Apply the adaptation/coarsening as often as possible */
        while (ref_lvl > 0 && num_elems_former_forest != t8_forest_get_global_num_elements(forest))
        {
            /* Update the adapt and interpolation data at the beginning of the iteration */
            cmc_t8_update_adapt_and_interpolation_data_beginning_of_iteration(adapt_data, interpolation_data, var_id);

            /* Update the number of elements of the former forest */
            num_elems_former_forest = t8_forest_get_global_num_elements(forest);

            /* Keep the 'forest' after the adaptation step */
            t8_forest_ref(forest);

            /* Create the adapted forest */
            coarsened_forest = t8_forest_new_adapt(forest, adapt_function, 0, 0, static_cast<void*>(&adapt_data));

            /* Allocate memory equal to the new elements */
            t8_data->vars[var_id]->var->data_new = new var_array_t(static_cast<size_t>(t8_forest_get_local_num_elements(coarsened_forest)), t8_data->vars[var_id]->get_type());

            /* Set the interpolation data accordingly */
            t8_forest_set_user_data(coarsened_forest, static_cast<void*>(&interpolation_data));
            
            /* Interpolate the element data onto the new coarsened forest */
            t8_forest_iterate_replace(coarsened_forest, forest, interpolation_function);

            /* Delete and Release the previous (fine/uncrompressed) data */
            if (t8_data->vars[var_id]->assets->initial_refinement_lvl != ref_lvl)
            {
                delete t8_data->vars[var_id]->var->data;
            }

            /* Assign the (coarsened/compressed) data */
            t8_data->vars[var_id]->var->data = t8_data->vars[var_id]->var->data_new;
            t8_data->vars[var_id]->var->data_new = nullptr;

            /* Free the former forest */
            t8_forest_unref(&forest);

            /* Switch to the coarsened forest */
            forest = coarsened_forest;

            /* Decrement the refinement level */
            --ref_lvl;

            /* In case of a parallel execution with distributed data, we need to repartition the forest, as well as the data at this stage */
            if (t8_data->use_distributed_data)
            {
                /* Repartition the variable's data */
                forest = cmc_t8_geo_data_repartition_during_compression(t8_data, coarsened_forest, var_id);
            }

            /* Save the adapated and/or partitioned forest temporarily in case the adapt-interpolate function needs it */
            adapt_data.forest_reference = forest;

            /* Update the adapt and interpolation struct at the end of the iteration */
            cmc_t8_update_adapt_and_interpolation_data_end_of_iteration(adapt_data, interpolation_data);
        }
    
        /* Free the former forest and Save the adapted forest */
        t8_forest_unref(&(t8_data->vars[var_id]->assets->forest));
        t8_data->vars[var_id]->assets->forest = forest;
        
        cmc_debug_msg("Adaptation/Compression of variable ", t8_data->vars[var_id]->var->name, " is finished.");
    }

    /* Delete allocations or perform any finalizing steps */
    cmc_t8_deconstruct_adapt_and_interpolate_data(adapt_data, interpolation_data);

    #endif
}

/**
 * @brief This functions calls the correspondning adapt_and_inerpoalte function concerning the used compression mode.
 * 
 * @param t8_data A pointer to the @struct cmc_t8_data holding the forest(s) and the variables
 * @param adapt_function The adapt_function to use
 * @param interpolation_function The interpolation_function to use
 */
static void
cmc_t8_adapt_interpolate_data_func(cmc_t8_data_t t8_data, t8_forest_adapt_t adapt_function, t8_forest_replace_t interpolation_function)
{
    #ifdef CMC_WITH_T8CODE
    cmc_assert(t8_data->compression_mode != CMC_T8_COMPRESSION_MODE::CMC_T8_COMPRESSION_UNDEFINED);
    
    /* Check if a 'one for all' or 'one for one' approach is considered */
    if (t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
        t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D)
    {
        cmc_t8_adapt_interpolate_data_func_one_for_all(t8_data, adapt_function, interpolation_function);
    } else if (t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ONE_2D ||
               t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ONE_3D)
    {
        cmc_t8_adapt_interpolate_data_func_one_for_one(t8_data, adapt_function, interpolation_function);
    } else
    {
        cmc_err_msg("An unknown compression mode was supplied.");
    }
    #endif
}
#endif
/** End STATIC Functions **/
/**************************/

/** Transform 3D variables into 2D variables by spliiting up a direction */
//TODO: update for parallelization
void
cmc_geo_data_transform_3d_var_to_2d(cmc_t8_data_t t8_data, const int var_id, const DATA_LAYOUT preferred_layout)
{
    #ifdef CMC_WITH_T8CODE
    /* Check if t8_data is supplied */
    cmc_assert(t8_data != nullptr);

    /* Save the ids of the variables which will be transformed */
    std::vector<int> ids_apply_transformation;

    if (var_id == CMC_APPLY_TRANSFORMATION_3D_TO_2D_TO_ALL_VARS)
    {
        /* Transform all 3D variables to 2D variables */
        ids_apply_transformation.reserve(t8_data->vars.size());
        for (size_t id{0}; id < t8_data->vars.size(); ++id)
        {
            /* Check if some vairables are eventually already two-dimensional */
            if (t8_data->vars[id]->var->num_dimensions != 2)
            {
                cmc_debug_msg("Variable (name: ", t8_data->vars[id]->var->name, ") will be transformed to corresponding 2D variables.");
                ids_apply_transformation.push_back(id);
            }
        }
        /* Return if all variables are already two-dimensional */
        if (ids_apply_transformation.size() == 0)
        {
            cmc_debug_msg("All variables are already two-dimensional.");
            return;
        }
    }
    else
    {
        /* Transform only the variable with the given id to 2D */
        cmc_assert(var_id >= 0 && var_id < static_cast<int>(t8_data->vars.size()));
        
        /* Save the id of the given variable (if it is a 3D variable) */
        if (t8_data->vars[var_id]->var->num_dimensions != 2)
        {
            cmc_debug_msg("Variable (name: ", t8_data->vars[var_id]->var->name, ") will be transformed to corresponding 2D variables.");
            ids_apply_transformation.push_back(var_id);
        } else
        {
            cmc_debug_msg("The variable is already two-dimensional.");
            return;
        }
    }

    /* Create a vetor holding all vars after the transformation */
    std::vector<cmc_t8_var_t> new_vars;
    /* Allocate a loose upper bound */
    //new_vars.reserve(t8_data->vars.size() * (1 << t8_data->vars[]>initial_refinement_lvl));
    new_vars.reserve(t8_data->vars.size() * std::max(std::initializer_list<size_t> {t8_data->geo_data->get_coord_length(CMC_COORD_IDS::CMC_LAT), t8_data->geo_data->get_coord_length(CMC_COORD_IDS::CMC_LON), t8_data->geo_data->get_coord_length(CMC_COORD_IDS::CMC_LEV)}));
    /* Copy unchanged variables */
    if (var_id == CMC_APPLY_TRANSFORMATION_3D_TO_2D_TO_ALL_VARS)
    {
        /* Copy all 2D variables */
        for (size_t i{0}; i < t8_data->vars.size(); ++i)
        {
            if (t8_data->vars[i]->var->num_dimensions == 2)
            {
                new_vars.push_back(t8_data->vars[i]);
            }
        }
    } else
    {
        /* Copy all variables except the one which is transferred to 2D variables */
        for (size_t i{0}; i < t8_data->vars.size(); ++i)
        {
            if (i != static_cast<size_t>(var_id))
            {
                new_vars.push_back(t8_data->vars[i]);
            }
        }
    }

    /* Transform the 3D variable(s) */
    for (size_t id{0}; id < ids_apply_transformation.size(); ++id)
    {
        cmc_assert(t8_data->vars[ids_apply_transformation[id]]->var->num_dimensions == 3);
        switch(preferred_layout)
        {
            case CMC_2D_LAT_LON:
                {
                    const size_t num_elems{t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LON]};
                    for (size_t lev{0}; lev < t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                    {
                        /* Create a new 2D variable with a corresponding name */
                        new_vars.emplace_back(new cmc_t8_var(t8_data->vars[ids_apply_transformation[id]]->var->name + "_lev" + std::to_string(lev)));
                        /* Allocate memory for the 2D data */
                        new_vars.back()->var->data = new var_array_t(num_elems, t8_data->vars[ids_apply_transformation[id]]->get_type());
                        /* Save the data layout */
                        new_vars.back()->var->data_layout = DATA_LAYOUT::CMC_2D_LAT_LON;
                        new_vars.back()->var->axis_ordering = std::vector<int>{CMC_COORD_IDS::CMC_LAT, CMC_COORD_IDS::CMC_LON};
                        /* Save the ordering scheme */
                        new_vars.back()->var->data_scheme = t8_data->vars[ids_apply_transformation[id]]->var->data_scheme;
                        /* Save information regarding the dimensionality */
                        new_vars.back()->var->num_dimensions = 2;
                        new_vars.back()->var->dim_lengths = t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths;
                        new_vars.back()->var->dim_lengths[CMC_COORD_IDS::CMC_LEV] = 1;
                        new_vars.back()->var->dimension_ids = t8_data->vars[ids_apply_transformation[id]]->var->dimension_ids;
                        new_vars.back()->var->dimension_ids[CMC_COORD_IDS::CMC_LEV] = CMC_COORDINATE_NOT_CONSIDERED;
                        /* Save information about meta data */
                        new_vars.back()->var->missing_value_present = t8_data->vars[ids_apply_transformation[id]]->var->missing_value_present;
                        new_vars.back()->var->applied_offset_scaling = t8_data->vars[ids_apply_transformation[id]]->var->applied_offset_scaling;
                        new_vars.back()->var->missing_value = t8_data->vars[ids_apply_transformation[id]]->var->missing_value;
                        new_vars.back()->var->add_offset = t8_data->vars[ids_apply_transformation[id]]->var->add_offset;
                        new_vars.back()->var->scale_factor = t8_data->vars[ids_apply_transformation[id]]->var->scale_factor;
                        new_vars.back()->var->var_type = t8_data->vars[ids_apply_transformation[id]]->var->var_type;
                        /* Set the data of the new 2D variable */
                        cmc_geo_data_fetch_2d_data(*(t8_data->vars[ids_apply_transformation[id]]), *(new_vars.back()), lev); 
                    }
                }
            break;
            case CMC_2D_LON_LAT:
                {
                    const size_t num_elems{t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LON]};
                    for (size_t lev{0}; lev < t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LEV]; ++lev)
                    {
                        /* Create a new 2D variable with a corresponding name */
                        new_vars.emplace_back(new cmc_t8_var(t8_data->vars[ids_apply_transformation[id]]->var->name + "_lev" + std::to_string(lev)));
                        /* Allocate memory for the 2D data */
                        new_vars.back()->var->data = new var_array_t(num_elems, t8_data->vars[ids_apply_transformation[id]]->get_type());
                        /* Save the data layout */
                        new_vars.back()->var->data_layout = DATA_LAYOUT::CMC_2D_LON_LAT;
                        new_vars.back()->var->axis_ordering = std::vector<int>{CMC_COORD_IDS::CMC_LON, CMC_COORD_IDS::CMC_LAT};
                        /* Save the ordering scheme */
                        new_vars.back()->var->data_scheme = t8_data->vars[ids_apply_transformation[id]]->var->data_scheme;
                        /* Save information reagarding the dimensionality */
                        new_vars.back()->var->num_dimensions = 2;
                        new_vars.back()->var->dim_lengths = t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths;
                        new_vars.back()->var->dim_lengths[CMC_COORD_IDS::CMC_LEV] = 1;
                        new_vars.back()->var->dimension_ids = t8_data->vars[ids_apply_transformation[id]]->var->dimension_ids;
                        new_vars.back()->var->dimension_ids[CMC_COORD_IDS::CMC_LEV] = CMC_COORDINATE_NOT_CONSIDERED;
                        /* Save information about meta data */
                        new_vars.back()->var->missing_value_present = t8_data->vars[ids_apply_transformation[id]]->var->missing_value_present;
                        new_vars.back()->var->applied_offset_scaling = t8_data->vars[ids_apply_transformation[id]]->var->applied_offset_scaling;
                        new_vars.back()->var->missing_value = t8_data->vars[ids_apply_transformation[id]]->var->missing_value;
                        new_vars.back()->var->add_offset = t8_data->vars[ids_apply_transformation[id]]->var->add_offset;
                        new_vars.back()->var->scale_factor = t8_data->vars[ids_apply_transformation[id]]->var->scale_factor;
                        new_vars.back()->var->var_type = t8_data->vars[ids_apply_transformation[id]]->var->var_type;
                        /* Set the data of the new 2D variable */
                        cmc_geo_data_fetch_2d_data(*(t8_data->vars[ids_apply_transformation[id]]), *(new_vars.back()), lev); 
                    }
                }
            break;
            case CMC_2D_LAT_LEV:
                {
                    const size_t num_elems{t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LEV]};
                    for (size_t lon{0}; lon < t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                    {
                        /* Create a new 2D variable with a corresponding name */
                        new_vars.emplace_back(new cmc_t8_var(t8_data->vars[ids_apply_transformation[id]]->var->name + "_lon" + std::to_string(lon)));
                        /* Allocate memory for the 2D data */
                        new_vars.back()->var->data = new var_array_t(num_elems, t8_data->vars[ids_apply_transformation[id]]->get_type());
                        /* Save the data layout */
                        new_vars.back()->var->data_layout = DATA_LAYOUT::CMC_2D_LAT_LEV;
                        new_vars.back()->var->axis_ordering = std::vector<int>{CMC_COORD_IDS::CMC_LAT, CMC_COORD_IDS::CMC_LEV};
                        /* Save the ordering scheme */
                        new_vars.back()->var->data_scheme = t8_data->vars[ids_apply_transformation[id]]->var->data_scheme;
                        /* Save information reagarding the dimensionality */
                        new_vars.back()->var->num_dimensions = 2;
                        new_vars.back()->var->dim_lengths = t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths;
                        new_vars.back()->var->dim_lengths[CMC_COORD_IDS::CMC_LON] = 1;
                        new_vars.back()->var->dimension_ids = t8_data->vars[ids_apply_transformation[id]]->var->dimension_ids;
                        new_vars.back()->var->dimension_ids[CMC_COORD_IDS::CMC_LON] = CMC_COORDINATE_NOT_CONSIDERED;
                        /* Save information about meta data */
                        new_vars.back()->var->missing_value_present = t8_data->vars[ids_apply_transformation[id]]->var->missing_value_present;
                        new_vars.back()->var->applied_offset_scaling = t8_data->vars[ids_apply_transformation[id]]->var->applied_offset_scaling;
                        new_vars.back()->var->missing_value = t8_data->vars[ids_apply_transformation[id]]->var->missing_value;
                        new_vars.back()->var->add_offset = t8_data->vars[ids_apply_transformation[id]]->var->add_offset;
                        new_vars.back()->var->scale_factor = t8_data->vars[ids_apply_transformation[id]]->var->scale_factor;
                        new_vars.back()->var->var_type = t8_data->vars[ids_apply_transformation[id]]->var->var_type;
                        /* Set the data of the new 2D variable */
                        cmc_geo_data_fetch_2d_data(*(t8_data->vars[ids_apply_transformation[id]]), *(new_vars.back()), lon); 
                    }
                }
            break;
            case CMC_2D_LEV_LAT:
                {
                    const size_t num_elems{t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LAT] * t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LEV]};
                    for (size_t lon{0}; lon < t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LON]; ++lon)
                    {
                        /* Create a new 2D variable with a corresponding name */
                        new_vars.emplace_back(new cmc_t8_var(t8_data->vars[ids_apply_transformation[id]]->var->name + "_lon" + std::to_string(lon)));
                        /* Allocate memory for the 2D data */
                        new_vars.back()->var->data = new var_array_t(num_elems, t8_data->vars[ids_apply_transformation[id]]->get_type());
                        /* Save the data layout */
                        new_vars.back()->var->data_layout = DATA_LAYOUT::CMC_2D_LEV_LAT;
                        new_vars.back()->var->axis_ordering = std::vector<int>{CMC_COORD_IDS::CMC_LEV, CMC_COORD_IDS::CMC_LAT};
                        /* Save the ordering scheme */
                        new_vars.back()->var->data_scheme = t8_data->vars[ids_apply_transformation[id]]->var->data_scheme;
                        /* Save information reagarding the dimensionality */
                        new_vars.back()->var->num_dimensions = 2;
                        new_vars.back()->var->dim_lengths = t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths;
                        new_vars.back()->var->dim_lengths[CMC_COORD_IDS::CMC_LON] = 1;
                        new_vars.back()->var->dimension_ids = t8_data->vars[ids_apply_transformation[id]]->var->dimension_ids;
                        new_vars.back()->var->dimension_ids[CMC_COORD_IDS::CMC_LON] = CMC_COORDINATE_NOT_CONSIDERED;
                        /* Save information about meta data */
                        new_vars.back()->var->missing_value_present = t8_data->vars[ids_apply_transformation[id]]->var->missing_value_present;
                        new_vars.back()->var->applied_offset_scaling = t8_data->vars[ids_apply_transformation[id]]->var->applied_offset_scaling;
                        new_vars.back()->var->missing_value = t8_data->vars[ids_apply_transformation[id]]->var->missing_value;
                        new_vars.back()->var->add_offset = t8_data->vars[ids_apply_transformation[id]]->var->add_offset;
                        new_vars.back()->var->scale_factor = t8_data->vars[ids_apply_transformation[id]]->var->scale_factor;
                        new_vars.back()->var->var_type = t8_data->vars[ids_apply_transformation[id]]->var->var_type;
                        /* Set the data of the new 2D variable */
                        cmc_geo_data_fetch_2d_data(*(t8_data->vars[ids_apply_transformation[id]]), *(new_vars.back()), lon); 
                    }
                }
            break;
            case CMC_2D_LON_LEV:
                {
                    const size_t num_elems{t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LON] * t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LEV]};
                    for (size_t lat{0}; lat < t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                    {
                        /* Create a new 2D variable with a corresponding name */
                        new_vars.emplace_back(new cmc_t8_var(t8_data->vars[ids_apply_transformation[id]]->var->name + "_lat" + std::to_string(lat)));
                        /* Allocate memory for the 2D data */
                        new_vars.back()->var->data = new var_array_t(num_elems, t8_data->vars[ids_apply_transformation[id]]->get_type());
                        /* Save the data layout */
                        new_vars.back()->var->data_layout = DATA_LAYOUT::CMC_2D_LON_LEV;
                        new_vars.back()->var->axis_ordering = std::vector<int>{CMC_COORD_IDS::CMC_LON, CMC_COORD_IDS::CMC_LEV};
                        /* Save the ordering scheme */
                        new_vars.back()->var->data_scheme = t8_data->vars[ids_apply_transformation[id]]->var->data_scheme;
                        /* Save information reagarding the dimensionality */
                        new_vars.back()->var->num_dimensions = 2;
                        new_vars.back()->var->dim_lengths = t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths;
                        new_vars.back()->var->dim_lengths[CMC_COORD_IDS::CMC_LAT] = 1;
                        new_vars.back()->var->dimension_ids = t8_data->vars[ids_apply_transformation[id]]->var->dimension_ids;
                        new_vars.back()->var->dimension_ids[CMC_COORD_IDS::CMC_LAT] = CMC_COORDINATE_NOT_CONSIDERED;
                        /* Save information about meta data */
                        new_vars.back()->var->missing_value_present = t8_data->vars[ids_apply_transformation[id]]->var->missing_value_present;
                        new_vars.back()->var->applied_offset_scaling = t8_data->vars[ids_apply_transformation[id]]->var->applied_offset_scaling;
                        new_vars.back()->var->missing_value = t8_data->vars[ids_apply_transformation[id]]->var->missing_value;
                        new_vars.back()->var->add_offset = t8_data->vars[ids_apply_transformation[id]]->var->add_offset;
                        new_vars.back()->var->scale_factor = t8_data->vars[ids_apply_transformation[id]]->var->scale_factor;
                        new_vars.back()->var->var_type = t8_data->vars[ids_apply_transformation[id]]->var->var_type;
                        /* Set the data of the new 2D variable */
                        cmc_geo_data_fetch_2d_data(*(t8_data->vars[ids_apply_transformation[id]]), *(new_vars.back()), lat); 
                    }
                }
            break;
            case CMC_2D_LEV_LON:
                {
                    const size_t num_elems{t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LON] * t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LEV]};
                    for (size_t lat{0}; lat < t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths[CMC_COORD_IDS::CMC_LAT]; ++lat)
                    {
                        /* Create a new 2D variable with a corresponding name */
                        new_vars.emplace_back(new cmc_t8_var(t8_data->vars[ids_apply_transformation[id]]->var->name + "_lat" + std::to_string(lat)));
                        /* Allocate memory for the 2D data */
                        new_vars.back()->var->data = new var_array_t(num_elems, t8_data->vars[ids_apply_transformation[id]]->get_type());
                        /* Save the data layout */
                        new_vars.back()->var->data_layout = DATA_LAYOUT::CMC_2D_LEV_LON;
                        new_vars.back()->var->axis_ordering = std::vector<int>{CMC_COORD_IDS::CMC_LEV, CMC_COORD_IDS::CMC_LON};
                        /* Save the ordering scheme */
                        new_vars.back()->var->data_scheme = t8_data->vars[ids_apply_transformation[id]]->var->data_scheme;
                        /* Save information reagarding the dimensionality */
                        new_vars.back()->var->num_dimensions = 2;
                        new_vars.back()->var->dim_lengths = t8_data->vars[ids_apply_transformation[id]]->var->dim_lengths;
                        new_vars.back()->var->dim_lengths[CMC_COORD_IDS::CMC_LAT] = 1;
                        new_vars.back()->var->dimension_ids = t8_data->vars[ids_apply_transformation[id]]->var->dimension_ids;
                        new_vars.back()->var->dimension_ids[CMC_COORD_IDS::CMC_LAT] = CMC_COORDINATE_NOT_CONSIDERED;
                        /* Save information about meta data */
                        new_vars.back()->var->missing_value_present = t8_data->vars[ids_apply_transformation[id]]->var->missing_value_present;
                        new_vars.back()->var->applied_offset_scaling = t8_data->vars[ids_apply_transformation[id]]->var->applied_offset_scaling;
                        new_vars.back()->var->missing_value = t8_data->vars[ids_apply_transformation[id]]->var->missing_value;
                        new_vars.back()->var->add_offset = t8_data->vars[ids_apply_transformation[id]]->var->add_offset;
                        new_vars.back()->var->scale_factor = t8_data->vars[ids_apply_transformation[id]]->var->scale_factor;
                        new_vars.back()->var->var_type = t8_data->vars[ids_apply_transformation[id]]->var->var_type;
                        /* Set the data of the new 2D variable */
                        cmc_geo_data_fetch_2d_data(*(t8_data->vars[ids_apply_transformation[id]]), *(new_vars.back()), lat); 
                    }
                }
            break;
            default:
                cmc_err_msg("There was no valid 2D data layout supplied.");
        }
    }

    /* Free/Deallocate the 3D variables */
    if (var_id == CMC_APPLY_TRANSFORMATION_3D_TO_2D_TO_ALL_VARS)
    {
        for (size_t i{0}; i < t8_data->vars.size(); ++i)
        {
            if (t8_data->vars[i]->var->num_dimensions == 3)
            {
                delete t8_data->vars[i];
            }
        }
    } else
    {
        delete t8_data->vars[var_id];
    }

    /* Switch the old variables vector with the new variables vector */
    std::swap(t8_data->vars, new_vars);

    /* Check if all variables are now defined on the same domain and the maximim dimensionality */
    /* Save the data layout of the first variable */
    DATA_LAYOUT layout{t8_data->vars[0]->var->data_layout};
    int max_dim{0};
    t8_data->variables_are_defined_on_the_same_domain = true;
    /* Check if all variables are defined in the same domain now */
    for (size_t id{1}; id < t8_data->vars.size(); ++id)
    {
        cmc_assert(t8_data->vars[id]->var->num_dimensions == 2 || t8_data->vars[id]->var->num_dimensions == 3);
        /* If the variables have a different data layout, several meshes will be needed in order to perform the compresseion */
        if (!compare_geo_domain_equality_of_data_layouts(layout, t8_data->vars[id]->var->data_layout))
        {
            /* Switch the flag indicating if all variables are defined on the same domain to false */
            t8_data->variables_are_defined_on_the_same_domain = false;
        }
        if (max_dim < t8_data->vars[id]->var->num_dimensions)
        {
            max_dim = t8_data->vars[id]->var->num_dimensions;
        }
    }
    /* Update the geo_data dimensionality */
    t8_data->geo_data->dim = max_dim;
    #endif
}

void
cmc_t8_apply_zcurve_ordering(cmc_t8_data& t8_data, const int var_id)
{
    #ifdef CMC_WITH_T8CODE
    /* Vector storing the ids of the variables which will be reordered */
    std::vector<int> ids_apply_reordering;

    /* Check whether all variables should be reordered or just a single one */
    if (var_id != CMC_APPLY_ZCURVE_TO_ALL_VARS)
    {
        /* In case only one variable should be reordered */
        cmc_assert(var_id >= 0 && var_id >= static_cast<int>(t8_data.vars.size()));
        /* Check if the data already is ordered compliant to the z-curve ordering */
        if (t8_data.vars[var_id]->var->data_scheme == CMC_DATA_ORDERING_SCHEME::CMC_GEO_DATA_Z_CURVE)
        {
            cmc_debug_msg("Data of variable ", t8_data.vars[var_id]->var->name, " is already in z-curve order.");
            return; 
        } else {
            /* Save the id in the vector of ids to which corresponding variables the reordering have to be applied */
            cmc_debug_msg("Variable (name: ", t8_data.vars[var_id]->var->name, ") will be reorderd compliant to the z-curve ordering.");
            ids_apply_reordering.push_back(var_id);
        }
    } else {
        /* If all variables should be considered */
        ids_apply_reordering.reserve(t8_data.vars.size());
        for (size_t id{0}; id < t8_data.vars.size(); ++id)
        {
            /* Check if some vairables are eventually already in z-curve ordering */
            if (t8_data.vars[id]->var->data_scheme != CMC_DATA_ORDERING_SCHEME::CMC_GEO_DATA_Z_CURVE)
            {
                cmc_debug_msg("Variable (name: ", t8_data.vars[id]->var->name, ") will be reorderd compliant to the z-curve ordering.");
                ids_apply_reordering.push_back(id);
            }
        }
        /* Return if all variables are already in z-curve order */
        if (ids_apply_reordering.size() == 0)
        {
            cmc_debug_msg("All variables are already ordered compliant to the z-curve order.");
            return;
        }
    }

    /* Check which compression_mode applies */
    switch (t8_data.compression_mode)
    {
        case CMC_T8_COMPRESSION_MODE::CMC_T8_COMPRESSION_UNDEFINED:
            cmc_err_msg("No compression mode was supplied. Therefore, the elements cannot be reordered by 'cmc_t8_apply_zcurve_ordering'()");
            break;
        case CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D:
            cmc_debug_msg("The variables will be reorderd (compliant to the z-curve order). The selected compression mode was 'one for all'");
            cmc_t8_apply_zcurve_ordering_one_for_all(t8_data, ids_apply_reordering, 2);
            break;
        case CMC_T8_COMPRESSION_MODE::ONE_FOR_ONE_2D:
            cmc_debug_msg("The variables will be reorderd (compliant to the z-curve order). The selected compression mode was 'one for one'");
            cmc_t8_apply_zcurve_ordering_one_for_one(t8_data, ids_apply_reordering, 2);
            break;
        case CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D:
            cmc_debug_msg("The variables will be reorderd (compliant to the z-curve order). The selected compression mode was 'one for all'");
            cmc_t8_apply_zcurve_ordering_one_for_all(t8_data, ids_apply_reordering, 3);
            break;
        case CMC_T8_COMPRESSION_MODE::ONE_FOR_ONE_3D:
            cmc_debug_msg("The variables will be reorderd (compliant to the z-curve order). The selected compression mode was 'one for one'");
            cmc_t8_apply_zcurve_ordering_one_for_one(t8_data, ids_apply_reordering, 3);
            break;
        case CMC_T8_COMPRESSION_MODE::GROUPED_2D:
            cmc_err_msg("Not implemented yet");
            break;
        case CMC_T8_COMPRESSION_MODE::GROUPED_3D:
            cmc_err_msg("Not implemented yet");
            break;
        
        default:
            cmc_err_msg("An unknown compression mode was supplied. Therefore, the elements cannot be reordered by 'cmc_t8_apply_zcurve_ordering'()");
    }
    #endif
}

void
cmc_t8_coarsen_data(cmc_t8_data_t t8_data, t8_forest_adapt_t adapt_function, t8_forest_replace_t interpolation_function)
{
    #ifdef CMC_WITH_T8CODE
    /** All variable data has to be in z-curve ordering.
     *  If the data is already properly ordered, this functions returns immediately.
    **/
    cmc_t8_apply_zcurve_ordering(*t8_data, CMC_APPLY_ZCURVE_TO_ALL_VARS);

    /** Interpolate and adapt the variables' data corresponding to the
     *  supplied 'adapt'- and 'interpolation'-function
    **/
    cmc_t8_adapt_interpolate_data_func(t8_data, adapt_function, interpolation_function);

    #endif
}


void
cmc_t8_apply_offset_and_scaling(cmc_t8_data_t t8_data, const int var_id)
{
    #ifdef CMC_WITH_T8CODE

    /* Vector storing the ids of the variables which will be scaled */
    std::vector<int> ids_apply_scaling;
    /* Check whether all variables should be scaled or just a single one */
    if (var_id != CMC_APPLY_OFFSET_AND_SCALING_TO_ALL_VARS)
    {
        /* In case only one variable should be scaled */
        cmc_assert(var_id >= 0 && var_id < static_cast<int>(t8_data->vars.size()));
        /* Check if the data of the variable is already scaled */
        if (t8_data->vars[var_id]->var->applied_offset_scaling)
        {
            cmc_debug_msg("Data of variable ", t8_data->vars[var_id]->var->name, " is already scaled.");
            return; 
        } else {
            /* Save the id in the vector of ids to which corresponding variables the scaling and offset will be applied */
            ids_apply_scaling.push_back(var_id);
        }
    } else {
        /* If all variables should be considered */
        ids_apply_scaling.reserve(t8_data->vars.size());
        for (size_t id{0}; id < t8_data->vars.size(); ++id)
        {
            /* Check if some vairables are eventually already scaled */
            if (!t8_data->vars[id]->var->applied_offset_scaling)
            {
                ids_apply_scaling.push_back(id);
            }
        }
        /* Return if all variables are already in z-curve order */
        if (ids_apply_scaling.size() == 0)
        {
            cmc_debug_msg("All variables are already scaled.");
            return;
        }
    }

    /* Define standard offset and scaling variables */
    double add_offset{0.0};
    double scale_factor{1.0};

    /* Iterate over all variables which should be scaled */
    for (size_t var_id{0}; var_id < ids_apply_scaling.size(); ++var_id)
    {
        /* Get the scale_factor and the offset of the for the data of the variable with the id 'ids_apply_scaling[var_id]' */
        scale_factor = cmc_get_universal_data<double>(t8_data->vars[ids_apply_scaling[var_id]]->var->scale_factor);
        add_offset = cmc_get_universal_data<double>(t8_data->vars[ids_apply_scaling[var_id]]->var->add_offset);

        /* Check if the meta data of the variable holds offset and scaling attributes and if so, if they already have been applied */
        if (!(cmc_approx_cmp(add_offset, 0.0) && cmc_approx_cmp(scale_factor, 1.0)))
        {
            cmc_debug_msg("Scaling and offset are applied to variable ", t8_data->vars[ids_apply_scaling[var_id]]->var->name);

            /* Check if missing values are present */
            if (t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value_present)
            {
                /** If missing values are present within the data **/
                /* Check if only scaling needs to be applied */
                if ((cmc_approx_cmp(add_offset, 0.0) && (!cmc_approx_cmp(scale_factor, 1.0))))
                {
                    /* Apply only scaling */
                    t8_data->vars[ids_apply_scaling[var_id]]->var->data->scale_with_missing_vals(scale_factor, t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value);
                    /* Check if data type has changed */
                    if (t8_data->vars[ids_apply_scaling[var_id]]->get_type() != static_cast<cmc_type>(t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value.index()))
                    {
                        cmc_debug_msg("The misssing value will be casted, since the scaling changed the data type.");
                        t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value = std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value);
                    }

                }
                /* Check if only an offset needs to be applied */
                else if ((!cmc_approx_cmp(add_offset, 0.0) && (cmc_approx_cmp(scale_factor, 1.0))))
                {
                    /* Apply only offset */
                    t8_data->vars[ids_apply_scaling[var_id]]->var->data->add_const_with_missing_vals(add_offset, t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value);
                    /* Check if data type has changed */
                    if (t8_data->vars[ids_apply_scaling[var_id]]->get_type() != static_cast<cmc_type>(t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value.index()))
                    {
                        cmc_debug_msg("The misssing value will be casted, since the scaling changed the data type.");
                        t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value = std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value);
                    }
                }
                /* Otherwise a scaling as well as an offset needs to be applied */
                else
                {
                    /* Apply offset and scaling */
                    t8_data->vars[ids_apply_scaling[var_id]]->var->data->axpy_scalar_with_missing_vals(scale_factor, add_offset, t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value);
                    /* Check if data type has changed */
                    if (t8_data->vars[ids_apply_scaling[var_id]]->get_type() != static_cast<cmc_type>(t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value.index()))
                    {
                        cmc_debug_msg("The misssing value will be casted, since the scaling changed the data type.");
                        t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value = std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value);
                    }
                } 
            }
            /* In case NO missing values are present within the data */
            else {
                /* Check if only scaling needs to be applied */
                if ((cmc_approx_cmp(add_offset, 0.0) && (!cmc_approx_cmp(scale_factor, 1.0))))
                {
                    /* Apply only scaling */
                    t8_data->vars[ids_apply_scaling[var_id]]->var->data->scale(t8_data->vars[ids_apply_scaling[var_id]]->var->scale_factor);
                    /* Check if data type has changed */
                    if (t8_data->vars[ids_apply_scaling[var_id]]->get_type() != static_cast<cmc_type>(t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value.index()))
                    {
                        cmc_debug_msg("The misssing value will be casted, since the scaling changed the data type.");
                        t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value = std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value);
                    }
                }
                /* Check of only an offset needs to be applied */
                else if ((!cmc_approx_cmp(add_offset, 0.0) && (cmc_approx_cmp(scale_factor, 1.0))))
                {
                    /* Apply only offset */
                    t8_data->vars[ids_apply_scaling[var_id]]->var->data->add_const(t8_data->vars[ids_apply_scaling[var_id]]->var->add_offset);
                    /* Check if data type has changed */
                    if (t8_data->vars[ids_apply_scaling[var_id]]->get_type() != static_cast<cmc_type>(t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value.index()))
                    {
                        cmc_debug_msg("The misssing value will be casted, since the scaling changed the data type.");
                        t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value = std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value);
                    }
                }
                /* In case sclaing as well as an offset needs to be applied */
                else
                {
                    /* Apply offset and scaling */
                    t8_data->vars[ids_apply_scaling[var_id]]->var->data->axpy_scalar(t8_data->vars[ids_apply_scaling[var_id]]->var->scale_factor, t8_data->vars[ids_apply_scaling[var_id]]->var->add_offset);
                    /* Check if data type has changed */
                    if (t8_data->vars[ids_apply_scaling[var_id]]->get_type() != static_cast<cmc_type>(t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value.index()))
                    {
                        cmc_debug_msg("The misssing value will be casted, since the scaling changed the data type.");

                        t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value = std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, t8_data->vars[ids_apply_scaling[var_id]]->var->missing_value);
                    }
                }
            }
        }
        /* Set the flag */
        t8_data->vars[ids_apply_scaling[var_id]]->var->applied_offset_scaling = true;
    }
    #endif
}

void
cmc_t8_refine_to_initial_level(cmc_t8_data_t t8_data)
{
    #ifdef CMC_WITH_T8CODE
    t8_gloidx_t num_elems_former_forest{0};

    t8_forest_t forest;
    t8_forest_t adapted_forest;

    cmc_t8_adapt_data adapt_data{t8_data};
    cmc_t8_interpolation_data interpolation_data{t8_data};

    if (t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
        t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D)
    {
        forest = t8_data->assets->forest;
        cmc_debug_msg("De-Compression of variables starts...");
        /* Start the decompression */
        while (num_elems_former_forest < t8_forest_get_global_num_elements(forest))
        {
            /* Keep the 'forest' after the adaptation step */
            t8_forest_ref(forest);

            /* Update the number of elements of the former forest */
            num_elems_former_forest = t8_forest_get_global_num_elements(forest);

            /* Create the adapted forest */
            adapted_forest = t8_forest_new_adapt(forest, cmc_t8_adapt_callback_refine_to_initial_lvl, 0, 0, static_cast<void*>(&adapt_data));

            /** Interpolation process starts **/
            for (size_t var_id{0}; var_id < t8_data->vars.size(); ++var_id)
            {
                /* Allocate memory equal to the new elements */
                t8_data->vars[var_id]->var->data_new = new var_array_t(static_cast<size_t>(t8_forest_get_local_num_elements(adapted_forest)), t8_data->vars[var_id]->get_type());   
            }

            /* Set the interpolation data accordingly */
            t8_forest_set_user_data(adapted_forest, static_cast<void*>(&interpolation_data));
            
            /* Interpolate the element data onto the new coarsened forest */
            t8_forest_iterate_replace(adapted_forest, forest, cmc_t8_geo_data_interpolate_plain_copy_values);

            /* Iterate over all variables */
            for (size_t var_id{0}; var_id < t8_data->vars.size(); ++var_id)
            {
                /* Deallocate previous (fine/uncompressed) data */
                delete t8_data->vars[var_id]->var->data;
                /* Assign the (coarsened/compressed) data */
                t8_data->vars[var_id]->var->data = t8_data->vars[var_id]->var->data_new;
            }
            #if 0
            /* Iterate over all variables */
            for (size_t var_id{0}; var_id < t8_data->vars.size(); ++var_id)
            {

                /* Allocate memory equal to the new elements */
                t8_data->vars[var_id]->var->data_new = new var_array_t(static_cast<size_t>(t8_forest_get_local_num_elements(adapted_forest)), t8_data->vars[var_id]->get_type());
                
                /* Set the interpolation data accordingly */
                t8_forest_set_user_data(adapted_forest, static_cast<void*>(&interpolation_data));
                /* Interpolate the element data onto the new coarsened forest */
                t8_forest_iterate_replace(adapted_forest, forest, cmc_t8_geo_data_interpolate_plain_copy_values);

                /* Deallocate previous (fine/uncompressed) data */
                delete t8_data->vars[var_id]->var->data;
                /* Assign the (coarsened/compressed) data */
                t8_data->vars[var_id]->var->data = t8_data->vars[var_id]->var->data_new;
            }
            #endif
            /** Interpolation process ends **/

            /* Free the former forest */
            t8_forest_unref(&forest);

            /* Switch to coarsened forest */
            forest = adapted_forest;
        }

        /* Free the former forest and Save the adapted forest */
        t8_data->assets->forest = forest;
        
        cmc_debug_msg("Adaptation/De-Compression is finished.");

    } 
    else if (t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ONE_2D ||
             t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ONE_3D)
    {
        /* Since every variable define its own mesh, we iterate over all variables */
        for (size_t var_id{0}; var_id < t8_data->vars.size(); ++var_id)
        {
            /* Save the current variable ID */
            adapt_data.current_var_id = var_id;
            interpolation_data.current_var_id = var_id;

            /* Get a pointer to the forest of the variable */
            forest = t8_data->vars[var_id]->assets->forest;

            /* Reset the number of former elements in the forest */
            num_elems_former_forest = 0;

            /* Apply the adaptation/coarsening as often as possible */
            while (num_elems_former_forest < t8_forest_get_global_num_elements(forest))
            {
                /* Update the number of elements of the former forest */
                num_elems_former_forest = t8_forest_get_global_num_elements(forest);

                /* Keep the 'forest' after the adaptation step */
                t8_forest_ref(forest);

                /* Create the adapted forest */
                adapted_forest = t8_forest_new_adapt(forest, cmc_t8_adapt_callback_refine_to_initial_lvl, 0, 0, static_cast<void*>(&adapt_data));

                /* Allocate memory equal to the new elements */
                t8_data->vars[var_id]->var->data_new = new var_array_t(static_cast<size_t>(t8_forest_get_local_num_elements(adapted_forest)), t8_data->vars[var_id]->get_type());

                /* Set the interpolation data accordingly */
                t8_forest_set_user_data(adapted_forest, static_cast<void*>(&interpolation_data));
                /* Interpolate the element data onto the new coarsened forest */
                t8_forest_iterate_replace(adapted_forest, forest, cmc_t8_geo_data_interpolate_plain_copy_values);

                /* Delete and Release the previous (fine/uncrompressed) data */
                delete t8_data->vars[var_id]->var->data;

                /* Assign the (coarsened/compressed) data */
                t8_data->vars[var_id]->var->data = t8_data->vars[var_id]->var->data_new;

                /* Free the former forest */
                t8_forest_unref(&forest);

                /* Switch to coarsened forest */
                forest = adapted_forest;
            }

            /* Free the former forest and Save the adapted forest */
            t8_data->vars[var_id]->assets->forest = forest;

            cmc_debug_msg("Adaptation/De-Compression of variable ", t8_data->vars[var_id]->var->name, " is finished.");
        }
    } else
    {
        cmc_err_msg("cmc_t8_data holds an unsupported compression mode");
    }

    #endif
}

void
cmc_t8_geo_data_set_error_criterium(cmc_t8_data_t t8_data, const double maximum_error_tolerance)
{
    #ifdef CMC_WITH_T8CODE
    /* Save the error tolerance */
    t8_data->settings.max_err = maximum_error_tolerance;

    /* Set a flag that the exclude criterion will be applied */
    if (t8_data->settings.compression_criterium == CMC_T8_COMPRESSION_CRITERIUM::CMC_CRITERIUM_UNDEFINED)
    {
        /* If no criterion has been specified, set the relative error criterion */
        t8_data->settings.compression_criterium = CMC_T8_COMPRESSION_CRITERIUM::CMC_REL_ERROR_THRESHOLD;
    } else if (t8_data->settings.compression_criterium != CMC_T8_COMPRESSION_CRITERIUM::CMC_REL_ERROR_THRESHOLD)
    {
        /* If another criterion already has been specified, set the combined flag */
        t8_data->settings.compression_criterium = CMC_T8_COMPRESSION_CRITERIUM::CMC_COMBINED_CRITERION;
    }
    #endif
}

template<typename T>
static std::pair<int,int>
cmc_t8_geo_data_get_area_threshold(const var_array_t& coordinate, T start_val, T end_val)
{
    #ifdef CMC_WITH_T8CODE
    /* Obtain the smaller and the bigger value of the threshold */
    const T smaller_val = start_val <= end_val ? start_val : end_val;
    const T bigger_val = start_val > end_val ? start_val : end_val;

    /* Check if there is a 'real' threshold supplied */
    if(bigger_val - smaller_val <= 0)
    {
        /* In this case there will be no threshold applied for the given coordinate dimension */
        cmc_warn_msg("The supplied exclude-area threshold will have no effect (please check the given min/max values for the desired area to be excluded from the compression).");
        return std::make_pair<int, int>(-1, INT_MAX);
    }

    bool flag_descending_order{false};

    /* Get the pointer to the coordinate data */
    T* coord_ptr = static_cast<T*>(coordinate.get_initial_data_ptr());

    /* Declare a return value */
    std::pair<int, int> ret_val{-1, INT_MAX};

    /* Check the ordering of the coordinate data */
    size_t iter{0};
    while (iter + 1 < coordinate.size())
    {
        if (*(coord_ptr + iter) == *(coord_ptr + iter + 1))
        {
            continue;
        } else if (*(coord_ptr + iter) > *(coord_ptr + iter + 1))
        {
            /* If the order is descending */
            flag_descending_order = true;
            break;
        } else
        {
            /* If anascending order is found, we just break since the flag's default is false */
            break;
        }
        
        /* Increment the iterator variable */
        ++iter;
    }

    /* Find the integer values resembling the actual domain of the threshold */
    if (flag_descending_order)
    {
        /* In case of a descending order */
        size_t i{0};
        for (; i < coordinate.size(); ++i)
        {
            if (*(coord_ptr + i) <= bigger_val)
            {
                /* We have found the start index */
                ret_val.second = static_cast<int>(coordinate.size() - i);
                break;
            }
        }
        /* Check for the end index */
        for (; i < coordinate.size(); ++i)
        {
            if (*(coord_ptr + i) <= smaller_val)
            {
                /* We have found the end index */
                ret_val.first = static_cast<int>(coordinate.size() - i);
                break;
            }
        }
    }
    /* In case of an ascending order */
    else
    {
        size_t i{0};
        for (; i < coordinate.size(); ++i)
        {
            if (*(coord_ptr + i) >= smaller_val)
            {
                /* We have found the start index */
                ret_val.first = static_cast<int>(i);
                break;
            }
        }
        /* Check for the end index */
        for (; i < coordinate.size(); ++i)
        {
            if (*(coord_ptr + i) >= bigger_val)
            {
                /* We have found the end index */
                ret_val.second = static_cast<int>(i);
                break;
            }
        }
    }

    return ret_val;
    #endif
}

// Currently, we assume that coordinates are ordered incrementally, e.g. longitude: -90, -85, ...,-5, 0, 5, ..., 85, 90 
void
cmc_t8_geo_data_set_exclude_area(cmc_t8_data_t t8_data, const CMC_COORD_IDS coord_id, const cmc_universal_type_t& starting_value, const cmc_universal_type_t& end_value)
{
    #ifdef CMC_WITH_T8CODE
    /* Currently, only possible for longitude, latitiude and elevation */
    cmc_assert(coord_id == CMC_COORD_IDS::CMC_LON || coord_id == CMC_COORD_IDS::CMC_LAT || coord_id == CMC_COORD_IDS::CMC_LEV);
    
    /* Get the coordinate array */
    var_array_t& coordinate = t8_data->geo_data->coords->operator[](coord_id);

    /* Check the data type */
    switch (coordinate.get_data_type())
    {
        case CMC_INT32_T:
        {
            /* Get the integer indices corresponding to the supplied threshold */
            std::pair<int, int> threshold_indices = cmc_t8_geo_data_get_area_threshold(coordinate, std::get<int32_t>(starting_value), std::get<int32_t>(end_value));

            /* Save the indices in the compression settings */
            t8_data->settings.exclude_area_start_indices[coord_id] = threshold_indices.first;
            t8_data->settings.exclude_area_end_indices[coord_id] = threshold_indices.second;
        }
        break;
        case CMC_FLOAT:
        {
            /* Get the integer indices corresponding to the supplied threshold */
            std::pair<int, int> threshold_indices = cmc_t8_geo_data_get_area_threshold(coordinate, std::get<float>(starting_value), std::get<float>(end_value));

            /* Save the indices in the compression settings */
            t8_data->settings.exclude_area_start_indices[coord_id] = threshold_indices.first;
            t8_data->settings.exclude_area_end_indices[coord_id] = threshold_indices.second;
        }
        break;
        case CMC_DOUBLE:
        {
            /* Get the integer indices corresponding to the supplied threshold */
            std::pair<int, int> threshold_indices = cmc_t8_geo_data_get_area_threshold(coordinate, std::get<double>(starting_value), std::get<double>(end_value));

            /* Save the indices in the compression settings */
            t8_data->settings.exclude_area_start_indices[coord_id] = threshold_indices.first;
            t8_data->settings.exclude_area_end_indices[coord_id] = threshold_indices.second;
        }
        break;
        case CMC_INT16_T:
        {
            /* Get the integer indices corresponding to the supplied threshold */
            std::pair<int, int> threshold_indices = cmc_t8_geo_data_get_area_threshold(coordinate, std::get<int16_t>(starting_value), std::get<int16_t>(end_value));

            /* Save the indices in the compression settings */
            t8_data->settings.exclude_area_start_indices[coord_id] = threshold_indices.first;
            t8_data->settings.exclude_area_end_indices[coord_id] = threshold_indices.second;
        }
        break;
        case CMC_INT64_T:
        {
            /* Get the integer indices corresponding to the supplied threshold */
            std::pair<int, int> threshold_indices = cmc_t8_geo_data_get_area_threshold(coordinate, std::get<int64_t>(starting_value), std::get<int64_t>(end_value));

            /* Save the indices in the compression settings */
            t8_data->settings.exclude_area_start_indices[coord_id] = threshold_indices.first;
            t8_data->settings.exclude_area_end_indices[coord_id] = threshold_indices.second;
        }
        break;
        case CMC_UINT64_T:
        {
            /* Get the integer indices corresponding to the supplied threshold */
            std::pair<int, int> threshold_indices = cmc_t8_geo_data_get_area_threshold(coordinate, std::get<uint64_t>(starting_value), std::get<uint64_t>(end_value));

            /* Save the indices in the compression settings */
            t8_data->settings.exclude_area_start_indices[coord_id] = threshold_indices.first;
            t8_data->settings.exclude_area_end_indices[coord_id] = threshold_indices.second;
        }
        break;
        case CMC_UINT32_T:
        {
            /* Get the integer indices corresponding to the supplied threshold */
            std::pair<int, int> threshold_indices = cmc_t8_geo_data_get_area_threshold(coordinate, std::get<uint32_t>(starting_value), std::get<uint32_t>(end_value));

            /* Save the indices in the compression settings */
            t8_data->settings.exclude_area_start_indices[coord_id] = threshold_indices.first;
            t8_data->settings.exclude_area_end_indices[coord_id] = threshold_indices.second;
        }
        break;
        case CMC_INT8_T:
        {
            /* Get the integer indices corresponding to the supplied threshold */
            std::pair<int, int> threshold_indices = cmc_t8_geo_data_get_area_threshold(coordinate, std::get<int8_t>(starting_value), std::get<int8_t>(end_value));

            /* Save the indices in the compression settings */
            t8_data->settings.exclude_area_start_indices[coord_id] = threshold_indices.first;
            t8_data->settings.exclude_area_end_indices[coord_id] = threshold_indices.second;
        }
        break;
        case CMC_UINT8_T:
        {
            /* Get the integer indices corresponding to the supplied threshold */
            std::pair<int, int> threshold_indices = cmc_t8_geo_data_get_area_threshold(coordinate, std::get<uint8_t>(starting_value), std::get<uint8_t>(end_value));

            /* Save the indices in the compression settings */
            t8_data->settings.exclude_area_start_indices[coord_id] = threshold_indices.first;
            t8_data->settings.exclude_area_end_indices[coord_id] = threshold_indices.second;
        }
        break;
        case CMC_UINT16_T:
        {
            /* Get the integer indices corresponding to the supplied threshold */
            std::pair<int, int> threshold_indices = cmc_t8_geo_data_get_area_threshold(coordinate, std::get<uint16_t>(starting_value), std::get<uint16_t>(end_value));

            /* Save the indices in the compression settings */
            t8_data->settings.exclude_area_start_indices[coord_id] = threshold_indices.first;
            t8_data->settings.exclude_area_end_indices[coord_id] = threshold_indices.second;
        }
        break;
        case CMC_BYTE:
            cmc_err_msg("Coordinate values of type byte are not supported.");
        break;
        case CMC_CHAR:
        {
            /* Get the integer indices corresponding to the supplied threshold */
            std::pair<int, int> threshold_indices = cmc_t8_geo_data_get_area_threshold(coordinate, std::get<char>(starting_value), std::get<char>(end_value));

            /* Save the indices in the compression settings */
            t8_data->settings.exclude_area_start_indices[coord_id] = threshold_indices.first;
            t8_data->settings.exclude_area_end_indices[coord_id] = threshold_indices.second;
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }
    
    /* Set a flag that the exclude criterion will be applied */
    if (t8_data->settings.compression_criterium == CMC_T8_COMPRESSION_CRITERIUM::CMC_CRITERIUM_UNDEFINED)
    {
        /* If no criterion has been specified, set the excldue area */
        t8_data->settings.compression_criterium = CMC_T8_COMPRESSION_CRITERIUM::CMC_EXCLUDE_AREA;
    } else if (t8_data->settings.compression_criterium != CMC_T8_COMPRESSION_CRITERIUM::CMC_EXCLUDE_AREA)
    {
        /* If another criterion already has been specified, set the combined flag */
        t8_data->settings.compression_criterium = CMC_T8_COMPRESSION_CRITERIUM::CMC_COMBINED_CRITERION;
    }

    std::cout << "Settings for start is now: " << t8_data->settings.exclude_area_start_indices[0] << ", " << t8_data->settings.exclude_area_start_indices[1] << ", " << t8_data->settings.exclude_area_start_indices[2] << ", " << t8_data->settings.exclude_area_start_indices[3] << ", " << std::endl;
    std::cout << "Settings for  end  is now: " << t8_data->settings.exclude_area_end_indices[0] << ", " << t8_data->settings.exclude_area_end_indices[1] << ", " << t8_data->settings.exclude_area_end_indices[2] << ", " << t8_data->settings.exclude_area_end_indices[3] << ", " << std::endl;
    #endif
}

/**
 * @brief This functions calculates a Morton-Index from a 2D or 3D cartesian coordinate.
 * 
 * @param lin_coords The cartesian coordinate whose Morton-Index will be calculated
 * @param dim The dimension of the coordinate point (i.e. 2D or 3D)
 * @return uint64_t The Morton-Index
 * @note This calculation is based on the bit-interleaving of the single coordinate dimensions
 */
static
uint64_t
cmc_get_morton_index(std::tuple<uint64_t, uint64_t, uint64_t> lin_coords, const int dim)
{
    cmc_assert(dim == 2 || dim == 3);

    if (dim == 2)
    {
        /* In a 2D case */
        uint64_t x{std::get<0>(lin_coords)};
        uint64_t y{std::get<1>(lin_coords)};

        cmc_assert(x < (static_cast<uint64_t>(1) << 32) && y < (static_cast<uint64_t>(1) << 32));

        /** Prepare all coordinates for interleaving by fragmenting the bits, like:
         * (...b4 b3 b2 b1 b0) -> (...0 b4 0 b3 0 b2 0 b1 0 b0) */
        x = (x | (x << 16)) & 0x0000FFFF0000FFFF;
        x = (x | (x << 8))  & 0x00FF00FF00FF00FF;
        x = (x | (x << 4))  & 0x0F0F0F0F0F0F0F0F;
        x = (x | (x << 2))  & 0x3333333333333333;
        x = (x | (x << 1))  & 0x5555555555555555;
        
        y = (y | (y << 16)) & 0x0000FFFF0000FFFF;
        y = (y | (y << 8))  & 0x00FF00FF00FF00FF;
        y = (y | (y << 4))  & 0x0F0F0F0F0F0F0F0F;
        y = (y | (y << 2))  & 0x3333333333333333;
        y = (y | (y << 1))  & 0x5555555555555555;

        /* Interleave both coordinates */
        return x | (y << 1);
    } else
    {
        /* In a 3D case */
        uint64_t x{std::get<0>(lin_coords)};
        uint64_t y{std::get<1>(lin_coords)};
        uint64_t z{std::get<2>(lin_coords)};
        cmc_assert(x < (1 << 21) && y < (1 << 21) && z < (1 << 21));

        /** Prepare all coordinates for interleaving by fragmenting the bits, like:
         * (...b4 b3 b2 b1 b0) -> (...0 0 b4 0 0 b3 0 0 b2 0 0 b1 0 0 b0) */
        x = (x | (x << 14)  | (x << 28)) & 0x0001FC000FE0007F;
        x = (x | (x << 10)) & 0x06007C3003E1801F;
        x = (x | ((x << 6)  & ~(0x00000C0000600000))) & 0x06181C30C0E18607;
        x = (x | (x << 2))  & 0x1248649243249219;
        x = (x | (x << 2))  & 0x1249249249249249;

        y = (y | (y << 14)  | (y << 28)) & 0x0001FC000FE0007F;
        y = (y | (y << 10)) & 0x06007C3003E1801F;
        y = (y | ((y << 6)  & ~(0x00000C0000600000))) & 0x06181C30C0E18607;
        y = (y | (y << 2))  & 0x1248649243249219;
        y = (y | (y << 2))  & 0x1249249249249249;

        z = (z | (z << 14)  | (z << 28)) & 0x0001FC000FE0007F;
        z = (z | (z << 10)) & 0x06007C3003E1801F;
        z = (z | ((z << 6)  & ~(0x00000C0000600000))) & 0x06181C30C0E18607;
        z = (z | (z << 2))  & 0x1248649243249219;
        z = (z | (z << 2))  & 0x1249249249249249;

        /* Interleave all three coordinates */
        return x | (y << 1) | (z << 2);
    }
}

/**
 * @brief This functions calculates which process need to receive local data in order to set up a partition compliant to t8code.
 *        Therefore this functions is used for example when 'outside' data (for example in a blockwise-partitined manner) is transferred
 *        to us and we need to redistribute it compliant to the Morton order and t8code's partitioning scheme
 * 
 * @param t8_data A pointer to the @struct cmc_t8_data holding the incorrectly ordered data and the initial forest (onto which we want to distribute the data)
 * @param send_list A reference to the send_list. This list will be filled with the communications which need to take place
 * @param offsets An array holding the partiotion offsets of the forest, we want to distribute the data onto. (The size of the array is @var comm_size)
 * @param comm_size The size of the communicator
 */
static
void
cmc_t8_data_fill_sender_list_for_initial_distribution(cmc_t8_data_t t8_data, std::map<int, cmc_mpi_t8_send_data>& send_list, uint64_t* offsets, const int comm_size)
{
    #ifdef CMC_WITH_T8CODE
    cmc_assert(t8_data->variables_are_defined_on_the_same_domain);

    int err, rank;

    /* Get the rank of the process */
    err = MPI_Comm_rank(t8_data->comm, &rank);
    cmc_mpi_check_err(err);

    /* Therefore, we are examining which variables are defined on the exact same domain with the same layout, since they can share the reference coordinates */
    std::map<DATA_LAYOUT, std::vector<size_t>> reference_coordinate_variable_groups;

    /* Create variable groups based on their data layout */
    for (size_t id{0}; id < t8_data->vars.size(); ++id)
    {
        cmc_assert(t8_data->vars[id]->var->num_dimensions == 2 || t8_data->vars[id]->var->num_dimensions == 3);
        
        /* Check if varibale(s) with the same data layout already have been grouped */
        auto iter = reference_coordinate_variable_groups.find(t8_data->vars[id]->var->data_layout);

        if (iter != reference_coordinate_variable_groups.end())
        {
            /* If the data layout already has been added to the groups */
            /* Save the id of the current variable */
            reference_coordinate_variable_groups[t8_data->vars[id]->var->data_layout].push_back(id);
        }
        else
        {
            /* If the data layout is not equal to any previous variables */
            /* Create a new group of variables and assign the current variable id */
            reference_coordinate_variable_groups[t8_data->vars[id]->var->data_layout] = std::vector{id};
            reference_coordinate_variable_groups[t8_data->vars[id]->var->data_layout].reserve(t8_data->vars.size());
        }
    }

    std::cout << "Die Map hat groeße: " << reference_coordinate_variable_groups.size() << std::endl;

    /* Iterate over all variable groups and create their reference coordinates */
    /* Iterate over all variable groups with the same data layout */
    for (auto group_iter{reference_coordinate_variable_groups.begin()}; group_iter != reference_coordinate_variable_groups.end(); ++group_iter)
    {
        std::cout << "Group hat num variables: " << (group_iter->second).size() << std::endl;

        /* Iterate over all coordinates in order to calculate their Morton index and the correct process ownership */
        /* Since the group is defined on the exact same domain, we choose the first variable of the group in order to access the reference coordinates */

        /* Based on the reference coordinate representation, the data will be evaluated */
        switch (t8_data->vars[(group_iter->second)[0]]->var->ref_coordinates->coordinate_representation)
        {
            case cmc_coordinate_type::CMC_CARTESIAN_COORDINATES:
            {
                /* Swap contents with a temporary vector (in order for automatic deallocation), because the cartesian coordinates will not be needed hereafter */
                std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> coordinates;
                std::swap(coordinates, t8_data->vars[(group_iter->second)[0]]->var->ref_coordinates->cartesian_coordinates);

                /* Get the reference of the morton indices vector */
                std::vector<uint64_t>& morton_indices = t8_data->vars[(group_iter->second)[0]]->var->ref_coordinates->morton_indices;
                morton_indices.reserve(coordinates.size());

                /* Dimensionality of the data */
                const int dim = t8_data->geo_data->dim;

                /* Iterate over all coordinates */
                for (auto coord_iter{coordinates.begin()}; coord_iter != coordinates.end(); ++coord_iter)
                {
                    /* Transform the cartesian coordinate into a morton index (since this accounts for the initial disitribution, the cartesian coordnates are in range of the initial refienement level) */
                    morton_indices.push_back(cmc_get_morton_index(*coord_iter, dim));
                }

                size_t index_counter{0};
                for (auto morton_iter{morton_indices.begin()}; morton_iter != morton_indices.end(); ++morton_iter, ++index_counter)
                {
                    //TODO: make this a binary search with lower_bound
                    int recv_rank = comm_size -1;
                    for (int ir{0}; ir < comm_size -1; ++ir)
                    {
                        if (*morton_iter >= offsets[ir] && *morton_iter < offsets[ir+1])
                        {
                            recv_rank = ir;
                            break;
                        }
                    }
                    //auto process_rank_iter = std::lower_bound(offsets, offsets + comm_size, *morton_iter);
                    //int recv_rank = (process_rank_iter  != offsets + comm_size ? std::distance(offsets, process_rank_iter) : comm_size -1);

                    /* Check for the receiving rank in the 'send_list' */
                    auto receiver = send_list.find(recv_rank);

                    /* Check if the rank is alraedy a receiving rank in the 'send_list' */
                    if (receiver != send_list.end())
                    {
                        /* If the receiving rank is already listed within the 'send_list', we add the data to it */
                        receiver->second.morton_indices.push_back(*morton_iter);

                        size_t var_id_within_group{0};

                        /* Itearte over all varibales within the group */
                        for (auto var_iter{group_iter->second.begin()}; var_iter != group_iter->second.end(); ++var_iter, ++var_id_within_group)
                        {
                            /* Push back the variable's data to the dynamic array */
                            receiver->second.data[var_id_within_group]->push_back(t8_data->vars[*var_iter]->var->data->operator[](index_counter));  
                        }
                    } else
                    {
                        /* If the rank is not already present as receiver */
                        /* Create a new 'send_list-entry' */
                        send_list[recv_rank] = cmc_mpi_t8_data();

                        /* Reserve some memory for the morton indices */
                        send_list[recv_rank].morton_indices.reserve(morton_indices.size() / comm_size + 1);
                        
                        /* Save the first morton indiex for this receiving rank */
                        send_list[recv_rank].morton_indices.push_back(*morton_iter);

                        /* Allocate memory for all variables within the reference group */
                        send_list[recv_rank].data.reserve(group_iter->second.size());

                        /* Push the data for each variable to the dynamic_array */
                        /* Create a dynamic array for all variables */
                        size_t var_id_within_group{0};
                        for (auto var_iter{group_iter->second.begin()}; var_iter != group_iter->second.end(); ++var_iter, ++var_id_within_group)
                        {
                            /* Create a new dynamic array */
                            std::cout << "wie groß ist: " << morton_indices.size() / comm_size + 1 << " und wie groß ist der type: " << cmc_type_to_bytes[t8_data->vars[*var_iter]->var->data->get_data_type()] << std::endl;
                            send_list[recv_rank].data.push_back(new var_dynamic_array_t{morton_indices.size() / comm_size + 1, t8_data->vars[*var_iter]->var->data->get_data_type()});
                            /* Push back the first variable's data to the dynamic array */
                            send_list[recv_rank].data[var_id_within_group]->push_back(t8_data->vars[*var_iter]->var->data->operator[](index_counter));
                            /* Save the (global) variable_id of the variable concerning the t8_data struct */
                            send_list[recv_rank].variable_indices.push_back(*var_iter);/////TODO DIEse abhaengigkeit aufloesen und den vector komplett loswerden
                            std::cout << "Hallo2" << std::endl;
                        }
                    }
                }

                /* Check if some elements retain on this process */
                auto self_receiving = send_list.find(rank);
                if (self_receiving != send_list.end())
                {
                    /* We save the amount of Morton indices which will retain on this process manually, since this will not be communicated otherwise */
                    send_list[rank].received_morton_indices = send_list[rank].morton_indices.size();
                }
            }
            break;
            case cmc_coordinate_type::CMC_MORTON_INDEX:
                cmc_err_msg("Not yet implemented");
            break;
            default:
                cmc_err_msg("An undefined coordinate type is set for the reference coordinate. Therefore, a parallel distribution over the domain cannot be defined.");
        }
    }
    #endif
}


static
std::tuple<uint64_t, uint64_t, uint64_t>
cmc_t8_calc_cart_coord_from_lin_idx_lon_lat(const std::vector<size_t>& dim_lengths, size_t linear_index)
{
    /* Get the current position of the latitude */
    const size_t lat_pos_in_current_lon = linear_index % (dim_lengths[CMC_COORD_IDS::CMC_LAT]);
    /* Increase the linear index to the next full hyperslab */
    linear_index += (dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - lat_pos_in_current_lon);
    /* Get the longitude coordinate */
    const uint64_t lon = static_cast<uint64_t>((linear_index / dim_lengths[CMC_COORD_IDS::CMC_LAT]));
    /* Return a tuple with the cartesian coordinates */
    return std::make_tuple(lon, static_cast<uint64_t>(dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - lat_pos_in_current_lon), 0);
}

static
std::tuple<uint64_t, uint64_t, uint64_t>
cmc_t8_calc_cart_coord_from_lin_idx_lat_lon(const std::vector<size_t>& dim_lengths, size_t linear_index)
{
    /* Get the current position of the longitude */
    const size_t lon_pos_in_current_lat = linear_index % (dim_lengths[CMC_COORD_IDS::CMC_LON]);
    /* Increase the linear index to the next full hyperslab */
    linear_index += (dim_lengths[CMC_COORD_IDS::CMC_LON] - 1 - lon_pos_in_current_lat);
    /* Get the latitude coordinate */
    const size_t lat = linear_index / dim_lengths[CMC_COORD_IDS::CMC_LON];
    /* Return a tuple with the cartesian coordinates */
    //cmc_debug_msg("Cartesian coordinate von linear index: ", linear_index, " ist lon: ", lon_pos_in_current_lat, ", lat: ",dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - lat);
    return std::make_tuple(static_cast<uint64_t>(lon_pos_in_current_lat), static_cast<uint64_t>(dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - lat), 0);
}

static
std::tuple<uint64_t, uint64_t, uint64_t>
cmc_t8_calc_cart_coord_from_lin_idx_lat_lev(const std::vector<size_t>& dim_lengths, size_t linear_index)
{
    /* Get the current position of the elevation */
    const size_t lev_pos_in_current_lat = linear_index % (dim_lengths[CMC_COORD_IDS::CMC_LEV]);
    /* Increase the linear index to the next full hyperslab */
    linear_index += (dim_lengths[CMC_COORD_IDS::CMC_LEV] - 1 - lev_pos_in_current_lat);
    /* Get the latitude coordinate */
    const size_t lat = linear_index / dim_lengths[CMC_COORD_IDS::CMC_LEV];
    /* Return a tuple with the cartesian coordinates */
    return std::make_tuple(static_cast<uint64_t>(dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - lat), lev_pos_in_current_lat , 0);
}

static
std::tuple<uint64_t, uint64_t, uint64_t>
cmc_t8_calc_cart_coord_from_lin_idx_lev_lat(const std::vector<size_t>& dim_lengths, size_t linear_index)
{
    /* Get the current position of the latitude */
    const size_t lat_pos_in_current_lev = linear_index % (dim_lengths[CMC_COORD_IDS::CMC_LAT]);
    /* Increase the linear index to the next full hyperslab */
    linear_index += (dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - lat_pos_in_current_lev);
    /* Get the elevation coordinate */
    const uint64_t lev = static_cast<uint64_t>(linear_index / dim_lengths[CMC_COORD_IDS::CMC_LAT]);
    /* Return a tuple with the cartesian coordinates */
    return std::make_tuple(static_cast<uint64_t>(dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - lat_pos_in_current_lev), lev , 0);
}

static
std::tuple<uint64_t, uint64_t, uint64_t>
cmc_t8_calc_cart_coord_from_lin_idx_lon_lev(const std::vector<size_t>& dim_lengths, size_t linear_index)
{
    /* Get the current position of the elevation */
    const size_t lev_pos_in_current_lon = linear_index % (dim_lengths[CMC_COORD_IDS::CMC_LEV]);
    /* Increase the linear index to the next full hyperslab */
    linear_index += (dim_lengths[CMC_COORD_IDS::CMC_LEV] - 1 - lev_pos_in_current_lon);
    /* Get the longitude coordinate */
    const uint64_t lon = static_cast<uint64_t>(linear_index / dim_lengths[CMC_COORD_IDS::CMC_LEV]);
    /* Return a tuple with the cartesian coordinates */
    return std::make_tuple(lon, static_cast<uint64_t>(lev_pos_in_current_lon), 0);
}

static
std::tuple<uint64_t, uint64_t, uint64_t>
cmc_t8_calc_cart_coord_from_lin_idx_lev_lon(const std::vector<size_t>& dim_lengths, size_t linear_index)
{
    /* Get the current position of the longitude */
    const size_t lon_pos_in_current_lev = linear_index % (dim_lengths[CMC_COORD_IDS::CMC_LON]);
    /* Increase the linear index to the next full hyperslab */
    linear_index += (dim_lengths[CMC_COORD_IDS::CMC_LON] - 1 - lon_pos_in_current_lev);
    /* Get the elevation coordinate */
    const uint64_t lev = static_cast<uint64_t>(linear_index / dim_lengths[CMC_COORD_IDS::CMC_LON]);
    /* Return a tuple with the cartesian coordinates */
    return std::make_tuple(static_cast<uint64_t>(lon_pos_in_current_lev), lev, 0);
}

static
std::tuple<uint64_t, uint64_t, uint64_t>
cmc_t8_calc_cart_coord_from_lin_idx_lon_lat_lev(const std::vector<size_t>& dim_lengths, size_t linear_index)
{
    /* Get the longitude position */
    const size_t lon_pos = linear_index / (dim_lengths[CMC_COORD_IDS::CMC_LAT] * dim_lengths[CMC_COORD_IDS::CMC_LEV]);
    /* Transform the linear index of the current slice (shift into the range of a single 'incomplete' slice) */
    linear_index -= (lon_pos * dim_lengths[CMC_COORD_IDS::CMC_LAT] * dim_lengths[CMC_COORD_IDS::CMC_LEV]);
    /* Get the elevation position */
    const size_t lev_pos = linear_index % (dim_lengths[CMC_COORD_IDS::CMC_LEV]);
    /* Increase the linear index to the next full slice */
    linear_index += (dim_lengths[CMC_COORD_IDS::CMC_LEV] - 1 - lev_pos);
    /* Get the latitude coordinate */
    const size_t lat_pos = linear_index / dim_lengths[CMC_COORD_IDS::CMC_LEV];
    /* Return a tuple with the cartesian coordinates */
    return std::make_tuple(static_cast<uint64_t>(lon_pos), static_cast<uint64_t>(dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - lat_pos), static_cast<uint64_t>(lev_pos));
}

static
std::tuple<uint64_t, uint64_t, uint64_t>
cmc_t8_calc_cart_coord_from_lin_idx_lev_lon_lat(const std::vector<size_t>& dim_lengths, size_t linear_index)
{
    /* Get the elevation position */
    const size_t lev_pos = linear_index / (dim_lengths[CMC_COORD_IDS::CMC_LON] * dim_lengths[CMC_COORD_IDS::CMC_LAT]);
    /* Transform the linear index of the current slice (shift into the range of a single 'incomplete' slice) */
    linear_index -= (lev_pos * dim_lengths[CMC_COORD_IDS::CMC_LON] * dim_lengths[CMC_COORD_IDS::CMC_LAT]);
    /* Get the latitude position */
    const size_t lat_pos = linear_index % (dim_lengths[CMC_COORD_IDS::CMC_LAT]);
    /* Increase the linear index to the next full slice */
    linear_index += (dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - lat_pos);
    /* Get the longitude coordinate */
    const uint64_t lon_pos = static_cast<uint64_t>(linear_index / dim_lengths[CMC_COORD_IDS::CMC_LAT]);
    /* Return a tuple with the cartesian coordinates */
    return std::make_tuple(lon_pos, static_cast<uint64_t>(dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - lat_pos), static_cast<uint64_t>(lev_pos));
}

static
std::tuple<uint64_t, uint64_t, uint64_t>
cmc_t8_calc_cart_coord_from_lin_idx_lon_lev_lat(const std::vector<size_t>& dim_lengths, size_t linear_index)
{
    /* Get the longitude position */
    const size_t lon_pos = linear_index / (dim_lengths[CMC_COORD_IDS::CMC_LAT] * dim_lengths[CMC_COORD_IDS::CMC_LEV]);
    /* Transform the linear index of the current slice (shift into the range of a single 'incomplete' slice) */
    linear_index -= (lon_pos * dim_lengths[CMC_COORD_IDS::CMC_LAT] * dim_lengths[CMC_COORD_IDS::CMC_LEV]);
    /* Get the latitude position */
    const size_t lat_pos = linear_index % (dim_lengths[CMC_COORD_IDS::CMC_LAT]);
    /* Increase the linear index to the next full slice */
    linear_index += (dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - lat_pos);
    /* Get the elevation coordinate */
    const uint64_t lev_pos = static_cast<uint64_t>(linear_index / dim_lengths[CMC_COORD_IDS::CMC_LAT]);
    /* Return a tuple with the cartesian coordinates */
    return std::make_tuple(static_cast<uint64_t>(lon_pos), static_cast<uint64_t>(dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - lat_pos), lev_pos);
}

static
std::tuple<uint64_t, uint64_t, uint64_t>
cmc_t8_calc_cart_coord_from_lin_idx_lev_lat_lon(const std::vector<size_t>& dim_lengths, size_t linear_index)
{
    /* Get the elevation position */
    const size_t lev_pos = linear_index / (dim_lengths[CMC_COORD_IDS::CMC_LON] * dim_lengths[CMC_COORD_IDS::CMC_LAT]);
    /* Transform the linear index of the current slice (shift into the range of a single 'incomplete' slice) */
    linear_index -= (lev_pos * dim_lengths[CMC_COORD_IDS::CMC_LON] * dim_lengths[CMC_COORD_IDS::CMC_LAT]);
    /* Get the longitude position */
    const size_t lon_pos = linear_index % (dim_lengths[CMC_COORD_IDS::CMC_LON]);
    /* Increase the linear index to the next full slice */
    linear_index += (dim_lengths[CMC_COORD_IDS::CMC_LON] - 1 - lon_pos);
    /* Get the latitude coordinate */
    const size_t lat_pos = linear_index / dim_lengths[CMC_COORD_IDS::CMC_LON];
    /* Return a tuple with the cartesian coordinates */
    return std::make_tuple(static_cast<uint64_t>(lon_pos), static_cast<uint64_t>(dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - lat_pos), static_cast<uint64_t>(lev_pos));
}

static
std::tuple<uint64_t, uint64_t, uint64_t>
cmc_t8_calc_cart_coord_from_lin_idx_lat_lev_lon(const std::vector<size_t>& dim_lengths, size_t linear_index)
{
    /* Get the latitude position */
    const size_t lat_pos = linear_index / (dim_lengths[CMC_COORD_IDS::CMC_LON] * dim_lengths[CMC_COORD_IDS::CMC_LEV]);
    /* Transform the linear index of the current slice (shift into the range of a single 'incomplete' slice) */
    linear_index -= (lat_pos * dim_lengths[CMC_COORD_IDS::CMC_LON] * dim_lengths[CMC_COORD_IDS::CMC_LEV]);
    /* Get the longitude position */
    const size_t lon_pos = linear_index % (dim_lengths[CMC_COORD_IDS::CMC_LON]);
    /* Increase the linear index to the next full slice */
    linear_index += (dim_lengths[CMC_COORD_IDS::CMC_LON] - 1 - lon_pos);
    /* Get the elevation coordinate */
    const uint64_t lev_pos = static_cast<uint64_t>(linear_index / dim_lengths[CMC_COORD_IDS::CMC_LON]);
    /* Return a tuple with the cartesian coordinates */
    return std::make_tuple(static_cast<uint64_t>(lon_pos), static_cast<uint64_t>(dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - lat_pos), lev_pos);
}

static
std::tuple<uint64_t, uint64_t, uint64_t>
cmc_t8_calc_cart_coord_from_lin_idx_lat_lon_lev(const std::vector<size_t>& dim_lengths, size_t linear_index)
{
    /* Get the latitude position */
    const size_t lat_pos = linear_index / (dim_lengths[CMC_COORD_IDS::CMC_LON] * dim_lengths[CMC_COORD_IDS::CMC_LEV]);
    /* Transform the linear index of the current slice (shift into the range of a single 'incomplete' slice) */
    linear_index -= (lat_pos * dim_lengths[CMC_COORD_IDS::CMC_LON] * dim_lengths[CMC_COORD_IDS::CMC_LEV]);
    /* Get the elevation position */
    const size_t lev_pos = linear_index % (dim_lengths[CMC_COORD_IDS::CMC_LEV]);
    /* Increase the linear index to the next full slice */
    linear_index += (dim_lengths[CMC_COORD_IDS::CMC_LEV] - 1 - lev_pos);
    /* Get the longitude coordinate */
    const uint64_t lon_pos = static_cast<uint64_t>(linear_index / dim_lengths[CMC_COORD_IDS::CMC_LEV]);
    /* Return a tuple with the cartesian coordinates */
    return std::make_tuple(lon_pos, static_cast<uint64_t>(dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - lat_pos), static_cast<uint64_t>(lev_pos));
}

/* A generic error function if the deduction of the correct 'calc_cart_coord'-function has failed */
static
std::tuple<uint64_t, uint64_t, uint64_t>
cmc_t8_calc_cart_coord_from_lin_idx_err(const std::vector<size_t>& dim_lengths, size_t linear_index)
{
    cmc_err_msg("The deduction of the correct function for calculating a cartesian coordinate from an linear index has failed.");
    return std::make_tuple(0UL, 0UL, 0UL);
}

/* Array holding all functions for calculating a cartesian coordinate from a linear index for each data layout */
const std::array<std::function<std::tuple<uint64_t, uint64_t, uint64_t>(const std::vector<size_t>&, size_t)>, 13> _offset_functions_lin_idx_to_cart
{
/* 2D Offset funcions */
cmc_t8_calc_cart_coord_from_lin_idx_lon_lat,
cmc_t8_calc_cart_coord_from_lin_idx_lat_lon,
cmc_t8_calc_cart_coord_from_lin_idx_lat_lev,
cmc_t8_calc_cart_coord_from_lin_idx_lev_lat,
cmc_t8_calc_cart_coord_from_lin_idx_lon_lev,
cmc_t8_calc_cart_coord_from_lin_idx_lev_lon,
/* 3D offset functions */
cmc_t8_calc_cart_coord_from_lin_idx_lon_lat_lev,
cmc_t8_calc_cart_coord_from_lin_idx_lon_lev_lat,
cmc_t8_calc_cart_coord_from_lin_idx_lev_lon_lat,
cmc_t8_calc_cart_coord_from_lin_idx_lev_lat_lon,
cmc_t8_calc_cart_coord_from_lin_idx_lat_lev_lon,
cmc_t8_calc_cart_coord_from_lin_idx_lat_lon_lev,
/* Error function */
cmc_t8_calc_cart_coord_from_lin_idx_err};


/* Return a function based on the data layout of the variable corresponding to the given id for calculating a cartesian coordinate from a linear index and the data_layout */
static
std::function<std::tuple<uint64_t, uint64_t, uint64_t>(const std::vector<size_t>&, size_t)>
cmc_t8_get_offset_function_lin_idx_to_cart_coord_based_on_data_layout(const DATA_LAYOUT data_layout)
{
    #ifdef CMC_WITH_T8CODE
    switch (data_layout)
    {
        case DATA_LAYOUT::CMC_2D_LON_LAT:
            return _offset_functions_lin_idx_to_cart[0];
        break;
        case DATA_LAYOUT::CMC_2D_LAT_LON:
            return _offset_functions_lin_idx_to_cart[1];
        break;
        case DATA_LAYOUT::CMC_2D_LAT_LEV:
            return _offset_functions_lin_idx_to_cart[2];
        break;
        case DATA_LAYOUT::CMC_2D_LEV_LAT:
            return _offset_functions_lin_idx_to_cart[3];
        break;
        case DATA_LAYOUT::CMC_2D_LON_LEV:
            return _offset_functions_lin_idx_to_cart[4];
        break;
        case DATA_LAYOUT::CMC_2D_LEV_LON:
            return _offset_functions_lin_idx_to_cart[5];
        break;
        case DATA_LAYOUT::CMC_3D_LON_LAT_LEV:
            return _offset_functions_lin_idx_to_cart[6];
        break;
        case DATA_LAYOUT::CMC_3D_LON_LEV_LAT:
            return _offset_functions_lin_idx_to_cart[7];
        break;
        case DATA_LAYOUT::CMC_3D_LEV_LON_LAT:
            return _offset_functions_lin_idx_to_cart[8];
        break;
        case DATA_LAYOUT::CMC_3D_LEV_LAT_LON:
            return _offset_functions_lin_idx_to_cart[9];
        break;
        case DATA_LAYOUT::CMC_3D_LAT_LEV_LON:
            return _offset_functions_lin_idx_to_cart[10];
        break;
        case DATA_LAYOUT::CMC_3D_LAT_LON_LEV:
            return _offset_functions_lin_idx_to_cart[11];
        break;
        default :
            cmc_err_msg("The variable contains an undefined data layout.");
            return _offset_functions_lin_idx_to_cart[12];
    }
    #endif
}

/**
 * @brief Adds a global offset to a local cartesian coordinate.
 * 
 * @param cart_coord The local cartesian coordinate
 * @param global_dim_lengths The global dimension lengths
 * @param start_offset The global offset this process has
 * @param local_dim_lengths The process-local dimension lengths
 * @param data_layout The data layout 
 */
static
void
add_start_offset_to_cart_coord(std::tuple<uint64_t, uint64_t, uint64_t>& cart_coord, const std::vector<size_t>& global_dim_lengths, const std::vector<size_t>& start_offset, const std::vector<size_t>& local_dim_lengths, const DATA_LAYOUT data_layout)
{
    cmc_assert(start_offset.size() >= 2);
    /* Obtain references to the elements of the tuple */
    auto& [x, y, z] = cart_coord;
    switch (data_layout)
    {
        case DATA_LAYOUT::CMC_2D_LON_LAT:
            x += start_offset[0];
            y = global_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - ((local_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - y) + start_offset[1]);
        break;
        case DATA_LAYOUT::CMC_2D_LAT_LON:
            x += start_offset[1];
            y = global_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - ((local_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - y) + start_offset[0]);
        break;
        case DATA_LAYOUT::CMC_2D_LAT_LEV:
            x = global_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - ((local_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - x) + start_offset[0]);
            y += start_offset[1];
        break;
        case DATA_LAYOUT::CMC_2D_LEV_LAT:
            x = global_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - ((local_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - x) + start_offset[1]);
            y += start_offset[0];
        break;
        case DATA_LAYOUT::CMC_2D_LON_LEV:
            x += start_offset[0];
            y += start_offset[1];
        break;
        case DATA_LAYOUT::CMC_2D_LEV_LON:
            x += start_offset[1];
            y += start_offset[0];
        break;
        case DATA_LAYOUT::CMC_3D_LON_LAT_LEV:
            x += start_offset[0];
            y = global_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - ((local_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - y) + start_offset[1]);
            z += start_offset[2];
        break;
        case DATA_LAYOUT::CMC_3D_LON_LEV_LAT:
            x += start_offset[0];
            y = global_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - ((local_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - y) + start_offset[2]);
            z += start_offset[1];
        break;
        case DATA_LAYOUT::CMC_3D_LEV_LON_LAT:
            x += start_offset[1];
            y = global_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - ((local_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - y) + start_offset[2]);
            z += start_offset[0];
        break;
        case DATA_LAYOUT::CMC_3D_LEV_LAT_LON:
            x += start_offset[2];
            y = global_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - ((local_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - y) + start_offset[1]);
            z += start_offset[0];
        break;
        case DATA_LAYOUT::CMC_3D_LAT_LEV_LON:
            x += start_offset[2];
            y = global_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - ((local_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - y) + start_offset[0]);
            z += start_offset[1];
        break;
        case DATA_LAYOUT::CMC_3D_LAT_LON_LEV:
            x += start_offset[1];
            y = global_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - ((local_dim_lengths[CMC_COORD_IDS::CMC_LAT] - 1 - y) + start_offset[0]);
            z += start_offset[2];
        break;
        default :
            cmc_err_msg("The variable contains an undefined data layout.");
    }
}

/**
 * @brief This functions creates linear integer coordinates for the variables regarding a global coordinate system
 * 
 * @param global_system The global coordinate system
 * @param var The local variable to consider
 * @param coord_type The current type of the coordinates
 * @return cmc_ref_coordinates_t The calculated integer reference coordinates for the variable @var var 
 */
static
cmc_ref_coordinates_t
cmc_create_ref_coordinates(cmc_global_coordinate_system_t global_system, cmc_var_t var, const cmc_coordinate_type coord_type)
{
    cmc_assert(global_system != nullptr);
    /* Create new reference coordinates */
    cmc_ref_coordinates_t ref_coords = new cmc_ref_coordinates();

    /* Craete avector with global dimension length*/
    std::vector<size_t> global_dim_lengths;
    global_dim_lengths.reserve(CMC_NUM_COORD_IDS);

    /* Assign the global dimension lengths */
    for (int coord_dim{0}; coord_dim < CMC_NUM_COORD_IDS; ++coord_dim)
    {
        global_dim_lengths.push_back(global_system->coords.operator[](static_cast<CMC_COORD_IDS>(coord_dim)).size());
        std::cout << "Global dim_length fuer coord_dim: " << global_dim_lengths.back() << std::endl;
    }
    /* Create the reference coordinates based on the coordinate type the data currently has */
    switch(coord_type)
    {
        case cmc_coordinate_type::CMC_CARTESIAN_COORDINATES:
        {
            cmc_err_msg("Not yet implemented");
        }
        break;
        case cmc_coordinate_type::CMC_MORTON_INDEX:
        {
            cmc_err_msg("Not yet implemented");
        }
        break;
        case cmc_coordinate_type::CMC_GEO_DOMAIN_DEFINED_BY_BOX:
        {
            // In this case, we assume the data to be sorted corectly compliant to a default ascending ordering scheme in each coordinate diemnsion 
            cmc_assert(var->start_ptr.size() > 0 && var->count_ptr.size() > 0);
            //TODO: Maybe check for sorting order of the values?

            /* Obtain an offset function for converting the boxed domain layout to cartesian reference coordinates */
            std::function<std::tuple<uint64_t, uint64_t, uint64_t>(const std::vector<size_t>&, size_t)> lin_idx_to_cart_coord = cmc_t8_get_offset_function_lin_idx_to_cart_coord_based_on_data_layout(var->data_layout);
            
            /* We 'linearize' over all dimensions (using a single loop to cover all dimensions) */
            const size_t total_length = std::reduce(var->count_ptr.begin(), var->count_ptr.end(), 1, std::multiplies<>());

            /* Allocate memory for all reference coordinates */
            ref_coords->cartesian_coordinates.reserve(total_length);

            /* Iterate linearily over all coordinates of the variable's data */
            for (size_t linear_id{0}; linear_id < total_length; ++linear_id)
            {
                /* Get a (local) cartesian coordinate description of the linear index */
                std::tuple<uint64_t, uint64_t, uint64_t> cart_coord = lin_idx_to_cart_coord(var->dim_lengths, linear_id);
                /* Apply the start vector offset for the (local) domain reference coordinates (in order to obtain global reference coordinates) */
                add_start_offset_to_cart_coord(cart_coord, global_dim_lengths, var->start_ptr, var->dim_lengths, var->data_layout);
                /* Save the cartesian reference coordinate */
                ref_coords->cartesian_coordinates.push_back(cart_coord);
            }

            /* Set the flag which coordinate type is represent in the reference coordinates */
            ref_coords->coordinate_representation = cmc_coordinate_type::CMC_CARTESIAN_COORDINATES;
        }
        break;
        default:
            cmc_err_msg("An unknwon coordinate type was supplied from which no refernece cooridnates are retrieveable.");
    }

    return ref_coords;
}

/**
 * @brief This functions calculates integer reference coordinates for all variables within the @var t8_data
 * 
 * @param t8_data The @struct cmc_t8_data holding the process-local data of the variables for which the reference coordinates will be calculated
 */
static
void
cmc_t8_data_create_cartesian_reference_coords_for_all_vars(cmc_t8_data_t t8_data)
{
    /* Construct the reference coordinates for the variables */
    /* Therefore, we are examining which variables are defined on the exact same domain with the same layout, since they can share the reference coordinates */
    std::map<DATA_LAYOUT, std::vector<size_t>> reference_coordinate_variable_groups;

    /* Create variable groups based on their data layout */
    for (size_t id{0}; id < t8_data->vars.size(); ++id)
    {
        cmc_assert(t8_data->vars[id]->var->num_dimensions == 2 || t8_data->vars[id]->var->num_dimensions == 3);
        
        /* Check if varibale(s) with the same data layout already have been grouped */
        auto iter = reference_coordinate_variable_groups.find(t8_data->vars[id]->var->data_layout);

        if (iter != reference_coordinate_variable_groups.end())
        {
            /* If the data layout already has been added to the groups */
            /* Save the id of the current variable */
            reference_coordinate_variable_groups[t8_data->vars[id]->var->data_layout].push_back(id);
        }
        else
        {
            /* If the data layout is not equal to any previous variables */
            /* Create a new group of variables and assign the current variable id */
            reference_coordinate_variable_groups[t8_data->vars[id]->var->data_layout] = std::vector{id};
            reference_coordinate_variable_groups[t8_data->vars[id]->var->data_layout].reserve(t8_data->vars.size());
        }
    }

    /* Iterate over all variable groups and create their reference coordinates */
    /* The references are saved generally in t8_geo_data in order to have a central ownership */
    
    /* Reserve memory for all reference coordinates */
    t8_data->geo_data->ref_coordinates.reserve(reference_coordinate_variable_groups.size());

    /* Iterate over all variable groups with the same data layout */
    for (auto iter{reference_coordinate_variable_groups.begin()}; iter != reference_coordinate_variable_groups.end(); ++iter)
    {
        /* Create the reference coordinates for the first variable of the group */
        t8_data->geo_data->ref_coordinates.push_back(cmc_create_ref_coordinates(t8_data->geo_data->coordinates, t8_data->vars[(iter->second)[0]]->var, cmc_coordinate_type::CMC_GEO_DOMAIN_DEFINED_BY_BOX));        

        /* Assign those reference coordinates to all variables from the group */
        for (size_t var_id{0}; var_id < (iter->second).size(); ++var_id)
        {
            /* Assign the newly created reference coordinates */
            t8_data->vars[(iter->second)[var_id]]->var->ref_coordinates = t8_data->geo_data->ref_coordinates.back();
        }
    }
}

#ifdef CMC_ENABLE_MPI
inline constexpr std::array<MPI_Datatype, cmc_type::CMC_NUM_TYPES> cmc_type_to_mpi_type{MPI_BYTE, MPI_INT8_T, MPI_CHAR, MPI_INT16_T, MPI_INT32_T, MPI_FLOAT, MPI_DOUBLE, MPI_UINT8_T, MPI_UINT16_T, MPI_UINT32_T, MPI_INT64_T, MPI_UINT64_T};
#endif


/**
 * @brief This funciton performs the initial communication which has been calculated previously with the fucntion cmc_t8_data_fill_sender_list_for_initial_distribution(...)
 * 
 * @param t8_data A pointer to the @struct cmc_t8_data holding the (process-local) data which will be distributed
 * @param send_list The list containg the communications which have to be performed
 * @param offsets An array holding the partition offsets of the forest compliant to which the data will be distributed (the array has the size @var comm_size)
 * @param comm_size The size of the communicator
 * @return std::map<int, cmc_mpi_t8_recv_data> A map holding the data which was received from the other processes
 */
static
std::map<int, cmc_mpi_t8_recv_data>
cmc_t8_data_mpi_initial_communication(cmc_t8_data_t t8_data, std::map<int, cmc_mpi_t8_send_data>& send_list, uint64_t* offsets, const int comm_size)
{
    #ifdef CMC_ENABLE_MPI

    cmc_debug_msg("Initial parallel data distribution begins.");

    int rank, err;

    /* Get the process-local rank of the communicator */
    err = MPI_Comm_rank(t8_data->comm, &rank);
    cmc_mpi_check_err(err);

    /* TODO: Big loop over these send and receive calls for different data distributions */
    
    /* Create a map saving the received data */
    std::map<int, cmc_mpi_t8_recv_data> recv_list;

    /* Define an arbitrarly tag for the Morton indices */
    const int tag_morton_indices = static_cast<int>(t8_data->vars.size());

    /* Create a vector of MPI requests */
    std::vector<MPI_Request> send_requests;
    send_requests.reserve(send_list.size());

    /* Counter for the amount of data which will reside on this process and does not have to be sent to other ranks */
    uint64_t count_of_objs_retained = 0;

    /* A flag indicating whether some objects are already correctly distributed or if all elements have been sent away to different ranks */
    bool flag_some_data_elements_are_kept_local = false;

    /* Iterate over the 'send_list' and send the messages to all processes */
    size_t send_req_idx = 0;

    /* Iterate over the send_list */
    for (auto iter = send_list.begin(); iter != send_list.end(); ++iter, ++send_req_idx)
    {
        /* Do not let the process send data to itself */
        if (iter->first != rank)
        {
            /* Distribute the Morton indices */
            err = MPI_Isend(iter->second.morton_indices.data(), iter->second.morton_indices.size(), MPI_UINT64_T, iter->first, tag_morton_indices, t8_data->comm, &send_requests[send_req_idx]);
            cmc_mpi_check_err(err);

            cmc_debug_msg("Sending ", iter->second.morton_indices.size(), " data elements to rank ", iter->first, " with tag ", tag_morton_indices);
        } else
        {
            /* Set the flag that some elements remain on this process */
            flag_some_data_elements_are_kept_local = true;
            /* Save the amount of elements which are kept local */
            count_of_objs_retained = static_cast<uint64_t>(iter->second.morton_indices.size());
        }
    }
    
    cmc_debug_msg("The amount of data elements which will reside on this processs are ", count_of_objs_retained, " data elements per variable.");
    
    /* Wait until all messages are staged */
    err = MPI_Barrier(t8_data->comm);
    cmc_mpi_check_err(err);

    /* A vector collecting thr IDs of ranks from which data is received */
    std::vector<int> partner_ranks;
    partner_ranks.reserve(send_list.size()); 

    /* Counter variable fot the while loop */
    size_t curr_recv_id = 0;

    /* Variable which will hold the amount of elements within the message */
    int recveied_objs = 0;

    /* MPI Status objects for the send commands */
    MPI_Status recv_stat, actual_recv_stat;

    /* Flags for working through the received messages */
    bool flag_messages_present = true;
    int mpi_msg_flag = 0;
    
    /* Receive the messages concerning the Morton indices */
    /* Since all messages have been staged before, the MPI_Iprobe call finds a matching message in each iteration */
    while (flag_messages_present)
    {
        /* Check for a message containing the morton indices from another rank */
        err = MPI_Iprobe(MPI_ANY_SOURCE, tag_morton_indices, t8_data->comm, &mpi_msg_flag, &recv_stat);
        cmc_mpi_check_err(err);

        /* Check if a message from any source with any tag is available */
        if (mpi_msg_flag)
        {
            /* Get the size of the message, respectively the amount of Morton indices which will be received */
            /* This number equals the amount of data which will be receied from the partner rank for each variable */
            err = MPI_Get_count(&recv_stat, MPI_UINT64_T, &recveied_objs);
            cmc_mpi_check_err(err);

            /* Check if the rank which has sent the data is already in the 'recv_list' */
            /* Check for the sending rank in the 'recv_list' */
            auto sender = recv_list.find(recv_stat.MPI_SOURCE);

            /* Check if the rank is not yet a sending rank in the 'recv_list' */
            if (sender == recv_list.end())
            {
                /* Create a new cmc_mpi_t8_recv_data in the recv_list */
                recv_list[recv_stat.MPI_SOURCE] = cmc_mpi_t8_data();
                /* Reserve memory for the incoming morton indices */
                recv_list[recv_stat.MPI_SOURCE].morton_indices.reserve(recveied_objs);
                /* Reserve memory for the vector containing the variable data and fill it with nullptrs */
                recv_list[recv_stat.MPI_SOURCE].data = std::vector<var_dynamic_array_t*>(t8_data->vars.size(), nullptr);
            }

            /* If a message is available, receive the Morton indices */
            err = MPI_Recv(recv_list[recv_stat.MPI_SOURCE].morton_indices.data(), recveied_objs, MPI_UINT64_T, recv_stat.MPI_SOURCE, recv_stat.MPI_TAG, t8_data->comm, &actual_recv_stat);
            cmc_mpi_check_err(err);
            
            /* Store the received morton_indices */
            recv_list[recv_stat.MPI_SOURCE].received_morton_indices = recveied_objs;

            cmc_debug_msg("Received ", recveied_objs, " data elements from rank ", actual_recv_stat.MPI_SOURCE, " with tag ", actual_recv_stat.MPI_TAG);

            /* Save the IDs of the ranks which have sent data */
            partner_ranks.push_back(actual_recv_stat.MPI_SOURCE);

            ++curr_recv_id;
        } else
        {
            /* If none messages are left */
            flag_messages_present = false;
        }
    }

    /** Send the variable's data to the other ranks **/

    /* A vector of requests for the send-operations */
    std::vector<MPI_Request> sending_requests(send_list.size() * t8_data->vars.size());

    /* An additional counter fot the loop */
    size_t rank_counter = 0;

    /* Iterate over the send_list */
    for (auto iter = send_list.begin(); iter != send_list.end(); ++iter, ++rank_counter)
    {
        /* Do not let the process send data to itself */
        if (iter->first != rank)
        {
            /* Iterate over all variables */
            for (size_t var_iter{0}; var_iter < iter->second.data.size(); ++var_iter)
            {
                cmc_debug_msg("Sending ", iter->second.data[var_iter]->size(), " data elements of variable ", t8_data->vars[iter->second.variable_indices[var_iter]]->var->name, " to rank ", iter->first);
                
                /* Send the data of the current variable */
                err = MPI_Isend(iter->second.data[var_iter]->get_initial_data_ptr(), iter->second.data[var_iter]->size(), MPI_DOUBLE, iter->first, iter->second.variable_indices[var_iter], t8_data->comm, &sending_requests[rank_counter * iter->second.data.size() + var_iter]);
                cmc_mpi_check_err(err);
            }
        }
    }

    /* Create a vector for the receiving statuses */
    MPI_Status receiving_stat;

    /* Receive all messages concerning the variable's data (one message for each partner rank and each variable) */
    for (size_t iter{0}; iter < partner_ranks.size() * t8_data->vars.size(); ++iter)
    {
        /* Check (and wait) for an incoming message */
        err = MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, t8_data->comm, &receiving_stat);
        cmc_mpi_check_err(err);

        /* Get the data type of the variable from from the t8_data->vars[tag]->var->get_type() or so */
        const cmc_type variable_type = t8_data->vars[receiving_stat.MPI_TAG]->get_type();

        /* Define a variable containing the size of the message */
        int receiving_objs = 0;

        /* Get the size of the receiving message */
        err = MPI_Get_count(&receiving_stat, cmc_type_to_mpi_type[variable_type], &receiving_objs);
        cmc_mpi_check_err(err);

        /* Allocate an array for receiving the data */
        recv_list[receiving_stat.MPI_SOURCE].data[receiving_stat.MPI_TAG] = new var_dynamic_array_t{static_cast<size_t>(receiving_objs), t8_data->vars[receiving_stat.MPI_TAG]->get_type()};

        /* Receive the message */
        err = MPI_Recv(recv_list[receiving_stat.MPI_SOURCE].data[receiving_stat.MPI_TAG]->get_initial_data_ptr(), receiving_objs, cmc_type_to_mpi_type[variable_type], receiving_stat.MPI_SOURCE, receiving_stat.MPI_TAG, t8_data->comm, MPI_STATUS_IGNORE);
        cmc_mpi_check_err(err);

        cmc_debug_msg("Received ", receiving_objs, " data elements of variable ", t8_data->vars[receiving_stat.MPI_TAG]->var->name, " from rank ", receiving_stat.MPI_SOURCE);
    }

    /* If some elements did not have to be communicated because they are already correctly distributed and local to this rank, we move them over from the send_list to the recv_list */
    if(flag_some_data_elements_are_kept_local)
    {
        /* Just move the local elements from the send_list to the recv_list */
        recv_list[rank] = std::move(send_list[rank]);
    }
    cmc_debug_msg("Halllooo");
    
    return recv_list;
    #else
    return std::map<int, cmc_mpi_t8_recv_data>();
    #endif
}

/**
 * @brief This function redistributes parallel data and orders is compliant to t8code's ordering scheme (i.e. Morton-Order)
 * 
 * @param t8_data A pointer to the @struct cmc_t8_data holding the distributed variable_s data which will be redistributed
 */
void
cmc_t8_geo_data_distribute_and_apply_ordering(cmc_t8_data_t t8_data)
{
    #ifdef CMC_WITH_T8CODE
    #ifdef CMC_ENABLE_MPI

    cmc_debug_msg("Redistributing data for t8code.");
    /* Only redistribute the data initially if the flag is set */
    if(t8_data->use_distributed_data)
    {
        //Currently, it is only possible to distribute data if all variables are defined on the same domain 
        //cmc_assert(t8_data->variables_are_defined_on_the_same_domain);
        const t8_forest_t forest = t8_data->initial_forest;
        int err, rank, size;

        /* Get the size of the communicator */
        err = MPI_Comm_size(t8_data->comm, &size);
        cmc_mpi_check_err(err);
        
        /* Get the (local) rank id */
        err = MPI_Comm_rank(t8_data->comm, &rank);
        cmc_mpi_check_err(err);

        /* Schemes for quadrilaterals and hexahedrons */
        t8_default_scheme_quad_c scheme_quad;
        t8_default_scheme_hex_c scheme_hex;
        t8_eclass_scheme_c* scheme_eclass;

        if (t8_data->geo_data->dim == 2)
        {
            scheme_eclass = static_cast<t8_eclass_scheme_c*>(&scheme_quad);
        } else
        {
            cmc_assert(t8_data->geo_data->dim == 3);
            scheme_eclass = static_cast<t8_eclass_scheme_c*>(&scheme_hex);
        }

        /* Get the offest of the first local element as a SFC index (at the initial refienment level) */
        uint64_t elem_offset = static_cast<uint64_t>(t8_element_get_linear_id (scheme_eclass,
                                                                               t8_forest_get_element_in_tree(forest, 0, 0),
                                                                               t8_data->geo_data->initial_refinement_lvl));

        //TODO: The offets do not have to be gathered, they might be accessible directly via t8code's internal data structure
        
        /* Allocate an offset array holding the SFC index (at the intitial refinement level) of each starting process-local element from each process */
        uint64_t* offsets = (uint64_t*) malloc(sizeof(uint64_t) * size);

        /* Distribute the offsets between all processes */
        err = MPI_Allgather(&elem_offset, 1, MPI_UINT64_T, offsets, 1, MPI_UINT64_T, t8_data->comm);
        cmc_mpi_check_err(err);

        /* Eventually communicate here the global coordinate domain (if it is not present on all processes) */

        /* Create reference coordinates in order to calculate the distribution later */
        cmc_t8_data_create_cartesian_reference_coords_for_all_vars(t8_data);

        /* A map collecting the receiver rank id as a key and the associated data which will be sent to it later */
        std::map<int, cmc_mpi_t8_send_data> send_list;

        /* Calculate the data which will be sent */
        cmc_t8_data_fill_sender_list_for_initial_distribution(t8_data, send_list, offsets, size);

        /* A map collecting the correctly distributed data after the communication */
        std::map<int, cmc_mpi_t8_send_data> recv_list = cmc_t8_data_mpi_initial_communication(t8_data, send_list, offsets, size);

        /** Now that the data is redistributed and residing in the recv_list, we have to sort/combine the data from the different ranks compliant to the z-curve order **/

        /** Since all indices are already in order morton order and the distribution in t8code is linearily, we can just fill up the missing morton indices (which are not present; in the range of the offsets of this rank) with the missing values */
        /** Die Luecken die nicht besetzt sind muessen gefuellt werden, das Problem: wir wissen nciht wie viele Elemente dort aufgefuellt werden muessen */
        /** Deswegen, muss ab dem ersten Element der Luecke durch den forest itertiert werden, bis wir auf eine Element treffen, das auf dem initial refinement ist und dessen element_id der des naechsten Elements direkt nach der Luecke entspricht */
        /** Diese Elemente muessen gezaehlt werden und anschließend dazwischenkopiert werden */
        cmc_debug_msg("Initial communication is over");
        /* Define a vector holding the sorted Morton indices from each rank */
        std::vector<std::vector<uint64_t>> sorted_morton_indices;
        sorted_morton_indices.reserve(recv_list.size());

        /* Count all elements */
        size_t num_elements = 0;

        /* Iterate over the recv_list */
        for (auto iter{recv_list.begin()}; iter != recv_list.end(); ++iter)
        {
            /* Save the amount of Morton indices */
            num_elements += iter->second.received_morton_indices;

            /* Create a new vector */
            sorted_morton_indices.emplace_back(std::vector<uint64_t>());

            /* Copy the unsorted Morton indices */
            std::copy(iter->second.morton_indices.data(), iter->second.morton_indices.data() + iter->second.received_morton_indices , std::back_inserter(sorted_morton_indices.back()));

            /* Sort the indices */
            std::sort(sorted_morton_indices.back().begin(), sorted_morton_indices.back().end());

            cmc_debug_msg("Size of sorted morton indices: ", sorted_morton_indices.back().size());
        }

        cmc_debug_msg("Size of sorted_morton_indices is: ", sorted_morton_indices.size());

        /* Define a vector which will hold the merged sequence of all Morton indices */
        std::vector<uint64_t> morton_indices;

        if(sorted_morton_indices.size() > 1)
        {
            /* Merge all sequences of sorted Morton indices */
            for (size_t vec_iter{0}; vec_iter < sorted_morton_indices.size() - 1; ++vec_iter)
            {
                /* Allocate a new vector capable of holding the merged sequence */
                std::vector<uint64_t> merged_sequence;

                /* Allocate memory for the vector */
                if (vec_iter < sorted_morton_indices.size() -2)
                {
                    /* Allocate just enough memeory for the merged sequence */
                    merged_sequence.reserve(sorted_morton_indices[vec_iter].size() + sorted_morton_indices[vec_iter + 1].size());
                } else
                {
                    /* In the last iteration, we are going to allocate as much elements as the forest holds locally, since we need to fill up some 'holes' within the data */
                    merged_sequence.reserve(t8_forest_get_local_num_elements(forest));
                }

                /* Merge the two sorted sequences of Morton indices */
                std::merge(sorted_morton_indices[vec_iter].begin(), sorted_morton_indices[vec_iter].end(), sorted_morton_indices[vec_iter +1].begin(), sorted_morton_indices[vec_iter +1].end(), std::back_inserter(merged_sequence));

                /* Swap the merged sequence */
                if (vec_iter < sorted_morton_indices.size() -2)
                {
                    /* Swap the merged sequence with the latter sorted_morton_indices-array */
                    std::swap(sorted_morton_indices[vec_iter + 1], merged_sequence);
                } else
                {
                    /* In the last iteration, swap the current merged-sequence with the outside vector, in order to access the fully sorted Morton indices array later */
                    std::swap(morton_indices, merged_sequence);
                }
            }
        } else if (sorted_morton_indices.size() == 1) 
        {
            /* If we received data only from one process, we do not have to merge */
            morton_indices = std::move(sorted_morton_indices.front());
        }

        /* Since the forest on which the data is mapped may eventually be greater than the geographical domain on which the variable's data is defined. There may be some holes within the Morton indices, which we need to fill missing values */
        std::vector<std::pair<size_t, size_t>> gaps_to_fill;

        /* We need to iterate later through the Morton indices and obtain the corresponding data and store it within the correct variables */
        /* Since we may have more than one variable we are going to create a mapping to use for each variables which will be filled later */
        std::vector<std::vector<uint64_t>> sequence_mappings;
        sequence_mappings.reserve(recv_list.size());

        /* Check if there is some real data on this process */
        if (morton_indices.size() > 0)
        {
            /* Check for gaps in the data at the beginning  */
            if (morton_indices[0] != offsets[rank])
            {
                /* In this case our first local Morton index is not equal to the first local element. Therefore, we have to insert some missing_values */

                int counter_skipped_elems = 0;
                t8_element_t* element;
                uint64_t morton_offset_calc = 0;

                while (morton_offset_calc < morton_indices[0] -1 && static_cast<size_t>(counter_skipped_elems + 1) < static_cast<size_t>(t8_forest_get_local_num_elements(forest) -1))
                {
                    ++counter_skipped_elems;
                    element = t8_forest_get_element_in_tree(forest, 0, counter_skipped_elems);
                    morton_offset_calc += std::pow(t8_element_num_children(scheme_eclass, element), t8_data->geo_data->initial_refinement_lvl - t8_element_level(scheme_eclass, element));
                }

                /* Save the position and length of the gap */
                gaps_to_fill.emplace_back(std::make_pair(0, counter_skipped_elems));

                /* Insert 'counter_skipped_elems' at this position, between the Morton indices we have checked. This has to be done in order to skip the prefilled gapsof missing values */
                morton_indices.insert(morton_indices.begin(), counter_skipped_elems, (offsets[rank] - 1 >= 0 ? offsets[rank] -1 : 0));
            }
            cmc_debug_msg("All gaps in morton indices have been found at start");

            /* We cannot start at morton_indices.begin() if we have added some missing_values before the first (data) Morton index */
            int offset_for_iteration = 0;
            if (gaps_to_fill.size() > 0)
            {
                /* Get the amount of elements to skip */
                offset_for_iteration = gaps_to_fill.back().second;
            }

            /* Check for gaps in the data somewhere in between */
            size_t current_length = morton_indices.size() - 1;

            for (size_t iter{static_cast<size_t>(offset_for_iteration)}; iter < current_length; ++iter)
            {
                /* Check if neighboring Morton_indices do not vary by exactly one */
                if (morton_indices[iter] + 1 < morton_indices[iter + 1])
                {
                    /** If not, we have a hole needed to be filled by missing values **/
                    int counter_skipped_elems = 0;
                    t8_element_t* element;
                    uint64_t morton_offset_calc = morton_indices[iter];

                    while (morton_offset_calc < morton_indices[iter + 1] -1 && iter + counter_skipped_elems + 1 < static_cast<size_t>(t8_forest_get_local_num_elements(forest) -1))
                    {
                        ++counter_skipped_elems;
                        element = t8_forest_get_element_in_tree(forest, 0, iter + counter_skipped_elems);
                        morton_offset_calc += std::pow(t8_element_num_children(scheme_eclass, element), t8_data->geo_data->initial_refinement_lvl - t8_element_level(scheme_eclass, element));
                    }

                    /* Save the position and length of the gap */
                    gaps_to_fill.emplace_back(std::make_pair(iter + 1, counter_skipped_elems));

                    /* Insert 'counter_skipped_elems' at this position, between the Morton indices we have checked. This has to be done in order to skip the prefilled gapsof missing values */
                    morton_indices.insert(std::next(morton_indices.begin(), iter + 1), counter_skipped_elems, morton_indices[iter]);

                    /* Advance the iterator to the correct position */
                    iter += counter_skipped_elems;
                    current_length += counter_skipped_elems;
                }
            }

            /* Check if there are missing values at the end of our data */
            /* If there are some Morton indices missing at the end, we can just add the missing values equal to the amount of (local forest elements - amount of local Morton indices) */
            if (morton_indices.size() < static_cast<size_t>(t8_forest_get_local_num_elements(forest)))
            {
                /* If this is case, the last Morton index should not be equal to the partition boundary */
                cmc_assert(((rank == size -1 ? t8_forest_get_global_num_elements(forest) - 1 : offsets[rank + 1] -1) - morton_indices.back()) > 0);

                /* Save the position and length of the gap */
                gaps_to_fill.emplace_back(std::make_pair(morton_indices.size(), t8_forest_get_local_num_elements(forest) - static_cast<int>(morton_indices.size())));

                /* Insert the missing values at the end */
                morton_indices.insert(morton_indices.end(), static_cast<size_t>(t8_forest_get_local_num_elements(forest)) - morton_indices.size(), morton_indices.back());//(rank == size -1 ? t8_forest_get_local_num_elements(forest) - 1 : offsets[rank + 1] - 1));
            }

            /* Iterate over the recv_list */
            for (auto iter{recv_list.begin()}; iter != recv_list.end(); ++iter)
            {
                /* Create a new vector for the current receiver and allocate memory for it */
                sequence_mappings.emplace_back(std::vector<uint64_t>());
                sequence_mappings.back().reserve(iter->second.received_morton_indices);

                /* Check the received data */
                for (int morton_iter{0}; morton_iter < iter->second.received_morton_indices; ++morton_iter)
                {
                    /* Perform a binary search for of the received Morton index */
                    auto bin_search_result = std::lower_bound(morton_indices.begin(), morton_indices.end(), iter->second.morton_indices[morton_iter]);
                    /* The element should be found in the fully sorted array of Morton indices. Therefore, the assertion */
                    cmc_assert(bin_search_result != morton_indices.end());
                    /* Now we store the position of the Morton-index (in the fully sorted array) within our mapping */
                    sequence_mappings.back().emplace_back(std::distance(morton_indices.begin(), bin_search_result));
                }
            }

        } else
        {
            /* Since no morton indices were sent to this process, we will fill everything with missing values */
            gaps_to_fill.emplace_back(std::make_pair(0, t8_forest_get_local_num_elements(forest)));
        }
        
        /* Allocate memory for all variables */
        /* This approach leads to the case that each rank holds data from each variable. Therefore, all variables are distributed */
        for (size_t var_iter{0}; var_iter < t8_data->vars.size(); ++var_iter)
        {
            /* Create a new var_array for the variable */
            t8_data->vars[var_iter]->var->data_new = new var_array_t(static_cast<size_t>(t8_forest_get_local_num_elements(forest)), t8_data->vars[var_iter]->get_type());  
            cmc_debug_msg("missing value is: ", std::get<double>(t8_data->vars[var_iter]->var->missing_value));
            /* Fill in the gaps with missing values for each variable */
            for (auto gap_iter{gaps_to_fill.begin()}; gap_iter != gaps_to_fill.end(); ++gap_iter)
            {
                /* Copy the amount of missing values to each gap */
                for (size_t num_elems_to_insert{0}; num_elems_to_insert < gap_iter->second; ++num_elems_to_insert)
                {
                    /* Copy-assign the variable's missing value */
                    t8_data->vars[var_iter]->var->data_new->assign(static_cast<size_t>(gap_iter->first + num_elems_to_insert), t8_data->vars[var_iter]->var->missing_value);
                }
            }
        }
        
        /* Define a counter for accessing at the current recv_list's position */
        int recv_list_counter = 0;

        /** Copy the data for each of the received variables from each rank **/

        /* Iterate over the recv_list */
        for (auto iter{recv_list.begin()}; iter != recv_list.end(); ++iter, ++recv_list_counter)
        {
            /* Iterate over the received variables */
            for (size_t var_iter{0}; var_iter < iter->second.data.size(); ++var_iter)
            {
                size_t size_of_data = t8_data->vars[var_iter]->get_data_size();

                std::byte* initial_data_ptr = static_cast<std::byte*>(iter->second.data[var_iter]->get_initial_data_ptr());
                std::byte* initial_data_new_ptr = static_cast<std::byte*>(t8_data->vars[var_iter]->get_initial_data_new_ptr());

                /* Fill in the actual variable's data */
                for (int morton_id{0}; morton_id < iter->second.received_morton_indices; ++morton_id)
                {
                    /* Copy the data value from the 'old' array (with the wrong ordering) to the new Morton-Curve-compliant ordering array */
                    memcpy(static_cast<void*>(initial_data_new_ptr + size_of_data * ((sequence_mappings[recv_list_counter])[morton_id])), static_cast<void*>(initial_data_ptr + size_of_data * morton_id), size_of_data);
                }
            }
        }

        /* Exchange the 'old' unordered data with the newly ordered array (compliant to the Morton curve) */
        for (size_t var_iter{0}; var_iter < t8_data->vars.size(); ++var_iter)
        {
            cmc_debug_msg("        bis hier-...", var_iter);
            /* Free the 'unordered' data and assign the 'newly ordered' data for each variable */
            t8_data->vars[var_iter]->var->switch_data();
            /* Set the Morton Curve ordering flag */
            t8_data->vars[var_iter]->var->data_scheme = CMC_DATA_ORDERING_SCHEME::CMC_GEO_DATA_Z_CURVE;
        }
    
        /* Free the allocated data */
        free(offsets);
        
        cmc_debug_msg("The data has been initially redistributed compliant to t8code.");

        err = MPI_Barrier(t8_data->comm);
        cmc_mpi_check_err(err);
    }

    #endif
    #endif
}

#if 0
static
cmc_var_vector_t
cmc_get_default_sorted_coords(cmc_var_vector_t initial_coords)
{
    /* Default sorting means ascending longitude. latitude and elevation coordinates */
    /* E.g. longitude = 0°, 2°, 4°, ..., 356°, 358° or longitude = -180°, -178°, ..., 0, 2*, ..., 178° */
    /*      latitude  = -90°, -88°, ..., 0°, 2°, ..., 90° */
    /*      elevation = 0km, 1km, 2km, ... */
    
    /* First, we need to check wich coordinate dimension are not sorted compliant to the default scheme */

    /* Sort the coordinates */

    /* Return the sorted coordinates */

}
#endif
