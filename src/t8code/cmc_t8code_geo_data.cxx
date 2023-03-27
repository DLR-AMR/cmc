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
    return x_coord + (dim_lengths[CMC_COORD_IDS::CMC_LAT] -1 -y_coord) * dim_lengths[CMC_COORD_IDS::CMC_LON];
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
    return y_coord + x_coord * dim_lengths[CMC_COORD_IDS::CMC_LEV];
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
    return x_coord + y_coord * dim_lengths[CMC_COORD_IDS::CMC_LAT];
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
            [[fallthrough]];
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_ERROR_THRESHOLD:
                adapt_data.initial_ref_lvl_ids.emplace_back(std::unordered_map<t8_locidx_t, t8_locidx_t>());
                #if 0
                adapt_data._counter.push_back(0);
                adapt_data._counter_nxt_lvl.push_back(0);
                #endif
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
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_ERROR_THRESHOLD:
            #if 0
                /* Allocate space for all unordered maps (for each varibale) */
                adapt_data.initial_ref_lvl_ids.reserve(t8_data->vars.size());
                adapt_data._counter.reserve(t8_data->vars.size());
                adapt_data._counter_nxt_lvl.reserve(t8_data->vars.size());
                /* Create a var_vector holding data which may be used by the interpolation */
                adapt_data.adapted_data = new var_vector_t();
                adapt_data.adapted_data->reserve(1);
            #endif
                /* Reset the adapt_data class */
                adapt_data.adapt_step = 0;
                /* Save the current variable ID */
                adapt_data.current_var_id = var_id;
                interpolation_data.current_var_id = var_id;
                /* Save the initial data */
                adapt_data.initial_ref_lvl_ids.emplace_back(std::unordered_map<t8_locidx_t, t8_locidx_t>());
                //t8_data->initial_data.push_back(t8_data->vars[var_id]->var->data);
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
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_ERROR_THRESHOLD:
                /* Reset the counters */
                #if 0
                adapt_data._counter[0] = 0;
                adapt_data._counter_nxt_lvl[0] = 0;
                #endif
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
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_ERROR_THRESHOLD:
                /* Reset the adapt counters */
                #if 0
                adapt_data._counter[var_id] = 0;
                adapt_data._counter_nxt_lvl[var_id] = 0;
                #endif
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
    //adapt_data._counter.reserve(adapt_data.t8_data->vars.size());
    //adapt_data._counter_nxt_lvl.reserve(adapt_data.t8_data->vars.size());
    /* Create a var_vector holding data which may be used by the interpolation */
    adapt_data.adapted_data = new var_vector_t();
    adapt_data.adapted_data->reserve(1);
    #endif
}

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
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_ERROR_THRESHOLD:
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
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_ERROR_THRESHOLD:
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
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_ERROR_THRESHOLD:
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
            case CMC_T8_COMPRESSION_CRITERIUM::CMC_ERROR_THRESHOLD:
                /* Delete the allocation of the var_vector */
                delete adapt_data.adapted_data;
            break;
            default:
                cmc_err_msg("The supplied lossy compression criterium is not yet implemented.");
        }
    }
    #endif
}

/** We are passing 't8_data' as a member of the amr_adapt_data to the forest, in order to have the maximum access to the current data of all variables
 *  This even allows us to check for example if certain thresholds are satisfied on/for all levels/layers/variables 
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

        /* Update the data striucts at the end of a iteration */
        cmc_t8_update_adapt_and_interpolation_data_end_of_iteration(adapt_data, interpolation_data);
    }

    /* Delete allocations or perform any finalizing steps */
    cmc_t8_deconstruct_adapt_and_interpolate_data(adapt_data, interpolation_data);

    /* Free the former forest and Save the adapted forest */
    t8_data->assets->forest = forest;
    
    cmc_debug_msg("Adaptation/Compression is finished.");

    #endif
}

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

            /* Switch to coarsened forest */
            forest = coarsened_forest;

            /* Decrement the refinement level */
            --ref_lvl;

            /* Update the adapt and inerpolation struct at the end of the iteration */
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

    /* Vector storing the ids of the variables which will be reordered */
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
cmc_t8_geo_data_set_error_criterium(cmc_t8_data_t t8_data, const double maximim_error_tolerance)
{
    #ifdef CMC_WITH_T8CODE
    /* Save the error tolerance */
    t8_data->settings.max_err = maximim_error_tolerance;
    t8_data->settings.compression_criterium = CMC_T8_COMPRESSION_CRITERIUM::CMC_CRITERIUM_UNDEFINED;
    #endif
}
