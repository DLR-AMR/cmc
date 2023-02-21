#include "cmc_t8code_geo_mesh.h"
#include "cmc_t8code_data.hxx"
#include "cmc_t8_adapt_callbacks.h"
#include "utilities/cmc_log_functions.h"
/////new try

/** Begin STATIC Functions **/
/** Calculate the maximum refinement level needed per direction/dimension in order to build an enclosing t8code mesh */
#ifdef CMC_WITH_T8CODE

static int
cmc_t8_calc_geo_refinement_level(const cmc_t8_data& t8_data, const int var_id)
{
    cmc_assert(var_id >= 0 && var_id < static_cast<int>(t8_data.vars.size()));
    const size_t max_elem_per_direction{*std::max_element(t8_data.vars[var_id]->var->dim_lengths.begin(), t8_data.vars[var_id]->var->dim_lengths.end())};
    /* Calculate the induced initial refinement level needed in order to build an enclosing mesh */
    return static_cast<int>(std::ceil(std::log2(max_elem_per_direction) + std::numeric_limits<double>::epsilon()));
}

static std::string
data_layout_to_string(const DATA_LAYOUT layout)
{
    switch (layout) {
        case DATA_LAYOUT::CMC_2D_LON_LAT:
            return std::string("CMC_2D_LON_LAT");
        break;
        case DATA_LAYOUT::CMC_2D_LAT_LON:
            return std::string("CMC_2D_LAT_LON");
        break;
        case DATA_LAYOUT::CMC_2D_LAT_LEV:
            return std::string("CMC_2D_LAT_LEV");
        break;
        case DATA_LAYOUT::CMC_2D_LEV_LAT:
            return std::string("CMC_2D_LEV_LAT");
        break;
        case DATA_LAYOUT::CMC_2D_LON_LEV:
            return std::string("CMC_2D_LON_LEV");
        break;
        case DATA_LAYOUT::CMC_2D_LEV_LON:
            return std::string("CMC_2D_LEV_LON");
        break;
        case DATA_LAYOUT::CMC_3D_LON_LAT_LEV:
            return std::string("CMC_3D_LON_LAT_LEV");
        break;
        case DATA_LAYOUT::CMC_3D_LON_LEV_LAT:
            return std::string("CMC_3D_LON_LEV_LAT");
        break;
        case DATA_LAYOUT::CMC_3D_LEV_LON_LAT:
            return std::string("CMC_3D_LEV_LON_LAT");
        break;
        case DATA_LAYOUT::CMC_3D_LEV_LAT_LON:
            return std::string("CMC_3D_LEV_LAT_LON");
        break;
        case DATA_LAYOUT::CMC_3D_LAT_LEV_LON:
            return std::string("CMC_3D_LAT_LEV_LON");
        break;
        case DATA_LAYOUT::CMC_3D_LAT_LON_LEV:
            return std::string("CMC_3D_LAT_LON_LEV");
        break;
        default :
            cmc_err_msg("The variable contains an undefined data layout.");
            return std::string("undefined data layout");
    }
}

/** Coarsen the initial uniform mesh which encloses the geo-spatial data */
static t8_forest_t
cmc_t8_coarsen_geo_mesh(cmc_t8_data& t8_data, t8_forest_t initial_forest, const int initial_refinement_lvl, const int var_id = -1)
{
    t8_forest_t forest{initial_forest};
    t8_forest_t forest_adapt;

    int refinement_step{0};
    t8_gloidx_t num_elems_former_forest{0};
    /* Create an 'adaption struct' */
    cmc_t8_adapt_data_t adapt_data{new cmc_t8_adapt_data(&t8_data)};
    /* Assign the current var id (the underlying geo-spatial data of the variable corresponding to the given id will be regarded during the coarsening) */
    /* If no id is explicitly given, it will be assumed that each variable is defined on the same geo-spatial domain. Therfore, the first variable id is (arbitrarily) assigned */
    adapt_data->current_var_id = (var_id <= -1 ? 0 : var_id);
    /* Set the partition for coarsening flag for the 'partition'-routine */
    const int partition_for_coarsening = 0;
    /* Adapt the forest as much as possible by coarsen the 'dummy' elements which do not resemble a geo-spatial data point */
    while (refinement_step < initial_refinement_lvl && num_elems_former_forest != t8_forest_get_global_num_elements(forest))
    {
        /* Initialize the new forest */
        t8_forest_init(&forest_adapt);
        /* Save the number of elements of the former forest */
        num_elems_former_forest = t8_forest_get_global_num_elements(forest);
        /* Set the user data (needed for the adaption step) */
        t8_forest_set_user_data(forest_adapt, static_cast<void*>(adapt_data));
        /* Adapt the forest accordingly to the callback function */
        t8_forest_set_adapt(forest_adapt, forest, cmc_t8_adapt_coarsen_geo_mesh_callback, 0);
        /* Partition the forest */
        t8_forest_set_partition(forest_adapt, forest, partition_for_coarsening);
        /* Commit the adapted forest (perform the adaption step) */
        t8_forest_commit(forest_adapt);
        /* Save the coarsened forest */
        forest = forest_adapt;

        /* Update iteration variable */
        ++refinement_step;
    }
    /* Deallocate the 'adapt data' */
    delete adapt_data;
    /* Print some information in Debug-Mode */
    cmc_debug_msg("The coarsened forest contains ", t8_forest_get_global_num_elements(forest_adapt)," elements.");
    cmc_debug_msg("Data dimensions are:"),
    cmc_debug_msg("#Latitude: ", t8_data.geo_data->get_coord_length(CMC_COORD_IDS::CMC_LAT), ", #Longitude: ", t8_data.geo_data->get_coord_length(CMC_COORD_IDS::CMC_LON), ", #Elevation: ", t8_data.geo_data->get_coord_length(CMC_COORD_IDS::CMC_LEV));

    /* Return the coarsened forest */
    return forest;
}
#endif
/** END STATIC Functions **/

bool
compare_geo_domain_equality_of_data_layouts(const DATA_LAYOUT first, const DATA_LAYOUT second)
{
    switch (first)
    {
        case DATA_LAYOUT::CMC_2D_LON_LAT:
            if (second == DATA_LAYOUT::CMC_2D_LON_LAT || second == DATA_LAYOUT::CMC_2D_LAT_LON)
            {
                return true;
            } else
            {
                return false;
            }
        break;
        case DATA_LAYOUT::CMC_2D_LAT_LON:
            if (second == DATA_LAYOUT::CMC_2D_LAT_LON || second == DATA_LAYOUT::CMC_2D_LON_LAT)
            {
                return true;
            } else
            {
                return false;
            }
        break;
        case DATA_LAYOUT::CMC_2D_LAT_LEV:
            if (second == DATA_LAYOUT::CMC_2D_LAT_LEV || second == DATA_LAYOUT::CMC_2D_LEV_LAT)
            {
                return true;
            } else
            {
                return false;
            }
        break;
        case DATA_LAYOUT::CMC_2D_LEV_LAT:
            if (second == DATA_LAYOUT::CMC_2D_LEV_LAT || second == DATA_LAYOUT::CMC_2D_LAT_LEV)
            {
                return true;
            } else
            {
                return false;
            }
        break;
        case DATA_LAYOUT::CMC_2D_LON_LEV:
            if (second == DATA_LAYOUT::CMC_2D_LON_LEV || second == DATA_LAYOUT::CMC_2D_LEV_LON)
            {
                return true;
            } else
            {
                return false;
            }
        break;
        case DATA_LAYOUT::CMC_2D_LEV_LON:
            if (second == DATA_LAYOUT::CMC_2D_LEV_LON || second == DATA_LAYOUT::CMC_2D_LON_LEV)
            {
                return true;
            } else
            {
                return false;
            }
        break;
        case DATA_LAYOUT::CMC_3D_LON_LAT_LEV:
            [[fallthrough]];
        case DATA_LAYOUT::CMC_3D_LON_LEV_LAT:
            [[fallthrough]];
        case DATA_LAYOUT::CMC_3D_LEV_LON_LAT:
            [[fallthrough]];
        case DATA_LAYOUT::CMC_3D_LEV_LAT_LON:
            [[fallthrough]];
        case DATA_LAYOUT::CMC_3D_LAT_LEV_LON:
            [[fallthrough]];
        case DATA_LAYOUT::CMC_3D_LAT_LON_LEV:
            /* If the second layout is a 3D layout too, the domain has to coincide with the first's domain */
            if (second > _INTERN_ID_END_2D_START_3D)
            {
                return true;
            } else
            {
                return false;
            }
        break;
        default :
            cmc_warn_msg("The data layouts are not comparable.");
            return false;
    }
}

/** Create a forest mesh which includes/encloses the given (geo-spatial) data */
void
cmc_t8_create_enclosing_geo_mesh(cmc_t8_data& t8_data)
{
    #ifdef CMC_WITH_T8CODE
    cmc_assert(t8_data.vars.size() > 0);
    cmc_assert(t8_data.vars[0]->var->num_dimensions == 2 || t8_data.vars[0]->var->num_dimensions == 3);
    t8_forest_t initial_forest{nullptr};
    /* Save the data layout of the first variable */
    DATA_LAYOUT layout{t8_data.vars[0]->var->data_layout};
    /* Check if all variables are defined in the same domain */
    for (size_t var_id{1}; var_id < t8_data.vars.size(); ++var_id)
    {
        cmc_assert(t8_data.vars[var_id]->var->num_dimensions == 2 || t8_data.vars[var_id]->var->num_dimensions == 3);
        /* If the variables have a different data layout, several meshes are needed in order to perform the compresseion */
        if (!compare_geo_domain_equality_of_data_layouts(layout, t8_data.vars[var_id]->var->data_layout))
        {
            /* Switch the flag indicating if all variables are defined on the same domain to false */
            t8_data.variables_are_defined_on_the_same_domain = false;
        }
    }

    /* Check if the compression mode is compliant with the variable's data domain */
    if (((t8_data.compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D) ||
        (t8_data.compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D)) &&
        !t8_data.variables_are_defined_on_the_same_domain)
    {
        cmc_err_msg("The choosen compression mode is of type 'One for All', but the geo-spatial data domain does not coincide with all supplied variables.\n",
                    "Please, either supply variables which are defined on the same coordinate dimensions or choose a different compression mode.");
        return;
    }

    if (t8_data.variables_are_defined_on_the_same_domain)
    {
        /* Create a new cmesh according to the dimension of the coordinate dimensions */
        t8_cmesh_t cmesh = t8_cmesh_new_hypercube((t8_data.geo_data->dim == 2 ? T8_ECLASS_QUAD : T8_ECLASS_HEX), t8_data.comm, 0, 0, 1);

        cmc_debug_msg("Built hypercube cmesh of dimension ", t8_data.geo_data->dim);
        
        /* Calculate the initial refinement level needed in order to circumvent the geo-grid with a t8code mesh */
        //t8_data.geo_data->initial_refinement_lvl = cmc_t8_calc_geo_refinement_level(t8_data.geo_data->coords);
        t8_data.geo_data->initial_refinement_lvl = cmc_t8_calc_geo_refinement_level(t8_data, 0);
        #if CMC_ENABLE_DEBUG
        /* Assertions in Debug-Mode */
        if (t8_data.geo_data->dim == 2)
        {
            cmc_assert(t8_data.geo_data->initial_refinement_lvl < P4EST_MAXLEVEL);
        } else if (t8_data.geo_data->dim == 3)
        {
            cmc_assert(t8_data.geo_data->initial_refinement_lvl < P8EST_MAXLEVEL);
        }
        #endif

        /* Create a new uniform forest with the former calculated initial refinement level */
        initial_forest = t8_forest_new_uniform(cmesh, t8_scheme_new_default_cxx(), t8_data.geo_data->initial_refinement_lvl, 0, t8_data.comm);

        /* Allocate assets (based on the compression mode) */
        if (t8_data.compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
            t8_data.compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D)
        {
            t8_data.assets = new cmc_t8_assets(t8_data.geo_data->initial_refinement_lvl);
            t8_data.assets_allocated = true;
            /* Save a reference for each variable to this asset */
            for (size_t var_id{0}; var_id < t8_data.vars.size(); ++var_id)
            {
                /* Save a pointer to the initial forest for each variable */
                t8_data.vars[var_id]->assets = t8_data.assets;
            }
        } else
        {
            /* Currently, only a compression mode of type 'One for One' is possible */
            for (size_t var_id{0}; var_id < t8_data.vars.size(); ++var_id)
            {
                /* Save a reference of the initial forest in each variable */ 
                t8_data.vars[var_id]->assets = new cmc_t8_assets(t8_data.geo_data->initial_refinement_lvl);
                t8_data.vars[var_id]->assets_allocated = true;
            }
        }

        cmc_debug_msg("Built forest (dim: ", t8_data.geo_data->dim, "D) with an intitial refinement level of ", t8_data.geo_data->initial_refinement_lvl);

        /* Coarsen the elements of the forest mesh which are not part of the simulation geo-grid */
        t8_data.initial_forest = cmc_t8_coarsen_geo_mesh(t8_data, initial_forest, t8_data.geo_data->initial_refinement_lvl);

        /* Save the initial forest and pointer to it accordingly for each variable (based on the compression mode) */
        if (t8_data.compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
            t8_data.compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D)
        {
            t8_forest_ref(t8_data.initial_forest);
            t8_data.assets->forest = t8_data.initial_forest;
        } else
        {
            /* Currently, only a compression mode of type 'One for One' is possible */
            for (size_t var_id{0}; var_id < t8_data.vars.size(); ++var_id)
            {
                /* Increase the reference counter for each variable holding a reference to it */
                t8_forest_ref(t8_data.initial_forest);
                /* Save a reference of the initial forest in each variable */ 
                t8_data.vars[var_id]->assets->forest = t8_data.initial_forest;
            }
        }
    } else
    {
        /* Currently, in this case the compression mode can only be of type 'One for One' */
        cmc_assert(t8_data.compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ONE_2D || t8_data.compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ONE_3D);

        /* Group variables with the same geo-spatial domain */
        std::vector<std::vector<int>> var_ids;

        /* Assign the other variables to groups */
        /* Flag indicating if a group has been found */
        bool assigned_to_group{false};
        for (size_t id{0}; id < t8_data.vars.size(); ++id)
        {
            /* Reset flag */
            assigned_to_group = false;
            /* Iterate over all previous groups */
            for (size_t group_id{0}; group_id < var_ids.size(); ++group_id)
            {
                /* Check if the data domain coincides with a group's geo-spatial data domain */
                if (compare_geo_domain_equality_of_data_layouts(t8_data.vars[var_ids[group_id].front()]->var->data_layout, t8_data.vars[id]->var->data_layout))
                {
                    assigned_to_group = true;
                    var_ids[group_id].push_back(id);
                    break;
                }
            }

            /* If no common group was found, create a new one */
            if (!assigned_to_group)
            {
                var_ids.push_back(std::vector<int>{});
                var_ids.back().reserve(t8_data.vars.size());
                var_ids.back().push_back(id);
            } 
        }

        /* Create an initial forest for each group */
        t8_cmesh_t cmesh{nullptr};
        t8_forest_t forest{nullptr};
        int initial_ref_lvl{0};

        cmc_debug_msg("Several uniform forests will be constructed. One for each present geo-spatial domain (e.g. Domain: lon x lat, lon x lev, etc.)");

        for (size_t group_id{0}; group_id < var_ids.size(); ++group_id)
        {
            /* Create a new cmesh according to the dimension of the coordinate dimensions */
            cmesh = t8_cmesh_new_hypercube((t8_data.vars[var_ids[group_id].front()]->var->num_dimensions == 2 ? T8_ECLASS_QUAD : T8_ECLASS_HEX), t8_data.comm, 0, 0, 1);

            /* Calculate the initial refinement level for this geo-spatial domain */
            initial_ref_lvl = cmc_t8_calc_geo_refinement_level(t8_data, var_ids[group_id].front());

            /* Allocate assets */
            for (size_t sub_id{0}; sub_id < var_ids[group_id].size(); ++sub_id)
            {
                /* Save a the initial refienment level for each variable */ 
                t8_data.vars[(var_ids[group_id])[sub_id]]->assets = new cmc_t8_assets(initial_ref_lvl);
                t8_data.vars[(var_ids[group_id])[sub_id]]->assets_allocated = true;
            }

            #if CMC_ENABLE_DEBUG
            /* Assertions in Debug-Mode */
            if (t8_data.vars[var_ids[group_id].front()]->var->num_dimensions == 2)
            {
                cmc_assert(initial_ref_lvl < P4EST_MAXLEVEL);
            } else if (t8_data.vars[var_ids[group_id].front()]->var->num_dimensions == 3)
            {
                cmc_assert(initial_ref_lvl < P8EST_MAXLEVEL);
            }
            #endif

            /* Create a new uniform forest with the former calculated initial refinement level */
            forest = t8_forest_new_uniform(cmesh, t8_scheme_new_default_cxx(), initial_ref_lvl, 0, t8_data.comm);

            cmc_debug_msg("Built forest (dim: ", t8_data.vars[var_ids[group_id].front()]->var->num_dimensions, "D) with an intitial refinement level of ", initial_ref_lvl, " for the Data Layout ", data_layout_to_string(t8_data.vars[var_ids[group_id].front()]->var->data_layout));

            /* Coarsen the elements of the forest mesh which are not part of the geo mesh */
            initial_forest = cmc_t8_coarsen_geo_mesh(t8_data, forest, initial_ref_lvl, var_ids[group_id].front());

            /* Save the initial forest and pointer to it accordingly */
            /* Currently, in this case the compression mode can only be of type 'One for One' (see assertion above) */
            for (size_t sub_id{0}; sub_id < var_ids[group_id].size(); ++sub_id)
            {
                /* Increase the reference counter for each variable holding a reference to it */
                t8_forest_ref(initial_forest);
                /* Save a reference of the initial forest in each variable */ 
                t8_data.vars[(var_ids[group_id])[sub_id]]->assets->forest = initial_forest;
            }
            /* Unref the forest once, sonce this forest is not saved anymore (except for one pointer in each variable) */
            t8_forest_unref(&initial_forest);
        }
    }
    #else
    cmc_err_msg("CMC is not compiled with t8code, please reconfigure with t8code linkage in order to use this function.");
    #endif
}

/* Check if a given element is inside the specified geo-spatial data domain */
//TODO: make not based on the "id of the coordiante" but on the actual value (e.g. 135.5 degree)
int
cmc_t8_elem_inside_geo_domain(const t8_element_t* element, t8_eclass_scheme_c* ts, const cmc_t8_data& t8_data, const int var_id,
                              const int lat_min, const int lat_max, const int lon_min, const int lon_max, const int lev_min, const int lev_max)
{
    #ifdef CMC_WITH_T8CODE
    int element_anchor[3];
    t8_default_scheme_common_c* ts_c = static_cast<t8_default_scheme_common_c*>(ts);
    #if CMC_ENABLE_DEBUG
    /* Some assertions in Debug-Mode */
    cmc_assert(t8_data.vars[var_id]->var->num_dimensions == 2 || t8_data.vars[var_id]->var->num_dimensions == 3);
    if (t8_data.vars[var_id]->var->num_dimensions == 2)
    {
        /* 2D case */
        cmc_assert(static_cast<t8_default_scheme_quad_c*>(ts));
    } else if (t8_data.vars[var_id]->var->num_dimensions == 3) {
        /* 3D case */
        cmc_assert(static_cast<t8_default_scheme_hex_c*>(ts));
    }
    #endif

    /* Maximum refinement level depending on the dimension of the data */
    int element_anchor_max_lvl{(t8_data.vars[var_id]->var->num_dimensions == 2 ? P4EST_MAXLEVEL : P8EST_MAXLEVEL)};
    /* Receive the integer anchor coodinates of the element */ 
    ts_c->t8_element_anchor (element, element_anchor);

    /* Transform coordinates into the range of the initial-refinement-level coordinate values */
    for (int index{0}; index < t8_data.vars[var_id]->var->num_dimensions; ++index)
    {
        element_anchor[index] >>= (element_anchor_max_lvl - t8_data.vars[var_id]->assets->initial_refinement_lvl);
    }
    /* Check if the anchor coordinates of the element lie within the "lat x lon x lev" mesh */
    /** \note: Different geo-spatial coordinate combinations are corresponding to certain x- and y- coordinates. This labeling is consistent througout the 'cmc_t8_...'-functions. 
     *         Data only concerning latitude and longitude will be ordered, s.t. longitude equals x, latitude equals y
     *         Data only concerning latitude and elevation will be ordered, s.t. latitude equals x, elevation equals y
     *         Data only concerning longitude and elevation will be ordered, s.t. longitude equals x, elevation equals y
    **/
    if (t8_data.vars[var_id]->var->num_dimensions == 2)
    {
        /* 2D case */
        switch (t8_data.vars[var_id]->var->data_layout)
        {
            case CMC_2D_LAT_LON:
                [[fallthrough]];
            case CMC_2D_LON_LAT:
                if (element_anchor[0] >= lon_min && element_anchor[0] < lon_max &&
                    element_anchor[1] >= lat_min && element_anchor[1] < lat_max)
                {
                    /* The 2D element is inside the "lon x lat" mesh */
                    return 1;
                }
            break;
            case CMC_2D_LAT_LEV:
                [[fallthrough]];
            case CMC_2D_LEV_LAT:
                if (element_anchor[0] >= lat_min && element_anchor[0] < lat_max &&
                    element_anchor[1] >= lev_min && element_anchor[1] < lev_max)
                {
                    /* The 2D element is inside the "lev x lat" mesh */
                    return 1;
                }
            break;
            case CMC_2D_LON_LEV:
                [[fallthrough]];
            case CMC_2D_LEV_LON:
                if (element_anchor[0] >= lon_min && element_anchor[0] < lon_max &&
                    element_anchor[1] >= lev_min && element_anchor[1] < lev_max)
                {
                    /* The 2D element is inside the "lev x lon" mesh */
                    return 1;
                }
            break;
            default:
                cmc_err_msg("There was no valid 2D data layout supplied.");
        }

    } else {
        /* 3D case */
        if (element_anchor[0] >= lon_min && element_anchor[0] < lon_max &&
            element_anchor[1] >= lat_min && element_anchor[1] < lat_max &&
            element_anchor[2] >= lev_min && element_anchor[2] < lev_max)
        {
            /* The 3D is inside the "lat x lon x lev" mesh */
            return 1;
        }
    }
    /* The element is not inside the "lat x lon x lev" mesh */
    return 0;
    #else
    return CMC_ERR;
    #endif
}

/* Check if a given element is inside the geo-spatial data domain */
int
cmc_t8_elem_inside_geo_mesh(const t8_element_t* element, t8_eclass_scheme_c* ts, const cmc_t8_data& t8_data, const int var_id)
{
    #ifdef CMC_WITH_T8CODE
    cmc_assert(t8_data.geo_data->coords->size() >= CMC_NUM_COORD_IDS);
    
    /* Check whether the element is inside the "lat x lon x lev" mesh or not */
    return cmc_t8_elem_inside_geo_domain(element, ts, t8_data, var_id, 0, t8_data.geo_data->get_coord_length(CMC_COORD_IDS::CMC_LAT), 0, t8_data.geo_data->get_coord_length(CMC_COORD_IDS::CMC_LON), 0, t8_data.geo_data->get_coord_length(CMC_COORD_IDS::CMC_LEV));
    #else
    return CMC_ERR;
    #endif
}
