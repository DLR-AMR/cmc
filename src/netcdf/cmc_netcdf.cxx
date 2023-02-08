#include "cmc_netcdf.h"
#include "utilities/cmc_container.hxx"
#include "utilities/cmc_geo_util.h"
#include "utilities/cmc_log_functions.h"

enum CMC_NC_STATUS {STATUS_UNDEFINED = 0, CMC_NC_INQ_COORDS, CMC_NC_ADD_VAR, CMC_NC_INQ_VAR, CMC_NC_READY_TO_COMPRESS};

struct cmc_nc_data
{
private:
    const int ncid;
public:
    
    cmc_nc_data(const int _ncid)
    : ncid{_ncid}{};
    ~cmc_nc_data(){};
    /* These variables only save cooridnate relative values for (latitude, longitude, leverage, time) */
    std::array<int, CMC_NUM_COORD_IDS> coord_dim_ids{CMC_COORDINATE_NOT_CONSIDERED, CMC_COORDINATE_NOT_CONSIDERED, CMC_COORDINATE_NOT_CONSIDERED, CMC_COORDINATE_NOT_CONSIDERED};
    std::array<int, CMC_NUM_COORD_IDS> coord_var_ids{CMC_COORDINATE_NOT_CONSIDERED, CMC_COORDINATE_NOT_CONSIDERED, CMC_COORDINATE_NOT_CONSIDERED, CMC_COORDINATE_NOT_CONSIDERED};
    std::array<size_t, CMC_NUM_COORD_IDS> coord_lengths{0, 0, 0, 0};

    /* Number of dimensions of the netCDF file */
    int num_dimensions{0};
    int num_global_atts{0}; //might be useful later
    int id_unlimited_dim{-1}; //might be useful later for time series data

    /* Vector for storing the coodinate variables */
    var_vector_t coords;

    /* Saves correponding data to each inquired data variable */
    std::vector<cmc_var_t> vars;

    /* These variables save all dimension names and sizes */
    std::vector<size_t> dimension_sizes{};
    std::vector<std::string> dimension_names;

    /* A status flag of the current mode the struct is in */
    CMC_NC_STATUS status{STATUS_UNDEFINED};

    MPI_Comm comm{MPI_COMM_WORLD}; //!< The communicator to use in a parallel environment

    int get_ncid() const;
};

/* Create a cmc_nc_data struct */
cmc_nc_data_t
cmc_nc_create(const int _ncid)
{
    return new cmc_nc_data{_ncid};
}

/* Set a communicator to use in a parallel environment */
void
cmc_nc_set_mpi_communicator(cmc_nc_data_t nc_data, MPI_Comm comm)
{
    nc_data->comm = comm;
}

/* Destroy/deallocate a cmc_nc_data struct */
void
cmc_nc_destroy(const cmc_nc_data_t nc_data)
{
    if (nc_data != nullptr)
    {
        delete nc_data;
    }
}

int
cmc_nc_data::get_ncid() const
{
    return this->ncid;
}

/************************************/
/****** BEGIN STATIC FUNCTIONS ******/
#ifdef CMC_WITH_NETCDF
/** Inquire the dimension lenghts and ids of the (geo-spatial) coordinate dimensions as well as their variables ids */ 
static void
cmc_nc_inquire_coordinate_dims(cmc_nc_data& nc_data)
{
    #ifdef CMC_WITH_NETCDF
    std::basic_string<char> nc_name;
    nc_name.reserve(NC_MAX_NAME +1);
    size_t dim_length;

    /* Inquire the most general information about the supplied netCDF file (like number of diemnsion, number of global attributes, number of variables, etc.)*/
    int err{nc_inq(nc_data.get_ncid(), &(nc_data.num_dimensions), NULL, &(nc_data.num_global_atts), &(nc_data.id_unlimited_dim))};
    cmc_nc_check_err(err);

    /* Reserve memory for all 'dimension_lenghts' and 'dimension_names' */
    nc_data.dimension_sizes.reserve(nc_data.num_dimensions);
    nc_data.dimension_names.reserve(nc_data.num_dimensions);
    /* Assign a default string as dimension name for each dimension */
    for (int dim_index{0}; dim_index < nc_data.num_dimensions; ++dim_index)
    {
        nc_data.dimension_names.push_back(std::string());
    }

    /* Loop over all dimension within the netCDF file */
    for (int dim_index{0}; dim_index < nc_data.num_dimensions; ++dim_index) 
    {
        /* Get the name and the length of the dimension which corresponds to the current id */
        err = nc_inq_dim(nc_data.get_ncid(), dim_index, &nc_name[0], &dim_length);
        cmc_nc_check_err(err);

        /* Check for coordinate dimension if a hint to the id was given */
        if (nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_LON] == dim_index)
        {
            nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_LON] = dim_index;
            nc_data.coord_lengths[CMC_COORD_IDS::CMC_LON] = dim_length;
        } else if (nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_LAT] == dim_index)
        {
            nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_LAT] = dim_index;
            nc_data.coord_lengths[CMC_COORD_IDS::CMC_LAT] = dim_length;   
        } else if (nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_LEV] == dim_index)
        {
            nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_LEV] = dim_index;
            nc_data.coord_lengths[CMC_COORD_IDS::CMC_LEV] = dim_length;
        } else if (nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_TIME] == dim_index)
        {
            nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_TIME] = dim_index;
            nc_data.coord_lengths[CMC_COORD_IDS::CMC_TIME] = dim_length;
        }
        /* In case no hints concerning the coordinate dimensions were given */
        /* Check if the name is equal to any of the considered (geo-spatial) coordinate variables */
        else if (nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_LON] == CMC_COORDINATE_NOT_CONSIDERED && 
                   (strcmp(nc_name.c_str(), "lon") == 0 || strcmp(nc_name.c_str(), "longitude") == 0))
        {
            nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_LON] = dim_index;
            nc_data.coord_lengths[CMC_COORD_IDS::CMC_LON] = dim_length;
        } else if (nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_TIME] == CMC_COORDINATE_NOT_CONSIDERED &&
                   (strcmp(nc_name.c_str(), "lat") == 0 || strcmp(nc_name.c_str(), "latitude") == 0))
        {
            nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_LAT] = dim_index;
            nc_data.coord_lengths[CMC_COORD_IDS::CMC_LAT] = dim_length;
        }
        else if (nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_TIME] == CMC_COORDINATE_NOT_CONSIDERED &&
                 (strcmp(nc_name.c_str(), "lev") == 0 || strcmp(nc_name.c_str(), "leverage") == 0 || strcmp(nc_name.c_str(), "height") == 0))
        {
            nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_LEV] = dim_index;
            nc_data.coord_lengths[CMC_COORD_IDS::CMC_LEV] = dim_length;
        } else if (nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_TIME] == CMC_COORDINATE_NOT_CONSIDERED && strcmp(nc_name.c_str(), "time") == 0)
        {
            nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_TIME] = dim_index;
            nc_data.coord_lengths[CMC_COORD_IDS::CMC_TIME] = dim_length;
        }

        /* Save the name and the length of each dimension within the netCDF file */
        nc_data.dimension_names[dim_index].assign(nc_name.data());
        nc_data.dimension_sizes[dim_index] = dim_length;
    }
    cmc_debug_msg("The following coordinate dimensions were retrieved:");
    cmc_debug_msg("LAT: ", nc_data.coord_lengths[CMC_COORD_IDS::CMC_LAT], ", LON: ", nc_data.coord_lengths[CMC_COORD_IDS::CMC_LON], ", LEV: ", nc_data.coord_lengths[CMC_COORD_IDS::CMC_LEV], ", T: ", nc_data.coord_lengths[CMC_COORD_IDS::CMC_TIME]);
    #endif
}

/** Inquire the ids of the coordinate variables and allocate memory in order to store their values */
static void
cmc_nc_inquire_coordiante_vars(cmc_nc_data& nc_data)
{
    #ifdef CMC_WITH_NETCDF
    int err, var_type;

    /* Loop over all possibly-considered coordinate variables (latitude, longitude, leverage, time) */
    for (int coord_ids{0}; coord_ids < CMC_NUM_COORD_IDS; ++coord_ids)
    {
        /* Check if the considered coordinate variable/dimension is supplied within the netCDF file */
        if (nc_data.coord_dim_ids[coord_ids] != CMC_COORDINATE_NOT_CONSIDERED)
        {
            /* Coordinate variables have a concerning coordinate dimension; the name of the dimension coincides with the name of the variable and is only dependent on "its own dimension" */
            /* Get the id of the coordinate variable */
            err = nc_inq_varid(nc_data.get_ncid(), nc_data.dimension_names[nc_data.coord_dim_ids[coord_ids]].c_str(), &(nc_data.coord_var_ids[coord_ids]));
            cmc_nc_check_err(err);

            /* Get the type of coordinate variables (e.g. float) */
            err = nc_inq_vartype(nc_data.get_ncid(), nc_data.coord_var_ids[coord_ids], &(var_type));
            cmc_nc_check_err(err);

            /* Allocate space for each (geo-spatial) coordinate variable */
            nc_data.coords.create_and_push_back(static_cast<size_t>(nc_data.coord_lengths[coord_ids]), static_cast<cmc_type>(var_type));
            cmc_debug_msg("Information about ", get_coord_name(static_cast<CMC_COORD_IDS>(coord_ids)), " was inquired.");
        } else
        {
            /* Push back a placeholder in order to retain the internal used coordinate order */
            nc_data.coords.push_back_placeholder();
        }
    }
    #endif
}

/** Store the geo-spatial coordinate data (latitude, longitude, leverage, time) */
static void
cmc_nc_inquire_coordinate_data(cmc_nc_data& nc_data)
{
    #ifdef CMC_WITH_NETCDF
    int err{0};

    /* Inquire the data of the latitude variable */
    if (nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_LAT] != CMC_COORDINATE_NOT_CONSIDERED)
    {
        err = nc_get_var(nc_data.get_ncid(), nc_data.coord_var_ids[CMC_COORD_IDS::CMC_LAT], nc_data.coords[CMC_COORD_IDS::CMC_LAT].get_initial_data_ptr());
        cmc_nc_check_err(err);
        cmc_debug_msg("Coordinate data of variable ", get_coord_name(CMC_LAT), " has been inquired.");
    }

    /* Inquire the data of the longitude variable */
    if (nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_LON] != CMC_COORDINATE_NOT_CONSIDERED)
    {
        err = nc_get_var(nc_data.get_ncid(), nc_data.coord_var_ids[CMC_COORD_IDS::CMC_LON], nc_data.coords[CMC_COORD_IDS::CMC_LON].get_initial_data_ptr());
        cmc_nc_check_err(err);
        cmc_debug_msg("Coordinate data of variable ", get_coord_name(CMC_LON), " has been inquired.");
    }

    /* Inquire the data of the leverage variable */
    if (nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_LEV] != CMC_COORDINATE_NOT_CONSIDERED)
    {
        err = nc_get_var(nc_data.get_ncid(), nc_data.coord_var_ids[CMC_COORD_IDS::CMC_LEV], nc_data.coords[CMC_COORD_IDS::CMC_LEV].get_initial_data_ptr());
        cmc_nc_check_err(err);
        cmc_debug_msg("Coordinate data of variable ", get_coord_name(CMC_LEV), " has been inquired.");
    }

    /* Inquire the data of the time variable */
    if (nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_TIME] != CMC_COORDINATE_NOT_CONSIDERED)
    {
        err = nc_get_var(nc_data.get_ncid(), nc_data.coord_var_ids[CMC_COORD_IDS::CMC_TIME], nc_data.coords[CMC_COORD_IDS::CMC_TIME].get_initial_data_ptr());
        cmc_nc_check_err(err);
        cmc_debug_msg("Coordinate data of variable ", get_coord_name(CMC_TIME), " has been inquired.");
    }
    #endif
}

/** Check if given attributes are stored alongside the inquired variables
 * \note Only if the data type of the attribute is (signed/unsigned) short, int, int64, as well as double or float, the data value will be stored (and applied later on) */  
static void
cmc_nc_check_for_attributes(cmc_nc_data& nc_data, const std::vector<const char*>& att_names)
{
    #ifdef CMC_WITH_NETCDF
    int err;
    int attribute_type;
    cmc_universal_type_t* att_ptr{nullptr};

    for (size_t i{0}; i < nc_data.vars.size(); ++i)
    {
        cmc_debug_msg("Checking attributes of variable ", nc_data.vars[i]->name, "...");
        for (auto iter_att{att_names.begin()}; iter_att != att_names.end(); ++iter_att)
        {
            /* Check which possible attribute is considered */
            if (strcmp(*iter_att, "add_offset") == 0)
            {
                att_ptr = &(nc_data.vars[i]->add_offset);
            } else if (strcmp(*iter_att, "scale_factor") == 0)
            {
                att_ptr = &(nc_data.vars[i]->scale_factor);
            } else if (strcmp(*iter_att, "missing_value") == 0)
            {
                att_ptr = &(nc_data.vars[i]->missing_value);
            } else
            {
                att_ptr = nullptr;
            }

            /* Check the data type of the attribute */
            err = nc_inq_atttype(nc_data.get_ncid(), nc_data.vars[i]->var_id, *iter_att, &attribute_type);

            if (err == NC_NOERR)
            {
                cmc_debug_msg("Found a/an '", *iter_att, "' attribute.");
                if (strcmp(*iter_att, "missing_value") == 0)
                {
                    nc_data.vars[i]->missing_value_present = true;
                }
                /* If the data type of the attribute has been successfully inquired, get the data value of the attribute */
                switch (attribute_type)
                {
                    case NC_SHORT:
                        int16_t short_att;
                        err = nc_get_att_short(nc_data.get_ncid(), nc_data.vars[i]->var_id, *iter_att, &short_att);
                        cmc_nc_check_err(err);
                        *att_ptr = short_att;
                        break;
                    case NC_INT:
                        int32_t int_att;
                        err = nc_get_att_int(nc_data.get_ncid(), nc_data.vars[i]->var_id, *iter_att, &int_att);
                        cmc_nc_check_err(err);
                        *att_ptr = int_att;
                        break;
                    case NC_FLOAT:
                        float float_att;
                        err = nc_get_att_float(nc_data.get_ncid(), nc_data.vars[i]->var_id, *iter_att, &float_att);
                        cmc_nc_check_err(err);
                        *att_ptr = float_att;
                        break;
                    case NC_DOUBLE:
                        double double_att;
                        err = nc_get_att_double(nc_data.get_ncid(), nc_data.vars[i]->var_id, *iter_att, &double_att);
                        cmc_nc_check_err(err);
                        *att_ptr = double_att;
                        break;
                    case NC_USHORT:
                        uint16_t ushort_att;
                        err = nc_get_att_ushort(nc_data.get_ncid(), nc_data.vars[i]->var_id, *iter_att, &ushort_att);
                        cmc_nc_check_err(err);
                        *att_ptr = ushort_att;
                        break;
                    case NC_UINT:
                        uint32_t uint_att;
                        err = nc_get_att_uint(nc_data.get_ncid(), nc_data.vars[i]->var_id, *iter_att, &uint_att);
                        cmc_nc_check_err(err);
                        *att_ptr = uint_att;
                        break;
                    case NC_INT64:
                        long long int longlong_att;
                        err = nc_get_att_longlong(nc_data.get_ncid(), nc_data.vars[i]->var_id, *iter_att, &longlong_att);
                        cmc_nc_check_err(err);
                        *att_ptr = static_cast<int64_t>(longlong_att);
                        break;
                    case NC_UINT64:
                        unsigned long long int ulonglong_att;
                        err = nc_get_att_ulonglong(nc_data.get_ncid(), nc_data.vars[i]->var_id, *iter_att, &ulonglong_att);
                        cmc_nc_check_err(err);
                        *att_ptr = static_cast<uint64_t>(ulonglong_att);
                        break;
                    default:
                        cmc_msg("The attribute ", *iter_att, " of variable ", nc_data.vars[i]->name, " could not be inquired (and will not be applied), since its data type is not supported.");
                }
            }
        }
    }
    #endif
};

/** Inquire meta data about the variables, for example the id within the netCDF file, 'offset'- or 'scale'-attributes, etc. */
static void
cmc_nc_inquire_var_meta_data(cmc_nc_data_t nc_data)
{
    #ifdef CMC_WITH_NETCDF
    int err{0};
    std::array<int, CMC_NC_MAX_VAR_DIMS> dim_ids;

    for (size_t i{0}; i < nc_data->vars.size(); ++i)
    {
        /* Inquire the variable id of the given variable's name */
        err = nc_inq_varid(nc_data->get_ncid(), nc_data->vars[i]->name.c_str(), &(nc_data->vars[i]->var_id));
        if (err == NC_NOERR)
        {
            /* If an ID is found, information like number of dimensions, data type, etc. will be inquired */
            err = nc_inq_var(nc_data->get_ncid(), nc_data->vars[i]->var_id, NULL, &(nc_data->vars[i]->var_type), &(nc_data->vars[i]->num_dimensions), &(dim_ids[0]), NULL);
            cmc_nc_check_err(err);
            
            /* Calculate the number of elements */
            nc_data->vars[i]->dimension_ids.reserve(nc_data->vars[i]->num_dimensions);
            for (int j{0}; j < nc_data->vars[i]->num_dimensions; ++j)
            {
                nc_data->vars[i]->dimension_ids.push_back(dim_ids[j]);
            }
            
        } else {
            /* Print an error message if no ID to the supplied name of the variable is found */
            cmc_msg("An error occured while inquiring the netCDf 'variable_id' corresponding to the supplied 'variable_name' (", nc_data->vars[i]->name, ").\n", "The netCDF error reads: \n");
            cmc_nc_check_err(err);
        }
    }
    
    /* Check if an 'offset'-attribute, 'missing_value'-attribute or 'scale_factor'-attribute is given for each variable */
    const std::vector<const char*> attributes{"add_offset", "scale_factor", "missing_value"};
    cmc_nc_check_for_attributes(*nc_data, attributes);
    #endif 
}

//TODO: Make the distribution more general. Currently, a somehow evenly hyperslab is assigned to each rank
static
void
cmc_nc_calc_offsets_for_parallel_reading(cmc_nc_data_t nc_data, const size_t current_var_id)
{
    #ifdef CMC_WITH_NETCDF_PAR
    //TODO: Currently, the data is divided into similar blocks (maybe this is an average approach between the amount of communications needed (in order to reorder the data) and the cache misses during the reading of the file) ?
    int err, size, rank;
    /* Get the size of the MPI communicator */
    err = MPI_Comm_size(nc_data->comm, &size);
    cmc_mpi_check_err(err);
    /* Get the rank id */
    err = MPI_Comm_rank(nc_data->comm, &rank);
    cmc_mpi_check_err(err);
    uint32_t leftover_dim_points{0}, evenly_split_dim_points{0};
    uint32_t global_offset{0};
    for (int dims{0}; dims < nc_data->vars[current_var_id]->num_dimensions; ++dims)
    {
        /* If a dimension is not considered in the data hyperslab, just skip it */
        if (nc_data->vars[current_var_id]->count_ptr[dims] <= 1)
        {
            continue;
        } else
        {
            cmc_assert(nc_data->vars[current_var_id]->count_ptr[dims] >= static_cast<size_t>(size));
            /* If the dimension is considered, calculate an offset for this dimension */
            //The assertion above indicates that at least one dimension point is assigned to each rank
            /* The leftover dimension points will be distributed */
            evenly_split_dim_points = static_cast<uint32_t>(nc_data->vars[current_var_id]->count_ptr[dims] / size);
            /* Calculate the offset for the start vector */
            global_offset = rank * evenly_split_dim_points;
            /* Calculate the points which cannot be evenly split between all processses */
            leftover_dim_points = nc_data->vars[current_var_id]->count_ptr[dims] % size;
            /* Adjust the offset for the leftover points */
            /* Rank 0 starts at position zero nevertheless, therefore, the offset cannot change */
            if (rank != 0)
            {
                if (rank < static_cast<int>(leftover_dim_points))
                {
                    /* The ranks (starting from the lowest to the highest id) obtain another additional dimension point when the diemnsion cannot be split evenly */
                    ++(evenly_split_dim_points);
                    /* Adjust the global offset */
                    global_offset += leftover_dim_points - rank;
                } else
                {
                    /* Adjust the global offset */
                    global_offset += leftover_dim_points;
                }
            }
            /* Save the start and count values for this dimension */
            nc_data->vars[current_var_id]->start_ptr[dims] += global_offset;
            nc_data->vars[current_var_id]->count_ptr[dims] = evenly_split_dim_points;
        }
    }
    #endif
}
#endif
/******* END STATIC FUNCTIONS *******/
/************************************/


/** Inquire all necessary information, especially about the geo-spatial coordinate dimension and variables (latitude, longitude, and leverage) */
void
cmc_inquire_coordinates(cmc_nc_data_t nc_data)
{
    #ifdef CMC_WITH_NETCDF
    /* Set status flag */
    nc_data->status = CMC_NC_STATUS::CMC_NC_INQ_COORDS;

    /* Inquire dimension names, sizes and ids */
    cmc_nc_inquire_coordinate_dims(*nc_data);

    /* Inquire the cooresponding (geo-spatial) coordinate variables (to the dimension) and their ids */
    cmc_nc_inquire_coordiante_vars(*nc_data);

    /* Retrieve the data of these coordinate variables */
    cmc_nc_inquire_coordinate_data(*nc_data);
    #endif
}


/** Inquire the data for each given variable (either the variable as a whole (start_ptr and count_ptr equal to nullptrs) or a specififed hyperslab) */
void
cmc_nc_inquire_var_data(cmc_nc_data_t nc_data, const size_t* start_ptr, const size_t* count_ptr)
{
    #ifdef CMC_WITH_NETCDF
    cmc_assert(nc_data->status < CMC_NC_STATUS::CMC_NC_INQ_VAR && nc_data->status > CMC_NC_STATUS::STATUS_UNDEFINED);
    cmc_assert(nc_data->vars.size() > 0);

    /* Set status flag */
    nc_data->status = CMC_NC_STATUS::CMC_NC_INQ_VAR;

    int err;
    size_t num_data_points{1};

    /* Vectors storing the start and count values for reading the data */
    std::vector<size_t> start_values;
    std::vector<size_t> count_values;

    #ifdef CMC_WITH_NETCDF_PAR
    int comm_size;
    #endif

    /* Read the meta data of the variable (for example 'data_type', 'missing_value', 'offset', 'scale_factor', ...) */
    cmc_nc_inquire_var_meta_data(nc_data);

    /* Preallocate data arrays */
    for (size_t i{0}; i < nc_data->vars.size(); ++i)
    {
        cmc_assert(nc_data->vars[i]->var_id != CMC_NC_VAR_NOT_CONSIDERED);

        /* Reset the number of data points */
        num_data_points = 1;

        /* Reserve space for storing the lengths of the data of each data dimension */
        nc_data->vars[i]->dim_lengths.reserve(nc_data->vars[i]->num_dimensions);
        /* Allocate memeory for the start and the count ptr */
        nc_data->vars[i]->start_ptr.reserve(nc_data->vars[i]->num_dimensions);
        nc_data->vars[i]->count_ptr.reserve(nc_data->vars[i]->num_dimensions);
        /* Save the start and count pointers associated with the hyperslab */
        if (start_ptr != nullptr && count_ptr != nullptr)
        {
            /* If only a hyperslab will be read */
            /* Store the data concerning the defined hyperslab */
            for (int dims{0}; dims < nc_data->vars[i]->num_dimensions; ++dims)
            {
                nc_data->vars[i]->start_ptr.push_back(static_cast<uint32_t>(start_ptr[dims]));
                nc_data->vars[i]->count_ptr.push_back(static_cast<uint32_t>(count_ptr[dims]));
            }
        } else
        {
            /* If the data of the whole varibale will be read */
            for (int dims{0}; dims < nc_data->vars[i]->num_dimensions; ++dims)
            {
                nc_data->vars[i]->start_ptr.push_back(0);
                nc_data->vars[i]->count_ptr.push_back(static_cast<uint32_t>(nc_data->dimension_sizes[nc_data->vars[i]->dimension_ids[dims]]));
            }
        }

        /* Adjust the start and count vectors in a parallel environment */
        #ifdef CMC_WITH_NETCDF_PAR
        /* Get the size of the communicator in order to check whether the data can be read in parallel or not */
        err = MPI_Comm_size(nc_data->comm, &comm_size);
        cmc_mpi_check_err(err);
        /* If there is more than one process active */
        if (comm_size > 1)
        {
            /* If the data may be read in distributed, calculate the offsets */
            /* Before this call the full dimension lengths of the variable are stored in the start_ptr and count_ptr vectors */
            cmc_nc_calc_offsets_for_parallel_reading(nc_data, i);
        }
        #endif

        /* Calculate the number of data points (neeeded for the array allocation) */
        for (int dims{0}; dims < nc_data->vars[i]->num_dimensions; ++dims)
        {
            (nc_data->vars[i]->dim_lengths).push_back(count_ptr[dims]);
            num_data_points *= nc_data->vars[i]->count_ptr[dims];
        }

        /* Allocate memory for an array able to hold the data of the variable */
        nc_data->vars[i]->data = new var_array_t{num_data_points, static_cast<cmc_type>(nc_data->vars[i]->var_type)};

        /* Allocate space for the vectors defining the data inquisition */
        start_values.reserve(nc_data->vars[i]->num_dimensions);
        count_values.reserve(nc_data->vars[i]->num_dimensions);

        /* Set the vectors for reading the data */
        for (int dim_id{0}; dim_id < nc_data->vars[i]->num_dimensions; ++dim_id)
        {
            start_values[dim_id] = static_cast<size_t>(nc_data->vars[i]->start_ptr[dim_id]);
            count_values[dim_id] = static_cast<size_t>(nc_data->vars[i]->count_ptr[dim_id]);
        }

        /* Decide whther the data will be read in parallel or serial */
        #ifdef CMC_WITH_NETCDF_PAR
        /* Get the size of the communicator in order to check whether the data can be read in parallel or not */
        err = MPI_Comm_size(nc_data->comm, &comm_size);
        cmc_mpi_check_err(err);
        /* If there is more than one process active */
        if (comm_size > 1)
        {
            /* If the data may be read in distributed, calculate the offsets */
            /* Before this call the full dimension lengths of the variable are stored in the start_ptr and count_ptr vectors */
            err = nc_get_vara(nc_data->get_ncid(), nc_data->vars[i]->var_id, &(start_values[0]), &(count_values[0]), nc_data->vars[i]->data->get_initial_data_ptr());
            cmc_nc_check_err(err);
        }
        else
        #endif
        {
        /* Inquire and store the variable's data */
        if (start_ptr == nullptr && count_ptr == nullptr)
        {
            /* Store data of the whole variable if no hyperslab is defined */
            err = nc_get_var(nc_data->get_ncid(), nc_data->vars[i]->var_id, nc_data->vars[i]->data->get_initial_data_ptr());
            cmc_nc_check_err(err);
        } else {
            /* Store only the data of the specified hyperslab */
            err = nc_get_vara(nc_data->get_ncid(), nc_data->vars[i]->var_id, &(start_values[0]), &(count_values[0]), nc_data->vars[i]->data->get_initial_data_ptr());
            cmc_nc_check_err(err);
        }
        }
        /* Save the data scheme */
        nc_data->vars[i]->data_scheme = CMC_DATA_ORDERING_SCHEME::CMC_GEO_DATA_LINEAR;

        cmc_debug_msg("The data of variable ", nc_data->vars[i]->name, " has been inquired.");
    }

    #endif
}

void
cmcc_nc_inquire_vars(cmc_nc_data_t nc_data, const size_t* start_ptr, const size_t* count_ptr, const int var_count, ...)
{
    #ifdef CMC_WITH_NETCDF

    va_list args;
    const char* current_var_name;

    /* Reserve space for all variables */
    nc_data->vars.reserve(var_count);

    /* Inquire information about dimensions and read coordinate variables */
    cmc_inquire_coordinates(nc_data);

    /* Initializing  variable argument list */
    va_start(args, var_count);

    /* Iterate over all variables */
    for (int _id{0}; _id < var_count; ++_id)
    {
        /* Get the next varibale name */
        current_var_name = va_arg(args, const char*);
        /* Create a variable based on its name */
        nc_data->vars.push_back(new cmc_var{std::string(current_var_name)});
    }

    /* End of varibale argument list */
    va_end(args);

    /* Inquire the data of these variables */
    cmc_nc_inquire_var_data(nc_data, start_ptr, count_ptr);

    /* Set status flag */
    nc_data->status = CMC_NC_STATUS::CMC_NC_READY_TO_COMPRESS;

    #else
    cmc_err_msg("CMC is not compiled with netCDF, please reconfigure with netCDF linkage in order to use this function.");
    #endif
}

/** Open the netCDF file */
int
cmc_nc_open(const char* path_to_file, MPI_Comm comm)
{
    #ifdef CMC_WITH_NETCDF
    int err{0};
    int ncid{0};

    #ifdef CMC_WITH_NETCDF_PAR
    int comm_size{0};
    err = MPI_Comm_size(comm, &comm_size);
    cmc_mpi_check_err(err);
    if (comm_size > 1)
    {
        MPI_Info info = MPI_INFO_NULL;
        err = nc_open_par(path_to_file, NC_NOWRITE, comm, info, &ncid);
        cmc_nc_check_err(err);
    } else
    #endif
    {
        err = nc__open(path_to_file, NC_NOWRITE, NULL, &ncid);
        cmc_nc_check_err(err);
    }

    return ncid;
    #else
    cmc_err_msg("CMC is not linked against netCDF. Please recompile the build with netCDF\n.");
    return CMC_ERR;
    #endif
}

/** Close the netCDF file */
void
cmc_nc_close(int ncid)
{
    #ifdef CMC_WITH_NETCDF
    int err{nc_close(ncid)};
    cmc_nc_check_err(err);
    #endif
}

void
cmc_nc_set_hint_latitude_dim(cmc_nc_data_t nc_data, const int lat_dim_id)
{
    #ifdef CMC_WITH_NETCDF
    cmc_assert(lat_dim_id >= 0);
    nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_LAT] = lat_dim_id;
    #endif
}

void
cmc_nc_set_hint_longitude_dim(cmc_nc_data_t nc_data, const int lon_dim_id)
{
    #ifdef CMC_WITH_NETCDF
    cmc_assert(lon_dim_id >= 0);
    nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_LON] = lon_dim_id;
    #endif
}

void
cmc_nc_set_hint_elevation_dim(cmc_nc_data_t nc_data, const int lev_dim_id)
{
    #ifdef CMC_WITH_NETCDF
    cmc_assert(lev_dim_id >= 0);
    nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_LEV] = lev_dim_id;
    #endif
}

void
cmc_nc_set_hint_time_dim(cmc_nc_data_t nc_data, const int time_dim_id)
{
    #ifdef CMC_WITH_NETCDF
    cmc_assert(time_dim_id >= 0);
    nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_TIME] = time_dim_id;
    #endif
}

void
cmc_nc_add_variable(cmc_nc_data_t nc_data, const char* var_name)
{
    #if CMC_WITH_NETCDF
    cmc_assert(nc_data->status <= CMC_NC_STATUS::CMC_NC_ADD_VAR);
    /* Set status flag */
    nc_data->status = CMC_NC_STATUS::CMC_NC_ADD_VAR;

    /* Allocate a new variable */
    nc_data->vars.push_back(new cmc_var{std::string(var_name)});
    #endif
}

void
cmcc_nc_inquire_added_vars(cmc_nc_data_t nc_data, const size_t* start_ptr, const size_t* count_ptr)
{
    #ifdef CMC_WITH_NETCDF
    cmc_assert(nc_data->vars.size() > 0 &&
               nc_data->status < CMC_NC_STATUS::CMC_NC_INQ_VAR);

    /* Inquire information about dimensions and read coordinate variables */
    cmc_inquire_coordinates(nc_data);

    /* Inquire the data of these variables */
    cmc_nc_inquire_var_data(nc_data, start_ptr, count_ptr);

    /* Set status flag */
    nc_data->status = CMC_NC_STATUS::CMC_NC_READY_TO_COMPRESS;

    #else
    cmc_err_msg("CMC is not compiled with netCDF, please reconfigure with netCDF linkage in order to use this function.");
    #endif
}

/** Create a 'cmc_nc_meta' class for each given variable name */
void
cmc_nc_preallocate_vars_by_name(cmc_nc_data_t nc_data, std::string&& var_name)
{
    #ifdef CMC_WITH_NETCDF
    /** \note This is the base function of the variadic template */
    /* Pre-allocate variables by constructing new meta classes and save the given names of these variables */
    nc_data->vars.push_back(new cmc_var(var_name));
    #endif
}

/* Internally used functions */
void
_cmc_nc_push_back_var(cmc_nc_data_t nc_data, std::string&& var_name)
{
    #ifdef CMC_WITH_NETCDF
    nc_data->vars.push_back(new cmc_var(std::move(var_name)));
    #endif
}
void
_cmc_nc_reserve_vars(cmc_nc_data_t nc_data, const size_t num_variables)
{
    #ifdef CMC_WITH_NETCDF
    nc_data->vars.reserve(num_variables);
    #endif
}

#ifdef CMC_WITH_T8CODE
void
_cmc_transform_nc_data_to_t8code_data(cmc_nc_data_t nc_data, cmc_t8_data_t t8_data)
{
    #ifdef CMC_WITH_NETCDF
    cmc_debug_msg("Setting up the compression with netCDF data and t8code.");

    /* Create a new geo_data struct */
    t8_data->geo_data = new cmc_t8_geo_data();
    /* Move the vector containing the coordinate variables */
    t8_data->geo_data->coords = new var_vector_t(std::move(nc_data->coords));
    
    /* Reserve memory for all variables */
    t8_data->vars.reserve(nc_data->vars.size());

    /* Retrieve the axis ordering of all variables */
    /* While doing this, save the maximum dimensionality of the variables */
    int max_dim{0};
    int update_num_dims;
    /* Create a new dimension length vector for this varibale. This member is used differently in nc_functions than in t8_functions... */
    std::vector<size_t> cmc_t8_dim_lengths(CMC_NUM_COORD_IDS);

    for (size_t var_id{0}; var_id < nc_data->vars.size(); ++var_id)
    {
        /* Reset update variable */
        update_num_dims = 0;

        /* Assign the netCDF variable to the t8_variable */
        t8_data->vars.push_back(new cmc_t8_var(nc_data->vars[var_id]));
        nc_data->vars[var_id] = nullptr;

        /* Check if the variable has a missing_value, otherwise we will assign an arbitrary missing_value */
        if (t8_data->vars[var_id]->var->missing_value.index() == 0)
        {
            //TODO: this is not super clean, because in the cmc_universal_type_t's first position is a std::byte which is theoretically normally usable, there should be a dummy used for this check */
            t8_data->vars[var_id]->var->assign_an_arbitrary_missing_value();
        }

        /* Nullify the new dim_lengths vector */
        std::fill(cmc_t8_dim_lengths.begin(), cmc_t8_dim_lengths.end(), 1);

        /* Reserve space for the data axis ordering */
        t8_data->vars[var_id]->var->axis_ordering.reserve(t8_data->vars[var_id]->var->num_dimensions);

        /* Retrieve the axis ordering of the variable */
        for(int i{0}; i < t8_data->vars[var_id]->var->num_dimensions; ++i)
        {
            if (t8_data->vars[var_id]->var->dimension_ids[i] == nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_LON])
            {
                if (t8_data->vars[var_id]->var->dim_lengths[i] > 1)
                {
                    t8_data->vars[var_id]->var->axis_ordering.emplace_back(CMC_COORD_IDS::CMC_LON);
                    cmc_t8_dim_lengths[CMC_COORD_IDS::CMC_LON] = t8_data->vars[var_id]->var->dim_lengths[i];
                } else
                {
                    ++update_num_dims;
                }
            } else if (t8_data->vars[var_id]->var->dimension_ids[i] == nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_LAT])
            {
                if (t8_data->vars[var_id]->var->dim_lengths[i] > 1)
                {
                    t8_data->vars[var_id]->var->axis_ordering.emplace_back(CMC_COORD_IDS::CMC_LAT);
                    cmc_t8_dim_lengths[CMC_COORD_IDS::CMC_LAT] = t8_data->vars[var_id]->var->dim_lengths[i];
                } else
                {
                    ++update_num_dims;
                }
            } else if (t8_data->vars[var_id]->var->dimension_ids[i] == nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_LEV])
            {
                if (t8_data->vars[var_id]->var->dim_lengths[i] > 1)
                {
                    t8_data->vars[var_id]->var->axis_ordering.emplace_back(CMC_COORD_IDS::CMC_LEV);
                    cmc_t8_dim_lengths[CMC_COORD_IDS::CMC_LEV] = t8_data->vars[var_id]->var->dim_lengths[i];
                } else
                {
                    ++update_num_dims;
                }
            } else if (t8_data->vars[var_id]->var->dimension_ids[i] == nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_TIME])
            {
                //Time series are currently not supported, therefore the time coordinates are skipped
                #if 0
                if (nc_data->coord_lengths[CMC_COORD_IDS::CMC_TIME] > 1)
                {
                    t8_data->vars[var_id]->var->axis_ordering.emplace_back(CMC_COORD_IDS::CMC_TIME);
                }
                #else
                ++update_num_dims;
                #endif
            } else
            {
                if (t8_data->vars[var_id]->var->dim_lengths[i] > 1)
                {
                    cmc_err_msg("The data (hyperslab) of variable (name: ", nc_data->vars[var_id]->name, ") is dependent on non-geo-spatial dimensions.");
                }
                ++update_num_dims;
            }
        }

        /* Update the dimension dependency of the variables */
        t8_data->vars[var_id]->var->num_dimensions -= update_num_dims;

        /* Save the maximum dimensionality */
        if (t8_data->vars[var_id]->var->num_dimensions > max_dim)
        {
            max_dim = t8_data->vars[var_id]->var->num_dimensions;
        }

        /* Set the new dimension lenghts */
        t8_data->vars[var_id]->var->dim_lengths = cmc_t8_dim_lengths;

        /* Retrieve the data layout of the variable */
        t8_data->vars[var_id]->var->get_data_layout_from_axis_ordering();
    }

    /* Save the maximum geo-spatial dimension of the data */
    t8_data->geo_data->dim = max_dim;

    #endif
}
#endif

#ifdef __cplusplus
[[noreturn]]
#else
_Noreturn
#endif
void
cmc_netcdf_exit(const int _err_code, const char* _location)
{
    std::cout << "CMC_NETCDF_EXIT is invoked..." << std::endl << "A netCDF-Error occured, Code: " << _err_code << std::endl << nc_strerror(_err_code) << std::endl << "In: " << _location << std::endl;
    std::exit(EXIT_FAILURE);
}
