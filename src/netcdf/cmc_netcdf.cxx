#include "netcdf/cmc_netcdf.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_geo_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"

#include <utility>

//TODO: remove
//#define CMC_WITH_NETCDF 1
//#define CMC_WITH_NETCDF_PAR 1

namespace cmc
{

[[noreturn]] void
NcExit(const int _err_code, const char* _location)
{
    std::cout << "CMC_NETCDF_EXIT is invoked..." << std::endl << "A netCDF-Error occured, Code: " << _err_code << std::endl << nc_strerror(_err_code) << std::endl << "In: " << _location << std::endl;
    std::exit(EXIT_FAILURE);
}

int
NcOpenSerial(const char* path_to_file)
{
    #ifdef CMC_WITH_NETCDF
    int err;
    int ncid;

    /* Open the file without explicit parallel access */
    err = nc__open(path_to_file, NC_NOWRITE, NULL, &ncid);
    NcCheckError(err);

    return ncid;
    #else
    cmc_err_msg("CMC is not linked against netCDF. Please recompile the build with netCDF\n.");
    return CMC_ERR;
    #endif
}

int
NcOpenParallel(const char* path_to_file, MPI_Comm comm)
{
    #ifdef CMC_WITH_NETCDF_PAR
    int err;
    int ncid;

    int comm_size;
    err = MPI_Comm_size(comm, &comm_size);
    MPICheckError(err);

    MPI_Info info = MPI_INFO_NULL;
    err = nc_open_par(path_to_file, NC_NOWRITE, comm, info, &ncid);
    NcCheckError(err);


    return ncid;
    #else
    cmc_err_msg("CMC is not linked against netCDF (at leat no parallel access functions are available). Please recompile the build with netCDF\n.");
    return CMC_ERR;
    #endif
}

void
NcData::NcOpen(const std::string& path_to_file, const NcOpeningMode mode, const MPI_Comm comm)
{
    #ifdef CMC_WITH_NETCDF
    #ifdef CMC_WITH_NETCDF_PAR
    if (mode == NcOpeningMode::Parallel)
    {
        ncid_ = NcOpenParallel(path_to_file.c_str(), comm);
    } else

    #endif /* CMC_WITH_NETCDF_PAR */
    {
        if (mode == NcOpeningMode::Parallel)
        {
            cmc_warn_msg("The netCDF file is ought to be opened for parallel access, although the parallel functionality is not accessible.");
            cmc_warn_msg("The file is opened for serial access.");
        }
        ncid_ = NcOpenSerial(path_to_file.c_str());
    }

    #else
    cmc_err_msg("CMC is not linked against netCDF. Please recompile the build with netCDF\n.");
    return CMC_ERR;
    #endif
}

void
NcData::SetHintLongitudeDimension(const int longitude_dimension_id)
{
    #ifdef CMC_WITH_NETCDF
    cmc_assert(longitude_dimension_id >= 0);
    coordinate_dimension_ids_[Dimension::Lon] = longitude_dimension_id;
    #endif
}

void
NcData::SetHintLatitudeDimension(const int latitude_dimension_id)
{
    #ifdef CMC_WITH_NETCDF
    cmc_assert(latitude_dimension_id >= 0);
    coordinate_dimension_ids_[Dimension::Lat] = latitude_dimension_id;
    #endif
}

void
NcData::SetHintHeightDimension(const int height_dimension_id)
{
    #ifdef CMC_WITH_NETCDF
    cmc_assert(height_dimension_id >= 0);
    coordinate_dimension_ids_[Dimension::Lev] = height_dimension_id;
    #endif
}

void
NcData::SetHintTimeDimension(const int time_dimension_id)
{
    #ifdef CMC_WITH_NETCDF
    cmc_assert(time_dimension_id >= 0);
    coordinate_dimension_ids_[Dimension::Time] = time_dimension_id;
    #endif
}

/** Inquire the dimension lenghts and ids of the (geo-spatial) coordinate dimensions as well as their variables ids */ 
void
NcData::InquireCoordinateDimensions()
{
    #ifdef CMC_WITH_NETCDF
    std::basic_string<char> nc_name;
    nc_name.reserve(NC_MAX_NAME +1);

    size_t dim_length;

    /* Reserve memory for all 'dimension_lengths_' and 'dimension_names_' */
    dimension_lengths_.reserve(num_dimensions_);
    dimension_names_.reserve(num_dimensions_);

    const int& number_dimensions = num_dimensions_;

    /* Loop over all dimension within the netCDF file */
    for (int dim_index = 0; dim_index < number_dimensions; ++dim_index) 
    {
        /* Get the name and the length of the dimension which corresponds to the current id */
        int err = nc_inq_dim(ncid_, dim_index, &nc_name[0], &dim_length);
        NcCheckError(err);

        /* Check if a hint to the dimension id was given (concerning the coordinate dimensions) */
        if (coordinate_dimension_ids_[Dimension::Lon] == dim_index)
        {
            coordinate_lengths_[Dimension::Lon] = static_cast<DomainIndex>(dim_length);
        } else if (coordinate_dimension_ids_[Dimension::Lat] == dim_index)
        {
            coordinate_lengths_[Dimension::Lat] = static_cast<DomainIndex>(dim_length);   
        } else if (coordinate_dimension_ids_[Dimension::Lev] == dim_index)
        {
            coordinate_lengths_[Dimension::Lev] = static_cast<DomainIndex>(dim_length);
        } else if (coordinate_dimension_ids_[Dimension::Time] == dim_index)
        {
            coordinate_lengths_[Dimension::Time] = static_cast<DomainIndex>(dim_length);
        }
        /* In case no hints concerning the coordinate dimensions were given, we test for some default names */
        /* Check if the name is equal to any of the considered (geo-spatial) coordinate variables */
        else if (coordinate_dimension_ids_[Dimension::Lon] == kDimensionNotConsiderdered && 
                   (strcmp(nc_name.c_str(), "lon") == 0 || strcmp(nc_name.c_str(), "longitude") == 0))
        {
            coordinate_dimension_ids_[Dimension::Lon] = dim_index;
            coordinate_lengths_[Dimension::Lon] = static_cast<DomainIndex>(dim_length);
        } else if (coordinate_dimension_ids_[Dimension::Lat] == kDimensionNotConsiderdered &&
                   (strcmp(nc_name.c_str(), "lat") == 0 || strcmp(nc_name.c_str(), "latitude") == 0))
        {
            coordinate_dimension_ids_[Dimension::Lat] = dim_index;
            coordinate_lengths_[Dimension::Lat] = static_cast<DomainIndex>(dim_length);
        }
        else if (coordinate_dimension_ids_[Dimension::Lev] == kDimensionNotConsiderdered &&
                 (strcmp(nc_name.c_str(), "lev") == 0 || strcmp(nc_name.c_str(), "level") == 0 || strcmp(nc_name.c_str(), "height") == 0 || strcmp(nc_name.c_str(), "elevation") == 0))
        {
            coordinate_dimension_ids_[Dimension::Lev] = dim_index;
            coordinate_lengths_[Dimension::Lev] = static_cast<DomainIndex>(dim_length);
        } else if (coordinate_dimension_ids_[Dimension::Time] == kDimensionNotConsiderdered && strcmp(nc_name.c_str(), "time") == 0)
        {
            coordinate_dimension_ids_[Dimension::Time] = dim_index;
            coordinate_lengths_[Dimension::Time] = static_cast<DomainIndex>(dim_length);
        }

        /* Save the name and the length of each dimension within the netCDF file */
        dimension_names_.emplace_back(nc_name.c_str());
        dimension_lengths_.push_back(static_cast<DomainIndex>(dim_length));
    }

    cmc_debug_msg("The following geo-spatial coordinate dimensions were retrieved from the netCDF file:");
    cmc_debug_msg("LON: ", coordinate_lengths_[Dimension::Lon], ", Lat: ", coordinate_lengths_[Dimension::Lat], ", LEV: ", coordinate_lengths_[Dimension::Lev], ", T: ", coordinate_lengths_[Dimension::Time]);
    
    #endif
}

#if 0
//Currently, there is nothing done with cooridnate variables.
//Only the data dimensions are interesting for the compression
//TODO: Add is anytime later
void
InquireCoordinateVariables()
{
    #ifdef CMC_WITH_NETCDF
    int err, var_type;

    /* Allocate a struct for the global coordinate system information */
    nc_data.coordinates = new cmc_global_coordinate_system();

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
            nc_data.coordinates->coords.create_and_push_back(static_cast<size_t>(nc_data.coord_lengths[coord_ids]), static_cast<cmc_type>(var_type));
            cmc_debug_msg("Information about ", get_coord_name(static_cast<CMC_COORD_IDS>(coord_ids)), " was inquired.");
        } else
        {
            /* Push back a placeholder in order to retain the internal used coordinate order */
            nc_data.coordinates->coords.push_back_placeholder();
        }
    }
    #endif
}

/** Store the geo-spatial coordinate data (latitude, longitude, leverage, time) */
static void
InquireCoordinateData(cmc_nc_data& nc_data)
{
    #ifdef CMC_WITH_NETCDF
    int err{0};

    /* Inquire the data of the latitude variable */
    if (nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_LAT] != CMC_COORDINATE_NOT_CONSIDERED)
    {
        err = nc_get_var(nc_data.get_ncid(), nc_data.coord_var_ids[CMC_COORD_IDS::CMC_LAT], nc_data.coordinates->coords[CMC_COORD_IDS::CMC_LAT].get_initial_data_ptr());
        cmc_nc_check_err(err);
        cmc_debug_msg("Coordinate data of variable ", get_coord_name(CMC_LAT), " has been inquired.");
    }

    /* Inquire the data of the longitude variable */
    if (nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_LON] != CMC_COORDINATE_NOT_CONSIDERED)
    {
        err = nc_get_var(nc_data.get_ncid(), nc_data.coord_var_ids[CMC_COORD_IDS::CMC_LON], nc_data.coordinates->coords[CMC_COORD_IDS::CMC_LON].get_initial_data_ptr());
        cmc_nc_check_err(err);
        cmc_debug_msg("Coordinate data of variable ", get_coord_name(CMC_LON), " has been inquired.");
    }

    /* Inquire the data of the leverage variable */
    if (nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_LEV] != CMC_COORDINATE_NOT_CONSIDERED)
    {
        err = nc_get_var(nc_data.get_ncid(), nc_data.coord_var_ids[CMC_COORD_IDS::CMC_LEV], nc_data.coordinates->coords[CMC_COORD_IDS::CMC_LEV].get_initial_data_ptr());
        cmc_nc_check_err(err);
        cmc_debug_msg("Coordinate data of variable ", get_coord_name(CMC_LEV), " has been inquired.");
    }

    /* Inquire the data of the time variable */
    if (nc_data.coord_dim_ids[CMC_COORD_IDS::CMC_TIME] != CMC_COORDINATE_NOT_CONSIDERED)
    {
        err = nc_get_var(nc_data.get_ncid(), nc_data.coord_var_ids[CMC_COORD_IDS::CMC_TIME], nc_data.coordinates->coords[CMC_COORD_IDS::CMC_TIME].get_initial_data_ptr());
        cmc_nc_check_err(err);
        cmc_debug_msg("Coordinate data of variable ", get_coord_name(CMC_TIME), " has been inquired.");
    }
    #endif
}
#endif

void
NcData::InquireCoordinates()
{
    #ifdef CMC_WITH_NETCDF
    /* Inquire the most general information about the supplied netCDF file (like number of diemnsion, number of global attributes, number of variables, etc.)*/
    int err = nc_inq(ncid_, &num_dimensions_, NULL, &num_global_attributes_, &id_unlimited_dimension_);
    NcCheckError(err);

    //TODO: When the coordinates are read in. Thex need to be adjusted to the given hyperslab later
    InquireCoordinateDimensions();
    //InquireCoordinateVariables();
    //InquireCoordinateData();

    #endif
}


static
int 
GetVariableId(const int ncid, const std::string& variable_name)
{
    int variable_id;

    /* Inquire the netCDF internal variable id by the given variable's name */
    const int err = nc_inq_varid(ncid, variable_name.c_str(), &variable_id);
    NcCheckError(err);

    return variable_id;
}

static
std::pair<bool, CmcUniversalType>
CheckVariableForAttribute(const int ncid, const int variable_id, const std::string& attribute_name)
{
    int attribute_type;

    int err = nc_inq_atttype(ncid, variable_id, attribute_name.c_str(), &attribute_type);

    CmcUniversalType attribute_value;

    if (err != NC_NOERR)
    {
        /* The specific attribute has not been found */
        return std::make_pair(false, attribute_value);
    }

    cmc_debug_msg("The attribute ", attribute_name, " has been found.");
    
    /* Otherwise, if the element has been found, we check it's type and inquire the attribute */
    switch (attribute_type)
    {
        case NC_DOUBLE:
            double double_att;
            err = nc_get_att_double(ncid, variable_id, attribute_name.c_str(), &double_att);
            NcCheckError(err);
            attribute_value = double_att;
            break;
        case NC_INT:
            int32_t int_att;
            err = nc_get_att_int(ncid, variable_id, attribute_name.c_str(), &int_att);
            NcCheckError(err);
            attribute_value = int_att;
            break;
        case NC_FLOAT:
            float float_att;
            err = nc_get_att_float(ncid, variable_id, attribute_name.c_str(), &float_att);
            NcCheckError(err);
            attribute_value = float_att;
            break;
        case NC_SHORT:
            int16_t short_att;
            err = nc_get_att_short(ncid, variable_id, attribute_name.c_str(), &short_att);
            NcCheckError(err);
            attribute_value = short_att;
            break;
        case NC_USHORT:
            uint16_t ushort_att;
            err = nc_get_att_ushort(ncid, variable_id, attribute_name.c_str(), &ushort_att);
            NcCheckError(err);
            attribute_value = ushort_att;
            break;
        case NC_UINT:
            uint32_t uint_att;
            err = nc_get_att_uint(ncid, variable_id, attribute_name.c_str(), &uint_att);
            NcCheckError(err);
            attribute_value = uint_att;
            break;
        case NC_INT64:
            long long int longlong_att;
            err = nc_get_att_longlong(ncid, variable_id, attribute_name.c_str(), &longlong_att);
            NcCheckError(err);
            attribute_value = static_cast<int64_t>(longlong_att);
            break;
        case NC_UINT64:
            unsigned long long int ulonglong_att;
            err = nc_get_att_ulonglong(ncid, variable_id, attribute_name.c_str(), &ulonglong_att);
            NcCheckError(err);
            attribute_value = static_cast<uint64_t>(ulonglong_att);
            break;
        default:
            cmc_debug_msg("The attribute ", attribute_name, " of variable with ID: ", variable_id, " could not be inquired (and will not be applied), since its data type is not supported.");
    }

    return std::make_pair(true, std::move(attribute_value));
}

static
DataLayout
GetLayoutFromDimensionOrdering(const std::vector<Dimension>& dimension_ordering)
{
    for (auto iter = dimension_ordering.begin(); iter != dimension_ordering.end(); ++iter)
    {
        cmc_debug_msg("dim ordering hat: ", *iter);
    }
    if (dimension_ordering.size() == 2)
    {
        /* Two dimensional layout */
        if (dimension_ordering[0] == Dimension::Lon && dimension_ordering[1] == Dimension::Lat)
            return DataLayout::Lon_Lat;
        else if (dimension_ordering[0] == Dimension::Lon && dimension_ordering[1] == Dimension::Lev)
            return DataLayout::Lon_Lev;
        else if (dimension_ordering[0] == Dimension::Lat && dimension_ordering[1] == Dimension::Lon)
            return DataLayout::Lat_Lon;
        else if (dimension_ordering[0] == Dimension::Lat && dimension_ordering[1] == Dimension::Lev)
            return DataLayout::Lat_Lev;
        else if (dimension_ordering[0] == Dimension::Lev && dimension_ordering[1] == Dimension::Lon)
            return DataLayout::Lev_Lon;
        else if (dimension_ordering[0] == Dimension::Lev && dimension_ordering[1] == Dimension::Lat)
            return DataLayout::Lev_Lat;
        else {
            cmc_err_msg("The two-dimensional data layout could not be deduced from the dimension ordering.");
            return DataLayout::LayoutUndefined;            
        } 
    } else if (dimension_ordering.size() == 3)
    {
        /* Three dimensional layout */
        if (dimension_ordering[0] == Dimension::Lev && dimension_ordering[1] == Dimension::Lon && dimension_ordering[2] == Dimension::Lat)
            return DataLayout::Lev_Lon_Lat;
        else if (dimension_ordering[0] == Dimension::Lev && dimension_ordering[1] == Dimension::Lat && dimension_ordering[2] == Dimension::Lon)
            return DataLayout::Lev_Lat_Lon;
        else if (dimension_ordering[0] == Dimension::Lat && dimension_ordering[1] == Dimension::Lev && dimension_ordering[2] == Dimension::Lon)
            return DataLayout::Lat_Lev_Lon;
        else if (dimension_ordering[0] == Dimension::Lat && dimension_ordering[1] == Dimension::Lon && dimension_ordering[2] == Dimension::Lev)
            return DataLayout::Lat_Lon_Lev;
        else if (dimension_ordering[0] == Dimension::Lon && dimension_ordering[1] == Dimension::Lat && dimension_ordering[2] == Dimension::Lev)
            return DataLayout::Lon_Lat_Lev;
        else if (dimension_ordering[0] == Dimension::Lon && dimension_ordering[1] == Dimension::Lev && dimension_ordering[2] == Dimension::Lat)
            return DataLayout::Lon_Lev_Lat;
        else {
            cmc_err_msg("The three-dimensional data layout could not be deduced from the dimension ordering.");
            return DataLayout::LayoutUndefined; 
        } 
    } else
    {
        cmc_err_msg("Only 2D and 3D variables are currently supported, therefore the DataLayout could not be deduced from the dimension ordering.");
        return DataLayout::LayoutUndefined;
    }
}

static
DataLayout
GetDataLayout(const int dimensionality, const Hyperslab& global_hyperslab, const CoordinateArray<int>& coordinate_dimension_ids_, const std::array<int, NC_MAX_VAR_DIMS>& variable_dimension_ids)
{
    std::vector<Dimension> data_dimensions;
    data_dimensions.reserve(dimensionality);

    int found_dimensions = 0;

    for (auto dim_iter = variable_dimension_ids.begin(); found_dimensions < dimensionality && dim_iter != variable_dimension_ids.end(); ++dim_iter)
    {
        cmc_debug_msg("dim iter hat: ", *dim_iter);
        if (const auto coordinate_dim_id_iter = std::find_if(coordinate_dimension_ids_.begin(), coordinate_dimension_ids_.end(), [&](const int& arg){return arg == *dim_iter;});
            coordinate_dim_id_iter != coordinate_dimension_ids_.end())
        {
            /* If the current dimension id corresponds to a coordinate dimension id (which is considered), we add the dimension to the vector */
            const Dimension dimension = static_cast<Dimension>(std::distance(coordinate_dimension_ids_.begin(), coordinate_dim_id_iter));
            if (global_hyperslab.GetDimensionLength(dimension) > 1)
            {
                data_dimensions.push_back(dimension);
                ++found_dimensions;
            }
        }
    }

    if (found_dimensions < dimensionality)
    {
        cmc_err_msg("The DataLayout could not be deduced. The variable is not only dependent on geo-spatial coordinate dimensions");
    }

    return GetLayoutFromDimensionOrdering(data_dimensions);
}


static
std::vector<Hyperslab>
DetermineDataDistribution(const Hyperslab& global_hyperslab)
{
    //TODO:: Update for the parallel case
    return std::vector<Hyperslab>{global_hyperslab};
}

static
std::pair<std::vector<size_t>, std::vector<size_t>>
GetStartAndCountValuesForVariable(const Hyperslab& hyperslab, const CoordinateArray<int>& coordinate_dimension_ids, int num_dimensions, const std::array<int, NC_MAX_VAR_DIMS>& dimension_ids)
{
    #if 0
    std::array<size_t, Dimension::NumCoordinates> start_vals{0,0,0,0};
    std::array<size_t, Dimension::NumCoordinates> count_vals{1,1,1,1};

    int array_accessor = 0;
    for (auto do_iter = dimension_ordering.begin(); do_iter != dimension_ordering.end(); ++do_iter, ++array_accessor)
    {
        start_vals[array_accessor] = static_cast<size_t>(hyperslab.GetDimensionStart(*do_iter));
        count_vals[array_accessor] = static_cast<size_t>(hyperslab.GetDimensionLength(*do_iter));
    }
    #endif

    std::vector<size_t> start_vals(num_dimensions, 0);
    std::vector<size_t> count_vals(num_dimensions, 1);

    for (int iter = 0; iter < num_dimensions; ++iter)
    {
        if (const auto coordinate_dim_id_iter = std::find_if(coordinate_dimension_ids.begin(), coordinate_dimension_ids.end(), [&](const int& arg){return arg == dimension_ids[iter];});
            coordinate_dim_id_iter != coordinate_dimension_ids.end())
        {
            const Dimension dimension = static_cast<Dimension>(std::distance(coordinate_dimension_ids.begin(), coordinate_dim_id_iter));
            start_vals[iter] = static_cast<size_t>(hyperslab.GetDimensionStart(dimension));
            count_vals[iter] = static_cast<size_t>(hyperslab.GetDimensionLength(dimension));
        }
        
    }

    return std::make_pair(std::move(start_vals), std::move(count_vals));
}

template<typename T>
static InputVar
SetUpInputVariable(const int ncid, const CoordinateArray<int>& coordinate_dimension_ids, const int variable_id, std::string&& variable_name, const DataLayout layout,
                   const DomainIndex num_values, std::vector<Hyperslab>&& hyperslabs, GeoDomain&& global_domain,
                   const int num_dimensions, const std::array<int, NC_MAX_VAR_DIMS>& dimension_ids)
{
    InputVariable<T> variable(std::move(variable_name), variable_id, layout);

    variable.SetGlobalDomain(std::move(global_domain));

    std::vector<T> local_data(num_values);
    DomainIndex offset = 0;
    for (auto hs_iter = hyperslabs.begin(); hs_iter != hyperslabs.end(); ++hs_iter)
    {
        //const std::vector<Dimension> dimension_ordering = GetDimensionVectorFromLayout(layout);
        const auto [start_vals, count_vals] = GetStartAndCountValuesForVariable(*hs_iter, coordinate_dimension_ids, num_dimensions, dimension_ids);

        const int err = nc_get_vara(ncid, variable_id, start_vals.data(), count_vals.data(), static_cast<void*>(&local_data[offset]));
        NcCheckError(err);

        offset += hs_iter->GetNumberCoordinates();
    }

    //std::vector<float> converted_data;
    //const int num_elems = 3600*1801;
    //converted_data.reserve(num_elems);
    //const float add_offset = 0.0970683052537568;
    //const float scale_factor = 2.96247040388686e-06;
    //for (int j = 0; j < num_elems; ++j)
    //{
    //    converted_data.push_back(scale_factor * static_cast<float>(local_data[j]) + add_offset);
    //}
    //FILE* file = fopen("era5_land_tp_lon_lat.bin", "wb");
    //fwrite(converted_data.data(), sizeof(float), num_elems, file);
    //fclose(file);
    //std::exit(1);
    
    variable.SetDataAndCoordinates(std::move(local_data), std::move(hyperslabs));

    return InputVar(std::move(variable));
}

InputVar
NcData::SetupVariableData(const int variable_nc_type, const int variable_id, std::string&& variable_name,
                          const DataLayout layout, const DomainIndex num_values, std::vector<Hyperslab>&& hyperslabs,
                          GeoDomain&& global_domain, const int num_dimensions, const std::array<int, NC_MAX_VAR_DIMS>& dimension_ids)
{
    switch (variable_nc_type)
    {
        case NC_DOUBLE:
            return SetUpInputVariable<double>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_FLOAT:
            return SetUpInputVariable<float>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_INT:
            return SetUpInputVariable<int32_t>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_UINT:
            return SetUpInputVariable<uint32_t>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_INT64:
            return SetUpInputVariable<int64_t>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_UINT64:
            return SetUpInputVariable<uint64_t>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_SHORT:
            return SetUpInputVariable<int16_t>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_USHORT:
            return SetUpInputVariable<uint16_t>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_BYTE:
            return SetUpInputVariable<int8_t>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_UBYTE:
            return SetUpInputVariable<uint8_t>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        case NC_CHAR:
            return SetUpInputVariable<char>(ncid_, coordinate_dimension_ids_, variable_id, std::move(variable_name), layout, num_values, std::move(hyperslabs), std::move(global_domain), num_dimensions, dimension_ids);
        break;
        default:
            cmc_err_msg("The variable of type (NC-Type: ", variable_nc_type, ") is not supported.");
            return InputVar();
    }
}

InputVar
NcData::InquireVariable(const Hyperslab& hyperslab, std::string&& variable_name)
{
    const int variable_id = GetVariableId(ncid_, variable_name);

    int variable_nc_type;
    int num_dimensions;
    std::array<int, NC_MAX_VAR_DIMS> dimension_ids;

    /* Inquire some information about the variable, e.g. the data type and its dimensions */
    int err = nc_inq_var(ncid_, variable_id, NULL, &variable_nc_type, &num_dimensions, dimension_ids.data(), NULL);
    NcCheckError(err);
    cmc_debug_msg("Num dims: ", num_dimensions, " IDs: ");
    for (auto iii = 0; iii < num_dimensions; ++iii)
    {
        cmc_debug_msg(dimension_ids[iii], ", ");
    }
    /* Check some attributes of the variable */
    const auto [is_missing_value_present, missing_value] = CheckVariableForAttribute(ncid_, variable_id, "missing_value");
    const auto [is_add_offset_present, add_offset] = CheckVariableForAttribute(ncid_, variable_id, "add_offset");
    const auto [is_scale_factor_present, scale_factor] = CheckVariableForAttribute(ncid_, variable_id, "scale_factor");

    const int data_dimensionality = hyperslab.GetDimensionality();

    /* Determine the data layout of the variable */
    const DataLayout data_layout = GetDataLayout(data_dimensionality, hyperslab, coordinate_dimension_ids_, dimension_ids);

    /* TODO: Add parallelization */
    //for now in serial only

    std::vector<Hyperslab> local_hyperslabs = DetermineDataDistribution(hyperslab);

    /* Get the number of values we will read from the file (for memory allocation) */
    DomainIndex num_local_data_values = 0;
    for (auto hs_iter = local_hyperslabs.begin(); hs_iter != local_hyperslabs.end(); ++hs_iter)
    {
        num_local_data_values += hs_iter->GetNumberCoordinates();
    }
    
    cmc_debug_msg("Num lcoal elems: ", num_local_data_values);

    /* Inquire the data of the variable and construct a InputVariable with it */
    InputVar variable = SetupVariableData(variable_nc_type, variable_id, std::move(variable_name), data_layout, num_local_data_values, std::move(local_hyperslabs), TransformHyperslabToGeoDomain(hyperslab), num_dimensions, dimension_ids);

    if (is_missing_value_present)
        variable.SetMissingValue(missing_value);
    
    if (is_add_offset_present)
        variable.SetAddOffset(add_offset);
    
    if (is_scale_factor_present)
        variable.SetScaleFactor(scale_factor);

    return variable;
}

[[nodiscard]]
std::vector<InputVar>&&
NcData::TransferData()
{
    if (_data_has_been_transfered_)
        cmc_err_msg("The data already has been transfered.");
    
    _data_has_been_transfered_ = true;
    
    return std::move(variables_);
}

}
