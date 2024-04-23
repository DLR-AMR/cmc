#include "netcdf/cmc_netcdf.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_geo_utilities.hxx"
#include "utilities/cmc_log_functions.h"

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

static
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

static
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
    
    return std::move(variables_);
}

}

#if 0

/*****************************/
/*****************************/
/***** OLD CODE **************/

#include "cmc_netcdf.h"
#include "utilities/cmc_container.h"
#include "utilities/cmc_geo_util.h"
#include "utilities/cmc_log_functions.h"
#include "mpi/cmc_mpi_io.h"

#define CMC_WITH_NETCDF 1
#define CMC_WITH_T8CODE 1

enum CMC_NC_STATUS {STATUS_UNDEFINED = 0, CMC_NC_INQ_COORDS, CMC_NC_ADD_VAR, CMC_NC_INQ_VAR, CMC_NC_READY_TO_COMPRESS};

struct cmc_nc_data
{
private:
    const int ncid;
public:
    
    cmc_nc_data(const int _ncid)
    : ncid{_ncid}{};
    cmc_nc_data(const int _ncid, const bool parallel_access)
    : ncid{_ncid}, use_distributed_data{parallel_access}{};
    ~cmc_nc_data(){
        if (coordinates != nullptr)
        {
            delete coordinates;
        }
    };
    /* These variables only save cooridnate relative values for (latitude, longitude, leverage, time) */
    std::array<int, CMC_NUM_COORD_IDS> coord_dim_ids{CMC_COORDINATE_NOT_CONSIDERED, CMC_COORDINATE_NOT_CONSIDERED, CMC_COORDINATE_NOT_CONSIDERED, CMC_COORDINATE_NOT_CONSIDERED};
    std::array<int, CMC_NUM_COORD_IDS> coord_var_ids{CMC_COORDINATE_NOT_CONSIDERED, CMC_COORDINATE_NOT_CONSIDERED, CMC_COORDINATE_NOT_CONSIDERED, CMC_COORDINATE_NOT_CONSIDERED};
    std::array<size_t, CMC_NUM_COORD_IDS> coord_lengths{0, 0, 0, 0};

    /* Number of dimensions of the netCDF file */
    int num_dimensions{0};
    int num_global_atts{0}; //might be useful later
    int id_unlimited_dim{-1}; //might be useful later for time series data

    /* Information about the global coordinate system */
    cmc_global_coordinate_system_t coordinates{nullptr};

    /* Saves correponding data to each inquired data variable */
    std::vector<cmc_var_t> vars;

    /* These variables save all dimension names and sizes */
    std::vector<size_t> dimension_sizes{};
    std::vector<std::string> dimension_names;

    /* A status flag of the current mode the struct is in */
    CMC_NC_STATUS status{STATUS_UNDEFINED};

    MPI_Comm comm{MPI_COMM_WORLD}; //!< The communicator to use in a parallel environment

    data_distribution_t data_distribution{DISTRIBUTION_UNDEFINED};
    std::vector<int> distribution_offsets;

    bool use_distributed_data{false};

    int get_ncid() const;
};

cmc_nc_data_t
cmc_nc_start(const char* path_to_file, const enum cmc_nc_opening_mode mode, const MPI_Comm comm)
{
    #ifdef CMC_WITH_NETCDF
    int err{0};
    int ncid;
    bool parallel_access{false};

    switch(mode)
    {
        case cmc_nc_opening_mode::CMC_NC_SERIAL:
            /* Open the file in a serial mode */
            ncid = cmc_nc_open_serial(path_to_file);
        break;
        case cmc_nc_opening_mode::CMC_NC_PARALLEL:
        {
            #ifdef CMC_WITH_NETCDF_PAR
            int comm_size{0};
            /* Get the size of the supplied communicator */
            err = MPI_Comm_size(comm, &comm_size);
            cmc_mpi_check_err(err);

            /* Check whether there are several processs in the communicator for actually reading in parallel from the file */
            if (comm_size > 1)
            {
                /* Open the file for parallel access */
                ncid = cmc_nc_open_parallel(path_to_file, comm);

                /* Set the flag for parallel access */
                parallel_access = true;
            } else
            {
                cmc_warn_msg("The file is ought to be opened in a parallel mode but there is not more than one process within the supplied communicator. Therfore, the fll be opened in a serial mode.");

                /* Open the file in a serial mode */
                ncid = cmc_nc_open_serial(path_to_file);
            }
            #else
            cmc_err_msg("NetCDF's parallel file access functions are not available. Please ensure that 'netcdf_par.h' is present when cmc is linked against netCDF.");
            #endif
        }
        break;
        default:
            cmc_err_msg("An unknown netCDF opening mode was supplied.");
    }

    /* Create a netCDF data struct and return it */
    return new cmc_nc_data{ncid, parallel_access};

    #else
    cmc_err_msg("CMC is not linked against netCDF. Please recompile the build with netCDF\n.");
    return CMC_ERR;
    #endif
}

void
cmc_nc_finish(cmc_nc_data_t nc_data)
{
    #ifdef CMC_WITH_NETCDF
    /* Close the netCDF File */
    cmc_nc_close(nc_data->get_ncid());

    /* Deallocate the nc_data */
    if (nc_data != nullptr)
    {
        delete nc_data;
    }
    #else
    cmc_err_msg("CMC is not linked against netCDF. Please recompile the build with netCDF\n.");
    #endif
}
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
    /* Save the supplied MPI communicator */
    nc_data->comm = comm;
}

/* Set a preferred (blockwise) paralled distribution for reading the data */
void
cmc_nc_set_blocked_reading(cmc_nc_data_t nc_data, const std::vector<int> blocked_domain_num_processes_per_dimension)
{
    /* Save the preferred reading distribution */
    nc_data->data_distribution = data_distribution_t::CMC_BLOCKED;
    /* Save the supplied offsets per dimension */
    nc_data->distribution_offsets = blocked_domain_num_processes_per_dimension;
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
cmc_nc_inquire_coordinate_data(cmc_nc_data& nc_data)
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

static
void
cmc_nc_update_global_coords(cmc_nc_data_t nc_data, const size_t* start_ptr, const size_t* count_ptr)
{
    #ifdef CMC_WITH_NETCDF

    /* The variables already have been pre-allocated and we take the axis ordering of a variable with the highest dimension */
    int num_dims = 0;
    int id_of_highest_var_ids = 0;
    int var_id = 0;

    /* Find the variable with the most dimensions */
    for (auto iter = nc_data->vars.begin(); iter != nc_data->vars.end(); ++iter, ++var_id)
    {
        int num_considered_dims = 0;
        for (int dims{0}; dims < (*iter)->num_dimensions; ++dims)
        {
            /* Check if the dimension is considered, and if so increment the dimensionality counter */
            if (count_ptr[dims] > 1)
            {
                ++num_considered_dims;
            }
        }
        /* If the considered dimensions are higher than the previous one, we save the amount and the id of the variable */
        if (num_considered_dims > num_dims)
        {
            num_dims = num_considered_dims;
            id_of_highest_var_ids = var_id;
        }
    }

    /* Eventually adjust the global coordinates */
    for (int dims{0}; dims < nc_data->vars[id_of_highest_var_ids]->num_dimensions; ++dims)
    {
        const int current_dim_id = nc_data->vars[id_of_highest_var_ids]->dimension_ids[dims];

        if (current_dim_id == nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_LAT])
        {
            if (count_ptr[dims] < nc_data->coordinates->coords.operator[](static_cast<size_t>(CMC_COORD_IDS::CMC_LAT)).size())
            {
                /* Only a part of the coordinate dimension is considered */
                nc_data->coordinates->coords.operator[](static_cast<size_t>(CMC_COORD_IDS::CMC_LAT)).crop_to(static_cast<size_t>(start_ptr[dims]), static_cast<size_t>(start_ptr[dims] + count_ptr[dims] - 1));
            }
        } else if (current_dim_id == nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_LON])
        {
            if (count_ptr[dims] < nc_data->coordinates->coords.operator[](static_cast<size_t>(CMC_COORD_IDS::CMC_LON)).size())
            {
                /* Only a part of the coordinate dimension is considered */
                nc_data->coordinates->coords.operator[](static_cast<size_t>(CMC_COORD_IDS::CMC_LON)).crop_to(static_cast<size_t>(start_ptr[dims]), static_cast<size_t>(start_ptr[dims] + count_ptr[dims] - 1));
            }
        } else if (current_dim_id == nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_LEV])
        {
            if (count_ptr[dims] < nc_data->coordinates->coords.operator[](static_cast<size_t>(CMC_COORD_IDS::CMC_LEV)).size())
            {
                /* Only a part of the coordinate dimension is considered */
                nc_data->coordinates->coords.operator[](static_cast<size_t>(CMC_COORD_IDS::CMC_LEV)).crop_to(static_cast<size_t>(start_ptr[dims]), static_cast<size_t>(start_ptr[dims] + count_ptr[dims] - 1));
            }
        } else if (current_dim_id == nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_TIME])
        {
            //Time-Series data is currently not supported
            continue;
        } else
        {
            continue;
        }
    }
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

    /* Read the meta data of the variable (for example 'data_type', 'missing_value', 'offset', 'scale_factor', ...) */
    cmc_nc_inquire_var_meta_data(nc_data);

    /* Beforehand, the whole coordinate dimensions have been read, but it is possible that only a certain view is considered. 
     * Therefore, we have to update the global coordinate system based on the start and count */
    cmc_nc_update_global_coords(nc_data, start_ptr, count_ptr);

    /* Preallocate data arrays */
    for (size_t i{0}; i < nc_data->vars.size(); ++i)
    {
        cmc_assert(nc_data->vars[i]->var_id != CMC_NC_VAR_NOT_CONSIDERED);

        /* Get the axis ordering of the data */
        nc_data->vars[i]->axis_ordering.reserve(nc_data->vars[i]->num_dimensions);

        /* Retrieve the axis ordering of the variable */
        for (int j{0}; j < nc_data->vars[i]->num_dimensions; ++j)
        {
            if (nc_data->vars[i]->dimension_ids[j] == nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_LON])
            {
                if (count_ptr[j] > 1)
                {
                    nc_data->vars[i]->axis_ordering.emplace_back(CMC_COORD_IDS::CMC_LON);
                    nc_data->vars[i]->dimension_considered[CMC_COORD_IDS::CMC_LON] = true;
                }
            } else if (nc_data->vars[i]->dimension_ids[j] == nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_LAT])
            {
                if (count_ptr[j] > 1)
                {
                    nc_data->vars[i]->axis_ordering.emplace_back(CMC_COORD_IDS::CMC_LAT);
                    nc_data->vars[i]->dimension_considered[CMC_COORD_IDS::CMC_LAT] = true;
                }
            } else if (nc_data->vars[i]->dimension_ids[j] == nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_LEV])
            {
                if (count_ptr[j] > 1)
                {
                    nc_data->vars[i]->axis_ordering.emplace_back(CMC_COORD_IDS::CMC_LEV);
                    nc_data->vars[i]->dimension_considered[CMC_COORD_IDS::CMC_LEV] = true;
                }
            } else if (nc_data->vars[i]->dimension_ids[j] == nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_TIME])
            {
                //Time series are currently not supported, therefore the time coordinates are skipped
                #if 0
                if (count_ptr[j] > 1)
                {
                    nc_data->vars[i]->axis_ordering.emplace_back(CMC_COORD_IDS::CMC_TIME);
                    nc_data->vars[i]->dimension_considered[CMC_COORD_IDS::CMC_TIME] = true;
                }
                #endif
            }
        }

        /* Reset the number of data points */
        num_data_points = 1;

        /* Reserve space for storing the lengths of the data of each data dimension */
        nc_data->vars[i]->dim_lengths.reserve(nc_data->vars[i]->num_dimensions);
        
        /* Reserve space for storing the start of the data of each data dimension */
        nc_data->vars[i]->dim_starts.reserve(nc_data->vars[i]->num_dimensions);

        /* Allocate memeory for the start and the count ptr */
        nc_data->vars[i]->start_ptr.reserve(nc_data->vars[i]->num_dimensions);
        nc_data->vars[i]->count_ptr.reserve(nc_data->vars[i]->num_dimensions);

        /* Save the global start and count pointers associated with the hyperslab */
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
            /* If the data of the whole varibale will be read, save the global start and count values */
            for (int dims{0}; dims < nc_data->vars[i]->num_dimensions; ++dims)
            {
                nc_data->vars[i]->start_ptr.push_back(0);
                nc_data->vars[i]->count_ptr.push_back(static_cast<uint32_t>(nc_data->dimension_sizes[nc_data->vars[i]->dimension_ids[dims]]));
            }
        }

        /* Check if the variables' data may be read in distributed */
        if (nc_data->use_distributed_data)
        {
            /* The data may be read in distributed */
            cmc_mpi_nc_data distributed_data;

            /* Calculate the offsets for each process based on the (preferred) data distribution */
            switch(nc_data->data_distribution)
            {
                case data_distribution_t::CMC_BLOCKED:
                    /* Calculate a blocked data dsitribution */
                    distributed_data = cmc_mpi_calculate_nc_reading_data_distribution_blocked(nc_data->vars[i]->start_ptr, nc_data->vars[i]->count_ptr, nc_data->comm, nc_data->distribution_offsets);
                break;
                case data_distribution_t::DISTRIBUTION_UNDEFINED:
                {
                    /* If no distribution has been specified, we default to a blocked distribution */
                    /* Size of the MPI communicator */
                    int comm_size{1};
                    #ifdef CMC_ENABLE_MPI
                    /* Get the actual size of the communicator */
                    int err = MPI_Comm_size(nc_data->comm, &comm_size);
                    cmc_mpi_check_err(err);
                    #endif

                    int num_considered_dims = 0;
                    /* Define a vector for a blocked distribution */
                    std::vector<int> p_distribution(nc_data->vars[i]->num_dimensions, 1);
                    /* Iterate over the count vector and inquire the general data dimensions (e.g. 2D or 3D) */
                    for (int dims{0}; dims < nc_data->vars[i]->num_dimensions; ++dims)
                    {
                        if (nc_data->vars[i]->count_ptr[dims] > 1)
                        {
                            ++num_considered_dims;
                        }
                    }

                    cmc_assert(num_considered_dims > 0);
                    
                    /* Calculate a n equal distribution per dimension */
                    const int equal_dim_distribution = static_cast<int>(std::pow(comm_size, 1.0/num_considered_dims));

                    /* Iterate again over the count vector in order to assign a blocked distribution at the correct dimensions */
                    for (int dims{0}; dims < nc_data->vars[i]->num_dimensions; ++dims)
                    {
                        if (nc_data->vars[i]->count_ptr[dims] > 1)
                        {
                            p_distribution[dims] = equal_dim_distribution;
                        }
                    }

                    /* Calculate a blocked data dsitribution */
                    distributed_data = cmc_mpi_calculate_nc_reading_data_distribution_blocked(nc_data->vars[i]->start_ptr, nc_data->vars[i]->count_ptr, nc_data->comm, p_distribution);
                    
                    /* Set that we have defaulted to a blcoekd scheme */
                    nc_data->data_distribution = data_distribution_t::CMC_BLOCKED;
                }
                break;
                default:
                    cmc_err_msg("An unknown data distribution for parallel reading was supplied.");
            }

            /* Update the start and count values for the variable to their local view */
            std::swap(nc_data->vars[i]->start_ptr, distributed_data.start_values);
            std::swap(nc_data->vars[i]->count_ptr, distributed_data.count_values);
        }

        /* Calculate the number of data points (neeeded for the array allocation) */
        for (int dims{0}; dims < nc_data->vars[i]->num_dimensions; ++dims)
        {
            /* Save the actual dimension length of the data */
            (nc_data->vars[i]->dim_lengths).push_back(nc_data->vars[i]->count_ptr[dims]);
            /* Save the dimension start */
            (nc_data->vars[i]->dim_starts).push_back(nc_data->vars[i]->start_ptr[dims]);
            /* Reduce the data points per dimension in order to obtain the overall amount of data points of the variable */
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
            start_values.push_back(static_cast<size_t>(nc_data->vars[i]->start_ptr[dim_id]));
            count_values.push_back(static_cast<size_t>(nc_data->vars[i]->count_ptr[dim_id]));
        }

        /* Differentiate between parallel and serial reading of the data */
        if (nc_data->use_distributed_data)
        {
            /** If the data is read in parallel **/
            /* Read in the process-local hyperslab of the data */
            err = nc_get_vara(nc_data->get_ncid(), nc_data->vars[i]->var_id, start_values.data(), count_values.data(), nc_data->vars[i]->data->get_initial_data_ptr());
            cmc_nc_check_err(err);
        }
        else
        {
            /** If the data is read in serial **/
            /* Distinguish between a serial read of a hyperslab or the serial data inquisition of the whole variable */
            if (start_ptr == nullptr && count_ptr == nullptr)
            {
                /* Store data of the whole variable if no hyperslab is defined */
                err = nc_get_var(nc_data->get_ncid(), nc_data->vars[i]->var_id, nc_data->vars[i]->data->get_initial_data_ptr());
                cmc_nc_check_err(err);
            }
            else
            {
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

    //TODO:  This has to be updated for parallel access as well

    /* Inquire the data of these variables */
    cmc_nc_inquire_var_data(nc_data, start_ptr, count_ptr);

    /* Set status flag */
    nc_data->status = CMC_NC_STATUS::CMC_NC_READY_TO_COMPRESS;

    #else
    cmc_err_msg("CMC is not compiled with netCDF, please reconfigure with netCDF linkage in order to use this function.");
    #endif
}

int
cmc_nc_open_serial(const char* path_to_file)
{
    #ifdef CMC_WITH_NETCDF
    int err{0};
    int ncid{0};

    /* Open the file without explicit parallel access */
    err = nc__open(path_to_file, NC_NOWRITE, NULL, &ncid);
    cmc_nc_check_err(err);

    return ncid;
    #else
    cmc_err_msg("CMC is not linked against netCDF. Please recompile the build with netCDF\n.");
    return CMC_ERR;
    #endif
}

int
cmc_nc_open_parallel(const char* path_to_file, MPI_Comm comm)
{
    #ifdef CMC_WITH_NETCDF_PAR
    int err{0};
    int ncid{0};

    int comm_size{0};
    err = MPI_Comm_size(comm, &comm_size);
    cmc_mpi_check_err(err);

    /* Is more than one process present */
    cmc_assert(comm_size > 1);

    MPI_Info info = MPI_INFO_NULL;
    err = nc_open_par(path_to_file, NC_NOWRITE, comm, info, &ncid);
    cmc_nc_check_err(err);


    return ncid;
    #else
    cmc_err_msg("CMC is not linked against netCDF (at leat no parallel access functions are available). Please recompile the build with netCDF\n.");
    return CMC_ERR;
    #endif
}

/** Open the netCDF file */
int
cmc_nc_open(const char* path_to_file, MPI_Comm comm)
{
    #ifdef CMC_WITH_NETCDF
    int err{0};

    #ifdef CMC_WITH_NETCDF_PAR
    int comm_size{0};
    err = MPI_Comm_size(comm, &comm_size);
    cmc_mpi_check_err(err);
    if (comm_size > 1)
    {
        return cmc_nc_open_parallel(path_to_file, comm);
    } else
    #endif
    {
        return cmc_nc_open_serial(path_to_file);
    }

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
    //t8_data->geo_data->coords = new var_vector_t(std::move(nc_data->coords));
    
    /* Save the global coordnate system information */
    t8_data->geo_data->coordinates = nc_data->coordinates;
    nc_data->coordinates = nullptr;

    /* Save the globale dimension lengths */
    t8_data->geo_data->global_dim_lengths.reserve(CMC_NUM_COORD_IDS);
    for (size_t i{0}; i < CMC_NUM_COORD_IDS; ++i)
    {
        t8_data->geo_data->global_dim_lengths.push_back(t8_data->geo_data->coordinates->coords.operator[](i).size());
    }

    //Time series are currently skipped
    t8_data->geo_data->global_dim_lengths[CMC_COORD_IDS::CMC_TIME] = 1;
    
    /* Reserve memory for all variables */
    t8_data->vars.reserve(nc_data->vars.size());

    /* Retrieve the axis ordering of all variables */
    /* While doing this, save the maximum dimensionality of the variables */
    int max_dim{0};
    int update_num_dims;
    /* Create a new dimension length vector for this varibale. This member is used differently in nc_functions than in t8_functions... */
    std::vector<size_t> cmc_t8_dim_lengths(CMC_NUM_COORD_IDS);
    /* Create a new dimension start vector for this variable. */
    std::vector<size_t> cmc_t8_dim_starts(CMC_NUM_COORD_IDS);

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

        /* 'Nullify' the new dim_lengths vector (fill it up with ones) */
        std::fill(cmc_t8_dim_lengths.begin(), cmc_t8_dim_lengths.end(), 0);

        /* 'Nullify' the new dim_lengths vector (fill it up with ones) */
        std::fill(cmc_t8_dim_starts.begin(), cmc_t8_dim_starts.end(), 0);

        #if 1
        /* Reserve space for the data axis ordering */
        t8_data->vars[var_id]->var->axis_ordering.reserve(t8_data->vars[var_id]->var->num_dimensions);

        /* Retrieve the axis ordering of the variable */
        for (int i{0}; i < t8_data->vars[var_id]->var->num_dimensions; ++i)
        {
            if (t8_data->vars[var_id]->var->dimension_ids[i] == nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_LON])
            {
                if (t8_data->vars[var_id]->var->dim_lengths[i] != 1 || t8_data->vars[var_id]->var->dimension_considered[CMC_COORD_IDS::CMC_LON])
                {
                    //t8_data->vars[var_id]->var->axis_ordering.emplace_back(CMC_COORD_IDS::CMC_LON);
                    cmc_t8_dim_lengths[CMC_COORD_IDS::CMC_LON] = t8_data->vars[var_id]->var->dim_lengths[i];
                    cmc_t8_dim_starts[CMC_COORD_IDS::CMC_LON] = t8_data->vars[var_id]->var->dim_starts[i];
                } else
                {
                    ++update_num_dims;
                }
            } else if (t8_data->vars[var_id]->var->dimension_ids[i] == nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_LAT])
            {
                if (t8_data->vars[var_id]->var->dim_lengths[i] != 1 || t8_data->vars[var_id]->var->dimension_considered[CMC_COORD_IDS::CMC_LAT])
                {
                    //t8_data->vars[var_id]->var->axis_ordering.emplace_back(CMC_COORD_IDS::CMC_LAT);
                    cmc_t8_dim_lengths[CMC_COORD_IDS::CMC_LAT] = t8_data->vars[var_id]->var->dim_lengths[i];
                    cmc_t8_dim_starts[CMC_COORD_IDS::CMC_LAT] = t8_data->vars[var_id]->var->dim_starts[i];
                } else
                {
                    ++update_num_dims;
                }
            } else if (t8_data->vars[var_id]->var->dimension_ids[i] == nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_LEV])
            {
                if (t8_data->vars[var_id]->var->dim_lengths[i] != 1 || t8_data->vars[var_id]->var->dimension_considered[CMC_COORD_IDS::CMC_LEV])
                {
                    //t8_data->vars[var_id]->var->axis_ordering.emplace_back(CMC_COORD_IDS::CMC_LEV);
                    cmc_t8_dim_lengths[CMC_COORD_IDS::CMC_LEV] = t8_data->vars[var_id]->var->dim_lengths[i];
                    cmc_t8_dim_starts[CMC_COORD_IDS::CMC_LEV] = t8_data->vars[var_id]->var->dim_starts[i];
                } else
                {
                    ++update_num_dims;
                }
            } else if (t8_data->vars[var_id]->var->dimension_ids[i] == nc_data->coord_dim_ids[CMC_COORD_IDS::CMC_TIME])
            {
                //Time series are currently not supported, therefore the time coordinates are skipped
                #if 0
                if (nc_data->coord_lengths[CMC_COORD_IDS::CMC_TIME] != 1 || t8_data->vars[var_id]->var->dimension_considered[CMC_COORD_IDS::CMC_TIME])
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
        t8_data->vars[var_id]->var->num_dimensions = t8_data->vars[var_id]->var->axis_ordering.size();
        #endif
        /* Save the maximum dimensionality */
        if (t8_data->vars[var_id]->var->num_dimensions > max_dim)
        {
            max_dim = t8_data->vars[var_id]->var->num_dimensions;
        }

        /* Update the start and count ptrs */
        size_t num_pts = 1;
        for (size_t j{0}; j < t8_data->vars[var_id]->var->count_ptr.size(); ++j)
        {
            /* Check if this dimension is considered in the data hyperslab */
            if (t8_data->vars[var_id]->var->count_ptr[j] == 1)
            {
                /* If the dimension is not considered, erase it */
                t8_data->vars[var_id]->var->count_ptr.erase(t8_data->vars[var_id]->var->count_ptr.begin() + j);
                t8_data->vars[var_id]->var->start_ptr.erase(t8_data->vars[var_id]->var->start_ptr.begin() + j);
                /* Decrement the counter j in order to potint to 'next' value (the one right after the erasure) */
                --j;
            } else
            {
                /* Multiply the considered dimension length */
                num_pts *= t8_data->vars[var_id]->var->count_ptr[j];
            }
        }

        /* Set the new dimension lenghts */
        t8_data->vars[var_id]->var->dim_lengths = cmc_t8_dim_lengths;

        /* Set the new dimension starts */
        t8_data->vars[var_id]->var->dim_starts = cmc_t8_dim_starts;

        /* Retrieve the data layout of the variable */
        t8_data->vars[var_id]->var->get_data_layout_from_axis_ordering();
    }

    /* Save the maximum geo-spatial dimension of the data */
    t8_data->geo_data->dim = max_dim;

    /* Save the flag whether the data is distributed among processes or not */
    t8_data->use_distributed_data = nc_data->use_distributed_data;

    /* Save distribution scheme of the data */
    t8_data->data_dist = nc_data->data_distribution;
    /* Set the information that the data source are netCDF files */
    t8_data->data_source = CMC_T8_DATA_INPUT::NETCDF_INPUT;
    
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

#endif