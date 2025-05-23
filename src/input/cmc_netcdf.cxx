#include "input/cmc_netcdf.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_geo_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "mpi/cmc_mpi_io.hxx"

namespace cmc::input::netcdf
{

void
Data::Open(const std::string& path_to_file, const nc::OpeningMode mode, const MPI_Comm comm)
{
    #ifdef CMC_WITH_NETCDF
    #ifdef CMC_WITH_NETCDF_PAR
    if (mode == nc::OpeningMode::Parallel)
    {
        ncid_ = nc::OpenParallel(path_to_file.c_str(), comm);
    } else

    #endif /* CMC_WITH_NETCDF_PAR */
    {
        if (mode == nc::OpeningMode::Parallel)
        {
            cmc_warn_msg("The netCDF file is ought to be opened for parallel access, although the parallel functionality is not accessible.");
            cmc_warn_msg("The file is opened for serial access.");
        }
        ncid_ = nc::OpenSerial(path_to_file.c_str());
    }

    #else
    cmc_err_msg("CMC is not linked against netCDF. Please recompile the build with netCDF\n.");
    return CMC_ERR;
    #endif
}

void
Data::SetHintLongitudeDimension(const int longitude_dimension_id)
{
    #ifdef CMC_WITH_NETCDF
    cmc_assert(longitude_dimension_id >= 0);
    coordinate_dimension_ids_[Dimension::Lon] = longitude_dimension_id;
    #endif
}

void
Data::SetHintLatitudeDimension(const int latitude_dimension_id)
{
    #ifdef CMC_WITH_NETCDF
    cmc_assert(latitude_dimension_id >= 0);
    coordinate_dimension_ids_[Dimension::Lat] = latitude_dimension_id;
    #endif
}

void
Data::SetHintHeightDimension(const int height_dimension_id)
{
    #ifdef CMC_WITH_NETCDF
    cmc_assert(height_dimension_id >= 0);
    coordinate_dimension_ids_[Dimension::Lev] = height_dimension_id;
    #endif
}

void
Data::SetHintTimeDimension(const int time_dimension_id)
{
    #ifdef CMC_WITH_NETCDF
    cmc_assert(time_dimension_id >= 0);
    coordinate_dimension_ids_[Dimension::Time] = time_dimension_id;
    #endif
}

/** Inquire the dimension lenghts and ids of the (geo-spatial) coordinate dimensions as well as their variables ids */ 
void
Data::InquireCoordinateDimensions()
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
        nc::CheckError(err);

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
    cmc_debug_msg("Lon: ", coordinate_lengths_[Dimension::Lon], ", Lat: ", coordinate_lengths_[Dimension::Lat], ", Lev: ", coordinate_lengths_[Dimension::Lev], ", Time: ", coordinate_lengths_[Dimension::Time]);
    
    #endif
}

void
Data::InquireCoordinates()
{
    #ifdef CMC_WITH_NETCDF
    /* Inquire the most general information about the supplied netCDF file (like number of diemnsion, number of global attributes, number of variables, etc.)*/
    int err = nc_inq(ncid_, &num_dimensions_, NULL, &num_global_attributes_, &id_unlimited_dimension_);
    nc::CheckError(err);

    InquireCoordinateDimensions();

    #endif
}

static
int 
GetVariableId(const int ncid, const std::string& variable_name)
{
    int variable_id;

    /* Inquire the netCDF internal variable id by the given variable's name */
    const int err = nc_inq_varid(ncid, variable_name.c_str(), &variable_id);
    nc::CheckError(err);

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
            nc::CheckError(err);
            attribute_value = double_att;
            break;
        case NC_INT:
            int32_t int_att;
            err = nc_get_att_int(ncid, variable_id, attribute_name.c_str(), &int_att);
            nc::CheckError(err);
            attribute_value = int_att;
            break;
        case NC_FLOAT:
            float float_att;
            err = nc_get_att_float(ncid, variable_id, attribute_name.c_str(), &float_att);
            nc::CheckError(err);
            attribute_value = float_att;
            break;
        case NC_SHORT:
            int16_t short_att;
            err = nc_get_att_short(ncid, variable_id, attribute_name.c_str(), &short_att);
            nc::CheckError(err);
            attribute_value = short_att;
            break;
        case NC_USHORT:
            uint16_t ushort_att;
            err = nc_get_att_ushort(ncid, variable_id, attribute_name.c_str(), &ushort_att);
            nc::CheckError(err);
            attribute_value = ushort_att;
            break;
        case NC_UINT:
            uint32_t uint_att;
            err = nc_get_att_uint(ncid, variable_id, attribute_name.c_str(), &uint_att);
            nc::CheckError(err);
            attribute_value = uint_att;
            break;
        case NC_INT64:
            long long int longlong_att;
            err = nc_get_att_longlong(ncid, variable_id, attribute_name.c_str(), &longlong_att);
            nc::CheckError(err);
            attribute_value = static_cast<int64_t>(longlong_att);
            break;
        case NC_UINT64:
            unsigned long long int ulonglong_att;
            err = nc_get_att_ulonglong(ncid, variable_id, attribute_name.c_str(), &ulonglong_att);
            nc::CheckError(err);
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
DetermineDataDistribution(const Hyperslab& global_hyperslab, const MPI_Comm comm)
{
    return DetermineSlicedDataDistribution(global_hyperslab, comm);
}

static
std::pair<std::vector<size_t>, std::vector<size_t>>
GetStartAndCountValuesForVariable(const Hyperslab& hyperslab, const CoordinateArray<int>& coordinate_dimension_ids, int num_dimensions, const std::array<int, NC_MAX_VAR_DIMS>& dimension_ids)
{
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
static input::Var
SetUpInputVariable(const int ncid, const CoordinateArray<int>& coordinate_dimension_ids, const int variable_id, std::string&& variable_name, const DataLayout layout,
                   const DomainIndex num_values, std::vector<Hyperslab>&& hyperslabs, GeoDomain&& global_domain,
                   const int num_dimensions, const std::array<int, NC_MAX_VAR_DIMS>& dimension_ids)
{
    input::Variable<T> variable(std::move(variable_name), variable_id, layout);

    variable.SetGlobalDomain(std::move(global_domain));

    std::vector<T> local_data(num_values);
    DomainIndex offset = 0;
    for (auto hs_iter = hyperslabs.begin(); hs_iter != hyperslabs.end(); ++hs_iter)
    {
        const auto [start_vals, count_vals] = GetStartAndCountValuesForVariable(*hs_iter, coordinate_dimension_ids, num_dimensions, dimension_ids);

        const int err = nc_get_vara(ncid, variable_id, start_vals.data(), count_vals.data(), static_cast<void*>(&local_data[offset]));
        nc::CheckError(err);

        offset += hs_iter->GetNumberCoordinates();
    }

    variable.SetDataAndCoordinates(std::move(local_data), std::move(hyperslabs));

    return input::Var(std::move(variable));
}

input::Var
Data::SetupVariableData(const int variable_nc_type, const int variable_id, std::string&& variable_name,
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
            return input::Var();
    }
}

input::Var
Data::InquireVariable(const Hyperslab& hyperslab, std::string&& variable_name)
{
    cmc_debug_msg("Variable ", variable_name, " will be inquired.");

    const int variable_id = GetVariableId(ncid_, variable_name);

    int variable_nc_type;
    int num_dimensions;
    std::array<int, NC_MAX_VAR_DIMS> dimension_ids;

    /* Inquire some information about the variable, e.g. the data type and its dimensions */
    int err = nc_inq_var(ncid_, variable_id, NULL, &variable_nc_type, &num_dimensions, dimension_ids.data(), NULL);
    nc::CheckError(err);
 
    /* Check some attributes of the variable */
    const auto [is_missing_value_present, missing_value] = CheckVariableForAttribute(ncid_, variable_id, "missing_value");
    const auto [is_add_offset_present, add_offset] = CheckVariableForAttribute(ncid_, variable_id, "add_offset");
    const auto [is_scale_factor_present, scale_factor] = CheckVariableForAttribute(ncid_, variable_id, "scale_factor");

    const int data_dimensionality = hyperslab.GetDimensionality();

    /* Determine the data layout of the variable */
    const DataLayout data_layout = GetDataLayout(data_dimensionality, hyperslab, coordinate_dimension_ids_, dimension_ids);

    /* Calculate a sliced distribution of the data */
    //TODO: Make a blocked distribution
    std::vector<Hyperslab> local_hyperslabs = DetermineDataDistribution(hyperslab, comm_);

    /* Get the number of values we will read from the file (for memory allocation) */
    DomainIndex num_local_data_values = 0;
    for (auto hs_iter = local_hyperslabs.begin(); hs_iter != local_hyperslabs.end(); ++hs_iter)
    {
        num_local_data_values += hs_iter->GetNumberCoordinates();
    }
    cmc_debug_msg("in inquire variables: Num lcoal values: ", num_local_data_values);
    /* Inquire the data of the variable and construct a input variable with it */
    input::Var variable = SetupVariableData(variable_nc_type, variable_id, std::move(variable_name), data_layout, num_local_data_values, std::move(local_hyperslabs), TransformHyperslabToGeoDomain(hyperslab), num_dimensions, dimension_ids);

    if (is_missing_value_present)
        variable.SetMissingValue(missing_value);
    
    if (is_add_offset_present)
        variable.SetAddOffset(add_offset);
    
    if (is_scale_factor_present)
        variable.SetScaleFactor(scale_factor);

    /* Set the MPI communcator */
    variable.SetMPIComm(comm_);
    cmc_debug_msg("Comm is: ", comm_);
    cmc_debug_msg("Variable ", variable.GetName(), " has been inquired.");
    
    return variable;
}

[[nodiscard]]
std::vector<input::Var>&&
Data::TransferData()
{
    if (_data_has_been_transfered_)
        cmc_err_msg("The data already has been transfered.");
    
    _data_has_been_transfered_ = true;

    return std::move(variables_);
}

}
