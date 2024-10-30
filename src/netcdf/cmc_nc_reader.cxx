#include "netcdf/cmc_nc_reader.hxx"

#include <algorithm>

namespace cmc
{

int
NcReader::NcOpen() 
{
    int ncid{-1};

#ifdef CMC_WITH_NETCDF_PAR
    if (comm_ != MPI_COMM_SELF)
    {
        /* Open for parallel access */
        const MPI_Info info = MPI_INFO_NULL;
        const int err = nc_open_par(file_name_.c_str(), NC_NOWRITE, comm_, info, &ncid);
        NcCheckError(err);
    } else
#endif
    {
        /* Otherwise open for serial access */
        const int err = nc__open(file_name_.c_str(), NC_NOWRITE, NULL, &ncid);
        NcCheckError(err);
    }

    is_file_opened_ = true;

    return ncid;
}

void
NcReader::NcClose(const int ncid)
{
    const int err = nc_close(ncid);
    NcCheckError(err);

    is_file_opened_ = false;
}

void
NcReader::InquireGeneralFileInformation(const int ncid)
{
    cmc_assert(is_file_opened_);

    if (has_general_information_been_inquired_) { return;}

    constexpr int kINvalidUnlimitedDimension = -1;

    int num_unlimited_dims{kINvalidUnlimitedDimension};

    /* Gather some information about the amount of dimensions, variables and attributes */
    int err = nc_inq(ncid, &num_dimensions_, &num_variables_, &num_global_attributes_, &num_unlimited_dims);
    NcCheckError(err);

    /* Get the netCDF format of this file */
    err = nc_inq_format(ncid, &netcdf_format_);
    NcCheckError(err);

    /* If there are unlimited dimensions, inquire their ids (In case it is not a netCDF-4 format, there is possibly only one unlimited dimension
     * and then it has to be the first declared dimension (=> id = 0), which is why this check is sufficient). In netCDF-4 format, the call nc_inq(...)
     * returns the number of unlimited diemnsion in the last parameter */
    if (num_unlimited_dims > 0)
    {
        /* Allocate the vector for the unlimited diemnsions */
        std::fill_n(std::back_inserter(unlimited_dimension_ids_), num_unlimited_dims, kINvalidUnlimitedDimension);

        /* Inquire the number of unlimited dimensions */
        err = nc_inq_unlimdims(ncid, NULL, unlimited_dimension_ids_.data());
        NcCheckError(err);
    } else
    {
        unlimited_dimension_ids_.push_back(num_unlimited_dims);
    }

    has_general_information_been_inquired_ = true;
}

void
NcReader::InquireGeneralFileInformation()
{
    if (has_general_information_been_inquired_) { return;}

    /* Open the file to be read */
    const int ncid = NcOpen();

    InquireGeneralFileInformation(ncid);

    /* Close the file after the reading process is finished */
    NcClose(ncid);
}

int
NcReader::FindVariableID(const int ncid, const std::string& variable_name)
{
    cmc_assert(is_file_opened_);

    InquireGeneralFileInformation(ncid);

    char var_name[NC_MAX_NAME];

    /* Iterate over all variables (their IDs, correspond to 0, ..., num_variables -1) */
    for (int var_id = 0; var_id < num_variables_; ++var_id)
    {
        /* Inquire the name of the variable */
        int err = nc_inq_varname(ncid, var_id, var_name);
        NcCheckError(err);

        /* Check if the variable name complies with the given one */
        if (std::strcmp(var_name, variable_name.c_str()) == 0)
        {
            /* We have found a variable with the given name */
            return var_id;
        }
    }

    cmc_warn_msg("No variable with the given name has been found in the netCDF file.");

    return NC_EBADTYPID;
}

GeneralHyperslab
NcReader::GetDataDomainAsGeneralHyperslab(const int ncid, const std::string& variable_name)
{
    cmc_assert(is_file_opened_);

    InquireGeneralFileInformation(ncid);

    /* Get the corresponding ID to the supplied variable name */
    const int var_id = FindVariableID(ncid, variable_name);

    /* Inquire the number of dimensions */
    int num_dims{0};
    int err = nc_inq_varndims(ncid, var_id, &num_dims);
    NcCheckError(err);

    /* Inquire the dimension ids */
    std::vector<int> dim_ids(num_dims, 0);
    err = nc_inq_vardimid(ncid, var_id, dim_ids.data());
    NcCheckError(err);

    /* Create vectors for the start and count values */
    std::vector<size_t> start_values(num_dims, 0);
    std::vector<size_t> count_values(num_dims, 0);

    /* Get the length of each dimension */
    int id{0};
    for (auto dim_id_iter = dim_ids.begin(); dim_id_iter != dim_ids.end(); ++dim_id_iter, ++id)
    {
        const int dim_err = nc_inq_dimlen(ncid, *dim_id_iter, &count_values[id]);
        NcCheckError(dim_err);
    }

    /* Create a GeneralHyperslab for the whole domain */
    return GeneralHyperslab(std::move(start_values), std::move(count_values));
}

GeneralHyperslab
NcReader::GetDataDomainAsGeneralHyperslab(const std::string& variable_name)
{
    /* Open the file to be read */
    const int ncid = NcOpen();

    /* Get the hyperslab of the variable from the opened file */
    GeneralHyperslab hs = GetDataDomainAsGeneralHyperslab(ncid, variable_name);

    /* Close the file after the reading process is finished */
    NcClose(ncid);

    return hs;
}

NcAttribute
ReadAttribute(const int ncid, const int var_id, const char* att_name, const nc_type type)
{
    switch (type)
    {
        case NC_BYTE:
        {
            signed char value{0};
            const int err = nc_get_att_schar(ncid, var_id, att_name, &value);
            NcCheckError(err);
            return NcAttribute(std::string(att_name), static_cast<char>(value));
        }
        break;
        case NC_CHAR:
        {
            signed char value{0};
            const int err = nc_get_att_schar(ncid, var_id, att_name, &value);
            NcCheckError(err);
            return NcAttribute(std::string(att_name), static_cast<char>(value));
        }
        break;
        case NC_SHORT:
        {
            short value{0};
            const int err = nc_get_att_short(ncid, var_id, att_name, &value);
            NcCheckError(err);
            return NcAttribute(std::string(att_name), static_cast<int16_t>(value));
        }
        break;
        case NC_INT:
        {
            int value{0};
            const int err = nc_get_att_int(ncid, var_id, att_name, &value);
            NcCheckError(err);
            return NcAttribute(std::string(att_name), static_cast<int32_t>(value));
        }
        break;
        case NC_FLOAT:
        {
            float value{static_cast<float>(0.0)};
            const int err = nc_get_att_float(ncid, var_id, att_name, &value);
            NcCheckError(err);
            return NcAttribute(std::string(att_name), value);
        }
        break;
        case NC_DOUBLE:
        {
            double value{0.0};
            const int err = nc_get_att_double(ncid, var_id, att_name, &value);
            NcCheckError(err);
            return NcAttribute(std::string(att_name), value);
        }
        break;
        case NC_USHORT:
        {
            unsigned short value{0};
            const int err = nc_get_att_ushort(ncid, var_id, att_name, &value);
            NcCheckError(err);
            return NcAttribute(std::string(att_name), static_cast<uint16_t>(value));
        }
        break;
        case NC_UINT:
        {
            unsigned int value{0};
            const int err = nc_get_att_uint(ncid, var_id, att_name, &value);
            NcCheckError(err);
            return NcAttribute(std::string(att_name), static_cast<uint32_t>(value));
        }
        break;
        case NC_INT64:
        {
            long long value{0};
            const int err = nc_get_att_longlong(ncid, var_id, att_name, &value);
            NcCheckError(err);
            return NcAttribute(std::string(att_name), static_cast<int64_t>(value));
        }
        break;
        case NC_UINT64:
        {
            unsigned long long value{0};
            const int err = nc_get_att_ulonglong(ncid, var_id, att_name, &value);
            NcCheckError(err);
            return NcAttribute(std::string(att_name), static_cast<uint64_t>(value));
        }
        break;
        default:
            cmc_err_msg("The netCDF attribute has an invalid data type. (String attributes cannot be processed yet).");
            return NcAttribute();
    }
}

std::vector<NcAttribute>
NcReader::InquireAttributes(const int ncid, const int var_id)
{
    /* Inquire the amount of attributes for the variable */
    int num_atts{0};
    if (var_id != NC_GLOBAL)
    {
        int err = nc_inq_varnatts(ncid, var_id, &num_atts);
        NcCheckError(err);
    } else
    {
        int err = nc_inq_natts(ncid, &num_atts);
        NcCheckError(err);
    }

    const int num_attributes = num_atts;


    std::vector<NcAttribute> attributes;

    /* Iterate over the amount of attributes; the attributes have the IDs 0, ..., (num_attributes -1) */
    char attribute_name[NC_MAX_NAME];
    for (int att_id = 0; att_id < num_attributes; ++att_id)
    {
        int err = nc_inq_attname(ncid, var_id, att_id, attribute_name);
        NcCheckError(err);

        /* Get the data type of the attribute */
        nc_type type{0};
        err = nc_inq_atttype(ncid, var_id, attribute_name, &type);
        NcCheckError(err);

        /* Get the number of values stroed within the attribute */
        size_t att_len{0};
        err = nc_inq_attlen(ncid, var_id, attribute_name, &att_len);
        NcCheckError(err);

        /* Only length-one attributes are supported; i.e. vectors/arrays are not compatible */
        if (att_len > 1)
        {
            cmc_warn_msg("Currently, there is only support for attributes of length one.");
            continue;
        }

        /* Read in the attribute */
        attributes.push_back(ReadAttribute(ncid, var_id, attribute_name, type));
    }

    return attributes;
}

std::vector<NcDimension>
NcReader::ConvertDimensionIDs(const int ncid, const std::vector<int>& dim_ids)
{
    cmc_assert(is_file_opened_);

    char dim_name[NC_MAX_NAME];

    std::vector<NcDimension> dimensions;
    dimensions.reserve(dim_ids.size());

    for (auto dim_iter = dim_ids.begin(); dim_iter != dim_ids.end(); ++dim_iter)
    {
        /* Get the name of the dimension */
        int err = nc_inq_dimname(ncid, *dim_iter, dim_name);
        NcCheckError(err);

        /* Get the length of the dimension */
        size_t dim_length{0};
        err = nc_inq_dimlen(ncid, *dim_iter, &dim_length);
        NcCheckError(err);

        /* Make an NcDimension out of the inquired information */
        dimensions.emplace_back(std::string(dim_name), dim_length, *dim_iter);
    }

    return dimensions;
}

std::vector<NcVariable>
NcReader::InquireVariableMetaData(const int ncid)
{
    char var_name[NC_MAX_NAME];

    std::vector<NcVariable> variables;
    variables.reserve(num_variables_);

    /* Iterate over all variables (their IDs, correspond to 0, ..., num_variables -1) */
    for (int var_id = 0; var_id < num_variables_; ++var_id)
    {
        /* Inquire the name of the variable */
        int err = nc_inq_varname(ncid, var_id, var_name);
        NcCheckError(err);

        /* Inquire the data type of the variable */
        nc_type type{NC_NAT};
        err = nc_inq_vartype(ncid, var_id, &type);
        NcCheckError(err);

        /* Inquire the number of dimensions */
        int num_dims{0};
        err = nc_inq_varndims(ncid, var_id, &num_dims);
        NcCheckError(err);

        /* Inquire the dimension ids */
        std::vector<int> dim_ids(num_dims, 0);
        err = nc_inq_vardimid(ncid, var_id, dim_ids.data());
        NcCheckError(err);

        /* Create an NcVariable out of the information */
        variables.emplace_back(InquireAttributes(ncid, var_id), ConvertDimensionIDs(ncid, dim_ids));
        variables.back().SetupSpecificVariable(var_name, ConvertNcTypeToCmcType(type));
    }

    return variables;
}

std::vector<NcAttribute>
NcReader::ReadGlobalAttrtibutes()
{
    /* Open the file to be read */
    const int ncid = NcOpen();

    /* Inquire some basic/meta information about the file and it's contents */
    InquireGeneralFileInformation(ncid);

    /* Inquire global attributes */
    const std::vector<NcAttribute> global_attributes = InquireAttributes(ncid, NC_GLOBAL);

    /* Close the file after the reading process is finished */
    NcClose(ncid);

    return global_attributes;
}

std::vector<NcVariable>
NcReader::ReadVariableMetaData()
{
    /* Open the file to be read */
    const int ncid = NcOpen();

    /* Inquire some basic/meta information about the file and it's contents */
    InquireGeneralFileInformation(ncid);

    /* Inquire the meta data of all variables */
    const std::vector<NcVariable> variables = InquireVariableMetaData(ncid);

    /* Close the file after the reading process is finished */
    NcClose(ncid);

    return variables;
}

std::pair<std::vector<NcVariable>, std::vector<NcAttribute>>
NcReader::ReadVariableMetaDataAndGlobalAttributes()
{
    /* Open the file to be read */
    const int ncid = NcOpen();

    /* Inquire some basic/meta information about the file and it's contents */
    InquireGeneralFileInformation(ncid);

    /* Inquire the meta data of all variables */
    const std::vector<NcVariable> variables = InquireVariableMetaData(ncid);

    /* Inquire global attributes */
    const std::vector<NcAttribute> global_attributes = InquireAttributes(ncid, NC_GLOBAL);

    /* Close the file after the reading process is finished */
    NcClose(ncid);

    return std::make_pair(std::move(variables), std::move(global_attributes));
}

void
NcReader::ReadVariableDataFromFile(const int ncid, const nc_type var_type, const std::string& var_name, const int var_id, const std::vector<GeneralHyperslab>& hyperslabs, NcVariable& variable)
{
    /* Iterate over all hyperslabs and count the amount of data which will be read */
    size_t num_data_values{0};
    for (auto hs_iter = hyperslabs.begin(); hs_iter != hyperslabs.end(); ++hs_iter)
    {
        num_data_values += hs_iter->GetNumberOfCoveredCoordinates();
    }

    /* Create a general variable capable of holding the data of the specified type */
    NcGeneralVariable general_variable = CreateSpecificVariable(var_type, var_name, var_id, num_data_values);

    /* Read in the data from the file */
    size_t values_read{0};
    for (auto hs_iter = hyperslabs.begin(); hs_iter != hyperslabs.end(); ++hs_iter)
    {
        /* Read the hyperslab data into the variable */
        std::visit([&](auto&& var){
            var.OverwriteDataFromFile(ncid, var_id, *hs_iter, values_read);
        }, general_variable);
        
        /* Add the read values to the offset variable */
        values_read += hs_iter->GetNumberOfCoveredCoordinates();
    }

    /* Set the variable data within the passed NcVariable which is ought to be filled */
    variable.SetSpecificVariable(std::move(general_variable));
}

std::pair<std::vector<NcVariable>, std::vector<NcAttribute>>
NcReader::ReadVariables()
{
    if (variable_stash_.empty())
    {
        cmc_warn_msg("No variables have been specified for reading.");
    }

    /* Open the file to be read */
    const int ncid = NcOpen();

    /* Inquire some basic/meta information about the file and it's contents */
    InquireGeneralFileInformation(ncid);

    std::vector<NcVariable> variables;
    variables.reserve(variable_stash_.size());

    /* Iterate over all stashed variables */
    for (auto stashed_var_iter = variable_stash_.begin(); stashed_var_iter != variable_stash_.end(); ++stashed_var_iter)
    {
        /* Inquire the ID of the variable */
        int var_id{-1};
        int err = nc_inq_varid(ncid, stashed_var_iter->name.c_str(), &var_id);
        NcCheckError(err);

        /* Inquire the data type of the variable */
        nc_type type{NC_NAT};
        err = nc_inq_vartype(ncid, var_id, &type);
        NcCheckError(err);

        /* Inquire the number of dimensions */
        int num_dims{0};
        err = nc_inq_varndims(ncid, var_id, &num_dims);
        NcCheckError(err);

        /* Inquire the dimension ids */
        std::vector<int> dim_ids(num_dims, 0);
        err = nc_inq_vardimid(ncid, var_id, dim_ids.data());
        NcCheckError(err);

        /* Create an NcVariable out of the information */
        variables.emplace_back(InquireAttributes(ncid, var_id), ConvertDimensionIDs(ncid, dim_ids));

        /* Inquire the data of the variable */
        ReadVariableDataFromFile(ncid, type, stashed_var_iter->name, var_id, stashed_var_iter->hyperslabs, variables.back());
    }

    /* Inquire global attributes */
    const std::vector<NcAttribute> global_attributes = InquireAttributes(ncid, NC_GLOBAL);

    /* Close the file after the reading process is finished */
    NcClose(ncid);

    return std::make_pair(std::move(variables), std::move(global_attributes));
}

NcVariable
NcReader::ReadVariable(const std::string& variable_name)
{
    /* Open the file to be read */
    const int ncid = NcOpen();

    /* Inquire some basic/meta information about the file and it's contents */
    InquireGeneralFileInformation(ncid);

    /* Get the corresponding variable ID */
    const int var_id = FindVariableID(ncid, variable_name);

    /* Inquire the data type of the variable */
    nc_type type{NC_NAT};
    int err = nc_inq_vartype(ncid, var_id, &type);
    NcCheckError(err);

    /* Inquire the number of dimensions */
    int num_dims{0};
    err = nc_inq_varndims(ncid, var_id, &num_dims);
    NcCheckError(err);

    /* Inquire the dimension ids */
    std::vector<int> dim_ids(num_dims, 0);
    err = nc_inq_vardimid(ncid, var_id, dim_ids.data());
    NcCheckError(err);

    /* Create an NcVariable out of the information */
    NcVariable nc_var(InquireAttributes(ncid, var_id), ConvertDimensionIDs(ncid, dim_ids));

    std::vector<GeneralHyperslab> hyperslab;
    hyperslab.push_back(GetDataDomainAsGeneralHyperslab(ncid, variable_name));

    /* Inquire the data of the variable */
    ReadVariableDataFromFile(ncid, type, variable_name, var_id, hyperslab, nc_var);

    /* Close the file after the reading process is finished */
    NcClose(ncid);

    return nc_var;
}

std::vector<NcAttribute>
NcReader::ReadVariableAttributes(const std::string& variable_name)
{
    /* Open the file to be read */
    const int ncid = NcOpen();

    /* Inquire some basic/meta information about the file and it's contents */
    InquireGeneralFileInformation(ncid);

    /* Get the corresponding variable ID */
    const int var_id = FindVariableID(ncid, variable_name);

    /* Inquire the attributes of the given variable */
    std::vector<NcAttribute> attributes = InquireAttributes(ncid, var_id);

    /* Close the file after the reading process is finished */
    NcClose(ncid);

    return attributes;
}

const std::string&
NcReader::GetFileName() const
{
    return file_name_;
}

int
NcReader::GetNetcdfFormat() const
{
    cmc_assert(has_general_information_been_inquired_);
    return netcdf_format_;
}

int
NcReader::GetNumberOfDimensions() const
{
    cmc_assert(has_general_information_been_inquired_);
    return num_dimensions_;
}

int
NcReader::GetNumberOfVariables() const
{
    cmc_assert(has_general_information_been_inquired_);
    return num_variables_;
}

int
NcReader::GetNumberOfGlobalAttributes() const
{
    cmc_assert(has_general_information_been_inquired_);
    return num_global_attributes_;
}

int
NcReader::GetNumberOfUnlimitedDimensions() const
{
    cmc_assert(has_general_information_been_inquired_);
    return unlimited_dimension_ids_.size();
}   

std::vector<int>
NcReader::GetUnlimitedDimensionIDs() const
{
    cmc_assert(has_general_information_been_inquired_);
    return unlimited_dimension_ids_;
}

CmcType
NcReader::GetTypeOfVariable(const std::string& variable_name)
{
    /* Open the file to be read */
    const int ncid = NcOpen();

    /* Inquire some basic/meta information about the file and it's contents */
    InquireGeneralFileInformation(ncid);

    /* Get the corresponding variable ID */
    const int var_id = FindVariableID(ncid, variable_name);

    nc_type var_type{NC_NAT};
    const int type_err = nc_inq_vartype(ncid, var_id, &var_type);
    NcCheckError(type_err);

    /* Close the file after the reading process is finished */
    NcClose(ncid);

    if (var_type == NC_NAT)
    {
        cmc_warn_msg("Either the variable type is invalid or no variable corresponds to the name ", variable_name, ". Therefore, no type information could be retrieved.");
    }

    return ConvertNcTypeToCmcType(var_type);
}

std::vector<NcDimension>
NcReader::GetVariableDimensions(const std::string& variable_name)
{
    /* Open the file to be read */
    const int ncid = NcOpen();

    /* Inquire some basic/meta information about the file and it's contents */
    InquireGeneralFileInformation(ncid);

    /* Get the corresponding variable ID */
    const int var_id = FindVariableID(ncid, variable_name);

    /* Inquire the number of dimensions */
    int num_dims{0};
    int err = nc_inq_varndims(ncid, var_id, &num_dims);
    NcCheckError(err);

    /* Inquire the dimension ids */
    std::vector<int> dim_ids(num_dims, 0);
    err = nc_inq_vardimid(ncid, var_id, dim_ids.data());
    NcCheckError(err);

    std::vector<NcDimension> dims = ConvertDimensionIDs(ncid, dim_ids);

    /* Close the file after the reading process is finished */
    NcClose(ncid);

    return dims;
}

void
NcReader::StashVariableForReading(const std::string& variable_name, const std::vector<GeneralHyperslab>& hyperslabs)
{
    variable_stash_.emplace_back(variable_name, hyperslabs);
}

void
NcReader::StashVariableForReading(const std::string& variable_name, std::vector<GeneralHyperslab>&& hyperslabs)
{
    variable_stash_.emplace_back(variable_name, std::move(hyperslabs));
}

void
NcReader::StashVariableForReading(const std::string& variable_name, const GeneralHyperslab& hyperslab)
{
    variable_stash_.emplace_back(variable_name, hyperslab);
}   

void
NcReader::StashVariableForReading(const std::string& variable_name, GeneralHyperslab&& hyperslab)
{
    variable_stash_.emplace_back(variable_name, std::move(hyperslab));
}

void
NcReader::ClearStashedVariables()
{
    variable_stash_.clear();
}

}
