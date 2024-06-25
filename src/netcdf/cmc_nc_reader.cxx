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
            cmc_err_msg("The netCDF attriute has an invalid data type.");
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

        /* Inquire the number of dimensions */
        int num_dims{0};
        err = nc_inq_varndims(ncid, var_id, &num_dims);
        NcCheckError(err);

        /* Inquire the data type of the variable */
        nc_type type{NC_NAT};
        err = nc_inq_vartype(ncid, var_id, &type);
        NcCheckError(err);

        /* Inquire the dimension ids */
        std::vector<int> dim_ids(num_dims, 0);
        err = nc_inq_vardimid(ncid, var_id, dim_ids.data());
        NcCheckError(err);

        /* Create an NcVariable out of the information */
        variables.emplace_back(InquireAttributes(ncid, var_id), ConvertDimensionIDs(ncid, dim_ids));
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

    const std::vector<NcVariable> variables = InquireVariableMetaData(ncid);

    /* Close the file after the reading process is finished */
    NcClose(ncid);

    return variables;
}

}
