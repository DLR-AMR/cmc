#include "netcdf/cmc_nc_writer.hxx"

namespace cmc::nc
{

int
Writer::Open() 
{
    if (!file_has_been_created_)
    {
        return Create();
    }

    int ncid{-1};

#ifdef CMC_WITH_NETCDF_PAR
    if (comm_ != MPI_COMM_SELF)
    {
        /* Open for parallel access */
        const MPI_Info info = MPI_INFO_NULL;
        const int err = nc_open_par(file_name_.c_str(), NC_WRITE, comm_, info, &ncid);
        CheckError(err);
    } else
#endif
    {
        /* Otherwise open for serial access */
        const int err = nc__open(file_name_.c_str(), NC_WRITE, NULL, &ncid);
        CheckError(err);
    }

    return ncid;
}

int
Writer::Create()
{
    int ncid{-1};

    if (netcdf_format_ != NC_64BIT_OFFSET && netcdf_format_ != NC_64BIT_DATA && netcdf_format_ != NC_NETCDF4 &&
        netcdf_format_ != NC_CLASSIC_MODEL && netcdf_format_ != NC_DISKLESS)
    {
        cmc_err_msg("The netCDF file could not be created since the supplied netcdf_format is invalid.");
    }

#ifdef CMC_WITH_NETCDF_PAR
    if (comm_ != MPI_COMM_SELF)
    {
        /* Open for parallel access */
        const MPI_Info info = MPI_INFO_NULL;
        const int err = nc_create_par(file_name_.c_str(), NC_CLOBBER | NC_WRITE | netcdf_format_, comm_, info, &ncid);
        CheckError(err);
    } else
#endif
    {
        /* Otherwise open for serial access */
        const int err = nc_create(file_name_.c_str(), NC_CLOBBER | NC_WRITE | netcdf_format_, &ncid);
        CheckError(err);
    }

    file_has_been_created_ = true;

    return ncid;
}

void
Writer::Close(const int ncid) const
{
    const int err = nc_close(ncid);
    CheckError(err);
}

void
Writer::ReserveVariables(const size_t num_variables)
{
    variables_.reserve(num_variables);
}

void
Writer::AddGlobalAttribute(const Attribute& attribute)
{
    global_attributes_.push_back(attribute);
}


void
Writer::AddGlobalAttribute(Attribute&& attribute)
{
    global_attributes_.push_back(std::move(attribute));
}

void
Writer::AddVariable(const Variable& variable)
{
    variables_.push_back(variable);
}

void
Writer::AddVariable(Variable&& variable)
{
    variables_.push_back(std::move(variable));
}   

void
Writer::ClearVariablesAndGlobalAttributes()
{
    variables_.clear();
    global_attributes_.clear();
}

std::vector<int>
Writer::DefineVariableDimensions(const int ncid, const std::vector<Dimension>& dimensions)
{
    std::vector<int> dim_ids(dimensions.size(), -1);
    int dim_count{0};

    for (auto dim_iter = dimensions.begin(); dim_iter != dimensions.end(); ++dim_iter, ++dim_count)
    {
        const int err = nc_def_dim(ncid, dim_iter->GetName().c_str(), dim_iter->GetLength(), &dim_ids[dim_count]);
        CheckError(err);
    }

    return dim_ids;
}

/* Overload construct which writes an attribute based on the parameter data type */
struct WriteAttribute
{
    WriteAttribute() = delete;
    WriteAttribute(const int ncid, const int var_id, const std::string name)
    : ncid_{ncid}, var_id_{var_id}, name_{name}{};
    
    void operator()(const int8_t val)
    {
        const signed char value = static_cast<signed char>(val);
        const int err = nc_put_att_schar(ncid_, var_id_, name_.c_str(), NC_BYTE, 1, &value);
        CheckError(err);
    }
    void operator()(const char val)
    {
        const signed char value = static_cast<signed char>(val);
        const int err = nc_put_att_schar(ncid_, var_id_, name_.c_str(), NC_CHAR, 1, &value);
        CheckError(err);
    }
    void operator()(const int16_t val)
    {
        const int err = nc_put_att_short(ncid_, var_id_, name_.c_str(), NC_SHORT, 1, &val);
        CheckError(err);
    }
    void operator()(const int32_t val)
    {
        const int err = nc_put_att_int(ncid_, var_id_, name_.c_str(), NC_INT, 1, &val);
        CheckError(err);
    }
    void operator()(const float val)
    {
        const int err = nc_put_att_float(ncid_, var_id_, name_.c_str(), NC_FLOAT, 1, &val);
        CheckError(err);
    }
    void operator()(const double val)
    {
        const int err = nc_put_att_double(ncid_, var_id_, name_.c_str(), NC_DOUBLE, 1, &val);
        CheckError(err);
    }
    void operator()(const uint8_t val)
    {
        const int err = nc_put_att_ubyte(ncid_, var_id_, name_.c_str(), NC_UBYTE, 1, &val);
        CheckError(err);
    }
    void operator()(const uint16_t val)
    {
        const int err = nc_put_att_ushort(ncid_, var_id_, name_.c_str(), NC_USHORT, 1, &val);
        CheckError(err);
    }
    void operator()(const uint32_t val)
    {
        const int err = nc_put_att_uint(ncid_, var_id_, name_.c_str(), NC_UINT, 1, &val);
        CheckError(err);
    }
    void operator()(const int64_t val)
    {
        const long long value = static_cast<long long>(val); 
        const int err = nc_put_att_longlong(ncid_, var_id_, name_.c_str(), NC_INT64, 1, &value);
        CheckError(err);
    }
    void operator()(const uint64_t val)
    {
        const unsigned long long value = static_cast<unsigned long long>(val); 
        const int err = nc_put_att_ulonglong(ncid_, var_id_, name_.c_str(), NC_UINT64, 1, &value);
        CheckError(err);
    }
private:
    const int ncid_;
    const int var_id_;
    const std::string name_;
};

void 
Writer::DefineAttributes(const int ncid, const int var_id, const std::vector<Attribute>& attributes)
{
    for (auto att_iter = attributes.begin(); att_iter != attributes.end(); ++att_iter)
    {
        std::visit(WriteAttribute(ncid, var_id, att_iter->GetName()), att_iter->GetValue());
    }
}

std::vector<int>
Writer::DefineVariables(const int ncid)
{
    std::vector<int> var_ids(variables_.size(), -1);
    int var_count{0};

    for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter, ++var_count)
    {
        /* Extract the dimensions from the variable and define them as dimension in the netCDF file */
        const std::vector<Dimension> var_dimensions = var_iter->GetDimensionsFromVariable();

        /* Define the dimensions of the variable */
        const std::vector<int> dimension_ids = DefineVariableDimensions(ncid, var_dimensions);

        /* Define the variable itself */
        const int err = nc_def_var(ncid, var_iter->GetName().c_str(), ConvertCmcTypeToNcType(var_iter->GetCmcType()), static_cast<int>(dimension_ids.size()), dimension_ids.data(), &var_ids[var_count]);
        CheckError(err);

        /* Create an attribute regarding the ID of the variable (since the dimensions are labelled accordingly) */
        var_iter->CreateIDAttribute();

        /* Define its attributes */
        DefineAttributes(ncid, var_ids[var_count], var_iter->GetAttributes());
    }

    return var_ids;
}

void
Writer::DefineGlobalAttributes(const int ncid)
{
    for (auto att_iter = global_attributes_.begin(); att_iter != global_attributes_.end(); ++att_iter)
    {
        std::visit(WriteAttribute(ncid, NC_GLOBAL, att_iter->GetName()), att_iter->GetValue());
    }
}

void 
Writer::WriteData(const int ncid, const std::vector<int>& variable_ids)
{
    int var_count = 0;
    for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter, ++var_count)
    {
        var_iter->WriteVariableData(ncid, variable_ids[var_count]);
    }
}

void
Writer::Write()
{
    const int ncid = Open();

    const std::vector<int> var_ids = DefineVariables(ncid);

    DefineGlobalAttributes(ncid);

    /* Switch from define to data mode */
    const int err = nc_enddef(ncid);
    CheckError(err);

    WriteData(ncid, var_ids);

    Close(ncid);
}
}
