#ifndef CMC_NC_IO_HXX
#define CMC_NC_IO_HXX

#include "cmc.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_hyperslab.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "netcdf/cmc_netcdf.hxx"
#include "netcdf/cmc_nc_io_conventions.hxx"

#include <string>
#include <vector>
#include <variant>

namespace cmc
{

struct NcDimension;
class NcAttribute;
template<typename T>
class NcSpecificVariable;
class NcVariable;

std::vector<NcAttribute>::const_iterator
FindAttribute(const std::vector<NcAttribute>& attributes, const std::string& attr_name);

class NcDimension
{
public:
    NcDimension(const std::string& name, const size_t length)
    : name_{name}, length_{length} {};
    NcDimension(std::string&& name, const size_t length)
    : name_{std::move(name)}, length_{length} {};
    NcDimension(const std::string& name, const size_t length, const int dim_id)
    : name_{name}, length_{length}, nc_dim_id{dim_id} {};
    NcDimension(std::string&& name, const size_t length, const int dim_id)
    : name_{std::move(name)}, length_{length}, nc_dim_id{dim_id} {};

    std::string GetName(){return name_;};
    const std::string& GetName() const {return name_;};
    size_t GetLength() const {return length_;};

    void AppendNumToName(const int number){name_.append(std::to_string(number));};

private:
    std::string name_;
    size_t length_;
    int nc_dim_id{-1};
}; 

class NcAttribute
{
public:
    NcAttribute() = default;
    NcAttribute(const std::string& name, const CmcUniversalType& value)
    : name_{name}, value_{value}{};
    NcAttribute(std::string&& name, const CmcUniversalType& value)
    : name_{std::move(name)}, value_{value}{};
    template<typename T> NcAttribute(const std::string& name, T value)
    : name_{name}, value_{value}{};
    template<typename T> NcAttribute(std::string&& name, T value)
    : name_{std::move(name)}, value_{value}{};
    ~NcAttribute() = default;

    NcAttribute(const NcAttribute& other) = default;
    NcAttribute& operator=(const NcAttribute& other) = default;
    NcAttribute(NcAttribute&& other) = default;
    NcAttribute& operator=(NcAttribute&& other) = default;

    const std::string& GetName() const;
    std::string GetName();
    
    const CmcUniversalType& GetValue() const;
    CmcUniversalType GetValue();

private:
    std::string name_;
    CmcUniversalType value_;
};


struct GeneralHyperslab
{
public:
    size_t GetNumberOfCoveredCoordinates() const
    {
        size_t num_coords{1};
        for (auto cv_iter = count_values.begin(); cv_iter != count_values.end(); ++cv_iter)
        {
            num_coords *= *cv_iter;
        }
        return num_coords;
    }

    std::vector<size_t> start_values;
    std::vector<size_t> count_values;
};


template<typename T>
class NcSpecificVariable
{
public:
    NcSpecificVariable() = default;
    NcSpecificVariable(const int id)
    : id_{id} {};
    NcSpecificVariable(const std::string& name, const int id)
    : name_{name}, id_{id} {};
    NcSpecificVariable(const std::string& name, const int id, const size_t size_hint)
    : name_{name}, id_{id} {
        data_ = std::vector<T>(size_hint);
    };
    ~NcSpecificVariable() = default;

    NcSpecificVariable(const NcSpecificVariable& other) = default;
    NcSpecificVariable& operator=(const NcSpecificVariable& other) = default;
    NcSpecificVariable(NcSpecificVariable&& other) = default;
    NcSpecificVariable& operator=(NcSpecificVariable&& other) = default;

    const std::string& GetName() const;
    std::vector<NcDimension> GetDimensionsFromVariable() const;
    CmcType GetCmcType() const;
    int GetID() const;

    const GeoDomain& GetGlobalDomain() const;
    void SetGlobalDomain(const GeoDomain& domain);
    void SetGlobalDomain(GeoDomain&& domain);

    void SetDataLayout(const DataLayout layout);

    void PushBack(const std::vector<T>& values, const Hyperslab& hyperslab);
    void PushBack(const std::vector<T>& values, Hyperslab&& hyperslab);
    void PushBack(std::vector<T>&& values, std::vector<Hyperslab>&& hyperslabs);
    void PushBack(const std::vector<T>& values, const std::vector<Hyperslab>& hyperslabs);
    void PushBack(const std::vector<T>& values, const int global_sfc_offset);
    void PushBack(std::vector<T>&& values, const int global_sfc_offset);

    void OverwriteDataFromFile(const int ncid, const int var_id, const GeneralHyperslab& hyperslab, const size_t start_offset = 0);

    const std::vector<T>& GetData() const;
    T GetMissingValue() const;
    void SetMissingValue(const CmcUniversalType& missing_value);
    void SetMissingValue(const T missing_value);

    void WriteVariableData(const int ncid, const int var_id) const;
private:
    std::string name_;
    int id_;
    std::vector<T> data_;
    T missing_value_;
    DataLayout layout_;
    DataFormat format_;
    std::vector<Hyperslab> hyperslabs_;
    LinearIndex global_sfc_offset_;
    GeoDomain global_domain_;
    bool is_data_stack_full_{false};
};

using NcGeneralVariable = std::variant<NcSpecificVariable<int8_t>, NcSpecificVariable<char>, NcSpecificVariable<int16_t>, NcSpecificVariable<int32_t>, NcSpecificVariable<float>, NcSpecificVariable<double>,
                                       NcSpecificVariable<uint8_t>, NcSpecificVariable<uint16_t>, NcSpecificVariable<uint32_t>, NcSpecificVariable<int64_t>, NcSpecificVariable<uint64_t>>;


class NcVariable
{
public:
    NcVariable() = default;
    ~NcVariable() = default;

    template<class T> NcVariable(const NcSpecificVariable<T>& variable)
    : variable_{variable} {};
    template<class T> NcVariable(NcSpecificVariable<T>&& variable)
    : variable_{std::move(variable)} {};
    template<class T> NcVariable(const NcSpecificVariable<T>& variable, const std::vector<NcAttribute>& attributes)
    : variable_{variable}, attributes_{attributes} {};
    template<class T> NcVariable(const NcSpecificVariable<T>& variable, std::vector<NcAttribute>&& attributes)
    : variable_{variable}, attributes_{std::move(attributes)} {};
    template<class T> NcVariable(NcSpecificVariable<T>&& variable, const std::vector<NcAttribute>& attributes)
    : variable_{std::move(variable)}, attributes_{attributes} {};
    template<class T> NcVariable(NcSpecificVariable<T>&& variable, std::vector<NcAttribute>&& attributes)
    : variable_{std::move(variable)}, attributes_{std::move(attributes)} {};
    template<class T> NcVariable(const NcSpecificVariable<T>& variable, const NcAttribute& attribute)
    : variable_{variable} {
        attributes_.push_back(attribute);
    };
    template<class T> NcVariable(const NcSpecificVariable<T>& variable, NcAttribute&& attribute)
    : variable_{variable} {
        attributes_.push_back(attribute);
    };
    template<class T> NcVariable(NcSpecificVariable<T>&& variable, const NcAttribute& attribute)
    : variable_{std::move(variable)} {
        attributes_.push_back(attribute);
    };
    template<class T> NcVariable(NcSpecificVariable<T>&& variable, NcAttribute&& attribute)
    : variable_{std::move(variable)} {
        attributes_.push_back(attribute);
    };

    NcVariable(std::vector<NcAttribute>&& attributes, std::vector<NcDimension>&& dimensions)
    : attributes_{std::move(attributes)}, dimensions_{std::move(dimensions)} {};

    NcVariable(const NcVariable& other) = default;
    NcVariable& operator=(const NcVariable& other) = default;
    NcVariable(NcVariable&& other) = default;
    NcVariable& operator=(NcVariable&& other) = default;
    
    template<class T> void SetSpecificVariable(const NcSpecificVariable<T>& variable);
    template<class T> void SetSpecificVariable(NcSpecificVariable<T>&& variable);
    void SetSpecificVariable(const NcGeneralVariable& variable);
    void SetSpecificVariable(NcGeneralVariable&& variable);

    std::vector<NcDimension> GetDimensionsFromVariable() const;
    const std::string& GetName() const;
    CmcType GetCmcType() const;
    void CreateIDAttribute();
    const std::vector<NcAttribute>& GetAttributes() const;
    void WriteVariableData(const int ncid, const int var_id) const;

    std::vector<NcDimension> GetDimensions() const;
    const NcGeneralVariable& GetVariable() const;
    template<typename T> NcSpecificVariable<T> DetachVariable();

private:
    NcGeneralVariable variable_;
    std::vector<NcAttribute> attributes_;
    std::vector<NcDimension> dimensions_;
    int nc_var_id{-1};
};

template<class T>
std::vector<NcDimension>
NcSpecificVariable<T>::GetDimensionsFromVariable() const 
{
    std::vector<NcDimension> nc_dims;

    switch(format_)
    {
        case DataFormat::LinearFormat:
            /* For example SFC indices */
            nc_dims.emplace_back("lin_index", global_domain_.GetNumberReferenceCoordsCovered());
            nc_dims.back().AppendNumToName(id_);
        break;
        default:
            /* Other representations */
            std::vector<Dimension> dimensions = GetDimensionVectorFromLayout(layout_);
            for (auto dim_iter = dimensions.begin(); dim_iter != dimensions.end(); ++dim_iter)
            {
                nc_dims.emplace_back(GetDimensionName(*dim_iter), global_domain_.GetDimensionLength(*dim_iter));
                nc_dims.back().AppendNumToName(id_);
            }
    }

    return nc_dims;
}

template<class T>
const std::string&
NcSpecificVariable<T>::GetName() const
{
    return name_;
}

template<class T>
int
NcSpecificVariable<T>::GetID() const
{
    return id_;
}

template<class T>
CmcType
NcSpecificVariable<T>::GetCmcType() const
{
    return ConvertToCmcType<T>();
}

template<class T>
const GeoDomain& 
NcSpecificVariable<T>::GetGlobalDomain() const
{
    return global_domain_;
}

template<class T>
void
NcSpecificVariable<T>::SetGlobalDomain(const GeoDomain& domain)
{
    global_domain_ = domain;
};

template<class T>
void
NcSpecificVariable<T>::SetGlobalDomain(GeoDomain&& domain)
{
    global_domain_ = std::move(domain);
};

template<class T>
void
NcSpecificVariable<T>::PushBack(const std::vector<T>& values, const Hyperslab& hyperslab)
{
    if (is_data_stack_full_)
    {
        cmc_err_msg("The data could not be push backed, since it would overwrite the previous set data.");
        return;
    }
    std::copy(values.begin(), values.end(), std::back_inserter(data_));
    hyperslabs_.push_back(hyperslab);

    format_ = DataFormat::HyperslabFormat;
};

template<class T>
void
NcSpecificVariable<T>::PushBack(const std::vector<T>& values, Hyperslab&& hyperslab)
{
    if (is_data_stack_full_)
    {
        cmc_err_msg("The data could not be push backed, since it would overwrite the previous set data.");
        return;
    }
    std::copy(values.begin(), values.end(), std::back_inserter(data_));
    hyperslabs_.push_back(std::move(hyperslab));

    format_ = DataFormat::HyperslabFormat;
};

template<class T>
void
NcSpecificVariable<T>::PushBack(std::vector<T>&& values, std::vector<Hyperslab>&& hyperslabs)
{
    if (is_data_stack_full_)
    {
        cmc_err_msg("The data could not be push backed, since it would overwrite the previous set data.");
        return;
    }

    data_ = std::move(values);
    hyperslabs_ = std::move(hyperslabs);

    format_ = DataFormat::HyperslabFormat;

    is_data_stack_full_ = true;
}


template<class T>
void
NcSpecificVariable<T>::PushBack(const std::vector<T>& values, const std::vector<Hyperslab>& hyperslabs)
{
    if (is_data_stack_full_)
    {
        cmc_err_msg("The data could not be push backed, since it would overwrite the previous set data.");
        return;
    }

    data_ = values;
    hyperslabs_ = hyperslabs;

    format_ = DataFormat::HyperslabFormat;

    is_data_stack_full_ = true;
}

template<class T>
void
NcSpecificVariable<T>::PushBack(std::vector<T>&& values, const int global_sfc_offset)
{
    if (is_data_stack_full_)
    {
        cmc_err_msg("The data could not be push backed, since it would overwrite the previous set data.");
        return;
    }

    data_ = std::move(values);
    global_sfc_offset_ = global_sfc_offset;
    format_ = DataFormat::LinearFormat;
    is_data_stack_full_ = true;   
}

template<class T>
void
NcSpecificVariable<T>::PushBack(const std::vector<T>& values, const int global_sfc_offset)
{
    if (is_data_stack_full_)
    {
        cmc_err_msg("The data could not be push backed, since it would overwrite the previous set data.");
        return;
    }

    data_ = values;
    global_sfc_offset_ = global_sfc_offset;
    format_ = DataFormat::LinearFormat;
    is_data_stack_full_ = true;   
}

template<class T>
void
NcSpecificVariable<T>::SetDataLayout(const DataLayout layout)
{
    layout_ = layout;
}

template<class T>
T
NcSpecificVariable<T>::GetMissingValue() const
{
    return missing_value_;
}

template<class T>
void
NcSpecificVariable<T>::SetMissingValue(const CmcUniversalType& missing_value)
{
    cmc_assert(std::holds_alternative<T>(missing_value));
    missing_value_ = std::get<T>(missing_value);
}

template<class T>
void
NcSpecificVariable<T>::SetMissingValue(const T missing_value)
{
    missing_value_ = missing_value;
}

template<class T>
void
NcSpecificVariable<T>::WriteVariableData(const int ncid, const int var_id) const
{
    switch(format_)
    {
        case DataFormat::LinearFormat:
            /* Emplace the Linear/SFC data */
            {
                const size_t start_val = global_sfc_offset_;
                const size_t count_val = data_.size();
                const int err = nc_put_vara(ncid, var_id, &start_val, & count_val, data_.data());
                NcCheckError(err);
            }
        break;
        case DataFormat::CartesianFormat:
            cmc_err_msg("The CartesianFormat is not yet supported.");
        break;
        case DataFormat::HyperslabFormat:
        {
            /* Iterate over all hyperslabs and emplace the data */
            const std::vector<Dimension> hs_dims = GetDimensionVectorFromLayout(layout_);
            const int dim = GetDimensionalityOfDataLayout(layout_);
            std::vector<size_t> start_ptr(dim, 0);
            std::vector<size_t> count_ptr(dim, 1);

            /* Keep track of the values that have already been written by previous hyperslabs */
            HyperslabIndex values_written = 0;

            /* Iterate over all local hyperslabs */
            for (auto hs_iter = hyperslabs_.begin(); hs_iter != hyperslabs_.end(); ++hs_iter)
            {
                /* Iterate over all diemnsions and fill the corresponding start and count ptr for writing the hyperslab data */
                int index = 0;
                for (auto hs_dim_iter = hs_dims.begin(); hs_dim_iter != hs_dims.end(); ++hs_dim_iter, ++index)
                {
                    start_ptr[index] = static_cast<size_t>(hs_iter->GetDimensionStart(*hs_dim_iter));
                    count_ptr[index] = static_cast<size_t>(hs_iter->GetDimensionLength(*hs_dim_iter));
                }

                /* Write the data corresponding to the hyperslab to the netCDF file */
                const int err = nc_put_vara(ncid, var_id, start_ptr.data(), count_ptr.data(), &data_[values_written]);
                NcCheckError(err);

                /* Add the number of written values to the offset variable */
                values_written += hs_iter->GetNumberCoordinates();
            }
        }
        break;
        default:
            cmc_err_msg("An unspecified or unknown DataFormat is stored within the variable.");
    }
}

template<class T>
void
NcSpecificVariable<T>::OverwriteDataFromFile(const int ncid, const int var_id, const GeneralHyperslab& hyperslab, const size_t start_offset)
{
    /* The size is already maxed out, since the data is not going to be push_backed, but just plain copied int the already present memory */
    cmc_assert(data_.size() >= hyperslab.GetNumberOfCoveredCoordinates() + start_offset);

    /* Overwrite the data with the variable_s data from the netCDF file */
    const int err = nc_get_vara(ncid, var_id, hyperslab.start_values.data(), hyperslab.count_values.data(), &data_[start_offset]);
    NcCheckError(err);
}

template<class T>
const std::vector<T>&
NcSpecificVariable<T>::GetData() const
{
    return data_;
}

nc_type
ConvertCmcTypeToNcType(const CmcType type);

NcGeneralVariable
CreateSpecificVariable(const nc_type type, const std::string& name, const int var_id, const size_t size_hint = 8);

template<class T>
void
NcVariable::SetSpecificVariable(const NcSpecificVariable<T>& variable)
{
    variable_ = variable;
}

template<class T> void
NcVariable::SetSpecificVariable(NcSpecificVariable<T>&& variable)
{
    variable_ = std::move(variable);
}

template<typename T>
NcSpecificVariable<T>
NcVariable::DetachVariable()
{
    NcSpecificVariable<T> var = std::get<NcSpecificVariable<T>>(std::move(variable_));

    variable_ = NcSpecificVariable<T>();

    return var;
}
}

#endif /* !CMC_NC_IO_HXX */
