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

namespace cmc::nc
{

struct Dimension;
class Attribute;
template<typename T>
class SpecificVariable;
class Variable;

std::vector<Attribute>::const_iterator
FindAttribute(const std::vector<Attribute>& attributes, const std::string& attr_name);

class Dimension
{
public:
    Dimension(const std::string& name, const size_t length)
    : name_{name}, length_{length} {};
    Dimension(std::string&& name, const size_t length)
    : name_{std::move(name)}, length_{length} {};
    Dimension(const std::string& name, const size_t length, const int dim_id)
    : name_{name}, length_{length}, nc_dim_id{dim_id} {};
    Dimension(std::string&& name, const size_t length, const int dim_id)
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

class Attribute
{
public:
    Attribute() = default;
    Attribute(const std::string& name, const CmcUniversalType& value)
    : name_{name}, value_{value}{};
    Attribute(std::string&& name, const CmcUniversalType& value)
    : name_{std::move(name)}, value_{value}{};
    template<typename T> Attribute(const std::string& name, T value)
    : name_{name}, value_{value}{};
    template<typename T> Attribute(std::string&& name, T value)
    : name_{std::move(name)}, value_{value}{};
    ~Attribute() = default;

    Attribute(const Attribute& other) = default;
    Attribute& operator=(const Attribute& other) = default;
    Attribute(Attribute&& other) = default;
    Attribute& operator=(Attribute&& other) = default;

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
    GeneralHyperslab(std::vector<size_t>&& start_vals, std::vector<size_t>&& count_vals)
    : start_values{std::move(start_vals)}, count_values{std::move(count_vals)}{};

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
class SpecificVariable
{
public:
    SpecificVariable() = default;
    SpecificVariable(const int id)
    : id_{id} {};
    SpecificVariable(const std::string& name)
    : name_{name} {};
    SpecificVariable(const std::string& name, const int id)
    : name_{name}, id_{id} {};
    SpecificVariable(const std::string& name, const int id, const size_t size_hint)
    : name_{name}, id_{id} {
        data_ = std::vector<T>(size_hint);
    };
    ~SpecificVariable() = default;

    SpecificVariable(const SpecificVariable& other) = default;
    SpecificVariable& operator=(const SpecificVariable& other) = default;
    SpecificVariable(SpecificVariable&& other) = default;
    SpecificVariable& operator=(SpecificVariable&& other) = default;

    const std::string& GetName() const;
    std::vector<Dimension> GetDimensionsFromVariable() const;
    CmcType GetCmcType() const;
    int GetID() const;

    void Reserve(const size_t num_elements);

    const GeoDomain& GetGlobalDomain() const;
    void SetGlobalDomain(const GeoDomain& domain);
    void SetGlobalDomain(GeoDomain&& domain);

    void SetDataLayout(const DataLayout layout);

    void SetData(const std::vector<T>& values, const Hyperslab& hyperslab);
    void SetData(const std::vector<T>& values, Hyperslab&& hyperslab);
    void SetData(std::vector<T>&& values, std::vector<Hyperslab>&& hyperslabs);
    void SetData(const std::vector<T>& values, const std::vector<Hyperslab>& hyperslabs);
    void SetData(const std::vector<T>& values, const int global_sfc_offset);
    void SetData(std::vector<T>&& values, const int global_sfc_offset);

    void PushBack(const std::vector<T>& values);

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
    DataFormat format_{DataFormat::LinearFormat};
    std::vector<Hyperslab> hyperslabs_;
    LinearIndex global_sfc_offset_{0};
    GeoDomain global_domain_;
    bool is_data_stack_full_{false};
};

using GeneralVariable = std::variant<SpecificVariable<int8_t>, SpecificVariable<char>, SpecificVariable<int16_t>, SpecificVariable<int32_t>, SpecificVariable<float>, SpecificVariable<double>,
                                       SpecificVariable<uint8_t>, SpecificVariable<uint16_t>, SpecificVariable<uint32_t>, SpecificVariable<int64_t>, SpecificVariable<uint64_t>>;


class Variable
{
public:
    Variable() = default;
    ~Variable() = default;

    template<class T> Variable(const SpecificVariable<T>& variable)
    : variable_{variable} {};
    template<class T> Variable(SpecificVariable<T>&& variable)
    : variable_{std::move(variable)} {};
    template<class T> Variable(const SpecificVariable<T>& variable, const std::vector<Attribute>& attributes)
    : variable_{variable}, attributes_{attributes} {};
    template<class T> Variable(const SpecificVariable<T>& variable, std::vector<Attribute>&& attributes)
    : variable_{variable}, attributes_{std::move(attributes)} {};
    template<class T> Variable(SpecificVariable<T>&& variable, const std::vector<Attribute>& attributes)
    : variable_{std::move(variable)}, attributes_{attributes} {};
    template<class T> Variable(SpecificVariable<T>&& variable, std::vector<Attribute>&& attributes)
    : variable_{std::move(variable)}, attributes_{std::move(attributes)} {};
    template<class T> Variable(const SpecificVariable<T>& variable, const Attribute& attribute)
    : variable_{variable} {
        attributes_.push_back(attribute);
    };
    template<class T> Variable(const SpecificVariable<T>& variable, Attribute&& attribute)
    : variable_{variable} {
        attributes_.push_back(attribute);
    };
    template<class T> Variable(SpecificVariable<T>&& variable, const Attribute& attribute)
    : variable_{std::move(variable)} {
        attributes_.push_back(attribute);
    };
    template<class T> Variable(SpecificVariable<T>&& variable, Attribute&& attribute)
    : variable_{std::move(variable)} {
        attributes_.push_back(attribute);
    };

    Variable(std::vector<Attribute>&& attributes, std::vector<Dimension>&& dimensions)
    : attributes_{std::move(attributes)}, dimensions_{std::move(dimensions)} {};

    Variable(const Variable& other) = default;
    Variable& operator=(const Variable& other) = default;
    Variable(Variable&& other) = default;
    Variable& operator=(Variable&& other) = default;
    
    template<class T> void SetSpecificVariable(const SpecificVariable<T>& variable);
    template<class T> void SetSpecificVariable(SpecificVariable<T>&& variable);
    void SetSpecificVariable(const GeneralVariable& variable);
    void SetSpecificVariable(GeneralVariable&& variable);

    void SetupSpecificVariable(const std::string& var_name, const CmcType type);

    std::vector<nc::Dimension> GetDimensionsFromVariable() const;
    const std::string& GetName() const;
    CmcType GetCmcType() const;
    void CreateIDAttribute();
    const std::vector<Attribute>& GetAttributes() const;
    void WriteVariableData(const int ncid, const int var_id) const;

    std::vector<Dimension> GetDimensions() const;
    const GeneralVariable& GetVariable() const;
    template<typename T> SpecificVariable<T> DetachVariable();

private:
    GeneralVariable variable_;
    std::vector<Attribute> attributes_;
    std::vector<Dimension> dimensions_;
    int nc_var_id{-1};
};

template<class T>
std::vector<nc::Dimension>
SpecificVariable<T>::GetDimensionsFromVariable() const 
{
    std::vector<nc::Dimension> nc_dims;

    switch(format_)
    {
        case DataFormat::LinearFormat:
            /* For example SFC indices */
            nc_dims.emplace_back(GetName() + "_lin_index", data_.size());
            nc_dims.back().AppendNumToName(id_);
        break;
        default:
            /* Other representations */
            std::vector<cmc::Dimension> dimensions = GetDimensionVectorFromLayout(layout_);
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
SpecificVariable<T>::GetName() const
{
    return name_;
}

template<class T>
int
SpecificVariable<T>::GetID() const
{
    return id_;
}

template<class T>
CmcType
SpecificVariable<T>::GetCmcType() const
{
    return ConvertToCmcType<T>();
}

template<class T>
void
SpecificVariable<T>::Reserve(const size_t num_elements)
{
    data_.reserve(num_elements);
};

template<class T>
const GeoDomain& 
SpecificVariable<T>::GetGlobalDomain() const
{
    return global_domain_;
}

template<class T>
void
SpecificVariable<T>::SetGlobalDomain(const GeoDomain& domain)
{
    global_domain_ = domain;
};

template<class T>
void
SpecificVariable<T>::SetGlobalDomain(GeoDomain&& domain)
{
    global_domain_ = std::move(domain);
};

template<class T>
void
SpecificVariable<T>::SetData(const std::vector<T>& values, const Hyperslab& hyperslab)
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
SpecificVariable<T>::SetData(const std::vector<T>& values, Hyperslab&& hyperslab)
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
SpecificVariable<T>::SetData(std::vector<T>&& values, std::vector<Hyperslab>&& hyperslabs)
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
SpecificVariable<T>::PushBack(const std::vector<T>& values)
{
    std::copy_n(values.begin(), values.size(), std::back_insert_iterator(data_));
}

template<class T>
void
SpecificVariable<T>::SetData(const std::vector<T>& values, const std::vector<Hyperslab>& hyperslabs)
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
SpecificVariable<T>::SetData(std::vector<T>&& values, const int global_sfc_offset)
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
SpecificVariable<T>::SetData(const std::vector<T>& values, const int global_sfc_offset)
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
SpecificVariable<T>::SetDataLayout(const DataLayout layout)
{
    layout_ = layout;
}

template<class T>
T
SpecificVariable<T>::GetMissingValue() const
{
    return missing_value_;
}

template<class T>
void
SpecificVariable<T>::SetMissingValue(const CmcUniversalType& missing_value)
{
    cmc_assert(std::holds_alternative<T>(missing_value));
    missing_value_ = std::get<T>(missing_value);
}

template<class T>
void
SpecificVariable<T>::SetMissingValue(const T missing_value)
{
    missing_value_ = missing_value;
}

template<class T>
void
SpecificVariable<T>::WriteVariableData(const int ncid, const int var_id) const
{
    switch(format_)
    {
        case DataFormat::LinearFormat:
            /* Emplace the Linear/SFC data */
            {
                const size_t start_val = global_sfc_offset_;
                const size_t count_val = data_.size();
                const int err = nc_put_vara(ncid, var_id, &start_val, & count_val, data_.data());
                CheckError(err);
            }
        break;
        case DataFormat::CartesianFormat:
            cmc_err_msg("The CartesianFormat is not yet supported.");
        break;
        case DataFormat::HyperslabFormat:
        {
            /* Iterate over all hyperslabs and emplace the data */
            const std::vector<cmc::Dimension> hs_dims = GetDimensionVectorFromLayout(layout_);
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
                CheckError(err);

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
SpecificVariable<T>::OverwriteDataFromFile(const int ncid, const int var_id, const GeneralHyperslab& hyperslab, const size_t start_offset)
{
    /* The size is already maxed out, since the data is not going to be push_backed, but just plain copied int the already present memory */
    cmc_assert(data_.size() >= hyperslab.GetNumberOfCoveredCoordinates() + start_offset);

    /* Overwrite the data with the variable_s data from the netCDF file */
    const int err = nc_get_vara(ncid, var_id, hyperslab.start_values.data(), hyperslab.count_values.data(), &data_[start_offset]);
    CheckError(err);
}

template<class T>
const std::vector<T>&
SpecificVariable<T>::GetData() const
{
    return data_;
}

nc_type
ConvertCmcTypeToNcType(const CmcType type);

CmcType
ConvertNcTypeToCmcType(const nc_type type);

GeneralVariable
CreateSpecificVariable(const nc_type type, const std::string& name, const int var_id, const size_t size_hint = 8);

template<class T>
void
Variable::SetSpecificVariable(const SpecificVariable<T>& variable)
{
    variable_ = variable;
}

template<class T> void
Variable::SetSpecificVariable(SpecificVariable<T>&& variable)
{
    variable_ = std::move(variable);
}

template<typename T>
SpecificVariable<T>
Variable::DetachVariable()
{
    SpecificVariable<T> var = std::get<SpecificVariable<T>>(std::move(variable_));

    variable_ = SpecificVariable<T>();

    return var;
}

}

#endif /* !CMC_NC_IO_HXX */
