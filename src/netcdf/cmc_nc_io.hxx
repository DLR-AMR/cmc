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
    Dimension() = default;
    Dimension(const std::string& name, const size_t length)
    : name_{name}, length_{length} {}
    Dimension(std::string&& name, const size_t length)
    : name_{std::move(name)}, length_{length} {}
    Dimension(const std::string& name, const size_t length, const int dim_id)
    : name_{name}, length_{length}, nc_dim_id{dim_id} {}
    Dimension(std::string&& name, const size_t length, const int dim_id)
    : name_{std::move(name)}, length_{length}, nc_dim_id{dim_id} {}

    std::string GetName(){return name_;}
    const std::string& GetName() const {return name_;}
    size_t GetLength() const {return length_;}

    void AppendNumToName(const int number){name_.append(std::to_string(number));}

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
    : name_{name}, value_{value}{}
    Attribute(std::string&& name, const CmcUniversalType& value)
    : name_{std::move(name)}, value_{value}{}
    template<typename T> Attribute(const std::string& name, T value)
    : name_{name}, value_{value}{}
    template<typename T> Attribute(std::string&& name, T value)
    : name_{std::move(name)}, value_{value}{}
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
    : start_values{std::move(start_vals)}, count_values{std::move(count_vals)}{}

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

template <typename T>
struct DataSlab
{
    DataSlab() = default;
    DataSlab(const std::vector<T>& data_, const uint64_t global_byte_offset_)
    : data(data_), global_byte_offset{global_byte_offset_} {};
    DataSlab(std::vector<T>&& data_, const uint64_t global_byte_offset_)
    : data(std::move(data_)), global_byte_offset{global_byte_offset_} {};

    std::vector<T> data;
    uint64_t global_byte_offset{0};
};

template<typename T>
class SpecificVariable
{
public:
    SpecificVariable() = default;
    SpecificVariable(const int id)
    : id_{id} {}
    SpecificVariable(const std::string& name)
    : name_{name} {}
    SpecificVariable(const std::string& name, const int id)
    : name_{name}, id_{id} {}
    SpecificVariable(const std::string& name, const int id, const size_t size_hint)
    : name_{name}, id_{id} {
        this->Reserve(size_hint);
    }
    ~SpecificVariable() = default;

    SpecificVariable(const SpecificVariable& other) = default;
    SpecificVariable& operator=(const SpecificVariable& other) = default;
    SpecificVariable(SpecificVariable&& other) = default;
    SpecificVariable& operator=(SpecificVariable&& other) = default;

    const std::string& GetName() const;

    void SetGlobalDimensionLength(const uint64_t global_dim_length);
    std::vector<Dimension> GetDimensionsFromVariable() const;

    CmcType GetCmcType() const;
    int GetID() const;

    void Reserve(const size_t num_data_slabs);

    void SetData(const std::vector<T>& values, const uint64_t global_byte_offset);
    void SetData(std::vector<T>&& values, const uint64_t global_byte_offset);

    //TODO:Rework
    void OverwriteDataFromFile(const int ncid, const int var_id, const GeneralHyperslab& hyperslab, const size_t start_offset = 0);

    const std::vector<DataSlab<T>>& GetData() const;

    void WriteVariableData(const int ncid, const int var_id) const;
private:
    std::string name_;
    int id_;
    Dimension dimension_;
    std::vector<DataSlab<T>> data_;
};

using GeneralVariable = std::variant<SpecificVariable<int8_t>, SpecificVariable<char>, SpecificVariable<int16_t>, SpecificVariable<int32_t>, SpecificVariable<float>, SpecificVariable<double>,
                                       SpecificVariable<uint8_t>, SpecificVariable<uint16_t>, SpecificVariable<uint32_t>, SpecificVariable<int64_t>, SpecificVariable<uint64_t>>;


class Variable
{
public:
    Variable() = default;
    ~Variable() = default;

    template<class T> Variable(const SpecificVariable<T>& variable)
    : variable_{variable} {}
    template<class T> Variable(SpecificVariable<T>&& variable)
    : variable_{std::move(variable)} {}
    template<class T> Variable(const SpecificVariable<T>& variable, const std::vector<Attribute>& attributes)
    : variable_{variable}, attributes_{attributes} {}
    template<class T> Variable(const SpecificVariable<T>& variable, std::vector<Attribute>&& attributes)
    : variable_{variable}, attributes_{std::move(attributes)} {}
    template<class T> Variable(SpecificVariable<T>&& variable, const std::vector<Attribute>& attributes)
    : variable_{std::move(variable)}, attributes_{attributes} {}
    template<class T> Variable(SpecificVariable<T>&& variable, std::vector<Attribute>&& attributes)
    : variable_{std::move(variable)}, attributes_{std::move(attributes)} {}
    template<class T> Variable(const SpecificVariable<T>& variable, const Attribute& attribute)
    : variable_{variable} {
        attributes_.push_back(attribute);
    }
    template<class T> Variable(const SpecificVariable<T>& variable, Attribute&& attribute)
    : variable_{variable} {
        attributes_.push_back(attribute);
    }
    template<class T> Variable(SpecificVariable<T>&& variable, const Attribute& attribute)
    : variable_{std::move(variable)} {
        attributes_.push_back(attribute);
    }
    template<class T> Variable(SpecificVariable<T>&& variable, Attribute&& attribute)
    : variable_{std::move(variable)} {
        attributes_.push_back(attribute);
    }

    Variable(std::vector<Attribute>&& attributes, std::vector<Dimension>&& dimensions)
    : attributes_{std::move(attributes)}, dimensions_{std::move(dimensions)} {}

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
void
SpecificVariable<T>::SetGlobalDimensionLength(const uint64_t global_dim_length)
{
    dimension_ = Dimension(GetName() + "_lin_index", global_dim_length);
}

template<class T>
std::vector<nc::Dimension>
SpecificVariable<T>::GetDimensionsFromVariable() const 
{
    std::vector<nc::Dimension> nc_dims{this->dimension_};

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
SpecificVariable<T>::Reserve(const size_t num_data_slabs)
{
    data_.reserve(num_data_slabs);
}

template<class T>
void
SpecificVariable<T>::SetData(std::vector<T>&& values, const uint64_t global_byte_offset)
{
    data_.emplace_back(std::move(values), global_byte_offset);
}

template<class T>
void
SpecificVariable<T>::SetData(const std::vector<T>& values, const uint64_t global_byte_offset)
{
    data_.emplace_back(values, global_byte_offset);
}

template<class T>
void
SpecificVariable<T>::WriteVariableData(const int ncid, const int var_id) const
{
    for (auto data_slab_iter = data_.begin(); data_slab_iter != data_.end(); ++data_slab_iter)
    {
        const size_t start_val = static_cast<size_t>(data_slab_iter->global_byte_offset);
        const size_t count_val = data_slab_iter->data.size();
        const int err = nc_put_vara(ncid, var_id, &start_val, &count_val, data_slab_iter->data.data());
        CheckError(err);
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
const std::vector<DataSlab<T>>&
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
