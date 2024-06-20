#ifndef CMC_NC_IO_HXX
#define CMC_NC_IO_HXX

#include "cmc.h"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_hyperslab.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_log_functions.h"
#include "netcdf/cmc_netcdf.hxx"

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

class NcDimension
{
public:
    NcDimension(const std::string& name, const size_t length)
    : name_{name}, length_{length} {};

    std::string GetName(){return name_;};
    const std::string& GetName() const {return name_;};
    size_t GetLength() const {return length_;};

    void AppendNumToName(const int number){name_.append(std::to_string(number));};

private:
    std::string name_;
    size_t length_;
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

template<typename T>
class NcSpecificVariable
{
public:
    NcSpecificVariable() = default;
    NcSpecificVariable(const int id)
    : id_{id} {};
    ~NcSpecificVariable() = default;

    NcSpecificVariable(const NcSpecificVariable& other) = default;
    NcSpecificVariable& operator=(const NcSpecificVariable& other) = default;
    NcSpecificVariable(NcSpecificVariable&& other) = default;
    NcSpecificVariable& operator=(NcSpecificVariable&& other) = default;

    const std::string& GetName() const;
    std::vector<NcDimension> GetDimensionsFromVariable() const;
    CmcType GetCmcType() const;
    int GetID() const;
    void WriteVariableData(const int ncid, const int var_id) const;
private:
    std::string name_;
    int id_;
    std::vector<T> data_;
    T missing_value_;
    DataLayout layout_;
    DataFormat format_;
    std::vector<Hyperslab> hyperslabs_;
    std::vector<LinearIndex> linear_indices_;
    GeoDomain global_domain_;
};

using NcGeneralVariable = std::variant<NcSpecificVariable<int8_t>, NcSpecificVariable<char>, NcSpecificVariable<int16_t>, NcSpecificVariable<int32_t>, NcSpecificVariable<float>, NcSpecificVariable<double>,
                                       NcSpecificVariable<uint8_t>, NcSpecificVariable<uint16_t>, NcSpecificVariable<uint32_t>, NcSpecificVariable<int64_t>, NcSpecificVariable<uint64_t>>;


class NcVariable
{
public:
    NcVariable() = default;
    ~NcVariable() = default;

    NcVariable(const NcVariable& other) = default;
    NcVariable& operator=(const NcVariable& other) = default;
    NcVariable(NcVariable&& other) = default;
    NcVariable& operator=(NcVariable&& other) = default;
    
    std::vector<NcDimension> GetDimensionsFromVariable() const;
    const std::string& GetName() const;
    CmcType GetCmcType() const;
    void CreateIDAttribute();
    const std::vector<NcAttribute>& GetAttributes() const;
    void WriteVariableData(const int ncid, const int var_id) const;
private:
    NcGeneralVariable variable_;
    std::vector<NcAttribute> attributes_;
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
void
NcSpecificVariable<T>::WriteVariableData(const int ncid, const int var_id) const
{
    switch(format_)
    {
        case DataFormat::LinearFormat:
            /* For example SFC indices */
            //TODO: Write the linearized data out, continue...
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

            /* Keep track of the values that hav ealready been written by previous hyperslabs */
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

constexpr nc_type
ConvertCmcTypeToNcType(const CmcType type);

}

#endif /* !CMC_NC_IO_HXX */
