#include "netcdf/cmc_nc_io.hxx"
#include "utilities/cmc_log_functions.h"

namespace cmc
{

const std::string&
NcAttribute::GetName() const
{
    return name_;
}

std::string
NcAttribute::GetName()
{
    return name_;
}

const CmcUniversalType&
NcAttribute::GetValue() const
{
    return value_;
}

CmcUniversalType
NcAttribute::GetValue()
{
    return value_;
}

const std::string&
NcVariable::GetName() const
{
    return std::visit([](auto&& var) -> const std::string& {
        return var.GetName();
    }, variable_);
}

std::vector<NcDimension>
NcVariable::GetDimensionsFromVariable() const
{
    return std::visit([](auto&& var) -> std::vector<NcDimension> {
        return var.GetDimensionsFromVariable();
    }, variable_);
}

CmcType
NcVariable::GetCmcType() const
{
    return std::visit([](auto&& var) -> CmcType {
        return var.GetCmcType();
    }, variable_);
}

const std::vector<NcAttribute>&
NcVariable::GetAttributes() const
{
    return attributes_;
}

void
NcVariable::CreateIDAttribute()
{
    const int id = std::visit([](auto&& var) -> int {
                        return var.GetID();
                    }, variable_);
    attributes_.emplace_back("id", id);
}

constexpr nc_type
ConvertCmcTypeToNcType(const CmcType type)
{
    enum CmcType {TypeUndefined = -1, Int8_t, Char, Int16_t, Int32_t, Float, Double, Uint8_t, Uint16_t, Uint32_t, Int64_t, Uint64_t, NumTypes};
    switch (type)
    {
        case CmcType::Int8_t:
            return NC_BYTE;
        break;
        case CmcType::Char:
            return NC_CHAR;
        break;
        case CmcType::Int16_t:
            return NC_SHORT;
        break;
        case CmcType::Int32_t:
            return NC_INT;
        break;
        case CmcType::Float:
            return NC_FLOAT;
        break;
        case CmcType::Double:
            return NC_DOUBLE;
        break;
        case CmcType::Uint8_t:
            return NC_UBYTE;
        break;
        case CmcType::Uint16_t:
            return NC_USHORT;
        break;
        case CmcType::Uint32_t:
            return NC_UINT;
        break;
        case CmcType::Int64_t:
            return NC_INT64;
        break;
        case CmcType::Uint64_t:
            return NC_UINT64;
        break;
        default:
            cmc_err_msg("A not supported CmcType has been supplied.");
            return NC_INT;
        
    }
}

}
