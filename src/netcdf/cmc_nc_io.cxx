#include "netcdf/cmc_nc_io.hxx"

#include <algorithm>

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

const NcGeneralVariable&
NcVariable::GetVariable() const
{
    return variable_;
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

void
NcVariable::WriteVariableData(const int ncid, const int var_id) const
{
    std::visit([&](auto&& var){
        var.WriteVariableData(ncid, var_id);
    }, variable_);
}

void
NcVariable::SetSpecificVariable(const NcGeneralVariable& variable)
{
    variable_ = variable;
}

void
NcVariable::SetSpecificVariable(NcGeneralVariable&& variable)
{
    variable_ = std::move(variable);
}

std::vector<NcDimension>
NcVariable::GetDimensions() const
{
    return dimensions_;
}

nc_type
ConvertCmcTypeToNcType(const CmcType type)
{
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

NcGeneralVariable
CreateSpecificVariable(const nc_type type, const std::string& name, const int var_id, const size_t size_hint)
{
    switch (type)
    {
        case NC_BYTE:
            return NcGeneralVariable{NcSpecificVariable<int8_t>(name, var_id, size_hint)};
        break;
        case NC_CHAR:
            return NcGeneralVariable{NcSpecificVariable<char>(name, var_id, size_hint)};
        break;
        case NC_SHORT:
            return NcGeneralVariable{NcSpecificVariable<int16_t>(name, var_id, size_hint)};
        break;
        case NC_INT:
            return NcGeneralVariable{NcSpecificVariable<int32_t>(name, var_id, size_hint)};
        break;
        case NC_FLOAT:
            return NcGeneralVariable{NcSpecificVariable<float>(name, var_id, size_hint)};
        break;
        case NC_DOUBLE:
            return NcGeneralVariable{NcSpecificVariable<double>(name, var_id, size_hint)};
        break;
        case NC_UBYTE:
            return NcGeneralVariable{NcSpecificVariable<uint8_t>(name, var_id, size_hint)};
        break;
        case NC_USHORT:
            return NcGeneralVariable{NcSpecificVariable<uint16_t>(name, var_id, size_hint)};
        break;
        case NC_UINT:
            return NcGeneralVariable{NcSpecificVariable<uint32_t>(name, var_id, size_hint)};
        break;
        case NC_INT64:
            return NcGeneralVariable{NcSpecificVariable<int64_t>(name, var_id, size_hint)};
        break;
        case NC_UINT64:
            return NcGeneralVariable{NcSpecificVariable<uint64_t>(name, var_id, size_hint)};
        break;
        default:
            cmc_err_msg("An unknown nc_type has been supplied.");
            return NcGeneralVariable{NcSpecificVariable<int32_t>(name, var_id, size_hint)};
    }
}

std::vector<NcAttribute>::const_iterator
FindAttribute(const std::vector<NcAttribute>& attributes, const std::string& attr_name)
{
    return std::find_if(attributes.begin(), attributes.end(), [&](const NcAttribute& attr){
        return !attr.GetName().compare("attr");
    });
}
}
