#include "netcdf/cmc_nc_io.hxx"

#include <algorithm>

namespace cmc::nc
{

const std::string&
Attribute::GetName() const
{
    return name_;
}

std::string
Attribute::GetName()
{
    return name_;
}

const CmcUniversalType&
Attribute::GetValue() const
{
    return value_;
}

CmcUniversalType
Attribute::GetValue()
{
    return value_;
}

const std::string&
Variable::GetName() const
{
    return std::visit([](auto&& var) -> const std::string& {
        return var.GetName();
    }, variable_);
}

std::vector<Dimension>
Variable::GetDimensionsFromVariable() const
{
    return std::visit([](auto&& var) -> std::vector<Dimension> {
        return var.GetDimensionsFromVariable();
    }, variable_);
}

const GeneralVariable&
Variable::GetVariable() const
{
    return variable_;
}

CmcType
Variable::GetCmcType() const
{
    return std::visit([](auto&& var) -> CmcType {
        return var.GetCmcType();
    }, variable_);
}

const std::vector<Attribute>&
Variable::GetAttributes() const
{
    return attributes_;
}

void
Variable::CreateIDAttribute()
{
    const int id = std::visit([](auto&& var) -> int {
                        return var.GetID();
                    }, variable_);
    attributes_.emplace_back("id", id);
}

void
Variable::WriteVariableData(const int ncid, const int var_id) const
{
    std::visit([&](auto&& var){
        var.WriteVariableData(ncid, var_id);
    }, variable_);
}

void
Variable::SetSpecificVariable(const GeneralVariable& variable)
{
    variable_ = variable;
}

void
Variable::SetSpecificVariable(GeneralVariable&& variable)
{
    variable_ = std::move(variable);
}

std::vector<Dimension>
Variable::GetDimensions() const
{
    return dimensions_;
}

void
Variable::SetupSpecificVariable(const std::string& var_name, const CmcType type)
{
    switch (type)
    {
        case CmcType::Int8_t:
            variable_ = SpecificVariable<int8_t>(var_name);
        break;
        case CmcType::Char:
            variable_ = SpecificVariable<char>(var_name);
        break;
        case CmcType::Int16_t:
            variable_ = SpecificVariable<int16_t>(var_name);
        break;
        case CmcType::Int32_t:
            variable_ = SpecificVariable<int32_t>(var_name);
        break;
        case CmcType::Float:
            variable_ = SpecificVariable<float>(var_name);
        break;
        case CmcType::Double:
            variable_ = SpecificVariable<double>(var_name);
        break;
        case CmcType::Uint8_t:
            variable_ = SpecificVariable<uint8_t>(var_name);
        break;
        case CmcType::Uint16_t:
            variable_ = SpecificVariable<uint16_t>(var_name);
        break;
        case CmcType::Uint32_t:
            variable_ = SpecificVariable<uint32_t>(var_name);
        break;
        case CmcType::Int64_t:
            variable_ = SpecificVariable<int64_t>(var_name);
        break;
        case CmcType::Uint64_t:
            variable_ = SpecificVariable<uint64_t>(var_name);
        break;
        default:
            cmc_err_msg("A not supported CmcType has been supplied.");
    }
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

CmcType
ConvertNcTypeToCmcType(const nc_type type)
{
    switch (type)
    {
        case NC_BYTE:
            return CmcType::Int8_t;
        break;
        case NC_CHAR:
            return CmcType::Char;
        break;
        case NC_SHORT:
            return CmcType::Int16_t;
        break;
        case NC_INT:
            return CmcType::Int32_t;
        break;
        case NC_FLOAT:
            return CmcType::Float;
        break;
        case NC_DOUBLE:
            return CmcType::Double;
        break;
        case NC_UBYTE:
            return CmcType::Uint8_t;
        break;
        case NC_USHORT:
            return CmcType::Uint16_t;
        break;
        case NC_UINT:
            return CmcType::Uint32_t;
        break;
        case NC_INT64:
            return CmcType::Int64_t;
        break;
        case NC_UINT64:
            return CmcType::Uint64_t;
        break;
        default:
            cmc_err_msg("A not supported nc_type has been supplied.");
            return CmcType::TypeUndefined;
    }
}

GeneralVariable
CreateSpecificVariable(const nc_type type, const std::string& name, const int var_id, const size_t size_hint)
{
    switch (type)
    {
        case NC_BYTE:
            return GeneralVariable{SpecificVariable<int8_t>(name, var_id, size_hint)};
        break;
        case NC_CHAR:
            return GeneralVariable{SpecificVariable<char>(name, var_id, size_hint)};
        break;
        case NC_SHORT:
            return GeneralVariable{SpecificVariable<int16_t>(name, var_id, size_hint)};
        break;
        case NC_INT:
            return GeneralVariable{SpecificVariable<int32_t>(name, var_id, size_hint)};
        break;
        case NC_FLOAT:
            return GeneralVariable{SpecificVariable<float>(name, var_id, size_hint)};
        break;
        case NC_DOUBLE:
            return GeneralVariable{SpecificVariable<double>(name, var_id, size_hint)};
        break;
        case NC_UBYTE:
            return GeneralVariable{SpecificVariable<uint8_t>(name, var_id, size_hint)};
        break;
        case NC_USHORT:
            return GeneralVariable{SpecificVariable<uint16_t>(name, var_id, size_hint)};
        break;
        case NC_UINT:
            return GeneralVariable{SpecificVariable<uint32_t>(name, var_id, size_hint)};
        break;
        case NC_INT64:
            return GeneralVariable{SpecificVariable<int64_t>(name, var_id, size_hint)};
        break;
        case NC_UINT64:
            return GeneralVariable{SpecificVariable<uint64_t>(name, var_id, size_hint)};
        break;
        default:
            cmc_err_msg("An unknown nc_type has been supplied.");
            return GeneralVariable{SpecificVariable<int32_t>(name, var_id, size_hint)};
    }
}

std::vector<Attribute>::const_iterator
FindAttribute(const std::vector<Attribute>& attributes, const std::string& attr_name)
{
    return std::find_if(attributes.begin(), attributes.end(), [&](const Attribute& attr){
        return !attr.GetName().compare(attr_name);
    });
}
}
