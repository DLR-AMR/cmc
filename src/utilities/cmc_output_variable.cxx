#include "utilities/cmc_output_variable.hxx"
#include "utilities/cmc_utilities.hxx"

namespace cmc
{

void
OutputVar::SetupOutputVar(const CmcType type, const int size_hint, const std::string& name, const DataLayout decompressed_layout, const GeoDomain& global_domain,const CmcUniversalType missing_value)
{
    switch(type)
    {
        case CmcType::Double:
        {
            OutputVariable<double> var(name, decompressed_layout, global_domain);
            var.AllocateData(size_hint);
            cmc_assert(std::holds_alternative<double>(missing_value));
            var.SetMissingValue(missing_value);
            var_ = CmcOutputVariable(std::move(var));
        }
        break;
        case CmcType::Float:
        {
            OutputVariable<float> var(name, decompressed_layout, global_domain);
            var.AllocateData(size_hint);
            cmc_assert(std::holds_alternative<float>(missing_value));
            var.SetMissingValue(missing_value);
            var_ = CmcOutputVariable(std::move(var));
        }
        break;
        case CmcType::Int32_t:
        {
            OutputVariable<int32_t> var(name, decompressed_layout, global_domain);
            var.AllocateData(size_hint);
            cmc_assert(std::holds_alternative<int32_t>(missing_value));
            var.SetMissingValue(missing_value);
            var_ = CmcOutputVariable(std::move(var));
        }
        break;
        case CmcType::Uint32_t:
        {
            OutputVariable<uint32_t> var(name, decompressed_layout, global_domain);
            var.AllocateData(size_hint);
            cmc_assert(std::holds_alternative<uint32_t>(missing_value));
            var.SetMissingValue(missing_value);
            var_ = CmcOutputVariable(std::move(var));
        }
        break;
        case CmcType::Int64_t:
        {
            OutputVariable<int64_t> var(name, decompressed_layout, global_domain);
            var.AllocateData(size_hint);
            cmc_assert(std::holds_alternative<int64_t>(missing_value));
            var.SetMissingValue(missing_value);
            var_ = CmcOutputVariable(std::move(var));
        }
        break;
        case CmcType::Uint64_t:
        {
            OutputVariable<uint64_t> var(name, decompressed_layout, global_domain);
            var.AllocateData(size_hint);
            cmc_assert(std::holds_alternative<uint64_t>(missing_value));
            var.SetMissingValue(missing_value);
            var_ = CmcOutputVariable(std::move(var));
        }
        break;
        case CmcType::Int16_t:
        {
            OutputVariable<int16_t> var(name, decompressed_layout, global_domain);
            var.AllocateData(size_hint);
            cmc_assert(std::holds_alternative<int16_t>(missing_value));
            var.SetMissingValue(missing_value);
            var_ = CmcOutputVariable(std::move(var));
        }
        break;
        case CmcType::Uint16_t:
        {
            OutputVariable<uint16_t> var(name, decompressed_layout, global_domain);
            var.AllocateData(size_hint);
            cmc_assert(std::holds_alternative<uint16_t>(missing_value));
            var.SetMissingValue(missing_value);
            var_ = CmcOutputVariable(std::move(var));
        }
        break;
        case CmcType::Int8_t:
        {
            OutputVariable<int8_t> var(name, decompressed_layout, global_domain);
            var.AllocateData(size_hint);
            cmc_assert(std::holds_alternative<int8_t>(missing_value));
            var.SetMissingValue(missing_value);
            var_ = CmcOutputVariable(std::move(var));
        }
        break;
        case CmcType::Uint8_t:
        {
            OutputVariable<uint8_t> var(name, decompressed_layout, global_domain);
            var.AllocateData(size_hint);
            cmc_assert(std::holds_alternative<uint8_t>(missing_value));
            var.SetMissingValue(missing_value);
            var_ = CmcOutputVariable(std::move(var));
        }
        break;
        default:
            cmc_err_msg("The supplied data type (", type, ") is not supported");
    }
}

}
