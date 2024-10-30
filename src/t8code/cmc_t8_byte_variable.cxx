
#include "t8code/cmc_t8_byte_variable.hxx"
#include "t8code/cmc_t8_adapt.hxx"
#include "t8code/cmc_t8_adapt_callbacks.h"

namespace cmc
{

void
ByteVar::SetUpByteVariable(const CmcType type, const std::string& name, const GeoDomain domain, const DataLayout layout, const DataLayout pre_compression_layout,
                           const int global_context_information, const CmcUniversalType missing_value)
{
    switch (type)
    {
        case CmcType::Int8_t:
            var_ = ByteVariable<int8_t>(name, domain, layout, pre_compression_layout, global_context_information, std::get<int8_t>(missing_value));
        break;
        case CmcType::Char:
            var_ = ByteVariable<char>(name, domain, layout, pre_compression_layout, global_context_information, std::get<char>(missing_value));
        break;
        case CmcType::Int16_t:
            var_ = ByteVariable<int16_t>(name, domain, layout, pre_compression_layout, global_context_information, std::get<int16_t>(missing_value));
        break;
        case CmcType::Int32_t:
            var_ = ByteVariable<int32_t>(name, domain, layout, pre_compression_layout, global_context_information, std::get<int32_t>(missing_value));
        break;
        case CmcType::Float:
            var_ = ByteVariable<float>(name, domain, layout, pre_compression_layout, global_context_information, std::get<float>(missing_value));
        break;
        case CmcType::Double:
            var_ = ByteVariable<double>(name, domain, layout, pre_compression_layout, global_context_information, std::get<double>(missing_value));
        break;
        case CmcType::Uint8_t:
            var_ = ByteVariable<uint8_t>(name, domain, layout, pre_compression_layout, global_context_information, std::get<uint8_t>(missing_value));
        break;
        case CmcType::Uint16_t:
            var_ = ByteVariable<uint16_t>(name, domain, layout, pre_compression_layout, global_context_information, std::get<uint16_t>(missing_value));
        break;
        case CmcType::Uint32_t:
            var_ = ByteVariable<uint32_t>(name, domain, layout, pre_compression_layout, global_context_information, std::get<uint32_t>(missing_value));
        break;
        case CmcType::Int64_t:
            var_ = ByteVariable<int64_t>(name, domain, layout, pre_compression_layout, global_context_information, std::get<int64_t>(missing_value));
        break;
        case CmcType::Uint64_t:
            var_ = ByteVariable<uint64_t>(name, domain, layout, pre_compression_layout, global_context_information, std::get<uint64_t>(missing_value));
        break;
        default:
            cmc_err_msg("An unknown data type has been supplied.");
    }
}

int
DetermineForestRefinementBits(std::vector<uint8_t>& serialized_variable, t8_forest_t forest)
{
    t8_forest_ref(forest);

    std::vector<std::vector<uint8_t>> serialized_forest_refinements;

    t8_forest_t forest_adapt = nullptr;

    int num_bytes = 0;

    while (t8_forest_get_local_num_elements(forest) > 1)
    {
        RefinementBits adapt_data(t8_forest_get_local_num_elements(forest));

        forest_adapt = t8_forest_new_adapt(forest, FindRefinementBits, 0, 0, static_cast<void*>(&adapt_data));

        forest = forest_adapt;

        serialized_forest_refinements.push_back(std::move(adapt_data.refinement_indicator));

        num_bytes += serialized_forest_refinements.back().size();
    }

    if (forest_adapt != nullptr)
    {
        t8_forest_unref(&forest_adapt);
    }

    for (auto lr_iter = serialized_forest_refinements.rbegin(); lr_iter != serialized_forest_refinements.rend(); ++lr_iter)
    {
        std::copy(lr_iter->begin(), lr_iter->end(), std::back_inserter(serialized_variable));
    }

    return num_bytes;
}

}
