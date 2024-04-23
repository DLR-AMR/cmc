
#include "t8code/cmc_t8_byte_variable.hxx"
#include "t8code/cmc_t8_adapt.hxx"
#include "t8code/cmc_t8_adapt_callbacks.h"

namespace cmc
{


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
