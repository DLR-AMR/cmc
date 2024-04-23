#include "t8code/cmc_t8_adapt.hxx"
#include "t8code/cmc_t8_adapt_callbacks.h"
#include "t8code/cmc_t8_mesh.hxx"

namespace cmc {

bool
AdaptData::IsCompressionProgressing() const
{
    /* If the the adpatation has not changed the number of elements, the compression has stagnated */
    return (previous_number_of_elements_ != new_number_of_elements_ ? true : false);
    //TODO: Remopve or let a flag control that when the compression should progress
    //return (previous_number_of_elements_ - new_number_of_elements_ <= 1000 ? false : true);
}

t8_forest_t
AdaptData::GetCurrentMesh() const
{
    cmc_assert(!variables_.empty());

    if (corresponding_variable_id_ == kMeshCorrespondsToAllVariables)
    {   
        return variables_.front().GetAmrMesh().GetMesh();
    } else
    {
        return variables_[corresponding_variable_id_].GetAmrMesh().GetMesh();
    }
}

int AdaptData::GetAdaptationStepCount() const
{
    return count_adaptation_step_;
}

int AdaptData::GetInitialRefinementLevelOfMesh() const
{
    cmc_assert(!variables_.empty());

    if (corresponding_variable_id_ == kMeshCorrespondsToAllVariables)
    {   
        return variables_.front().GetAmrMesh().GetInitialRefinementLevel();
    } else
    {
        return variables_[corresponding_variable_id_].GetAmrMesh().GetInitialRefinementLevel();
    }
}

void
AdaptData::SetCurrentMesh(t8_forest_t forest)
{
    cmc_assert(!variables_.empty());

    if (corresponding_variable_id_ == kMeshCorrespondsToAllVariables)
    {   
        for (auto iter = variables_.begin(); iter != variables_.end(); ++iter)
        {
            iter->GetAmrMesh().SetMesh(forest);
        }
    } else
    {
        variables_[corresponding_variable_id_].GetAmrMesh().SetMesh(forest);
    }
}

t8_forest_t
AdaptData::RepartitionData(t8_forest_t adapted_forest)
{
    //TODO:: implement

    /* Repartition the forest */
    
    /* Repartition the variables */

    /* Return the repartitioned forest */

    return adapted_forest;
}

t8_forest_adapt_t
AdaptData::GetAdaptationFunction() const
{
    /* Only return One FOr One or One for All function */
    if (mode_ == CompressionMode::OneForOne)
    {   
        return PerformAdaptiveCoarseningOneForOne;
    } else
    {
        cmc_err_msg("Not implemented yet");
        return PerformAdaptiveCoarseningOneForOne;
    }
}

}
