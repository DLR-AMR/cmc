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
    /* Keep the adapted forest for the partition step */
    t8_forest_ref(adapted_forest);

    /* Repartition the forest */
    t8_forest_t partitioned_forest;
    t8_forest_init(&partitioned_forest);
    const int partition_for_coarsening = 0; //TODO: change to one in the future
    t8_forest_set_partition(partitioned_forest, adapted_forest, partition_for_coarsening);
    t8_forest_commit(partitioned_forest);

    /* Repartition the variables */
    for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter)
    {
        var_iter->RepartitionData(adapted_forest, partitioned_forest);
        
    }

    /* Deallocate the addapted forest */
    t8_forest_unref(&adapted_forest);
    cmc_debug_msg("End of AdaptData->RepartitionData");
    /* Return the repartitioned forest */
    return partitioned_forest;
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
