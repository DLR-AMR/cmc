#include "t8code/cmc_t8_adapt.hxx"
#include "t8code/cmc_t8_adapt_callbacks.h"
#include "t8code/cmc_t8_mesh.hxx"

namespace cmc {

bool
AdaptData::IsCompressionProgressing() const
{
    /* If the the adpatation has not changed the number of elements, the compression has stagnated */
    return (previous_number_of_elements_ != new_number_of_elements_ ? true : false);
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
        return GetCurrentCompressionVariable().GetAmrMesh().GetMesh();
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
        return GetCurrentCompressionVariable().GetAmrMesh().GetInitialRefinementLevel();
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
        GetCurrentCompressionVariable().GetAmrMesh().SetMesh(forest);
    }
}

t8_forest_t
AdaptData::RepartitionData(t8_forest_t adapted_forest)
{
    /* Keep the adapted forest for the partition step */
    t8_forest_ref(adapted_forest);
    cmc_debug_msg("Forest in partition data hat num ocal elems: ", t8_forest_get_local_num_elements(adapted_forest));
    /* Repartition the forest */
    t8_forest_t partitioned_forest;
    t8_forest_init(&partitioned_forest);
    const int partition_for_coarsening = 0; //TODO: change to one in the future
    t8_forest_set_partition(partitioned_forest, adapted_forest, partition_for_coarsening);
    t8_forest_commit(partitioned_forest);

    if (corresponding_variable_id_ == kMeshCorrespondsToAllVariables)
    {   
        /* In a One For All compression mode */
        /* Repartition the variables */
        for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter)
        {
            var_iter->RepartitionData(adapted_forest, partitioned_forest);
        }
    } else
    {
        /* In a One For One compression mode */
        /* Repartition the variable */
        Var& var = GetCurrentCompressionVariable();
        var.RepartitionData(adapted_forest, partitioned_forest);
    }

    /* Deallocate the adapted forest */
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

void
AdaptData::UpdateCompressionData()
{
    Var& compression_variable = GetCurrentCompressionVariable();
    compression_variable.UpdateCompressionData();
}

void
AdaptData::FinalizeCompressionIteration()
{
    /* Store the number of elements obtained after the adaptation */
    new_number_of_elements_ = t8_forest_get_global_num_elements(GetCurrentMesh());
    ++count_adaptation_step_;
}

void
AdaptData::InitializeCompressionIteration()
{
    previous_number_of_elements_ = t8_forest_get_global_num_elements(GetCurrentMesh());
    /* Allocate the "..._new_" vectors which are used and filled during the compression */
    if (corresponding_variable_id_ == kMeshCorrespondsToAllVariables)
    {
        /* In a One For All compression mode */
        for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter)
        {
            var_iter->InitializeVariableForCompressionIteration();
        }
    } else
    {
        /* In a One For One compression mode */
        GetCurrentCompressionVariable().InitializeVariableForCompressionIteration();
    }
}

Var&
AdaptData::GetCurrentCompressionVariable()
{
    /* Search for the variable which is currently compressed */
    auto var_iter = std::find_if(variables_.begin(), variables_.end(), [this](const Var& var){return var.GetInternalID() == corresponding_variable_id_;});
    cmc_assert(var_iter != variables_.end());
    return *var_iter;
}

const Var&
AdaptData::GetCurrentCompressionVariable() const
{
    /* Search for the variable which is currently compressed */
    auto var_iter = std::find_if(variables_.begin(), variables_.end(), [this](const Var& var){return var.GetInternalID() == corresponding_variable_id_;});
    cmc_assert(var_iter != variables_.end());
    return *var_iter;
}

}
