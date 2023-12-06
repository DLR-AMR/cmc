#include "t8code/cmc_t8_adapt.hxx"
#include "t8code/cmc_t8_mesh.hxx"

namespace cmc {

inline bool
AdaptData::IsCompressionProgressing() const
{
    /* If the the adpatation has not changed the number of elements, the compression has stagnated */
    return (previous_number_of_elements_ != new_number_of_elements_ ? true : false);
}

inline t8_forest_t
AdaptData::GetCurrentMesh() const
{
    cmc_assert(!variables->empty());

    if (corresponding_variable_id_ == kMeshCorrespondsToAllVariables)
    {   
        return variables_->front().GetAmrMesh().GetMesh();
    } else
    {
        return variables_[corresponding_variable_id].GetAmrMesh().GetMesh();
    }
}

inline void
AdaptData::SetCurrentMesh(t8_forest_t forest)
{
    cmc_assert(!variables->empty());

    if (corresponding_variable_id_ == kMeshCorrespondsToAllVariables)
    {   
        for (auto iter = variables->begin(); iter != variables->end(); ++iter)
        {
            iter->GetAmrMesh().SetMesh(forest);
        }
    } else
    {
        variables_[corresponding_variable_id].GetAmrMesh().SetMesh(forest);
    }
}

void
AdaptData::InterpolateData(t8_forest_t coarsened_forest)
{
    cmc_assert(coarsened_forest != nullptr);

    /* We are storing the global number of elements (which is needed in order to determine whether the compression continues) */
    previous_number_of_elements_ = new_number_of_elements_;
    new_number_of_elements_ = t8_forest_get_global_num_elements(coarsened_forest);

    //...
}

t8_forest_t
AdaptData::RepartitionData(t8_forest_t adapted_forest)
{
    /* Repartition the forest */
    
    /* Repartition the variables */

    /* Return the repartitioned forest */
    
}

}