#ifndef CMC_ADAPTIVE_COARSENING_HXX
#define CMC_ADAPTIVE_COARSENING_HXX

#include "lossy/cmc_ac_compression_variable.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "utilities/cmc_compression_settings.hxx"
#include "utilities/cmc_log_functions.hxx"

#include <stdexcept>

namespace cmc::lossy
{

template<typename T>
class DefaultAdaptData : public ICompressionAdaptData
{
public:
    DefaultAdaptData() = delete;
    DefaultAdaptData(AbstractCompressionVariable<T>* variable, const CompressionSettings& settings)
    : base_variable_{variable}, compression_settings_{settings} {};

    bool IsCompressionProgressing() override;
    void InitializeCompressionIteration() override;
    void FinalizeCompressionIteration() override;
    void CompleteInterpolation(const t8_forest_t previous_forest, const t8_forest_t adapted_forest) override;
    void RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest) override;
    cmc::t8::AdaptationFn GetAdaptationFunction() override;
    int EvaluateCoarsening(const int tree_id, const int elem_id, const int num_elements, const std::vector<PermittedError> permitted_errors) override;
    int EvaluateCoarsening([[maybe_unused]] const int tree_id, [[maybe_unused]] const int elem_id, [[maybe_unused]] const int num_elements, [[maybe_unused]] const std::vector<PermittedError> permitted_errors, [[maybe_unused]] const CmcUniversalType missing_value) override {return CMC_ERR;};
    int LeaveElementUnchanged(const int tree_id, const int local_elem_id) override;

protected:
    AbstractCompressionVariable<T>* base_variable_{nullptr};
    const CompressionSettings& compression_settings_;
    t8_gloidx_t previous_number_of_elements_{-1};
    t8_gloidx_t new_number_of_elements_{0};
    int count_adaptation_step_{0};
};

template <typename T>
inline ICompressionAdaptData*
CreateAdaptationClass(AbstractCompressionVariable<T>* abstract_var, const CompressionSettings& settings)
{
    return new DefaultAdaptData<T>(abstract_var, settings);
}

inline void
DestroyAdaptationClass(ICompressionAdaptData* iadapt_data)
{
    delete iadapt_data;
}

template<typename T>
inline t8_locidx_t
DefaultAdaptiveCoarsening (t8_forest_t forest,
                           t8_forest_t forest_from,
                           t8_locidx_t which_tree,
                           t8_locidx_t lelement_id,
                           t8_eclass_scheme_c * ts,
                           const int is_family,
                           const int num_elements,
                           t8_element_t * elements[])
{
    /* Retrieve the adapt_data */
    ICompressionAdaptData* _adapt_data = static_cast<ICompressionAdaptData*>(t8_forest_get_user_data(forest));
    cmc_assert(_adapt_data != nullptr);

    /* Cast the base-class pointer back to the derived type in order to use DefaultAdaptData<T>'s member functions */
    DefaultAdaptData<T>* iadapt_data = dynamic_cast<DefaultAdaptData<T>*>(_adapt_data);

    /* Check if a family is supplied to the adaptation function */
    if (is_family == 0)
    {
        /* If there is no family, the element stays unchanged */
        const int ret_val = iadapt_data->LeaveElementUnchanged(which_tree, lelement_id);
        return ret_val;
    }

    /* Get the permitted errors for the element */
    const std::vector<PermittedError> permitted_errors;

    /* Check if the coarsening complies to the error bound */
    const int ret_val = iadapt_data->EvaluateCoarsening(which_tree, lelement_id, num_elements, permitted_errors);
    
    return ret_val;
}

template<class T>
class CompressionVariable : public AbstractCompressionVariable<T>
{
public:
    CompressionVariable(const std::string& name, t8_forest_t initial_mesh, const std::vector<T>& variable_data)
    : AbstractCompressionVariable<T>()
    {
        if (static_cast<size_t>(t8_forest_get_local_num_elements(initial_mesh)) != variable_data.size())
        {
            throw std::invalid_argument("The number of local mesh elements does not match the amount of data elements.");
        }

        AbstractCompressionVariable<T>::name_ = name;
        AbstractCompressionVariable<T>::mesh_ = AmrMesh(initial_mesh);
        AbstractCompressionVariable<T>::data_ = variable_data;
        AbstractCompressionVariable<T>::adaptation_creator_ = CreateAdaptationClass<T>;
        AbstractCompressionVariable<T>::adaptation_destructor_ = DestroyAdaptationClass;
    };

private:

};

template <typename T>
bool
DefaultAdaptData<T>::IsCompressionProgressing()
{
    return (previous_number_of_elements_ != new_number_of_elements_ ? true : false);
}

template <typename T>
void
DefaultAdaptData<T>::InitializeCompressionIteration()
{
    previous_number_of_elements_ = t8_forest_get_global_num_elements(base_variable_->GetAmrMesh().GetMesh());
}

template <typename T>
void
DefaultAdaptData<T>::FinalizeCompressionIteration()
{
    new_number_of_elements_ = t8_forest_get_global_num_elements(base_variable_->GetAmrMesh().GetMesh());
    ++count_adaptation_step_;
}

template <typename T>
void
DefaultAdaptData<T>::CompleteInterpolation(const t8_forest_t previous_forest, const t8_forest_t adapted_forest)
{

}

template <typename T>
void
DefaultAdaptData<T>::RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest)
{

}

template <typename T>
cmc::t8::AdaptationFn
DefaultAdaptData<T>::GetAdaptationFunction()
{
    return DefaultAdaptiveCoarsening<T>;
}


template <typename T>
int
DefaultAdaptData<T>::EvaluateCoarsening(const int tree_id, const int elem_id, const int num_elements, const std::vector<PermittedError> permitted_errors)
{
    return base_variable_->EvaluateCoarsening(tree_id, elem_id, num_elements, permitted_errors);
}

template <typename T>
int
DefaultAdaptData<T>::LeaveElementUnchanged(const int tree_id, const int elem_id)
{
    return base_variable_->LeaveElementUnchanged(tree_id, elem_id);
}

}


#endif /* !CMC_ADAPTIVE_COARSENING_HXX */
