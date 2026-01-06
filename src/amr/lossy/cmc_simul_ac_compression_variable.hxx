#ifndef CMC_SIMUL_AC_COMPRESSION_VARIABLE_HXX
#define CMC_SIMUL_AC_COMPRESSION_VARIABLE_HXX

#include "cmc_config.h"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_variable_utilities.hxx"
#include "utilities/cmc_compression_settings.hxx"
#include "utilities/cmc_iface_compression_adapt_data.hxx"

#ifdef CMC_WITH_T8CODE
#include <t8.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_scheme.hxx> 
#include <t8_forest/t8_forest_iterate.h> 
#endif

#include <vector>
#include <functional>
#include <string>

namespace cmc::lossy::simul
{

template <typename T>
class AbstractCompressionVariable;

template<typename T>
using AdaptCreator = std::function<ICompressionAdaptData*(AbstractCompressionVariable<T>*,const CompressionSettings&)>;

template<typename T>
using AdaptDestructor = std::function<void(ICompressionAdaptData*)>;


template <typename T>
struct Variable
{
    Variable(const std::string& name_, const std::vector<T>& data_)
    : name(name_), data(data_) {}

    std::string name;
    std::vector<T> data;
};

template <typename T>
class AbstractCompressionVariable
{
public:
    void Compress(const CompressionSettings& settings);

    void WriteVTKFile(const std::string& file_name);

    void SetVariable(const std::string& name, const std::vector<T>& data);
    void SetVariable(const Variable<T>& variable);

    size_t GetNumVariables() const {return variables_.size();}

    AmrMesh& GetAmrMesh() {return mesh_;};
    const AmrMesh& GetAmrMesh() const {return mesh_;};

    std::vector<VariableUtilities<T>>& GetVariableUtilities() {return utilities_;};
    const std::vector<VariableUtilities<T>>& GetVariableUtilities() const {return utilities_;};

    void SetInterpolation(Interpolate<T> interpolation_function);

    int EvaluateCoarsening(const int tree_id, const int elem_id, const int num_elements, const std::vector<PermittedError> permitted_errors);

    int LeaveElementUnchanged(const int tree_id, const int elem_id);

    virtual ~AbstractCompressionVariable(){};
protected:
    AbstractCompressionVariable() = default;

    std::vector<Variable<T>> variables_; //!< The actual data of the variable
    std::vector<std::vector<T>> data_new_; //!< A helper variable for the adaptation
    AmrMesh mesh_; //!< The mesh on which the variable is defined
    std::vector<VariableUtilities<T>> utilities_; //!< Utilities that are needed in order to track the errors
    AdaptCreator<T> adaptation_creator_; //!< A function pointer which is used to create the wished adaptation structure
    AdaptDestructor<T> adaptation_destructor_; //!< A function pointer which is used to destruct the adaptation structure

private:
    void AllocateCoarseningIteration();
    void SwitchToAdaptedData();
    ICompressionAdaptData* CreateAdaptData(const CompressionSettings& settings) {return adaptation_creator_(this, settings);};
    t8_forest_t RepartitionMesh(t8_forest_t adapted_forest);
};

template <typename T>
void
AbstractCompressionVariable<T>::AllocateCoarseningIteration()
{
    data_new_ = std::vector<std::vector<T>>(variables_.size());

    const size_t allocation_elems = mesh_.GetNumberLocalElements() / (2 << mesh_.GetDimensionality()) + 8;
    for (size_t var_id{0}; var_id < variables_.size(); ++var_id)
    {
        data_new_[var_id].reserve(allocation_elems);
    }
}

template <typename T>
void
AbstractCompressionVariable<T>::SwitchToAdaptedData()
{
    for (size_t var_id{0}; var_id < variables_.size(); ++var_id)
    {
        variables_[var_id].data.swap(data_new_[var_id]);
    }
    data_new_.clear();
}

template <typename T>
void
AbstractCompressionVariable<T>::SetInterpolation(Interpolate<T> interpolation_function)
{
    for (size_t var_id{0}; var_id < variables_.size(); ++var_id)
    {
        utilities_[var_id].SetInterpolation(interpolation_function);
    }
}

template <typename T>
void
AbstractCompressionVariable<T>::SetVariable(const std::string& name, const std::vector<T>& data)
{
    variables_.push_back(Variable<T>(name, data));
}

template <typename T>
void
AbstractCompressionVariable<T>::SetVariable(const Variable<T>& variable)
{
    variables_.push_back(variable);
}


template <typename T>
inline void
AbstractCompressionVariable<T>::Compress(const CompressionSettings& settings)
{
    cmc_debug_msg("Compression of ", variables_.size(), " variables (by the means of adaptive coarsening) starts...");
    cmc_debug_msg("Considered variables are: ");
    for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter)
    {
        cmc_debug_msg(var_iter->name);
    }

    /* Get the number of variables */
    const size_t num_variables = variables_.size();

    /* Set up the inaccuracy storage for all variables */
    utilities_ = std::vector<VariableUtilities<T>>(num_variables);
    for (size_t var_id{0}; var_id < num_variables; ++var_id)
    {
        utilities_[var_id].SetUpInaccuracyStorage(mesh_.GetNumberLocalElements());
    }

    /* We create the adapt data based on the compression settings, the forest and the variables to consider during the adaptation/coarsening */
    ICompressionAdaptData* adapt_data = this->CreateAdaptData(settings);

    while (adapt_data->IsCompressionProgressing())
    {
        cmc_debug_msg("A coarsening iteration is initialized.");

        /* Initialize/Allocate for a coarsening iteration*/
        this->AllocateCoarseningIteration();
        adapt_data->InitializeCompressionIteration();

        /* Get and indicate to keep the 'previous forest' after the adaptation step */
        t8_forest_t previous_forest = mesh_.GetMesh();
        t8_forest_ref(previous_forest);

        /* Perform a coarsening iteration */
        t8_forest_t adapted_forest = t8_forest_new_adapt(previous_forest, adapt_data->GetAdaptationFunction(), 0, 0, static_cast<void*>(adapt_data));
        cmc_debug_msg("The mesh adaptation step is finished; resulting in ", t8_forest_get_global_num_leaf_elements(adapted_forest), " global elements");

        /* Complete the interpolation step by storing the newly computed adapted data alongside its deviations */
        this->SwitchToAdaptedData();

        for (size_t var_id{0}; var_id < num_variables; ++var_id)
        {
            utilities_[var_id].SwitchDeviations();
        }

        adapt_data->CompleteInterpolation(previous_forest, adapted_forest);

        /* Free the former forest */
        t8_forest_unref(&previous_forest);

        /* Repartition the mesh */
        t8_forest_t partitioned_forest = RepartitionMesh(adapted_forest);

        /* Repartition the data */
        //this->RepartitionData(adapted_forest, partitioned_forest);
        //adapt_data->RepartitionData(adapted_forest, partitioned_forest);

        cmc_debug_msg("The mesh and the data has been re-partitioned.");

        /* Free the former forest and store the adapted/repartitioned mesh */
        t8_forest_unref(&adapted_forest);
        mesh_.SetMesh(partitioned_forest);

        /* Finalize the comrpession iteration */
        adapt_data->FinalizeCompressionIteration();

        cmc_debug_msg("The coarsening iteration is finished.");
    }

    /* Free the adapt data structure */
    this->adaptation_destructor_(adapt_data);
    cmc_debug_msg("Compression of variables has been finished.");
}


template <typename T>
inline t8_forest_t
AbstractCompressionVariable<T>::RepartitionMesh(t8_forest_t adapted_forest)
{
    /* Keep the not-partitioned forest */
    t8_forest_ref(adapted_forest);

    /* Allocate a forest */
    t8_forest_t partitioned_forest;
    t8_forest_init(&partitioned_forest);

    /* Partition the forest */
    const int partition_for_coarsening = 0; //TODO: change to 'one' when partition for coarsening is in t8code
    t8_forest_set_partition(partitioned_forest, adapted_forest, partition_for_coarsening);
    t8_forest_commit(partitioned_forest);

    return partitioned_forest;
}


template <typename T>
inline int
AbstractCompressionVariable<T>::EvaluateCoarsening(const int tree_id, const int local_elem_id, const int num_elements, const std::vector<PermittedError> permitted_errors)
{
    /* Compute the local element id in the data array */
    const int elem_id = t8_forest_get_tree_element_offset(mesh_.GetMesh(), tree_id) + local_elem_id;
    
    bool is_error_satisfied_for_all_variables = true;

    std::vector<T> interpolated_values;
    interpolated_values.reserve(variables_.size());

    std::vector<double> introduced_inaccuracy;
    introduced_inaccuracy.reserve(variables_.size());

    for (size_t var_id{0}; var_id != variables_.size(); ++var_id)
    {
        /* Obtain the values to be interpolated */
        const VectorView<T> values(variables_[var_id].data.data() + elem_id, num_elements);

        /* Interpolate the values */
        const T interpolated_value = utilities_[var_id].Interpolation(values, mesh_.GetMesh(), tree_id, elem_id, num_elements);

        /* Get the deviations of the elements */
        const std::vector<double> previous_deviations = utilities_[var_id].GetPreviousDeviations(elem_id, num_elements);

        /* Evaluate whether the interpolation is error-bound-conformal */
        const ErrorCompliance evaluation = utilities_[var_id].IsCoarseningErrorCompliant(permitted_errors, values, previous_deviations, interpolated_value);

        if (not evaluation.is_error_threshold_satisfied)
        {
            /* If the error is not satisfied, we cannot perform the coarsening */
            is_error_satisfied_for_all_variables = false;
            break;
        }

        interpolated_values.push_back(interpolated_value);
        introduced_inaccuracy.push_back(evaluation.max_introduced_error);
    }

    if (is_error_satisfied_for_all_variables)
    {
        /* If the error criterion is fullfilled for all variables */
        for (size_t var_id{0}; var_id < variables_.size(); ++var_id)
        {
            data_new_[var_id].push_back(interpolated_values[var_id]);
            utilities_[var_id].StoreInaccuracy(elem_id, introduced_inaccuracy[var_id]);
        }

        return t8::kCoarsenElements;
    } else
    {
        /* If the error criterion is NOT fullfilled for all variables */
        for (size_t var_id{0}; var_id < variables_.size(); ++var_id)
        {
            data_new_[var_id].push_back(variables_[var_id].data[elem_id]);
            utilities_[var_id].TransferPreviousDeviation(elem_id);
        }

        return t8::kLeaveElementUnchanged;
    }

}


template <typename T>
inline int
AbstractCompressionVariable<T>::LeaveElementUnchanged(const int tree_id, const int local_elem_id) 
{
    /* Compute the local element id in the data array */
    const int elem_id = t8_forest_get_tree_element_offset(mesh_.GetMesh(), tree_id) + local_elem_id;
    
    for (size_t var_id{0}; var_id < variables_.size(); ++var_id)
    {
        /* Get the value */
        const T value = variables_[var_id].data[elem_id];

        /* Store the unchanged value in the new data vector */
        data_new_[var_id].push_back(value);

        /* Transfer the previous inaccuracy */
        utilities_[var_id].TransferPreviousDeviation(elem_id);
    }

    return t8::kLeaveElementUnchanged;
}

template <typename T>
inline void
AbstractCompressionVariable<T>::WriteVTKFile(const std::string& file_name)
{
    const size_t num_variables = variables_.size();
    std::vector<std::vector<double>> double_data;
    double_data.reserve(num_variables);

    for (size_t var_id{0}; var_id < variables_.size(); ++var_id)
    {
        double_data.emplace_back();
        double_data.back().reserve(variables_[var_id].data.size());

        for (auto val_iter = variables_[var_id].data.begin(); val_iter != variables_[var_id].data.end(); ++val_iter)
        {
            double_data.back().push_back(static_cast<double>(*val_iter));
        }
    }

    t8_vtk_data_field_t vtk_data[num_variables];
    for (size_t var_id{0}; var_id < variables_.size(); ++var_id)
    {
        snprintf (vtk_data[var_id].description, BUFSIZ, variables_[var_id].name.c_str());
        vtk_data[var_id].type = T8_VTK_SCALAR;
        vtk_data[var_id].data = double_data[var_id].data();
    }

    t8_forest_write_vtk_ext (mesh_.GetMesh(), file_name.c_str(), 0, 0, 0, 0, 0, 0, 0, num_variables, vtk_data);

    cmc_debug_msg("The variables have been written to the file: ", file_name);
}


}


#endif /* !CMC_SIMUL_AC_COMPRESSION_VARIABLE_HXX */
