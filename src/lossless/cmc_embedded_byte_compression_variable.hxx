#ifndef LOSSLESS_CMC_EMBEDDED_BYTE_COMPRESSION_VARIABLE_HXX
#define LOSSLESS_CMC_EMBEDDED_BYTE_COMPRESSION_VARIABLE_HXX

#include "cmc_config.h"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_entropy_coder.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "utilities/cmc_compression_schema.hxx"
#include "utilities/cmc_embedded_variable_attributes.hxx"
#include "input/cmc_input_variable.hxx"
#include "mesh_compression/cmc_iface_embedded_mesh_encoder.hxx"
#include "utilities/cmc_embedded_mesh_utilities.hxx"

#include "mpi/cmc_mpi.hxx"
#include "t8code/cmc_t8_mpi.hxx"

#include <vector>
#include <functional>
#include <memory>
#include <algorithm>
#include <variant>

#include <type_traits>

#include <t8_forest/t8_forest_vtk.h>

namespace cmc::lossless::embedded
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

/**
 * @brief A strcut holding the data for an extraction process.
 * 
 * @tparam T The origianl data type of the underlying data (e.g. float)
 */
template<typename T>
struct ExtractionData
{
    ExtractionData(const CompressionValue<T>& coarse_val, std::vector<CompressionValue<T>>&& fine_vals)
    : coarse_value(coarse_val), fine_values(std::move(fine_vals)) {};
    ExtractionData(CompressionValue<T>&& coarse_val, std::vector<CompressionValue<T>>&& fine_vals)
    : coarse_value(std::move(coarse_val)), fine_values(std::move(fine_vals)) {};
    
    CompressionValue<T> coarse_value;
    std::vector<CompressionValue<T>> fine_values;
};

/**
 * @brief A strcut holding the data for an adapation process that leaves the element unchanged.
 * 
 * @tparam T The origianl data type of the underlying data (e.g. float)
 */
template<typename T>
struct UnchangedData
{
    UnchangedData(const CompressionValue<T>& coarse_val, const CompressionValue<T>& fine_val)
    : coarse_value(coarse_val), fine_value(fine_val) {};
    UnchangedData(CompressionValue<T>&& coarse_val, CompressionValue<T>&& fine_val)
    : coarse_value(std::move(coarse_val)), fine_value(std::move(fine_val)) {};
    
    CompressionValue<T> coarse_value;
    CompressionValue<T> fine_value;
};

/* Forward declarations */
template <typename T>
class AbstractEmbeddedByteCompressionVariable;
template <typename T>
class IEmbeddedCompressionAdaptData;

template<typename T>
using AdaptCreator = std::function<IEmbeddedCompressionAdaptData<T>*(AbstractEmbeddedByteCompressionVariable<T>*)>;

template<typename T>
using AdaptDestructor = std::function<void(IEmbeddedCompressionAdaptData<T>*)>;

/**
 * @brief The Interface/Template for a variable that performs lossless compression on the serialized data
 * in a byte-/bit-wise fashion. The compression algorithm is fixed and may be specialized in a derived class
 * with a derived adaptation data (\see IEmbeddedCompressionAdaptData) in order to fit the compression for the 
 * given needs.
 * 
 * @tparam T The origianl data type of the underlying data (e.g. float)
 */
template <typename T>
class AbstractEmbeddedByteCompressionVariable
{
public:
    void Compress();

    const std::string& GetName() const {return name_;};

    size_t Size() const {return data_.size();};

    const AmrMesh& GetAmrMesh() const {return mesh_;};
    const std::vector<CompressionValue<T>>& GetData() const {return data_;};

    virtual ~AbstractEmbeddedByteCompressionVariable(){};

    void MoveEncodedEntropyCodesInto(std::vector<std::vector<uint8_t>>& vec_to_hold_encoded_levelwise_entropy_codes)
    {
        vec_to_hold_encoded_levelwise_entropy_codes = std::move(buffered_entropy_codes_);
        buffered_entropy_codes_ = std::vector<std::vector<uint8_t>>();
    };

    void MoveEncodedDataInto(std::vector<std::vector<uint8_t>>& vec_to_hold_encoded_levelwise_data)
    {
        vec_to_hold_encoded_levelwise_data = std::move(buffered_encoded_data_);
        buffered_encoded_data_ = std::vector<std::vector<uint8_t>>();
    };

    void MoveEncodedMeshInto(std::vector<std::vector<uint8_t>>& vec_to_hold_encoded_levelwise_mesh)
    {
        vec_to_hold_encoded_levelwise_mesh = std::move(buffered_encoded_mesh_);
        buffered_encoded_mesh_ = std::vector<std::vector<uint8_t>>();
    };

    const std::vector<std::vector<uint8_t>>& GetEncodedEntropyCodes() const {return buffered_entropy_codes_;}
    const std::vector<std::vector<uint8_t>>& GetEncodedData() const {return buffered_encoded_data_;};
    const std::vector<std::vector<uint8_t>>& GetEncodedMesh() const {return buffered_encoded_mesh_;};

    bool AreMeshRefinementBitsStored() const {return store_refinement_indication_bits_;};
    MPI_Comm GetMPIComm() const {return comm_;};
    const VariableAttributes<T>& GetVariableAttributes() const {return attributes_;};

    virtual CompressionSchema GetCompressionSchema() const = 0;

    friend IEmbeddedCompressionAdaptData<T>;
protected:
    AbstractEmbeddedByteCompressionVariable() = delete;
    AbstractEmbeddedByteCompressionVariable(input::Var& input_variable)
    {
        this->SetupInputVariableForCompression(input_variable);
    };

    void SetName(const std::string& name) {name_ = name;};
    void SetAmrMesh(const AmrMesh& mesh) {mesh_ = mesh;};
    void SetAmrMesh(AmrMesh&& mesh) {mesh_ = std::move(mesh);};
    void SetData(const std::vector<T>& initial_data);
    void SetData(const std::vector<SerializedCompressionValue<sizeof(T)>>& initial_data);
    void SetData(std::vector<SerializedCompressionValue<sizeof(T)>>&& initial_data);
    void SetMPIComm(const MPI_Comm comm) {comm_ = comm;};
    void SetAttributes(const cmc::VariableAttributes<T>& attributes) {attributes_ = attributes;}
    void SetAttributes(cmc::VariableAttributes<T>&& attributes) {attributes_ = std::move(attributes);}
    void IndicateWhetherMeshRefinementBitsWillBeStored(const bool store_indication_bits) {store_refinement_indication_bits_ = true;}
    
    AdaptCreator<T> adaptation_creator_; //!< A function pointer which is used to create the wished adaptation structure
    AdaptDestructor<T> adaptation_destructor_; //!< A function pointer which is used to destruct the adaptation structure
    
    std::unique_ptr<mesh_compression::IEmbeddedMeshEncoder> mesh_encoder_{nullptr};
private:
    VectorView<CompressionValue<T>> GetView(const int start_index, const int count) const;
    VectorView<CompressionValue<T>> GetView(const int tree_id, const int lelement_index, const int count) const;

    void StoreExtractedValues(const int tree_id, const int lelement_id, const int num_elements, const ExtractionData<T>& extracted_values);
    void StoreUnchangedElement(const int tree_id, const int lelement_id, const UnchangedData<T>& extracted_value);

    void IndicateElementStaysUnchanged() {if (store_refinement_indication_bits_) {mesh_encoder_->IndicateElementStaysUnchanged();}};
    void IndicateCoarsening() {if (store_refinement_indication_bits_) {mesh_encoder_->IndicateCoarsening();}};

    bool IsValidForCompression() const;
    void AllocateExtractionIteration() {data_new_.reserve(mesh_.GetNumberLocalElements() / (2 << mesh_.GetDimensionality()) + 8);}
    void SwitchToExtractedData() {data_.swap(data_new_); data_new_.clear();};
    IEmbeddedCompressionAdaptData<T>* CreateAdaptData() {return adaptation_creator_(this);};
    t8_forest_t RepartitionMesh(t8_forest_t adapted_forest);
    void RepartitionData(t8_forest_t adapted_forest, t8_forest_t partitioned_forest);

    void SetupInputVariableForCompression(input::Var& input_variable);
    void DistributeDataOnInitialMesh(input::Variable<T>& input_variable);
    void SortLocalDataOnInitialMesh(input::Variable<T>& input_variable);
    std::vector<VariableRecvMessage> ReceiveInitialData(input::Variable<T>& input_variable);
    std::pair<std::vector<VariableSendMessage>, std::vector<MPI_Request>> SendInitialData(input::Variable<T>& input_variable);
    void SortInitialDataIntoVariables(input::Variable<T>& input_variable, const std::vector<VariableRecvMessage>& messages);
    std::vector<input::IndexReduction> UpdateLinearIndicesToTheInitialMesh();
    AmrMesh BuildInitialMesh(const input::Variable<T>& input_variable);


    void PerformTailTruncationOnInitialData(); //TODO: delete


    std::string name_; //!< The name of the variable
    AmrMesh mesh_; //!< The mesh on which the variable is defined
    std::vector<CompressionValue<T>> data_; //!< The current data of the variable 
    std::vector<CompressionValue<T>> data_new_; //!< A helper variable for the adaptation

    std::vector<std::vector<uint8_t>> buffered_entropy_codes_; //!< Level-wise storage of the entropy codes, e.g. LZC or prefix lengths
    std::vector<std::vector<uint8_t>> buffered_encoded_data_; //!< Level-wise storage for the encoded data
    std::vector<std::vector<uint8_t>> buffered_encoded_mesh_; //!< Level-wwise storage for the encoded mesh

    cmc::VariableAttributes<T> attributes_;
    bool store_refinement_indication_bits_{true};
    
    MPI_Comm comm_{MPI_COMM_NULL}; //!< The MPI communicator to use
};

/**
 * @brief Interface/Template for the adaptation data used within the lossless compression of the variable 
 * 
 * @tparam T The origianl data type of the underlying data (e.g. float)
 */
template <typename T>
class IEmbeddedCompressionAdaptData
{
public:
    IEmbeddedCompressionAdaptData() = delete;
    IEmbeddedCompressionAdaptData(AbstractEmbeddedByteCompressionVariable<T>* variable)
    : base_variable_{variable} {};

    bool IsCompressionProgressing() const;

    virtual void InitializeExtractionIteration() = 0;
    virtual void FinalizeExtractionIteration() = 0;
    virtual void CompleteExtractionIteration(const t8_forest_t previous_forest, const t8_forest_t adapted_forest) = 0;
    virtual void RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest) = 0;

    int ExtractValue(const int which_tree, const int lelement_id, const int num_elements);
    int LeaveElementUnchanged(const int which_tree, const int lelement_id);

    virtual std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeLevelData(const std::vector<CompressionValue<T>>& level_values) const = 0;
    virtual std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const = 0;
    
    MPI_Comm GetMPIComm() const {return base_variable_->GetMPIComm();};

    virtual ~IEmbeddedCompressionAdaptData(){};

    bool IsValidForCompression() const;
protected:
    virtual ExtractionData<T> PerformExtraction(const int which_tree, const int lelement_id, const int num_elements, const VectorView<CompressionValue<T>> values) = 0;
    virtual UnchangedData<T> ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value) = 0;

    std::unique_ptr<entropy_coding::IByteCompressionEntropyCoder> entropy_coder_{nullptr}; //!< The entropy coder to use in order to encode information
private:
    AbstractEmbeddedByteCompressionVariable<T>* const base_variable_{nullptr};
};

#if 1
//TODO: delete block
enum CompressionCriterion {CriterionUndefined, RelativeErrorThreshold, AbsoluteErrorThreshold};

struct PermittedError
{
    PermittedError() = delete;
    PermittedError(const CompressionCriterion etype, const double permitted_error)
    : criterion{etype}, error{permitted_error}{};

    const CompressionCriterion criterion{CompressionCriterion::CriterionUndefined};
    const double error{0.0};
};

struct ErrorCompliance
{
    ErrorCompliance() = delete;
    ErrorCompliance(const bool is_error_threshold_fulfilled, const double max_error)
    : is_error_threshold_satisfied{is_error_threshold_fulfilled}, max_introduced_error{max_error}{};

    const bool is_error_threshold_satisfied;
    const double max_introduced_error;
};

template<typename T>
auto
ComputeSingleRelativeDeviation(const T initial_value, const T nominal_value)
 -> std::enable_if_t<std::is_signed_v<T>, double>
{
    /* Calculate the relative deviation */
    return static_cast<double>(std::abs(initial_value - nominal_value)) / static_cast<double>(std::abs(initial_value));
}

template<typename T>
auto
ComputeSingleRelativeDeviation(const T initial_value, const T nominal_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>
{
    /* Calculate the relative deviation */
    return static_cast<double>((initial_value > nominal_value ? initial_value - nominal_value : nominal_value - initial_value)) / static_cast<double>(initial_value);
}

template<typename T>
auto
ComputeSingleAbsoluteDeviation(const T initial_value, const T nominal_value)
 -> std::enable_if_t<std::is_signed_v<T>, double>
{
    /* Calculate the absolute deviation */
    return static_cast<double>(std::abs(initial_value - nominal_value));
}

template<typename T>
auto
ComputeSingleAbsoluteDeviation(const T initial_value, const T nominal_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>
{
    /* Calculate the absolute deviation */
    return static_cast<double>((initial_value > nominal_value ? initial_value - nominal_value : nominal_value - initial_value));
}


template<typename T>
auto
ComputeSingleAbsoluteDeviationSkipMissingValues(const T initial_value, const T nominal_value, const T& missing_value)
 -> std::enable_if_t<std::is_signed_v<T>, double>
{
    if (!ApproxCompare(missing_value, initial_value))
    {
        /* Calculate the absolute deviation */
        return static_cast<double>(std::abs(initial_value - nominal_value));
    } else
    {
        return 0.0;
    }
}

template<typename T>
auto
ComputeSingleAbsoluteDeviationSkipMissingValues(const T initial_value, const T nominal_value, const T& missing_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>
{
    if (!ApproxCompare(missing_value, initial_value))
    {
        /* Calculate the absolute deviation */
        return static_cast<double>((initial_value > nominal_value ? initial_value - nominal_value : nominal_value - initial_value));
    } else
    {
        return 0.0;
    }
}

template<typename T>
auto
ComputeSingleRelativeDeviationSkipMissingValues(const T initial_value, const T nominal_value, const T& missing_value)
 -> std::enable_if_t<std::is_signed_v<T>, double>
{
    if (!ApproxCompare(missing_value, initial_value))
    {
        /* Calculate the relative deviation */
        return static_cast<double>(std::abs(initial_value - nominal_value)) / static_cast<double>(std::abs(initial_value));
    } else
    {
        return 0.0;
    }
}

template<typename T>
auto
ComputeSingleRelativeDeviationSkipMissingValues(const T initial_value, const T nominal_value, const T& missing_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>
{
    if (!ApproxCompare(missing_value, initial_value))
    {
        /* Calculate the relative deviation */
        return static_cast<double>((initial_value > nominal_value ? initial_value - nominal_value : nominal_value - initial_value)) / static_cast<double>(initial_value);
    } else
    {
        return 0.0;
    }
}

template<class T>
ErrorCompliance
IsValueErrorCompliant(const std::vector<PermittedError>& permitted_errors, const T initial_value, const T nominal_value, const T& missing_value)
{
    /* Get the absolute inaccuracy for all values (the absolute error is needed in any case) */
    const double abs_inaccuracy = ComputeSingleAbsoluteDeviationSkipMissingValues(initial_value, nominal_value, missing_value);

    for (auto pe_iter = permitted_errors.begin(); pe_iter != permitted_errors.end(); ++pe_iter)
    {
        switch (pe_iter->criterion)
        {
            case CompressionCriterion::AbsoluteErrorThreshold:
            {
                /* Check if it is compliant with the permitted error */
                if (abs_inaccuracy > pe_iter->error)
                {
                    return ErrorCompliance(false, 0.0);
                }
            }
            break;
            case CompressionCriterion::RelativeErrorThreshold:
            {
                /* Get the relative inaccuracy for all values */
                const double rel_inaccuracy = ComputeSingleRelativeDeviationSkipMissingValues(initial_value, nominal_value, missing_value);
                /* Check if it is compliant with the permitted error */
                if (rel_inaccuracy > pe_iter->error)
                {
                    return ErrorCompliance(false, 0.0);
                }
            }
            break;
                default:
                cmc_err_msg("The error specifications hold an unrecognized criterion.");
                return ErrorCompliance(false, 0.0);
        }
    }

    /* If the funciton reaches this point, the interpolated value complies with the permitted errors */
    return ErrorCompliance(true, abs_inaccuracy);
}

template<typename T>
CompressionValue<T>
GetMaximumTailToggledValue(const int index, const std::vector<PermittedError>& permitted_errors, const CompressionValue<T>& initial_serialized_value, const T& missing_value)
{
    bool is_toogling_progressing = true;
    CompressionValue<T> toggled_value = initial_serialized_value;

    const T initial_val = initial_serialized_value.template ReinterpretDataAs<T>();

    int iteration_count = 0;
    const int max_iteration_count = sizeof(T) * CHAR_BIT;

    while (is_toogling_progressing && iteration_count < max_iteration_count)
    {
        const CompressionValue<T> save_previous_value = toggled_value;

        /* Toggle all ones up until the next unset bit (inclusive) */
        toggled_value.ToggleTailUntilNextUnsetBit();
        const T reinterpreted_value = toggled_value.template ReinterpretDataAs<T>();

        /* Check if it is error compliant */
        const ErrorCompliance error_evaluation = IsValueErrorCompliant(permitted_errors, initial_val, reinterpreted_value, missing_value);

        if (!error_evaluation.is_error_threshold_satisfied)
        {
            /* Revert the changes to the value */
            toggled_value = save_previous_value;
            is_toogling_progressing = false;
        }

        ++iteration_count;
    }

    return toggled_value;
}


template<typename T>
CompressionValue<T>
GetMaximumTailClearedValue(const int index, const std::vector<PermittedError>& permitted_errors, const CompressionValue<T>& initial_serialized_value, const T& missing_value)
{
    bool is_clearing_progressing = true;
    CompressionValue<T> cleared_value = initial_serialized_value;

    const T initial_val = initial_serialized_value.template ReinterpretDataAs<T>();

    int iteration_count = 0;
    const int max_iteration_count = sizeof(T) * CHAR_BIT;

    while (is_clearing_progressing && iteration_count < max_iteration_count)
    {
        const CompressionValue<T> save_previous_value = cleared_value;

        /* Clear the next set bit from the tail */
        cleared_value.ClearNextSetBitFromTail();
        const T reinterpreted_value = cleared_value.template ReinterpretDataAs<T>();

        /* Check if it is error compliant */
        const ErrorCompliance error_evaluation = IsValueErrorCompliant(permitted_errors, initial_val, reinterpreted_value, missing_value);

        if (!error_evaluation.is_error_threshold_satisfied)
        {
            /* Revert the changes to the value */
            cleared_value = save_previous_value;
            is_clearing_progressing = false;
        }

        ++iteration_count;
    }

    return cleared_value;
}


template<typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::PerformTailTruncationOnInitialData()
{
    const T missing_value = attributes_.GetMissingValue();

    /* Iterate through the serialized values and try to emplace as many zeros at the tail as possible (compliant to the error threshold) */
    int index = 0;
    for (auto val_iter = data_.begin(); val_iter != data_.end(); ++val_iter, ++index)
    {
        const T re_val = val_iter->template ReinterpretDataAs<T>();

        if (!ApproxCompare(re_val, missing_value))
        {
            /* Get the permitted error for the current values */
            const std::vector<PermittedError> permitted_errors = std::vector<PermittedError>{PermittedError(CompressionCriterion::RelativeErrorThreshold, 0.01)};

            /* Get the value which has been transformed by toggling as many ones from the back while setting the succeeding 'zero' bits to one */
            const CompressionValue<T> toggled_value = GetMaximumTailToggledValue(index, permitted_errors, *val_iter, missing_value);

            /* Get the value which has been transformed by clearing as many of the last set bits as possible */
            const CompressionValue<T> cleared_value = GetMaximumTailClearedValue(index, permitted_errors, *val_iter, missing_value);

            /* Check which approach leads to more zero bits at the end */
            const int num_toogled_trailing_zeros = toggled_value.GetNumberTrailingZeros();
            const int num_cleared_trailing_zeros = cleared_value.GetNumberTrailingZeros();

            /* Replace the initial value with the transformed one */
            if (num_toogled_trailing_zeros >= num_cleared_trailing_zeros)
            {
                /* If the toggling approach has been more successfull */
                *val_iter = toggled_value;
            } else
            {
                /* If the clearing approach has been more successfull */
                *val_iter = cleared_value;
            }

            /* Update the trail bit count for the new value */
            val_iter->UpdateTailBitCount();
            //cmc_debug_msg("Trailing Zeros: Toggled: ", num_toogled_trailing_zeros, ", Cleared: ", num_cleared_trailing_zeros);
        } else
        {
            /* In order to not change missing values, we are just able to trim their trailing zeros */
            val_iter->UpdateTailBitCount();
        }
    }
}


#endif




template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::SetupInputVariableForCompression(input::Var& input_var)
{
    cmc_debug_msg("The InputVariable will be set up for compression.");
    //cmc_debug_msg("At the start: num coords:", input_variable.GetNumberCoordinates());
    /* get the MPI communicator from the variable */
    comm_ = input_var.GetMPIComm();

    /* Potentially, apply scaling values and offsets if defined (potentially, the datatype changes by this call) */
    input_var.ApplyScalingAndOffset();
    
    if (not std::holds_alternative<input::Variable<T>>(input_var.GetInternalVariant()))
    {
        cmc_err_msg("The data type of the input variable differs from the embedded byte comrpession variable.");
    }

    /* Get the actual variable from the input var wrapper */
    input::Variable<T> input_variable = std::get<input::Variable<T>>(input_var.GetInternalVariant());
    cmc_debug_msg("After sclaing and transform : num coords:", input_variable.GetNumberCoordinates());
    /* Generate the initial embedded mesh */
    mesh_ = this->BuildInitialMesh(input_variable);

    /* Distribute the data in the mesh */
    this->DistributeDataOnInitialMesh(input_variable);

    /* Now, we extract the data and the attributes */
    attributes_ = VariableAttributes<T>(input_variable.GetGlobalDomain(), input_variable.GetMissingValue(), input_variable.GetInitialDataLayout(), input_variable.GetPreCompressionDataLayout(),
                                        input_variable.GetGlobalContextInformation());


    #if 0
    /* Create a new vtk field holding the element data arrays */
    t8_vtk_data_field_t *vtk_data = new t8_vtk_data_field_t[1];

    cmc_debug_msg("Size of input data: ", input_variable.GetDataForReading().size());
    cmc_debug_msg("Local mesh elems: ", mesh_.GetNumberLocalElements(), ", init levl: ", mesh_.GetInitialRefinementLevel(), ", and dim: ", mesh_.GetDimensionality());
    /* Set the type of the data and pointer to the data */
    std::ignore = snprintf(vtk_data[0].description, 4, "t2m");
    vtk_data[0].type = T8_VTK_SCALAR;

    std::vector<float> init_data = input_variable.GetDataForReading();
    std::vector<double> converted_data;
    converted_data.reserve(init_data.size());

    for (auto val_iter = init_data.begin(); val_iter != init_data.end(); ++val_iter)
    {
        converted_data.push_back(static_cast<double>(*val_iter));
    }

    vtk_data[0].data = converted_data.data();
    t8_forest_t forest = mesh_.GetMesh();
    cmc_debug_msg("Vor vtk write: num local elems: ", t8_forest_get_local_num_elements(forest));
    const int vtk_err = t8_forest_vtk_write_file(forest, "example_t2m_input", 0, 1, 0, 0, 0, 1, vtk_data);
    
    if (vtk_err == 0)
        cmc_err_msg("An error occrued during the creation of the t8code-forest vtk file.");
    
    delete[] vtk_data;

    #endif





    /* Get the data from the variable and store it as CompressionValues */
    this->SetData(input_variable.GetDataForReading());
    
    cmc_debug_msg("The setup of the InputVariable for compression has been succcessfull.");

    /* Nur zum testen gerade; funktioniert nur für PrefixAMR */
    this->PerformTailTruncationOnInitialData();
}

template <typename T>
inline bool
IEmbeddedCompressionAdaptData<T>::IsValidForCompression() const
{
    if (base_variable_ == nullptr)
    {
        cmc_err_msg("The pointer to the base variable is not set. Therefore, no compression can be applied.");
        return false;
    }
    if (entropy_coder_ == nullptr)
    {
        cmc_err_msg("The entropy coder is not set. Therefore, no compression can be applied.");
        return false;
    }

    return true;
}

/**
 * @brief This function indicates whether the compression continues or not.
 * The compression is carried out, until only the root elements of the mesh
 * are left over and further coarsening can be applied.
 * 
 * @tparam T The data type of the underlying data (e.g. float)
 * @return true If the compression continues
 * @return false If the compression is finished
 */
template <typename T>
bool
IEmbeddedCompressionAdaptData<T>::IsCompressionProgressing() const
{
    /* Perform the extraction until only the root elements are left */
    return (base_variable_->GetAmrMesh().GetNumberGlobalElements() > base_variable_->GetAmrMesh().GetNumberGlobalTrees());
}

/**
 * @brief This funciton is called during the compression and extracts a value which will be stored on the coarser level
 * (i.e. after coarsening the family of elements). Moreover, the remaining values on the finer level may be adjusted
 * as well. This function calls the custom "PerformExtraction" which needs to be implemented by the derived class.
 * 
 * @tparam T The data type of the underlying data (e.g. float)
 * @param which_tree The lcoal tree id from which the elements are taken
 * @param lelement_id The tree-local start index of the family of elements
 * @param num_elements The number of elements corresponding to this family
 * @return int The return value indicates that this family of elements will be coarsened
 */
template <typename T>
int
IEmbeddedCompressionAdaptData<T>::ExtractValue(const int which_tree, const int lelement_id, const int num_elements)
{
    cmc_assert(which_tree >= 0 && lelement_id >= 0);
    cmc_assert(num_elements > 1);

    /* We indicate to that a coarsening is applied to the family of elements */
    base_variable_->IndicateCoarsening();

    /* Get the corresponding values */
    const VectorView<CompressionValue<T>> values = base_variable_->GetView(which_tree, lelement_id, num_elements);

    /* Extract the coarse values, and potentially alter the remaining fine values */
    const ExtractionData<T> extracted_values = PerformExtraction(which_tree, lelement_id, num_elements, values);

    /* Store the extracted values wihtin the variable */
    base_variable_->StoreExtractedValues(which_tree, lelement_id, num_elements, extracted_values);

    return cmc::t8::kCoarsenElements;
}

/**
 * @brief This function is called during the compression and supplies the (potential) altered value for the
 * element (which remains unchanged in the mmesh) after the adaptation which will be stored for the next adaptation
 * iteration. Moreover the left-over value remaining "in the old data vector" can be altered as well, if needed.
 * This function calls the custom "LeaveElementUnchanged" which needs to be implemented by the derived class.
 * 
 * @tparam T The data type of the underlying data
 * @param which_tree The lcoal tree id from which the element is taken
 * @param lelement_id The tree-local index of the element
 * @return int The return value indicates that this element will remain unchanged
 */
template <typename T>
int
IEmbeddedCompressionAdaptData<T>::LeaveElementUnchanged(const int which_tree, const int lelement_id)
{
    cmc_assert(which_tree >= 0 && lelement_id >= 0);

     /* We indicate to that the element stays unchanged */
    base_variable_->IndicateElementStaysUnchanged();

    /* Get the corresponding value */
    const VectorView<CompressionValue<T>> value = base_variable_->GetView(which_tree, lelement_id, 1);

    /* Leave the element unchanged */
    const UnchangedData<T> extracted_values = this->ElementStaysUnchanged(which_tree, lelement_id, value.front());

    /* Store the unchanged data */
    base_variable_->StoreUnchangedElement(which_tree, lelement_id, extracted_values);

    return cmc::t8::kLeaveElementUnchanged;
}

/**
 * @brief The adaptation function which is used for the lossless compression variables.
 * In case a family is passed to this callback, an extraction process is always performed.
 * 
 * @return t8_locidx_t Indicates whether the element stays unchanged or if the family of elements
 * will be coarsened
 */
template<typename T>
inline t8_locidx_t
LosslessByteCompression (t8_forest_t forest,
                         [[maybe_unused]] t8_forest_t forest_from,
                         t8_locidx_t which_tree,
                         t8_locidx_t lelement_id,
                         [[maybe_unused]] t8_eclass_scheme_c * ts,
                         const int is_family,
                         const int num_elements,
                         [[maybe_unused]] t8_element_t * elements[])
{
    /* Retrieve the adapt_data */
    IEmbeddedCompressionAdaptData<T>* adapt_data = static_cast<IEmbeddedCompressionAdaptData<T>*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    /* Check if a family is supplied to the adaptation function */
    if (is_family == 0)
    {
        /* If there is no family, the element stays unchanged */
        const int ret_val = adapt_data->LeaveElementUnchanged(which_tree, lelement_id);
        return ret_val;
    } else
    {
        /* Extract a value of the family and coarsen it */
        const int ret_val = adapt_data->ExtractValue(which_tree, lelement_id, num_elements);
        return ret_val;
    }
}

template <typename T>
inline void
AbstractEmbeddedByteCompressionVariable<T>::Compress()
{
    cmc_assert(this->IsValidForCompression());
    cmc_debug_msg("Lossless compression of variable ", this->name_, " starts...");

    /* We create the adapt data based on the compression settings, the forest and the variables to consider during the adaptation/coarsening */
    IEmbeddedCompressionAdaptData<T>* adapt_data = this->CreateAdaptData();

    cmc_assert(adapt_data->IsValidForCompression());

    while (adapt_data->IsCompressionProgressing())
    {
        cmc_debug_msg("A coarsening iteration is initialized.");

        /* Initialize/Allocate for a coarsening iteration */
        this->AllocateExtractionIteration();
        adapt_data->InitializeExtractionIteration();
        if (AreMeshRefinementBitsStored()){mesh_encoder_->IntializeCompressionIteration(mesh_.GetNumberLocalElements());}

        /* Get and indicate to keep the 'previous forest' after the adaptation step */
        t8_forest_t previous_forest = mesh_.GetMesh();
        t8_forest_ref(previous_forest);

        /* Perform a coarsening iteration */
        t8_forest_t adapted_forest = t8_forest_new_adapt(previous_forest, LosslessByteCompression<T>, 0, 0, static_cast<void*>(adapt_data));
        cmc_debug_msg("The mesh adaptation step is finished; resulting in ", t8_forest_get_global_num_elements(adapted_forest), " global elements");

        /* Complete the interpolation step by storing the newly computed adapted data alongside its deviations */
        adapt_data->CompleteExtractionIteration(previous_forest, adapted_forest);

        /* Free the former forest */
        t8_forest_unref(&previous_forest);

        /* Encode the data of this level and store it within a buffer */
        auto [encoded_entropy_codes, encoded_data] = adapt_data->EncodeLevelData(data_);
        buffered_entropy_codes_.push_back(std::move(encoded_entropy_codes));
        buffered_encoded_data_.push_back(std::move(encoded_data));

        /* Once the data is buffered, we can overwrite it with the adapted data */
        this->SwitchToExtractedData();

        /* Repartition the mesh */
        t8_forest_t partitioned_forest = RepartitionMesh(adapted_forest);

        /* Get the encoded the mesh adapatations */
        if (AreMeshRefinementBitsStored())
        {
            /* If the mesh refinement bits are stored */
            std::vector<uint8_t> encoded_mesh_data = mesh_encoder_->GetPartitionedEncodedLevelData(adapted_forest, partitioned_forest, this->GetMPIComm());
            buffered_encoded_mesh_.push_back(std::move(encoded_mesh_data));
        } else
        {
            /* To keep the encoded streams in sync, we emplace an empty vector even if no refinement bits are stored */
            buffered_encoded_mesh_.push_back(std::vector<uint8_t>());
        }

        /* Repartition the data */
        this->RepartitionData(adapted_forest, partitioned_forest);
        adapt_data->RepartitionData(adapted_forest, partitioned_forest);

        cmc_debug_msg("The mesh and the data has been re-partitioned.");

        /* Free the former forest and store the adapted/repartitioned mesh */
        t8_forest_unref(&adapted_forest);
        mesh_.SetMesh(partitioned_forest);

        /* Finalize the compression iteration */
        adapt_data->FinalizeExtractionIteration();
        if (AreMeshRefinementBitsStored()){mesh_encoder_->FinalizeCompressionIteration();}
        cmc_debug_msg("The coarsening iteration is finished.");
    }

    /* At last, we need to encode the root level data */
    auto [encoded_root_entropy_codes, encoded_root_data] = adapt_data->EncodeRootLevelData(data_);
    buffered_entropy_codes_.push_back(std::move(encoded_root_entropy_codes));
    buffered_encoded_data_.push_back(std::move(encoded_root_data));

    /* At last, we need to encode the root level of the mesh */
    std::vector<uint8_t> encoded_root_mesh = mesh_encoder_->EncodeRootLevelMesh(mesh_, attributes_.GetGlobalDomain());
    buffered_encoded_mesh_.push_back(std::move(encoded_root_mesh));

    /* Free the adapt data structure */
    this->adaptation_destructor_(adapt_data);
    cmc_debug_msg("Compression of variable ", this->name_, " is finished.");
}

template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::SetData(const std::vector<T>& initial_data)
{
    data_.reserve(initial_data.size());

    for (auto val_iter = initial_data.begin(); val_iter != initial_data.end(); ++val_iter)
    {
        data_.emplace_back(CompressionValue<T>(*val_iter));
    }
}

template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::SetData(const std::vector<SerializedCompressionValue<sizeof(T)>>& initial_data)
{
    data_ = initial_data;
}

template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::SetData(std::vector<SerializedCompressionValue<sizeof(T)>>&& initial_data)
{
    data_ = std::move(initial_data);
}

template <typename T>
VectorView<CompressionValue<T>>
AbstractEmbeddedByteCompressionVariable<T>::GetView(const int start_index, const int count) const
{
    cmc_assert(start_index >= 0 && count >= 0);
    cmc_assert(static_cast<size_t>(start_index + count) <= data_.size());

    return VectorView(&data_[start_index], count);
}


template <typename T>
VectorView<CompressionValue<T>>
AbstractEmbeddedByteCompressionVariable<T>::GetView(const int tree_id, const int lelement_index, const int count) const
{
    cmc_assert(tree_id >= 0 && tree_id < t8_forest_get_num_local_trees(mesh_.GetMesh()));

    /* Compute the start offset in the local contiguous array */
    const int start_index = t8_forest_get_tree_element_offset (mesh_.GetMesh(), tree_id) + lelement_index;

    cmc_assert(start_index >= 0 && count >= 0);
    cmc_assert(static_cast<size_t>(start_index + count) <= data_.size());

    return VectorView(&data_[start_index], count);

}

template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::StoreExtractedValues(const int tree_id, const int lelement_id, const int num_elements, const ExtractionData<T>& extracted_values)
{
    cmc_assert(tree_id >= 0 && tree_id < t8_forest_get_num_local_trees(mesh_.GetMesh()));
    cmc_assert(static_cast<size_t>(num_elements) == extracted_values.fine_values.size());

    /* Compute the start offset in the local contiguous array */
    const int elem_start_index = t8_forest_get_tree_element_offset (mesh_.GetMesh(), tree_id) + lelement_id;

    cmc_assert(elem_start_index >= 0 && num_elements >= 0);
    cmc_assert(static_cast<size_t>(elem_start_index + num_elements) <= data_.size());

    /* Store the adjusted fine values */
    std::copy_n(extracted_values.fine_values.begin(), num_elements, &data_[elem_start_index]);

    /* Store the extracted coarse value for the next coarsening iteration */
    data_new_.push_back(extracted_values.coarse_value);
}

template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::StoreUnchangedElement(const int tree_id, const int lelement_id, const UnchangedData<T>& extracted_value)
{
    cmc_assert(tree_id >= 0 && tree_id < t8_forest_get_num_local_trees(mesh_.GetMesh()));

    /* Get the element id in the contiguous array of all local elements */
    const int elem_id = t8_forest_get_tree_element_offset (mesh_.GetMesh(), tree_id) + lelement_id;
    
    cmc_assert(elem_id >= 0 && static_cast<size_t>(elem_id) < data_.size());

    /* Set the adjusted remaining element value */
    data_[elem_id] = extracted_value.fine_value;

    /* Set the new value for the next coarsening iteration */
    data_new_.push_back(extracted_value.coarse_value);
}


template <typename T>
inline t8_forest_t
AbstractEmbeddedByteCompressionVariable<T>::RepartitionMesh(t8_forest_t adapted_forest)
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
inline void
AbstractEmbeddedByteCompressionVariable<T>::RepartitionData(t8_forest_t adapted_forest, t8_forest_t partitioned_forest)
{
    /* Create an sc_array_t wrapper of the variable's data */
    sc_array_t* in_data = sc_array_new_data (static_cast<void*>(data_.data()), sizeof(CompressionValue<T>), data_.size());

    cmc_debug_msg("Number of local data elements before partitioning: ", data_.size());
    cmc_debug_msg("Number of local mesh elements before partitioning: ", t8_forest_get_local_num_elements(adapted_forest));
    cmc_debug_msg("Size of a single data element: ", in_data->elem_size);

    /* Allocate memory for the partitioned data */
    const t8_locidx_t new_num_elems = t8_forest_get_local_num_elements(partitioned_forest);
    data_new_ = std::vector<CompressionValue<T>>(new_num_elems);

    cmc_debug_msg("Number of local data elements after partitioning: ", data_new_.size());
    cmc_debug_msg("Number of local mesh elements after partitioning: ", new_num_elems);

    /* Create a wrapper for the freshly allocated partitioned data */
    sc_array_t* out_data = sc_array_new_data (static_cast<void*>(data_new_.data()), sizeof(CompressionValue<T>), data_new_.size());

    cmc_assert(static_cast<size_t>(t8_forest_get_local_num_elements(adapted_forest)) == data_.size());
    cmc_assert(static_cast<size_t>(t8_forest_get_local_num_elements(partitioned_forest)) == data_new_.size());

    /* Partition the variables data */
    t8_forest_partition_data(adapted_forest, partitioned_forest, in_data, out_data);

    /* Destroy the array wrappers */
    sc_array_destroy(in_data);
    sc_array_destroy(out_data);

    /* Set the variable's data to the newly partitioned data */
    SwitchToExtractedData();
    cmc_debug_msg("Partitioning of mesh and data elements has been finished.");
}

template <typename T>
inline bool
AbstractEmbeddedByteCompressionVariable<T>::IsValidForCompression() const 
{
    if (name_.empty())
    {
        cmc_err_msg("The variable needs a name in order to store the compressed output. Therefore, no compression can be applied.");
        return false;
    }
    if (data_.empty())
    {
        cmc_err_msg("There is no data attached to the variable. Therefore, no compression can be applied.");
        return false;
    }
    if (not mesh_.IsValid())
    {
        cmc_err_msg("The mesh is not valid. Therefore, no compression can be applied.");
        return false;
    }
    if (mesh_encoder_ == nullptr)
    {
        cmc_err_msg("The mesh encoder is not set. Therefore, no compression can be applied.");
        return false;
    }

#if CMC_ENABLE_MPI
    if (comm_ == MPI_COMM_NULL)
    {
        cmc_err_msg("The MPI-Communicator is not set. Therefore, no compression can be applied.");
        return false;
    }

    MPI_Comm mesh_comm = t8_forest_get_mpicomm(mesh_.GetMesh());
    int are_mpi_comms_equal{0};
    const int ret_val_comm_compare = MPI_Comm_compare(comm_, mesh_comm, &are_mpi_comms_equal);
    MPICheckError(ret_val_comm_compare);
    
    if (are_mpi_comms_equal != MPI_IDENT)
    {
        cmc_err_msg("The MPI-Communicator of the mesh and the variable differs. Therefore, no compression can be applied.");
        return false;
    }

#endif

    return true;
}












//////////////////////////////////
/////////////////////////////////
/////////////////////////////////









template <typename T>
AmrMesh
AbstractEmbeddedByteCompressionVariable<T>::BuildInitialMesh(const input::Variable<T>& input_variable)
{
    //cmc_assert(input_variable.IsValid());

    /* Since all varibales are defined on the same global domain, we are able to take the domain of the first variable */
    const GeoDomain& global_domain = input_variable.GetGlobalDomain();

    /* The mesh layout is equal to the layout of the domains (the exact ordering of the dimensions is allowed to vary) */
    const DataLayout initial_mesh_layout = input_variable.GetInitialDataLayout();

    /* Build the actual embedded mesh based on the given features */
    auto [initial_forest, initial_refinement_level, dimensionality] = BuildInitialEmbeddedMesh(global_domain, initial_mesh_layout, GetMPIComm());

    return AmrMesh(initial_forest, initial_refinement_level, dimensionality);
}

template <typename T>
std::vector<input::IndexReduction>
AbstractEmbeddedByteCompressionVariable<T>::UpdateLinearIndicesToTheInitialMesh()
{
    const t8_locidx_t num_local_elements = mesh_.GetNumberLocalElements();

    const int initial_refinement_level = mesh_.GetInitialRefinementLevel();

    const t8_eclass_t eclass = t8_forest_get_eclass(mesh_.GetMesh(), 0);

    /* Get the scheme of the forest's only tree */
    t8_eclass_scheme_c* ts =  t8_forest_get_eclass_scheme (mesh_.GetMesh(), eclass);

    /* Get the element scheme */
    t8_default_scheme_common_c* ts_c = static_cast<t8_default_scheme_common_c*>(ts);

    const t8_locidx_t first_ltree_id = 0;

    const int num_children = ts_c->t8_element_num_children(t8_forest_get_element_in_tree(mesh_.GetMesh(), first_ltree_id, 0));

    const MortonIndex linear_index_start_elem = GetMortonIndexOnLevel(t8_forest_get_element_in_tree(mesh_.GetMesh(), first_ltree_id, 0),
                                                   ts_c, t8_eclass_to_dimension[eclass], initial_refinement_level);
    
    /* Locally, this update function reduces the global_index to zero for the first local element */
    std::vector<input::IndexReduction> index_correction;
    index_correction.emplace_back(linear_index_start_elem, linear_index_start_elem);

    /* All elements that are not holding data need to be accumulated in order to be subtracted additionally for the local index correction */
    MortonIndex skipped_indices = linear_index_start_elem;

    bool coarse_element_streak = false;

    /* Iterate through all local elements and find how the global Morton indices need to be adjusted in order to comply the local data ordering */
    for (auto iter = 0; iter < num_local_elements; ++iter)
    {
        const t8_element_t* elem = t8_forest_get_element_in_tree(mesh_.GetMesh(), 0, iter);

        if (t8_element_level(ts, elem) != initial_refinement_level)
        {
            /* Get the number of uniform indices which were skipped by this element which lays outside of the domain */
            skipped_indices += std::pow(num_children, initial_refinement_level - t8_element_level(ts, elem)) - 1;
            coarse_element_streak = true;
        } else if (coarse_element_streak)
        {
            /* We accumulate the amount of skipped indices (with regard to the initial refinement level) and store the offset 
             * once we have reached again an element on the initial refinement level */
            const MortonIndex uniform_index_of_elem = GetMortonIndexOnLevel(elem, ts_c, t8_eclass_to_dimension[eclass], initial_refinement_level);
            index_correction.emplace_back(uniform_index_of_elem, skipped_indices);
            coarse_element_streak = false;
        }
    }

    /* Return the correction scheme for the global indices */
    return index_correction;
}

template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::SortInitialDataIntoVariables(input::Variable<T>& input_variable, const std::vector<VariableRecvMessage>& messages)
{
    cmc_debug_msg("\n\n\nIn Sort Initial Data into variables !!!!\n\n\n");
    /* Update the uniform morton indices to linear element IDs corresponding to the initial mesh */
    input::UpdateLinearIndices local_indices_update(UpdateLinearIndicesToTheInitialMesh());

    /* Setup the variable with the right amount of global elements */
    input_variable.SetUpFilledVariable(mesh_.GetNumberLocalElements(), input_variable.GetMissingValue());

    cmc_debug_msg("Size of messages: ", messages.size());
    /* Iterate over all messages and assign their data to the correct variables */
    for (auto msg_iter = messages.begin(); msg_iter != messages.end(); ++msg_iter)
    {
        cmc_debug_msg("input_variable.GetInternID(): ", input_variable.GetInternID(), ", msg_iter->GetVariableID(): ", msg_iter->GetVariableID());
        cmc_assert(input_variable.GetInternID() == msg_iter->GetVariableID());
        if (input_variable.GetInternID() == msg_iter->GetVariableID())
        {
            cmc_debug_msg("\n\n\nIs this happeing!!!!\n\n\n");
            /* Assign the data from the messag at the right position within the variable */
            input_variable.AssignDataAtLinearIndices(*msg_iter, local_indices_update);
        } else
        {
            cmc_warn_msg("A message has been received, but the variable ID does not correspond to any (local) variable.");
        }
    }
}

template <typename T>
[[nodiscard]]
std::pair<std::vector<VariableSendMessage>, std::vector<MPI_Request>>
AbstractEmbeddedByteCompressionVariable<T>::SendInitialData(input::Variable<T>& input_variable)
{
    cmc_debug_msg("Num coors in send data before transform: ", input_variable.GetNumberCoordinates());

    input_variable.TransformCoordinatesToMortonIndices();

    //cmc_assert(input_variable.GetActiveDataFormat());
    cmc_debug_msg("IN Send intiial data: active format: ", input_variable.GetActiveDataFormat());
    int comm_size{1};

    const int ret_val = MPI_Comm_size(GetMPIComm(), &comm_size);
    MPICheckError(ret_val);

    /* Gather the partitioning of the initial mesh */
    const DataOffsets offsets = GatherGlobalDataOffsets(mesh_, GetMPIComm());

    /* Inquire all data that has to be sent */
    std::vector<VariableSendMessage> send_messages;

    input::ReceiverMap<T> send_data = input_variable.GatherDataToBeDistributed(offsets);
    input::AppendSendData(send_messages, std::move(send_data));

    cmc_debug_msg("Messages to send: ", send_messages.size());

    /* A vector collecting all send requests */
    std::vector<MPI_Request> send_requests;
    send_requests.reserve(2 * send_messages.size());

    /* Send all messages */
    for (auto sm_iter = send_messages.begin(); sm_iter != send_messages.end(); ++sm_iter)
    {
        /* Send the message and receive the returned requests */
        auto [req_morton_ids, req_data] = sm_iter->Send(GetMPIComm());

        /* Store the requests */
        send_requests.push_back(std::move(req_morton_ids));
        send_requests.push_back(std::move(req_data));
    }

    return std::make_pair(std::move(send_messages), std::move(send_requests));
}

template <typename T>
std::vector<VariableRecvMessage>
AbstractEmbeddedByteCompressionVariable<T>::ReceiveInitialData(input::Variable<T>& input_variable)
{
    cmc_debug_msg("In receive initial data\n\n");
    /* Receive all messages */
    bool are_messages_incoming{true};

    std::vector<VariableRecvMessage> recv_messages;

    /* Message receiving loop */
    while(are_messages_incoming)
    {
        int message_flag{0};
        MPI_Status probe_status;

        /* Check if there are any messages incoming */
        int err = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, GetMPIComm(), &message_flag, &probe_status);
        MPICheckError(err);

        /* Get the message if there is one */
        if (message_flag != 0)
        {
            /* Get the necessary information for retrieving the message */
            const int source = probe_status.MPI_SOURCE;
            const int tag = probe_status.MPI_TAG;
            const int variable_id = GetIDFromTag(tag);
            const CmcType var_type = input_variable.GetType();
            const MPI_Datatype data_type = ConvertCmcTypeToMPIType(var_type);

            cmc_debug_msg("Message is received: ", "source: ", source, ", tag: ", tag, ", variable_id: ", variable_id);
            /* Get the number of elements from this message */
            int num_elems{0};
            if (IsMessageADataMessage(tag))
            {
                /* If the data contains actual data values */
                err = MPI_Get_count(&probe_status, data_type, &num_elems);
                MPICheckError(err);
            } else
            {
                /* If the data contains Morton indices */
                err = MPI_Get_count(&probe_status, MPI_MORTON_INDEX_T, &num_elems);
                MPICheckError(err);
            }

            /** We only do receive two messages from each process that sends data, since it is collected prior to the communication ,
              * Therefore, we need to check if we have already received elements from the rank this messag originated from. */
            std::vector<cmc::VariableRecvMessage>::iterator msg_iter = std::find_if(recv_messages.begin(), recv_messages.end(),
                                                                                    [&source, &variable_id](auto& msg){return (source == msg.GetSendingRank() && variable_id == msg.GetVariableID());});
            
            /* We receive (potentially) two messages from each process for each variable */

            /* Check if there has been a message found coming from this rank */
            if (msg_iter == recv_messages.end())
            {
                /* We have not yet received a message from this rank */
                /* Allocate a variable message */
                recv_messages.emplace_back(var_type, CreateVariableMessage(source, variable_id, var_type, num_elems));

                /* Assign the newly created message to the iterator */
                msg_iter = std::prev(recv_messages.end());
            }

            /* Check if it is a Morton indices message or the actual data message */
            if (IsMessageADataMessage(tag))
            {
                void* dataptr = msg_iter->GetInitialDataPtr();
                /* Receive the actual data from this process */
                err = MPI_Recv(dataptr, num_elems, data_type, source, tag, GetMPIComm(), MPI_STATUS_IGNORE);
                MPICheckError(err); 
            } else
            {
                void* miptr = msg_iter->GetInitialMortonIndicesPtr();
                /* Otherwise, we receive the sent Morton indices */
                err = MPI_Recv(miptr, num_elems, MPI_MORTON_INDEX_T, source, tag, GetMPIComm(), MPI_STATUS_IGNORE);
                MPICheckError(err);
            } 
        } else
        {
            /* If all messages have been processed, we break the receiving loop */
            are_messages_incoming = false;
        }
    }

    return recv_messages;
}

/* Sorting the initial data only locally */
template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::SortLocalDataOnInitialMesh(input::Variable<T>& input_variable)
{
    /* Update the uniform morton indices to linear element IDs corresponding to the initial mesh */
    input::UpdateLinearIndices local_indices_update(UpdateLinearIndicesToTheInitialMesh());

    input::Variable<T> sorted_input_variable = input::HollowCopy(input_variable);

    sorted_input_variable.SetUpFilledVariable(mesh_.GetNumberLocalElements(), input_variable.GetMissingValue());
    sorted_input_variable.AssignDataAtLinearIndices(input_variable, local_indices_update);

    /* Swap the sorted input variables with the initial input variables */
    std::swap(input_variable, sorted_input_variable);
}

template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::DistributeDataOnInitialMesh(input::Variable<T>& input_variable)
{
    #ifdef CMC_ENABLE_MPI
    cmc_debug_msg("The intial data will be distributed on the embedded mesh.");
    cmc_assert(mesh_.IsValid());

    /* Gather and send all the data that has to be communicated. 
     * The send messages has to be returned in order to no get deallocated before 
     * the actual sending has happend. Therefore, we keep them here 'alive' until 
     * all send_requests have been completed */
    auto [send_messages, send_requests] = SendInitialData(input_variable);

    /* Wait until all messages have been staged */
    int err = MPI_Barrier(comm_);
    MPICheckError(err);

    /* After all messages have been staged, we will receive them */
    const std::vector<VariableRecvMessage> received_messages = ReceiveInitialData(input_variable);

    /* Sort the data accordingly to the Morton indices */
    SortInitialDataIntoVariables(input_variable, received_messages);

    /* Wait until any send messages have been completed */
    err = MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUS_IGNORE);
    MPICheckError(err);

    cmc_debug_msg("The intial data has been distributed successfully on the embedded mesh.");

    #else
    /* Call a locally sorting function which setups the data compliant to the Morton order */
    cmc_assert(mesh.IsValid());
    SortLocalDataOnInitialMesh();
    #endif
}

}

#endif /* !LOSSLESS_CMC_EMBEDDED_BYTE_COMPRESSION_VARIABLE_HXX */
