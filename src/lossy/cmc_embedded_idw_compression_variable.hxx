#ifndef LOSSY_CMC_EMBEDDED_IDW_COMPRESSION_VARIABLE_HXX
#define LOSSY_CMC_EMBEDDED_IDW_COMPRESSION_VARIABLE_HXX

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
#include "utilities/cmc_iface_abstract_embedded_byte_compression_variable.hxx"
#include "utilities/cmc_compression_settings.hxx"

#include "mpi/cmc_mpi.hxx"
#include "t8code/cmc_t8_mpi.hxx"

#include <vector>
#include <functional>
#include <memory>
#include <algorithm>
#include <variant>
#include <type_traits>

namespace cmc::lossy::embedded::idw
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

/* Set the data type for inaccuracy tracking */
using err_type_t = double;

struct PointData
{
    PointData() = default;
    PointData(const float x_, const float y_, const float z_)
    : x{x_}, y{y_}, z{z_}{};

    float x{0.0f}, y{0.0f}, z{0.0f};
};

template<typename T>
struct ElementData
{
    ElementData() = default;

    T value;
    PointData coordinates;
};

inline
double
Dist(const PointData& coords1, const PointData& coords2)
{
    return std::sqrt((coords1.x - coords2.x) * (coords1.x - coords2.x) + (coords1.y - coords2.y) * (coords1.y - coords2.y) + (coords1.z - coords2.z) * (coords1.z - coords2.z));
}

inline
double
Dist(const PointData& coords1, const std::vector<double>& coords2)
{
    return std::sqrt((static_cast<double>(coords1.x) - coords2[0]) * (static_cast<double>(coords1.x) - coords2[0]) + (static_cast<double>(coords1.y) - coords2[1]) * (static_cast<double>(coords1.y) - coords2[1]) + (static_cast<double>(coords1.z) - coords2[2]) * (static_cast<double>(coords1.z) - coords2[2]));
}

/**
 * @brief A struct holding the data for an extraction process.
 * 
 * @tparam T The origianl data type of the underlying data (e.g. float)
 */
template<typename T>
struct ExtractionData
{
    ExtractionData(const ElementData<T>& element_data)
    : extracted_values(element_data) {};
    ExtractionData(ElementData<T>&& element_data)
    : extracted_values(std::move(element_data)) {};

    ElementData<T> extracted_values;
};

template<typename T>
struct ResidualData
{
    ResidualData() = default;
    ResidualData(const std::vector<CompressionValue<T>>& residuals)
    : values(residuals) {};
    ResidualData(std::vector<CompressionValue<T>>&& residuals)
    : values(std::move(residuals)) {};

    std::vector<CompressionValue<T>> values;
};

/**
 * @brief A strcut holding the data for an adapation process that leaves the element unchanged.
 * 
 * @tparam T The origianl data type of the underlying data (e.g. float)
 */
template<typename T>
struct UnchangedData
{
    UnchangedData(const ElementData<T>& element_data)
    : extracted_values(element_data) {};
    UnchangedData(ElementData<T>&& element_data)
    : extracted_values(std::move(element_data)) {};

    ElementData<T> extracted_values;
};


/* Forward declarations */
template <typename T>
class AbstractEmbeddedByteCompressionVariable;
template <typename T>
class IEmbeddedCompressionAdaptData;

template<typename T>
using AdaptCreator = std::function<IEmbeddedCompressionAdaptData<T>*(AbstractEmbeddedByteCompressionVariable<T>*, const CompressionSettings&)>;

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
class AbstractEmbeddedByteCompressionVariable : public IEmbeddedByteCompressionVariable<T>
{
public:
    void Compress() override;

    const std::string& GetName() const override {return name_;};

    size_t Size() const override {return data_.size();};

    const AmrMesh& GetAmrMesh() const override {return mesh_;};
    const std::vector<CompressionValue<T>>& GetData() const {return data_;};

    virtual ~AbstractEmbeddedByteCompressionVariable(){};

    void MoveEncodedEntropyCodesInto(std::vector<std::vector<uint8_t>>& vec_to_hold_encoded_levelwise_entropy_codes) override
    {
        vec_to_hold_encoded_levelwise_entropy_codes = std::move(buffered_entropy_codes_);
        buffered_entropy_codes_ = std::vector<std::vector<uint8_t>>();
    };

    void MoveEncodedDataInto(std::vector<std::vector<uint8_t>>& vec_to_hold_encoded_levelwise_data) override
    {
        vec_to_hold_encoded_levelwise_data = std::move(buffered_encoded_data_);
        buffered_encoded_data_ = std::vector<std::vector<uint8_t>>();
    };

    void MoveEncodedMeshInto(std::vector<std::vector<uint8_t>>& vec_to_hold_encoded_levelwise_mesh) override
    {
        vec_to_hold_encoded_levelwise_mesh = std::move(buffered_encoded_mesh_);
        buffered_encoded_mesh_ = std::vector<std::vector<uint8_t>>();
    };

    const std::vector<std::vector<uint8_t>>& GetEncodedEntropyCodes() const override {return buffered_entropy_codes_;}
    const std::vector<std::vector<uint8_t>>& GetEncodedData() const override {return buffered_encoded_data_;};
    const std::vector<std::vector<uint8_t>>& GetEncodedMesh() const override {return buffered_encoded_mesh_;};

    bool AreMeshRefinementBitsStored() const override {return store_refinement_indication_bits_;};
    MPI_Comm GetMPIComm() const override {return comm_;};
    const VariableAttributes<T>& GetVariableAttributes() const override {return attributes_;};

    virtual CompressionSchema GetCompressionSchema() const = 0;

    int GetInitialMaximumRefinementLevel() const {return max_initial_refinement_level_;}

    friend IEmbeddedCompressionAdaptData<T>;
protected:
    AbstractEmbeddedByteCompressionVariable() = delete;
    AbstractEmbeddedByteCompressionVariable(const CompressionSettings& settings, input::Var& input_variable)
    : settings_(settings) {
        this->SetupInputVariableForCompression(input_variable);
    };
    AbstractEmbeddedByteCompressionVariable(CompressionSettings&& settings, input::Var& input_variable)
    : settings_(std::move(settings)) {
        this->SetupInputVariableForCompression(input_variable);
    };

    virtual void PreCompressionProcessing([[maybe_unused]] std::vector<CompressionValue<T>>& initial_data){}

    void SetName(const std::string& name) {name_ = name;};
    void SetAmrMesh(const AmrMesh& mesh) {mesh_ = mesh;};
    void SetAmrMesh(AmrMesh&& mesh) {mesh_ = std::move(mesh);};
    void SetMPIComm(const MPI_Comm comm) {comm_ = comm;};
    void SetAttributes(const cmc::VariableAttributes<T>& attributes) {attributes_ = attributes;}
    void SetAttributes(cmc::VariableAttributes<T>&& attributes) {attributes_ = std::move(attributes);}
    void IndicateWhetherMeshRefinementBitsWillBeStored(const bool store_indication_bits) {store_refinement_indication_bits_ = true;}
    
    AdaptCreator<T> adaptation_creator_; //!< A function pointer which is used to create the wished adaptation structure
    AdaptDestructor<T> adaptation_destructor_; //!< A function pointer which is used to destruct the adaptation structure
    
    std::unique_ptr<mesh_compression::IEmbeddedMeshEncoder> mesh_encoder_{nullptr};
private:
    VectorView<T> GetViewOnInitialData(const int start_index, const int count) const;
    VectorView<ElementData<T>> GetView(const int tree_id, const int lelement_index, const int count) const;
    T GetInitialDataValueAtIndex(const int local_index) const;
    const ElementData<T>& GetDataValueAtIndex(const int tree_id, const int lelement_index) const;

    VectorView<ElementData<T>> GetViewOnAdaptedData(const int start_index, const int count) const;
    VectorView<ElementData<T>> GetViewOnAdaptedData(t8_forest_t adapted_forest, const int tree_id, const int lelement_index, const int count) const;
    const ElementData<T>& GetAdaptedDataValueAtIndex(const int local_index) const;

    VectorView<err_type_t> GetViewOnRemainingPermittedErrors(const int tree_id, const int lelement_index, const int count) const;
    err_type_t GetAdaptedRemainingPermittedError(t8_forest_t adapted_forest, const int tree_id, const int lelement_index) const;
    err_type_t GetRemainingPermittedError(const int tree_id, const int lelement_index) const;
    err_type_t GetNewRemainingPermittedError(t8_forest_t adapted_forest, const int tree_id, const int lelement_index) const;
    void StoreAdaptedRemainingPermittedError(const err_type_t remaining_error);
    void UpdateNewRemainingPermittedError(t8_forest_t adapted_forest, const int tree_id, const int lelement_index, const err_type_t new_remaining_error);
    void ClearAdaptedRemainingPermittedErrors();
    void ReserveAdaptedRemainingPermittedError(const int num_elements_coarsened_forest);
    void SwitchToAdaptedRemainingErrors();

    void StoreExtractedValues(const int tree_id, const int lelement_id, const int num_elements, const ExtractionData<T>& extracted_values);
    void StoreUnchangedElement(const int tree_id, const int lelement_id, const UnchangedData<T>& extracted_value);
    void StoreResiduals(const t8_forest_t forest, const int tree_id, const int lelement_id, const int num_elements, const ResidualData<T>& residuals);

    void IndicateElementStaysUnchanged() {if (store_refinement_indication_bits_) {mesh_encoder_->IndicateElementStaysUnchanged();}};
    void IndicateCoarsening() {if (store_refinement_indication_bits_) {mesh_encoder_->IndicateCoarsening();}};

    bool IsValidForCompression() const;
    void AllocateExtractionIteration() {data_new_.reserve(mesh_.GetNumberLocalElements() / (2 << mesh_.GetDimensionality()) + 8);}
    void SwitchToExtractedData() {data_.swap(data_new_); data_new_.clear();};
    IEmbeddedCompressionAdaptData<T>* CreateAdaptData() {return adaptation_creator_(this, settings_);};
    t8_forest_t RepartitionMesh(t8_forest_t adapted_forest);
    void RepartitionData(t8_forest_t adapted_forest, t8_forest_t partitioned_forest);
    void RepartitionRemainingPermittedErrors(t8_forest_t adapted_forest, t8_forest_t partitioned_forest);

    void SetupInputVariableForCompression(input::Var& input_variable);
    void DistributeDataOnInitialMesh(input::Variable<T>& input_variable);
    void SortLocalDataOnInitialMesh(input::Variable<T>& input_variable);
    std::vector<VariableRecvMessage> ReceiveInitialData(input::Variable<T>& input_variable);
    std::pair<std::vector<VariableSendMessage>, std::vector<MPI_Request>> SendInitialData(input::Variable<T>& input_variable);
    void SortInitialDataIntoVariables(input::Variable<T>& input_variable, const std::vector<VariableRecvMessage>& messages);
    std::vector<input::IndexReduction> UpdateLinearIndicesToTheInitialMesh();
    AmrMesh BuildInitialMesh(const input::Variable<T>& input_variable);
    void DetermineInitialMaximumRefinementLevelAndAbsoluteErrorBounds();

    std::string name_; //!< The name of the variable
    AmrMesh mesh_; //!< The mesh on which the variable is defined
    std::vector<T> initial_data_; //!< The current data of the variable 

    std::vector<ElementData<T>> data_;
    std::vector<ElementData<T>> data_new_; //!< A helper variable for the adaptation

    std::vector<CompressionValue<T>> residuals;

    std::vector<err_type_t> remaining_permitted_error_; //!< Tracking the left-over permitted inaccuracy per element 
    std::vector<err_type_t> remaining_permitted_error_new_; //!< A helper variable for the adaptation 

    std::vector<std::vector<uint8_t>> buffered_entropy_codes_; //!< Level-wise storage of the entropy codes, e.g. LZC or prefix lengths
    std::vector<std::vector<uint8_t>> buffered_encoded_data_; //!< Level-wise storage for the encoded data
    std::vector<std::vector<uint8_t>> buffered_encoded_mesh_; //!< Level-wwise storage for the encoded mesh

    cmc::VariableAttributes<T> attributes_;
    bool store_refinement_indication_bits_{true};
    
    const CompressionSettings settings_;
    MPI_Comm comm_{MPI_COMM_NULL}; //!< The MPI communicator to use
    int max_initial_refinement_level_{-1};

    //std::vector<CompressionValue<T>> test_data_;
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
    IEmbeddedCompressionAdaptData(AbstractEmbeddedByteCompressionVariable<T>* variable, const CompressionSettings& settings)
    : base_variable_{variable}, settings_{settings} {};

    bool IsCompressionProgressing() const;

    virtual void InitializeExtractionIteration() = 0;
    void FinalizeExtractionIteration();
    void CompleteExtractionIteration(const t8_forest_t previous_forest, const t8_forest_t adapted_forest);
    virtual void RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest) = 0;

    int ExtractValue(t8_forest_t forest, t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                     const t8_scheme_c * ts, const int num_elements, t8_element_t * elements[]);
    int LeaveElementUnchanged(t8_forest_t forest, t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
        const t8_scheme_c * ts, const int num_elements, t8_element_t * elements[]);

    virtual std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeLevelData(const std::vector<CompressionValue<T>>& level_values) const = 0;
    virtual std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const = 0;
    
    virtual void TestQuantizationEntropyCoding() = 0;
    
    MPI_Comm GetMPIComm() const {return base_variable_->GetMPIComm();};

    virtual ~IEmbeddedCompressionAdaptData(){};

    bool IsValidForCompression() const;
    const VariableAttributes<T>& GetVariableAttributes() const {cmc_assert(base_variable_ != nullptr); return base_variable_->GetVariableAttributes();}
    int GetInitialMaximumRefinementLevel() {return base_variable_->GetInitialMaximumRefinementLevel();}
    int GetCurrentCompressionStep() const {return compression_step_;}
    
    const AmrMesh& GetAmrMesh() const {return base_variable_->GetAmrMesh();};

    VectorView<T> GetViewOnInitialData(const int start_index, const int count) const {return base_variable_->GetViewOnInitialData(start_index, count);}
    VectorView<ElementData<T>> GetView(const int tree_id, const int lelement_index, const int count) const {return base_variable_->GetView(tree_id, lelement_index, count);}
    T GetInitialDataValueAtIndex(const int local_index) const {return base_variable_->GetInitialDataValueAtIndex(local_index);}
    
    VectorView<ElementData<T>> GetViewOnAdaptedData(const int start_index, const int count) const {return base_variable_->GetViewOnAdaptedData(start_index, count);}
    VectorView<ElementData<T>> GetViewOnAdaptedData(t8_forest_t adapted_forest, const int tree_id, const int lelement_index, const int count) const {return base_variable_->GetView(adapted_forest, tree_id, lelement_index, count);}
    const ElementData<T>& GetAdaptedDataValueAtIndex(const int local_index) const {return base_variable_->GetAdaptedDataValueAtIndex(local_index);}
    
    VectorView<err_type_t> GetViewOnRemainingPermittedErrors(const int tree_id, const int lelement_index, const int count) const {return base_variable_->GetViewOnRemainingPermittedErrors(tree_id, lelement_index, count);}
    err_type_t GetAdaptedRemainingPermittedError(t8_forest_t adapted_forest, const int tree_id, const int lelement_index) const {return base_variable_->GetAdaptedRemainingPermittedError(adapted_forest, tree_id, lelement_index);}
    err_type_t GetRemainingPermittedError(const int tree_id, const int lelement_index) const {return base_variable_->GetRemainingPermittedError(tree_id, lelement_index);}

    std::vector<err_type_t> GatherCoarseFaceRemainingPermittedErrors(t8_forest_t forest, const int tree_id, const t8_eclass_t tree_class, const t8_scheme_c *scheme, const t8_element_t* elem) const;

    void StoreAdaptedRemainingPermittedError(const err_type_t remaining_error) {return base_variable_->StoreAdaptedRemainingPermittedError(remaining_error);}
    void UpdateNewRemainingPermittedError(t8_forest_t adapted_forest, const int tree_id, const int lelement_index, const err_type_t new_remaining_error) {return base_variable_->UpdateNewRemainingPermittedError(adapted_forest, tree_id, lelement_index, new_remaining_error);}
    
    void StoreResiduals(t8_forest_t previous_forest, const int tree_id, const int lelement_id, const int num_elements, const ResidualData<T>& residuals) {return base_variable_->StoreResiduals(previous_forest, tree_id, lelement_id, num_elements, residuals);}
    
    virtual std::pair<ResidualData<T>, err_type_t> ComputeResiduals(t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
        const t8_eclass_t tree_class, const t8_scheme_c *scheme, const int refine,
        const int num_outgoing, const t8_locidx_t first_outgoing, const int num_incoming,
        const t8_locidx_t first_incoming) = 0;

    protected:
    virtual ExtractionData<T> PerformExtraction(t8_forest_t forest, t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                const t8_scheme_c * ts, const int num_elements, t8_element_t * elements[]) = 0;
    virtual UnchangedData<T> ElementStaysUnchanged(t8_forest_t forest, t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
        const t8_scheme_c * ts, const int num_elements, t8_element_t * elements[]) = 0;
    virtual void FinalizingExtractionIteration() = 0;
    
    std::unique_ptr<entropy_coding::IByteCompressionEntropyCoder> entropy_coder_{nullptr}; //!< The entropy coder to use in order to encode information
private:

    AbstractEmbeddedByteCompressionVariable<T>* const base_variable_{nullptr};
    const CompressionSettings& settings_;
    int compression_step_{0};
};

template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::SetupInputVariableForCompression(input::Var& input_var)
{
    cmc_debug_msg("The InputVariable will be set up for compression.");

    /* get the MPI communicator from the variable */
    comm_ = input_var.GetMPIComm();

    //const std::vector<T>& values = std::get<input::Variable<T>>(input_var.GetInternalVariant()).GetDataForReading();
    //FILE* file_out = fopen("direct_output_data_bin_reader.cmc", "wb");
    //fwrite(values.data(), sizeof(T), values.size(), file_out);
    //fclose(file_out);
    //cmc_err_msg("End here");

    /* Potentially, apply scaling values and offsets if defined (potentially, the datatype may be changed by this call) */
    input_var.ApplyScalingAndOffset();
    
    if (not std::holds_alternative<input::Variable<T>>(input_var.GetInternalVariant()))
    {
        cmc_err_msg("The data type of the input variable differs from the embedded byte comrpession variable.");
    }

    /* Get the actual variable from the input var wrapper */
    input::Variable<T> input_variable = std::get<input::Variable<T>>(input_var.GetInternalVariant());

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
    cmc_debug_msg("Vor vtk write: num local elems: ", t8_forest_get_local_num_leaf_elements(forest));
    const int vtk_err = t8_forest_vtk_write_file(forest, "example_t2m_input", 0, 1, 0, 0, 0, 1, vtk_data);
    
    if (vtk_err == 0)
        cmc_err_msg("An error occrued during the creation of the t8code-forest vtk file.");
    
    delete[] vtk_data;

    //FILE* file_out = fopen("initial_input_prefix_amr_data.cmc", "wb");
    //fwrite(input_variable.GetDataForReading().data(), sizeof(T), input_variable.GetDataForReading().size(), file_out);
    //fclose(file_out);

    #endif


    /* Get the data from the variable and store it as CompressionValues */
    initial_data_ = input_variable.GetDataForReading();
    
    cmc_debug_msg("The setup of the InputVariable for compression has been successfull.");
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
IEmbeddedCompressionAdaptData<T>::ExtractValue(t8_forest_t forest, t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                               const t8_scheme_c * ts, const int num_elements, t8_element_t * elements[])
{
    cmc_assert(which_tree >= 0 && lelement_id >= 0);
    cmc_assert(num_elements > 1);

    /* We indicate to that a coarsening is applied to the family of elements */
    base_variable_->IndicateCoarsening();

    /* Extract the coarse values, and potentially alter the remaining fine values */
    const ExtractionData<T> extracted_values = PerformExtraction(forest, which_tree, tree_class, lelement_id, ts, num_elements, elements);

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
IEmbeddedCompressionAdaptData<T>::LeaveElementUnchanged(t8_forest_t forest, t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
    const t8_scheme_c * ts, const int num_elements, t8_element_t * elements[])
{
    cmc_assert(which_tree >= 0 && lelement_id >= 0);

    /* We indicate to that the element stays unchanged */
    base_variable_->IndicateElementStaysUnchanged();

    /* Leave the element unchanged */
    const UnchangedData<T> extracted_values = this->ElementStaysUnchanged(forest, which_tree, tree_class, lelement_id, ts, num_elements, elements);

    /* Store the unchanged data */
    base_variable_->StoreUnchangedElement(which_tree, lelement_id, extracted_values);

    return cmc::t8::kLeaveElementUnchanged;
}

template <typename T>
void
IEmbeddedCompressionAdaptData<T>::FinalizeExtractionIteration()
{
    ++compression_step_;
    this->FinalizingExtractionIteration();
}

constexpr int kIndicationOfCoarseningDuringAdaptation = -1;

template <typename T>
void
ResidualComputation(t8_forest_t previous_forest, t8_forest_t coarsened_forest, t8_locidx_t tree_id,
                    const t8_eclass_t tree_class, const t8_scheme_c *scheme, const int refine,
                    const int num_outgoing, const t8_locidx_t first_outgoing, const int num_incoming,
                    const t8_locidx_t first_incoming)
{
    IEmbeddedCompressionAdaptData<T>* adapt_data = static_cast<IEmbeddedCompressionAdaptData<T>*>(t8_forest_get_user_data(coarsened_forest));
    cmc_assert(adapt_data != nullptr);

    /* Compute the residuals and the maximum introduced absolute error */
    auto [residuals, new_remaining_permitted_error] = adapt_data->ComputeResiduals(previous_forest, coarsened_forest, tree_id, tree_class, scheme, refine, num_outgoing,
                                                             first_outgoing, num_incoming, first_incoming);

    /* Update the values on the finer mesh to the just calculated residuals */
    adapt_data->StoreResiduals(previous_forest, tree_id, first_outgoing, num_outgoing, residuals);

    /* Store the currenty maximum introduced error; Update the remaining permitted error */
    adapt_data->StoreAdaptedRemainingPermittedError(new_remaining_permitted_error);
}


template <typename T>
std::vector<err_type_t>
IEmbeddedCompressionAdaptData<T>::GatherCoarseFaceRemainingPermittedErrors(t8_forest_t forest, const int tree_id, const t8_eclass_t tree_class, const t8_scheme_c *scheme, const t8_element_t* elem) const 
{

    /* Get the number of faces for this element */
    const int num_faces = scheme->element_get_num_faces(tree_class, elem);

    std::vector<err_type_t> remaining_errors(num_faces, std::numeric_limits<err_type_t>::max());

    for (int face_idx = 0; face_idx < num_faces; ++face_idx)
    {
        t8_element_t** neighbor_leaves;
        int* dual_faces;
        int num_neighbors{0};
        t8_locidx_t* neighbor_element_indices;
        t8_eclass_t neighbor_tree_class;

        /* Gather the face neighbor via this face */
        t8_forest_leaf_face_neighbors (forest, tree_id, elem, &neighbor_leaves, face_idx, &dual_faces, &num_neighbors,
                                       &neighbor_element_indices, &neighbor_tree_class, 1);

        /* Check if there is a neighboring element at the face */
        if (num_neighbors > 0)
        {
            /* Since the forest is balanced, there should be exactly one face neighbor element */
            cmc_assert(num_neighbors == 1);
            cmc_assert(neighbor_element_indices[0] < t8_forest_get_local_num_leaf_elements(forest));

            /* Get the remaining permitted error of this face */
            const err_type_t coarse_face_remaining_permitted_error = this->GetAdaptedRemainingPermittedError(forest, tree_id, neighbor_element_indices[0]);

            /* Store the gathered remaining permitted error */
            remaining_errors[face_idx] = coarse_face_remaining_permitted_error;

            /* Deallocate the memory for the face neighbor construction */
            scheme->element_destroy (neighbor_tree_class, num_neighbors, neighbor_leaves);
            T8_FREE (neighbor_leaves);
            T8_FREE (neighbor_element_indices);
            T8_FREE (dual_faces);
        }
    
    }
    /* In case there is no face neighbor at the requested face */
    return remaining_errors;
}

template <typename T>
void
ExchangeNewRemainingPermittedErrors(t8_forest_t previous_forest, t8_forest_t coarsened_forest, t8_locidx_t tree_id,
    const t8_eclass_t tree_class, const t8_scheme_c *scheme, const int refine,
    const int num_outgoing, const t8_locidx_t first_outgoing, const int num_incoming,
    const t8_locidx_t first_incoming)
{
    if (refine != cmc::lossy::embedded::idw::kIndicationOfCoarseningDuringAdaptation)
    {
        /* In case the elements does not have been refined, there is nothing to be done,
         * since the the maximum remaining permitted error for this elemeent has already been set by the residual computation */
        return;
    }

    /* In case a family of elements have been coarsened, we need to gather the minimum remaining permitted error from the faceneighbors and update it */
    IEmbeddedCompressionAdaptData<T>* adapt_data = static_cast<IEmbeddedCompressionAdaptData<T>*>(t8_forest_get_user_data(coarsened_forest));
    cmc_assert(adapt_data != nullptr);

    /* Get the coarse element that has been introduced in the forest */
    const t8_element_t* elem = t8_forest_get_leaf_element_in_tree(coarsened_forest, tree_id, first_incoming);

    /* Gather all remaining permitted errors from the face neighbors */
    std::vector<err_type_t> remaining_errors = adapt_data->GatherCoarseFaceRemainingPermittedErrors(coarsened_forest, tree_id, tree_class, scheme, elem);

    /* We need to check whether the permitted error from this coarse element is smnaller than those from the face neighbors as well */
    remaining_errors.push_back(adapt_data->GetAdaptedRemainingPermittedError(coarsened_forest, tree_id, first_incoming));

    /* Find the minimum remaining permitted error */
    auto min_iter = std::min_element(remaining_errors.begin(), remaining_errors.end());
    cmc_assert(min_iter != remaining_errors.end());

    const err_type_t min_remaining_permitted_error = *min_iter;

    /* Update the remaining permitted error for this element */
    adapt_data->UpdateNewRemainingPermittedError(coarsened_forest, tree_id, first_incoming, min_remaining_permitted_error);
}

template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::StoreResiduals(const t8_forest_t previous_forest, const int tree_id, const int lelement_id, const int num_elements, const ResidualData<T>& residuals)
{
    #if 0
    cmc_assert(tree_id >= 0 && tree_id < t8_forest_get_num_local_trees(previous_forest));

    cmc_assert(static_cast<size_t>(num_elements) == residuals.values.size());

    /* Compute the start offset in the local contiguous array */
    const int elem_start_index = t8_forest_get_tree_element_offset (previous_forest, tree_id) + lelement_id;

    cmc_assert(elem_start_index >= 0 && num_elements >= 0);
    cmc_assert(static_cast<size_t>(elem_start_index + num_elements) <= data_.size());

    /* Store the adjusted fine values */
    std::copy_n(residuals.values.begin(), num_elements, &residuals[elem_start_index]);
    #endif
    //std::copy_n(residuals.values.begin(), num_elements, std::back_inserter(test_data_));
    //test_data_.pop_back();
}


template <typename T>
void
IEmbeddedCompressionAdaptData<T>::CompleteExtractionIteration(const t8_forest_t previous_forest, const t8_forest_t adapted_forest)
{
    cmc_assert(previous_forest != nullptr && adapted_forest != nullptr);
    cmc_debug_msg("The residual computation starts...");

    /* Set the user data for the forests */
    t8_forest_set_user_data(previous_forest, static_cast<void*>(this));
    t8_forest_set_user_data(adapted_forest, static_cast<void*>(this));

    /* ALlocate memory for the new remaining permitted errors */
    this->base_variable_->ReserveAdaptedRemainingPermittedError(t8_forest_get_local_num_leaf_elements(adapted_forest));

    /* Perform the computation of the reisduals if a coarsening has been taken place */
    t8_forest_iterate_replace (adapted_forest, previous_forest, ResidualComputation<T>);

    this->TestQuantizationEntropyCoding();

    /* At this point, we do not need the remaining_permitted_error_ from the base variable anymore. Therefore, we will already fill it with the updated data */
    this->base_variable_->SwitchToAdaptedRemainingErrors();

    /* We need to perform an additional iterate replace in order to update the new remaining permitted errors 
     * via the face neighbors */
    t8_forest_iterate_replace (adapted_forest, previous_forest, ExchangeNewRemainingPermittedErrors<T>);

    /* After the remaining permitted errors have been updated via the face neighbors, we do not need the remaining_permitted_errors_new_ anymore. They will be filled again by the next adaptation step */
    this->base_variable_->ClearAdaptedRemainingPermittedErrors();

    cmc_debug_msg("The residual computation has been finished and the remaining permitted errors per element have been updated accordingly.");
}

inline
bool CheckIfElementIsReadyForCoarsening(const int initial_max_ref_level, const int elem_level, const int compression_step)
{
    if (elem_level == initial_max_ref_level - compression_step)
    {
        return true;
    } else
    {
        return false;
    }
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
LossyByteCompression ([[maybe_unused]] t8_forest_t forest,
                         t8_forest_t forest_from,
                         t8_locidx_t which_tree,
                         const t8_eclass_t tree_class,
                         t8_locidx_t lelement_id,
                         const t8_scheme_c * ts,
                         const int is_family,
                         const int num_elements,
                         t8_element_t * elements[])
{
    /* Retrieve the adapt_data */
    IEmbeddedCompressionAdaptData<T>* adapt_data = static_cast<IEmbeddedCompressionAdaptData<T>*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    /* Check if a family is supplied to the adaptation function */
    if (is_family == 0 || not CheckIfElementIsReadyForCoarsening(adapt_data->GetInitialMaximumRefinementLevel(), ts->element_get_level(tree_class, elements[0]), adapt_data->GetCurrentCompressionStep()))
    {
        /* If there is no family, the element stays unchanged */
        const int ret_val = adapt_data->LeaveElementUnchanged(forest_from, which_tree, tree_class, lelement_id,
            ts, num_elements, elements);
        return ret_val;
    } else
    {
        /* Extract a value of the family and coarsen it */
        const int ret_val = adapt_data->ExtractValue(forest_from, which_tree, tree_class, lelement_id,
                                                     ts, num_elements, elements);
        return ret_val;
    }
}

template <typename T>
inline void
AbstractEmbeddedByteCompressionVariable<T>::Compress()
{
    /* Potentially, create a pre-compression processing step */
    //this->PreCompressionProcessing(data_);

    //cmc_assert(this->IsValidForCompression());
    cmc_debug_msg("Lossy compression of variable ", this->name_, " starts...");

    /* Determine the maximum present refinement level and setup the initial remaining permitted errors */
    this->DetermineInitialMaximumRefinementLevelAndAbsoluteErrorBounds();
    cmc_assert(this->GetInitialMaximumRefinementLevel() > 0);
    cmc_assert(remaining_permitted_error_.size() == mesh_.GetNumberLocalElements());

    /* We create the adapt data based on the compression settings, the forest and the variables to consider during the adaptation/coarsening */
    IEmbeddedCompressionAdaptData<T>* adapt_data = this->CreateAdaptData();

    //cmc_assert(adapt_data->IsValidForCompression());

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

        //test_data_.clear();
        //test_data_.reserve(t8_forest_get_local_num_leaf_elements(previous_forest));

        /* Perform a coarsening iteration */
        t8_forest_t adapted_forest = t8_forest_new_adapt(previous_forest, LossyByteCompression<T>, 0, 1, static_cast<void*>(adapt_data)); //Adapt the forest accordingly and create a ghost layer
        cmc_debug_msg("The mesh adaptation step is finished; resulting in ", t8_forest_get_global_num_leaf_elements(adapted_forest), " global elements.");

        /* Complete the interpolation step by storing the newly computed adapted data alongside its deviations */
        adapt_data->CompleteExtractionIteration(previous_forest, adapted_forest);

        /* Free the former forest */
        t8_forest_unref(&previous_forest);
        
        cmc_err_msg("Currently stop here for testing prediction ");
        /* Encode the data of this level and store it within a buffer */
        //auto [encoded_entropy_codes, encoded_data] = adapt_data->EncodeLevelData(test_data_);
        #if 0
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
        #endif
    }

    #if 0
    /* At last, we need to encode the root level data */
    auto [encoded_root_entropy_codes, encoded_root_data] = adapt_data->EncodeRootLevelData(data_);
    buffered_entropy_codes_.push_back(std::move(encoded_root_entropy_codes));
    buffered_encoded_data_.push_back(std::move(encoded_root_data));

    /* At last, we need to encode the root level of the mesh */
    std::vector<uint8_t> encoded_root_mesh = mesh_encoder_->EncodeRootLevelMesh(mesh_, attributes_.GetGlobalDomain());
    buffered_encoded_mesh_.push_back(std::move(encoded_root_mesh));
    #endif 

    /* Free the adapt data structure */
    this->adaptation_destructor_(adapt_data);
    cmc_debug_msg("Compression of variable ", this->name_, " is finished.");
}


template <typename T>
VectorView<T>
AbstractEmbeddedByteCompressionVariable<T>::GetViewOnInitialData(const int start_index, const int count) const
{
    cmc_assert(start_index >= 0 && count >= 0);
    cmc_assert(static_cast<size_t>(start_index + count) <= initial_data_.size());

    return VectorView(&initial_data_[start_index], count);
}


template <typename T>
VectorView<ElementData<T>>
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
T
AbstractEmbeddedByteCompressionVariable<T>::GetInitialDataValueAtIndex(const int local_index) const
{
    cmc_assert(static_cast<size_t>(local_index) < initial_data_.size());
    return initial_data_[local_index];
}

template <typename T>
const ElementData<T>&
AbstractEmbeddedByteCompressionVariable<T>::GetDataValueAtIndex(const int tree_id, const int lelement_index) const
{
    cmc_assert(tree_id >= 0 && tree_id < t8_forest_get_num_local_trees(mesh_.GetMesh()));

    /* Compute the start offset in the local contiguous array */
    const int start_index = t8_forest_get_tree_element_offset (mesh_.GetMesh(), tree_id) + lelement_index;

    cmc_assert(start_index >= 0 && static_cast<size_t>(start_index) < data_.size());

    return data_[start_index];
}

template <typename T>
VectorView<ElementData<T>>
AbstractEmbeddedByteCompressionVariable<T>::GetViewOnAdaptedData(const int start_index, const int count) const
{
    cmc_assert(not data_new_.empty());
    cmc_assert(start_index >= 0 && count >= 0);
    cmc_assert(static_cast<size_t>(start_index + count) <= data_new_.size());

    return VectorView(&data_new_[start_index], count);
}


template <typename T>
VectorView<ElementData<T>>
AbstractEmbeddedByteCompressionVariable<T>::GetViewOnAdaptedData(t8_forest_t adapted_forest, const int tree_id, const int lelement_index, const int count) const
{
    cmc_assert(not data_new_.empty());
    cmc_assert(tree_id >= 0 && tree_id < t8_forest_get_num_local_trees(adapted_forest));

    /* Compute the start offset in the local contiguous array */
    const int start_index = t8_forest_get_tree_element_offset (adapted_forest, tree_id) + lelement_index;

    cmc_assert(start_index >= 0 && count >= 0);
    cmc_assert(static_cast<size_t>(start_index + count) <= data_new_.size());

    return VectorView(&data_new_[start_index], count);

}

template <typename T>
const ElementData<T>&
AbstractEmbeddedByteCompressionVariable<T>::GetAdaptedDataValueAtIndex(const int local_index) const
{
    cmc_assert(local_index < data_new_.size());
    return data_new_[local_index];
}

template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::StoreExtractedValues(const int tree_id, const int lelement_id, const int num_elements, const ExtractionData<T>& extracted_values)
{
    cmc_assert(tree_id >= 0 && tree_id < t8_forest_get_num_local_trees(mesh_.GetMesh()));

    #if 0

    cmc_assert(static_cast<size_t>(num_elements) == extracted_values.fine_values.size());

    /* Compute the start offset in the local contiguous array */
    const int elem_start_index = t8_forest_get_tree_element_offset (mesh_.GetMesh(), tree_id) + lelement_id;

    cmc_assert(elem_start_index >= 0 && num_elements >= 0);
    cmc_assert(static_cast<size_t>(elem_start_index + num_elements) <= data_.size());

    /* Store the adjusted fine values */
    std::copy_n(extracted_values.fine_values.begin(), num_elements, &data_[elem_start_index]);

    #endif

    /* Store the extracted coarse value for the next coarsening iteration */
    data_new_.push_back(extracted_values.extracted_values);
}

template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::StoreUnchangedElement(const int tree_id, const int lelement_id, const UnchangedData<T>& extracted_value)
{
    cmc_assert(tree_id >= 0 && tree_id < t8_forest_get_num_local_trees(mesh_.GetMesh()));

    #if 0

    /* Get the element id in the contiguous array of all local elements */
    const int elem_id = t8_forest_get_tree_element_offset (mesh_.GetMesh(), tree_id) + lelement_id;
    
    cmc_assert(elem_id >= 0 && static_cast<size_t>(elem_id) < data_.size());

    /* Set the adjusted remaining element value */
    data_[elem_id] = extracted_value.fine_value;
    
    #endif

    /* Set the new value for the next coarsening iteration */
    data_new_.push_back(extracted_value.extracted_values);
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

    /* Construct the ghost layer */
    t8_forest_set_ghost (partitioned_forest, 1, T8_GHOST_FACES);

    t8_forest_commit(partitioned_forest);

    return partitioned_forest;
}

template <typename T>
inline void
AbstractEmbeddedByteCompressionVariable<T>::RepartitionData(t8_forest_t adapted_forest, t8_forest_t partitioned_forest)
{
    cmc_assert(static_cast<size_t>(t8_forest_get_local_num_leaf_elements(adapted_forest)) == data_.size());

    /* Create an sc_array_t wrapper of the variable's data */
    sc_array_t* in_data = sc_array_new_data (static_cast<void*>(data_.data()), sizeof(CompressionValue<T>), data_.size());

    cmc_debug_msg("Number of local data elements before partitioning: ", data_.size());
    cmc_debug_msg("Number of local mesh elements before partitioning: ", t8_forest_get_local_num_leaf_elements(adapted_forest));
    cmc_debug_msg("Size of a single data element: ", in_data->elem_size);

    /* Allocate memory for the partitioned data */
    const t8_locidx_t new_num_elems = t8_forest_get_local_num_leaf_elements(partitioned_forest);
    data_new_ = std::vector<CompressionValue<T>>(new_num_elems);

    cmc_debug_msg("Number of local data elements after partitioning: ", data_new_.size());
    cmc_debug_msg("Number of local mesh elements after partitioning: ", new_num_elems);

    /* Create a wrapper for the freshly allocated partitioned data */
    sc_array_t* out_data = sc_array_new_data (static_cast<void*>(data_new_.data()), sizeof(CompressionValue<T>), data_new_.size());

    cmc_assert(static_cast<size_t>(t8_forest_get_local_num_leaf_elements(adapted_forest)) == data_.size());
    cmc_assert(static_cast<size_t>(t8_forest_get_local_num_leaf_elements(partitioned_forest)) == data_new_.size());

    /* Partition the variables data */
    t8_forest_partition_data(adapted_forest, partitioned_forest, in_data, out_data);

    /* Destroy the array wrappers */
    sc_array_destroy(in_data);
    sc_array_destroy(out_data);

    /* Set the variable's data to the newly partitioned data */
    SwitchToExtractedData();
    cmc_debug_msg("Partitioning of data elements has been finished.");
}

template <typename T>
inline void
AbstractEmbeddedByteCompressionVariable<T>::RepartitionRemainingPermittedErrors(t8_forest_t adapted_forest, t8_forest_t partitioned_forest)
{
    cmc_assert(static_cast<size_t>(t8_forest_get_local_num_leaf_elements(adapted_forest)) == remaining_permitted_error_.size());

    /* Create an sc_array_t wrapper of the variable's remaining errors */
    sc_array_t* in_data = sc_array_new_data (static_cast<void*>(remaining_permitted_error_.data()), sizeof(err_type_t), remaining_permitted_error_.size());

    cmc_debug_msg("Number of local data elements before partitioning: ", remaining_permitted_error_.size());
    cmc_debug_msg("Number of local mesh elements before partitioning: ", t8_forest_get_local_num_leaf_elements(adapted_forest));
    cmc_debug_msg("Size of a single data element: ", in_data->elem_size);

    /* Allocate memory for the partitioned remaining errors */
    const t8_locidx_t new_num_elems = t8_forest_get_local_num_leaf_elements(partitioned_forest);
    remaining_permitted_error_new_ = std::vector<err_type_t>(new_num_elems);

    cmc_debug_msg("Number of local data elements after partitioning: ", remaining_permitted_error_new_.size());
    cmc_debug_msg("Number of local mesh elements after partitioning: ", new_num_elems);

    /* Create a wrapper for the freshly allocated partitioned remaining errors */
    sc_array_t* out_data = sc_array_new_data (static_cast<void*>(remaining_permitted_error_new_.data()), sizeof(err_type_t), remaining_permitted_error_new_.size());

    cmc_assert(static_cast<size_t>(t8_forest_get_local_num_leaf_elements(adapted_forest)) == remaining_permitted_error_.size());
    cmc_assert(static_cast<size_t>(t8_forest_get_local_num_leaf_elements(partitioned_forest)) == remaining_permitted_error_new_.size());

    /* Partition the variables data */
    t8_forest_partition_data(adapted_forest, partitioned_forest, in_data, out_data);

    /* Destroy the array wrappers */
    sc_array_destroy(in_data);
    sc_array_destroy(out_data);

    /* Set the variable's remaining errors to the newly partitioned remaining errors */
    remaining_permitted_error_.swap(remaining_permitted_error_new_);
    remaining_permitted_error_new_.clear();

    cmc_debug_msg("Partitioning of remaining permitted errors per element has been finished.");
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

#ifdef CMC_ENABLE_MPI
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

    /* Balance the forest */
    t8_forest_t initial_forest_balanced;
    t8_forest_init(&initial_forest_balanced);
    t8_forest_set_balance (initial_forest_balanced, initial_forest, 0);
    t8_forest_set_ghost (initial_forest_balanced, 1, T8_GHOST_FACES);
    t8_forest_commit (initial_forest_balanced);

    return AmrMesh(initial_forest_balanced, initial_refinement_level, dimensionality);
}

template <typename T>
std::vector<input::IndexReduction>
AbstractEmbeddedByteCompressionVariable<T>::UpdateLinearIndicesToTheInitialMesh()
{
    const t8_locidx_t num_local_elements = mesh_.GetNumberLocalElements();

    const int initial_refinement_level = mesh_.GetInitialRefinementLevel();

    const t8_eclass_t eclass = t8_forest_get_eclass(mesh_.GetMesh(), 0);

    /* Get the scheme of the forest's only tree */
    const t8_scheme_c* scheme =  t8_forest_get_scheme(mesh_.GetMesh());

    const t8_locidx_t first_ltree_id = 0;

    const int num_children = scheme->element_get_num_children(eclass, t8_forest_get_leaf_element_in_tree(mesh_.GetMesh(), first_ltree_id, 0));

    const MortonIndex linear_index_start_elem = GetMortonIndexOnLevel(eclass, t8_forest_get_leaf_element_in_tree(mesh_.GetMesh(), first_ltree_id, 0),
                                                   scheme, t8_eclass_to_dimension[eclass], initial_refinement_level);
    
    /* Locally, this update function reduces the global_index to zero for the first local element */
    std::vector<input::IndexReduction> index_correction;
    index_correction.emplace_back(linear_index_start_elem, linear_index_start_elem);

    /* All elements that are not holding data need to be accumulated in order to be subtracted additionally for the local index correction */
    MortonIndex skipped_indices = linear_index_start_elem;

    bool coarse_element_streak = false;

    /* Iterate through all local elements and find how the global Morton indices need to be adjusted in order to comply the local data ordering */
    for (auto iter = 0; iter < num_local_elements; ++iter)
    {
        const t8_element_t* elem = t8_forest_get_leaf_element_in_tree(mesh_.GetMesh(), 0, iter);

        if (scheme->element_get_level(eclass, elem) != initial_refinement_level)
        {
            /* Get the number of uniform indices which were skipped by this element which lays outside of the domain */
            skipped_indices += std::pow(num_children, initial_refinement_level - scheme->element_get_level(eclass, elem)) - 1;
            coarse_element_streak = true;
        } else if (coarse_element_streak)
        {
            /* We accumulate the amount of skipped indices (with regard to the initial refinement level) and store the offset 
             * once we have reached again an element on the initial refinement level */
            const MortonIndex uniform_index_of_elem = GetMortonIndexOnLevel(eclass, elem, scheme, t8_eclass_to_dimension[eclass], initial_refinement_level);
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

template <typename T>
err_type_t
GetPermittedAbsoluteErrorThreshold(const std::vector<PermittedError>& permitted_errors, const T& initial_value)
{
    err_type_t current_abs_permitted_err{std::numeric_limits<err_type_t>::max()};

    for (auto err_iter = permitted_errors.begin(); err_iter != permitted_errors.end(); ++err_iter)
    {
        if (err_iter->criterion == CompressionCriterion::AbsoluteErrorThreshold)
        {
            /* If it is an absolute error criterion, we check whether the permitted error is lower than the current one */
            if (current_abs_permitted_err > err_iter->error)
            {
                current_abs_permitted_err = err_iter->error;
            }
        } else if (err_iter->criterion == CompressionCriterion::RelativeErrorThreshold)
        {
            /* In case it is a relative error threshold, we compute the permitted absolute deviation */
            err_type_t rel_abs_err;
            if constexpr (std::is_signed_v<T>)
            {
                rel_abs_err = err_iter->error * std::abs(initial_value);
            } else
            {
                rel_abs_err = err_iter->error * initial_value;
            }

            /* If the permitted error is smaller then the currently assigned one, we update the error bound */
            if (current_abs_permitted_err > rel_abs_err)
            {
                current_abs_permitted_err = rel_abs_err;
            }
        } else
        {
            cmc_err_msg("An undefined error criterion has been supplied.");
        }
    }

    return current_abs_permitted_err;
}

template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::DetermineInitialMaximumRefinementLevelAndAbsoluteErrorBounds()
{
    t8_forest_t mesh = this->GetAmrMesh().GetMesh();
    const t8_scheme_c* scheme =  t8_forest_get_scheme(mesh);

    const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees(mesh);
    int val_idx = 0;

    /* Allocate memory for the permitted errors */
    remaining_permitted_error_.reserve(this->GetAmrMesh().GetNumberLocalElements());

    /* Iterate over all elements in all trees */
    for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx)
    {
        const t8_eclass_t tree_class = t8_forest_get_tree_class (mesh, tree_idx);
        const t8_locidx_t  num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (mesh, tree_idx);
        for (t8_locidx_t elem_idx = 0; elem_idx < num_elements_in_tree; ++elem_idx, ++val_idx)
        {
            /* Get the current element */
            const t8_element_t* element = t8_forest_get_leaf_element_in_tree (mesh, tree_idx, elem_idx);

            /* Get the level of the element */
            const int elem_level = scheme->element_get_level(tree_class, element);

            /* Check whether it is greater than the currently highest level */
            if (max_initial_refinement_level_ < elem_level)
            {
                max_initial_refinement_level_ = elem_level;
            }

            /* Check the permitted error for this element */
            std::vector<PermittedError> permitted_errors = settings_.FindRestrictingErrors(mesh, tree_idx, elem_idx, scheme, 1, &element);

            /* Compute the start offset in the local contiguous array */
            const int start_index = t8_forest_get_tree_element_offset (mesh_.GetMesh(), tree_idx) + elem_idx;

            /* Get the initial value if the element */
            const T initial_value = GetInitialDataValueAtIndex(start_index);

            /* Check whether it is an absolute and or relative error threshold */
            const err_type_t permitted_abs_error = GetPermittedAbsoluteErrorThreshold<T>(permitted_errors, initial_value);

            /* Store the permitted absolute error */
            remaining_permitted_error_.push_back(permitted_abs_error);
        }
    }

    cmc_debug_msg("The maximum present element refinement level is ", max_initial_refinement_level_);
    cmc_debug_msg("The initially permitted errors per element have been set up.");
}

template <typename T>
VectorView<err_type_t>
AbstractEmbeddedByteCompressionVariable<T>::GetViewOnRemainingPermittedErrors(const int tree_id, const int lelement_index, const int count) const
{
    cmc_assert(not remaining_permitted_error_.empty());
    cmc_assert(tree_id >= 0 && tree_id < t8_forest_get_num_local_trees(mesh_.GetMesh()));

    /* Compute the start offset in the local contiguous array */
    const int start_index = t8_forest_get_tree_element_offset (mesh_.GetMesh(), tree_id) + lelement_index;

    cmc_assert(start_index >= 0 && count >= 0);
    cmc_assert(static_cast<size_t>(start_index + count) <= remaining_permitted_error_.size());

    return VectorView(&remaining_permitted_error_[start_index], count);
}

template <typename T>
err_type_t
AbstractEmbeddedByteCompressionVariable<T>::GetAdaptedRemainingPermittedError(t8_forest_t adapted_forest, const int tree_id, const int lelement_index) const
{
    cmc_assert(not remaining_permitted_error_new_.empty());
    cmc_assert(tree_id >= 0 && tree_id < t8_forest_get_num_local_trees(adapted_forest));

    /* Compute the start offset in the local contiguous array */
    const int start_index = t8_forest_get_tree_element_offset (adapted_forest, tree_id) + lelement_index;

    cmc_assert(start_index >= 0);
    cmc_assert(static_cast<size_t>(start_index) < remaining_permitted_error_new_.size());

    return remaining_permitted_error_new_[start_index];
}

template <typename T>
err_type_t
AbstractEmbeddedByteCompressionVariable<T>::GetRemainingPermittedError(const int tree_id, const int lelement_index) const
{
    cmc_assert(not remaining_permitted_error_.empty());
    cmc_assert(tree_id >= 0 && tree_id < t8_forest_get_num_local_trees(mesh_.GetMesh()));

    /* Compute the start offset in the local contiguous array */
    const int start_index = t8_forest_get_tree_element_offset (mesh_.GetMesh(), tree_id) + lelement_index;

    cmc_assert(start_index >= 0);
    cmc_assert(static_cast<size_t>(start_index) < remaining_permitted_error_.size());

    return remaining_permitted_error_[start_index];
}

template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::StoreAdaptedRemainingPermittedError(const err_type_t remaining_error)
{
    /* Store the adapted remaining error */
    remaining_permitted_error_new_.push_back(remaining_error);
}

/* This method is only intended to be called after the the remaining errors of the previous forest has been set up to resemble the remaining errors of the adapted forest */
template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::UpdateNewRemainingPermittedError(t8_forest_t adapted_forest, const int tree_id, const int lelement_index, const err_type_t new_remaining_error)
{
    cmc_assert(tree_id >= 0 && tree_id < t8_forest_get_num_local_trees(adapted_forest));
    cmc_assert(t8_forest_get_tree_element_offset (adapted_forest, tree_id) + lelement_index < t8_forest_get_local_num_leaf_elements(adapted_forest));

    /* Compute the start offset in the local contiguous array */
    const int start_index = t8_forest_get_tree_element_offset (adapted_forest, tree_id) + lelement_index;

    cmc_assert(start_index >= 0);
    cmc_assert(static_cast<size_t>(start_index) < remaining_permitted_error_.size());

    remaining_permitted_error_[start_index] = new_remaining_error;
}

/* This method is only intended to be called after the the remaining errors of the previous forest has been set up to resemble the remaining errors of the adapted forest */
template <typename T>
err_type_t
AbstractEmbeddedByteCompressionVariable<T>::GetNewRemainingPermittedError(t8_forest_t adapted_forest, const int tree_id, const int lelement_index) const
{
    cmc_assert(not remaining_permitted_error_.empty());
    cmc_assert(tree_id >= 0 && tree_id < t8_forest_get_num_local_trees(adapted_forest));

    /* Compute the start offset in the local contiguous array */
    const int start_index = t8_forest_get_tree_element_offset (adapted_forest, tree_id) + lelement_index;

    cmc_assert(start_index >= 0);
    cmc_assert(static_cast<size_t>(start_index) < remaining_permitted_error_.size());

    return remaining_permitted_error_[start_index];
}

template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::ClearAdaptedRemainingPermittedErrors()
{
    remaining_permitted_error_new_.clear();
}

template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::ReserveAdaptedRemainingPermittedError(const int num_elements_coarsened_forest)
{
    remaining_permitted_error_new_.reserve(num_elements_coarsened_forest);
}

template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::SwitchToAdaptedRemainingErrors()
{
    remaining_permitted_error_ = remaining_permitted_error_new_;
}

}

#endif /* !LOSSY_CMC_EMBEDDED_IDW_COMPRESSION_VARIABLE_HXX */
