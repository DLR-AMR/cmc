#ifndef CMC_T8_DATA_VARIABLES_HXX
#define CMC_T8_DATA_VARIABLES_HXX
/**
 * @file cmc_t8_data_variables.hxx
 */
#include "t8code/cmc_t8_data_variables_forward.hxx"
#include "utilities/cmc_variable_transformer_forward.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_interpolation.hxx"
#include "t8code/cmc_t8_adapt_track_inaccuracy.hxx"
#include "utilities/cmc_compression_settings.hxx"
#include "utilities/cmc_prefix.hxx"
#include "t8code/cmc_t8_variables.hxx"

#ifdef CMC_WITH_NETCDF
#include "netcdf/cmc_netcdf.hxx"
#endif

#include <cstdlib>
#include <memory>
#include <type_traits>
#include <vector>

namespace cmc
{

template<class T>
class Variable
{
public:
    Variable() = default;
    ~Variable(){};

    Variable(const Variable& other) = default;
    Variable& operator=(const Variable& other) = default;
    Variable(Variable&& other) = default;
    Variable& operator=(Variable&& other) = default;

    template<typename U> auto operator[](U index) const -> std::enable_if_t<std::is_integral_v<U>, T> {return data_[index];};
    template<typename U> auto operator[](U index) -> std::enable_if_t<std::is_integral_v<U>, T&> {return data_[index];};

    size_t size() const {return data_.size();};
    void push_back(const T& value) {data_.push_back(value);};
    void push_back(T&& value) {data_.push_back(std::move(value));};
    template<typename U> auto push_back(const U& value) const -> std::enable_if_t<!std::is_same_v<U,T>, void> {cmc_err_msg("The variable is of a different data type.");};

    void SetName(const std::string& name) {name_ = name;};
    const std::string& GetName() const {return name_;};

    bool IsValidForCompression() const;
    AmrMesh& GetAmrMesh() {return mesh_;};
    const AmrMesh& GetAmrMesh() const {return mesh_;};
    VariableAttributes<T>& GetVariableAttributes() {return attributes_;};
    const VariableAttributes<T>& GetVariableAttributes() const {return attributes_;};
    VariableUtilities<T>& GetVariableUtilities() {return utilities_;};
    const VariableUtilities<T>& GetVariableUtilities() const {return utilities_;};

    void SetInterpolation(Interpolate2<T> interpolation_function) {utilities_.SetInterpolation(interpolation_function);};
    void SetInaccuracyStorage(const TrackingOption& inaccuracy_tracking_option) {utilities_.SetInaccuracyStorage(inaccuracy_tracking_option);};
    void SetUpInaccuracyStorage(const size_t size_hint = kInvalidSizeHintForInaccuracyContainer) {utilities_.SetUpInaccuracyStorage(size_hint);};
    void SetUpCompressionCriteria(const CompressionSpecifications& variable_specifications);

    void ApplyInterpolationResult(const int index, const ErrorCompliance& evaluation) {utilities_.StoreInaccuracy(index, evaluation.max_introduced_error);};
    void PopInterpolationResult() {data_new_.pop_back();};

    void LeaveElementUnchanged(const int index_previous_values);

    void AllocateForDecompression(const size_t size_hint);
    void DecompressElementConstantly(const int index, const int num_insertions);

    ErrorCompliance EvaluateCoarsening(const std::vector<PermittedError>& permitted_errors, const t8_forest_t previous_mesh, const t8_locidx_t start_index, const int num_elements);

    void InitializeVariableForCompressionIteration();
    void UpdateCompressionData();
    void UpdateDecompressionData();

    const GeoDomain& GetGlobalDomain() const {return global_domain_;};
    CmcUniversalType GetMissingValue() const {return CmcUniversalType(attributes_.GetMissingValue());};

    DataLayout GetInitialDataLayout() const {return attributes_.GetInitialDataLayout();};
    Dimension HasSplitDimension() const {return attributes_.HasSplitDimension();};
    DataLayout GetPreCompressionLayout() const {return attributes_.GetPreCompressionLayout();};
    int GetGlobalContextInformation() const {return attributes_.GetGlobalContextInformation();};

    std::vector<PermittedError> GetPermittedError(const int num_elements, const t8_element_t* elements[], t8_eclass_scheme_c * ts) const;

    CmcUniversalType GetValueAt(const int index) const;
    std::vector<double> GetDataAsDoubleVector() const;
    std::vector<CompressionValue<sizeof(T)>> GetDataAsCompressionValues() const;

    std::vector<ErrorCompliance> CheckErrorBoundsForValues(const std::vector<T>& alternative_values, const int index) const;
    double  GetRemainingMaxAllowedAbsoluteError(const int index) const;
    void SetVariableDataInNCFile(const int ncid, const int var_id) const;

    void RepartitionData(t8_forest_t adapted_forest, t8_forest_t partitioned_forest);

    friend class TransformerInputToCompressionVariable;
    friend class TransformerCompressionToOutputVariable;
    friend class TransformerCompressionToByteVariable;
    friend class Var;

private:
    void SwitchToDataNew();
    void StoreInterpolationResult(const T interpolated_value);
    void SetMissingValueInNCFile(const int ncid, const int var_id, const int nc_type) const;


    std::string name_; //!< The name of the variable
    std::vector<T> data_; //!< The actual data of the variable
    std::vector<T> data_new_;
    AmrMesh mesh_; //!< The mesh on which the variable is defined
    VariableAttributes<T> attributes_;
    VariableUtilities<T> utilities_;
    GeoDomain global_domain_;
};


class Var
{
public:
    Var() = delete;
    ~Var() = default;

    Var(const int id, const CmcType type, CmcVariable&& variable)
    : id_{id}, type_{type}, var_{std::move(variable)}{};
    Var(const int id, const int internal_id, const CmcType type, CmcVariable&& variable)
    : id_{id}, internal_id_{internal_id}, type_{type}, var_{std::move(variable)}{};
    Var(const Var& other) = default;
    Var(Var&& other) = default;
    
    Var& operator=(const Var& other);
    Var& operator=(Var&& other) = default;

    int GetID() const;
    int GetInternalID() const;
    const std::string& GetName() const;
    const GeoDomain& GetGlobalDomain() const;
    DataLayout GetInitialDataLayout() const;
    CmcType GetType() const {return type_;};
    AmrMesh& GetAmrMesh();
    const AmrMesh& GetAmrMesh() const;
    Dimension HasSplitDimension() const;
    CmcUniversalType GetMissingValue() const;
    DataLayout GetPreCompressionLayout() const;
    int GetGlobalContextInformation() const;

    void SetAmrMesh(const AmrMesh& mesh);

    void SetUpCompressionCriteria(const CompressionSpecifications& variable_specifications);
    void SetUpInaccuracyStorage(const size_t size_hint = kInvalidSizeHintForInaccuracyContainer);

    std::vector<PermittedError> GetPermittedError(const int num_elements, const t8_element_t* elements[], t8_eclass_scheme_c * ts) const;
    ErrorCompliance EvaluateCoarsening(const std::vector<PermittedError>& permitted_errors, const t8_forest_t previous_mesh, const t8_locidx_t start_index, const int num_elements);

    void ApplyInterpolation(const int lelement_id, const ErrorCompliance& evaluation);
    void PopInterpolation();

    void LeaveElementUnchanged(const t8_locidx_t index_previous_values);

    void UpdateCompressionData();
    void UpdateDecompressionData();

    void InitializeVariableForCompressionIteration();
    void AllocateForDecompression(const size_t size_hint);

    void DecompressElementConstantly(const int index, const int num_insertions);

    size_t size() const;
    CmcUniversalType GetValueAt(const int index) const;
    std::vector<double> GetDataAsDoubleVector() const;

    bool IsValidForCompression() const;

    const CmcVariable& GetInternalVariable () const;
    void SetVariableDataInNCFile(const int ncid, const int var_id) const;
    void SetMissingValueInNCFile(const int ncid, const int var_id, const int nc_type) const;

    void RepartitionData(t8_forest_t adapted_forest, t8_forest_t partitioned_forest);

    friend class TransformerInputToCompressionVariable;
    friend class TransformerCompressionToOutputVariable;
    friend class TransformerCompressionToByteVariable;
private:

    int id_;
    int internal_id_;
    CmcType type_;
    CmcVariable var_;
};

/** VARIABLE<T> MEMBER FUNCITONS **/
template<typename T>
void
Variable<T>::SetUpCompressionCriteria(const CompressionSpecifications& variable_specifications)
{
    utilities_.SetUpVariableErrorDomains(variable_specifications, global_domain_);
}

template<typename T>
bool 
Variable<T>::IsValidForCompression() const
{
    const VariableUtilities<T>& util = GetVariableUtilities();
    const AmrMesh& mesh = GetAmrMesh();
    
    if (!util.IsInaccuracyStorageAllocated())
        return false;
    if (static_cast<size_t>(mesh.GetNumberLocalElements()) != data_.size())
        return false;
    if (!data_new_.empty())
        return false;
    if (!global_domain_.IsValid())
        return false;
    return true;
}


template<typename T>
std::vector<PermittedError>
Variable<T>::GetPermittedError(const int num_elements, const t8_element_t* elements[], t8_eclass_scheme_c * ts) const
{
    /* Iterate over all error domains for all elements and find the minimum error */
    bool is_abs_error_present{false};
    bool is_rel_error_present{false};
    double abs_err{std::numeric_limits<double>::max()};
    double rel_err{std::numeric_limits<double>::max()};

    /* The first error domain is always the general criterion on the whole domain */
    auto ed_iter = utilities_.GetErrorDomainsBegin();
    if (IsAnyElementWithinGlobalDomain(num_elements, elements, ts, ed_iter->GetDomain(), mesh_.GetInitialRefinementLevel(), attributes_.GetInitialDataLayout()))
    {
        if (ed_iter->GetCriterion() == CompressionCriterion::AbsoluteErrorThreshold)
        {
            is_abs_error_present = true;
            if (ed_iter->GetError() < abs_err)
            {
                abs_err = ed_iter->GetError();
            }
        } else if (ed_iter->GetCriterion() == CompressionCriterion::RelativeErrorThreshold)
        {
            is_rel_error_present = true;
            if (ed_iter->GetError() < rel_err)
            {
                rel_err = ed_iter->GetError();
            }
        }
    }

    /* Check all additional error domains */
    for (++ed_iter; ed_iter != utilities_.GetErrorDomainsEnd(); ++ed_iter)
    {
        if (IsAnyElementWithinGeoDomain(num_elements, elements, ts, ed_iter->GetDomain(), mesh_.GetInitialRefinementLevel(), attributes_.GetInitialDataLayout()))
        {
            if (ed_iter->GetCriterion() == CompressionCriterion::AbsoluteErrorThreshold)
            {
                is_abs_error_present = true;
                if (ed_iter->GetError() < abs_err)
                {
                    abs_err = ed_iter->GetError();
                }
            } else if (ed_iter->GetCriterion() == CompressionCriterion::RelativeErrorThreshold)
            {
                is_rel_error_present = true;
                if (ed_iter->GetError() < rel_err)
                {
                    rel_err = ed_iter->GetError();
                }
            }
        }
    }

    if (is_abs_error_present && is_rel_error_present)
    {
        return std::vector<PermittedError>{PermittedError{CompressionCriterion::AbsoluteErrorThreshold, abs_err}, PermittedError{CompressionCriterion::RelativeErrorThreshold, rel_err}};
    } else if (is_abs_error_present)
    {
        return std::vector<PermittedError>{PermittedError{CompressionCriterion::AbsoluteErrorThreshold, abs_err}};
    } else if (is_rel_error_present)
    {
        return std::vector<PermittedError>{PermittedError{CompressionCriterion::RelativeErrorThreshold, rel_err}};
    } else
    {
        cmc_err_msg("There could not be an error associated to the supplied elements.");
        return std::vector<PermittedError>();
    }
}

template<typename T>
void
Variable<T>::LeaveElementUnchanged(const int index_previous_values)
{
    data_new_.push_back(data_[index_previous_values]);
    utilities_.TransferPreviousDeviation(index_previous_values);
}

template<typename T>
void
Variable<T>::AllocateForDecompression(const size_t size_hint)
{
    data_new_.clear();
    data_new_.reserve(size_hint);
}

template<typename T>
void
Variable<T>::DecompressElementConstantly(const int index, const int num_insertions)
{
    std::fill_n(std::back_inserter(data_new_), num_insertions, data_[index]);
}

template<typename T>
ErrorCompliance
Variable<T>::EvaluateCoarsening(const std::vector<PermittedError>& permitted_errors, const t8_forest_t previous_mesh, const t8_locidx_t start_index, const int num_elements)
{
    const T missing_value = attributes_.GetMissingValue();
    
    /* Obtain the values to be interpolated */
    const VectorView<T> values(data_.data() + start_index, num_elements);

    const T interpolated_value = utilities_.Interpolate(values, previous_mesh, start_index, num_elements, missing_value);

    const std::vector<double> previous_deviations = utilities_.GetPreviousDeviations(start_index, num_elements);

    const ErrorCompliance error_evaluation = utilities_.IsCoarseningErrorCompliant(permitted_errors, values, previous_deviations, interpolated_value, missing_value);
    
    StoreInterpolationResult(interpolated_value);
    
    return error_evaluation;
}

template<typename T>
void
Variable<T>::UpdateCompressionData()
{
    SwitchToDataNew();
    utilities_.SwitchDeviations();
}

template<typename T>
void
Variable<T>::UpdateDecompressionData()
{
    SwitchToDataNew();
}

template<typename T>
void
Variable<T>::InitializeVariableForCompressionIteration()
{
    const t8_locidx_t num_local_elements = GetAmrMesh().GetNumberLocalElements();
    data_new_.reserve(num_local_elements);
    utilities_.AllocateDeviationStorage(num_local_elements);
}

template<typename T>
CmcUniversalType 
Variable<T>::GetValueAt(const int index) const
{
    cmc_assert(data_.size() > static_cast<size_t>(index));
    return CmcUniversalType(data_[index]);
}

template<typename T>
std::vector<double>
Variable<T>::GetDataAsDoubleVector() const
{
    if constexpr (std::is_same_v<T, double>)
    {
        return data_;
    } else
    {
        std::vector<double> double_data;
        double_data.reserve(data_.size());

        for (auto iter = data_.begin(); iter != data_.end(); ++iter)
        {
            if (!std::isnan(*iter))
            {
                double_data.push_back(static_cast<double>(*iter));
            } else
            {
                double_data.push_back(attributes_.GetMissingValue());
            }
        }

        return double_data;
    }
}

template<typename T>
void
Variable<T>::SetVariableDataInNCFile(const int ncid, const int var_id) const
{
    #ifdef CMC_WITH_NETCDF
    if (data_.size() > 0)
    {
        int err = nc_put_var(ncid, var_id, static_cast<const void*>(data_.data()));
        NcCheckError(err);
    } else
    {
        /* We set a single missing value within the file (in order to not define a unlimited dimension for the variable) */
        T missing_val = attributes_.GetMissingValue();
        int err = nc_put_var(ncid, var_id, static_cast<const void*>(&missing_val));
        NcCheckError(err);
    }
    #endif
}

template<typename T>
void
Variable<T>::SwitchToDataNew()
{
    data_.swap(data_new_);
    data_new_.clear();
}

template<typename T>
void
Variable<T>::StoreInterpolationResult(const T interpolated_value)
{
    data_new_.push_back(interpolated_value);
}

template<typename T>
void
Variable<T>::SetMissingValueInNCFile(const int ncid, const int var_id, const int nc_type) const
{
    #ifdef CMC_WITH_NETCDF
    attributes_.SetMissingValueInNCFile(ncid, var_id, nc_type);
    #endif
}

template<typename T>
std::vector<ErrorCompliance>
Variable<T>::CheckErrorBoundsForValues(const std::vector<T>& alternative_values, const int index) const
{
    const t8_eclass_t eclass = t8_forest_get_eclass(mesh_.GetMesh(), 0);

    const t8_element_t* element = t8_forest_get_element_in_tree(mesh_.GetMesh(), 0, index);
    
    t8_eclass_scheme_c* ts =  t8_forest_get_eclass_scheme (mesh_.GetMesh(), eclass);

    const std::vector<PermittedError> permitted_errors = GetPermittedError(1, &element, ts);

    const std::vector<double> previous_deviations = utilities_.GetPreviousDeviations(index, 1);

    return utilities_.AreAlternativeValuesErrorCompliant(alternative_values, permitted_errors, VectorView<T>(data_.data() + index, 1), previous_deviations, attributes_.GetMissingValue());
}

template<typename T>
double 
Variable<T>::GetRemainingMaxAllowedAbsoluteError(const int index) const
{
    const t8_eclass_t eclass = t8_forest_get_eclass(mesh_.GetMesh(), 0);

    const t8_element_t* element = t8_forest_get_element_in_tree(mesh_.GetMesh(), 0, index);
    
    t8_eclass_scheme_c* ts =  t8_forest_get_eclass_scheme (mesh_.GetMesh(), eclass);

    const std::vector<PermittedError> permitted_errors = GetPermittedError(1, &element, ts);

    const double current_abs_deviation = utilities_.GetPreviousDeviation(index);

    //TODO: Continue
    
    return 0.0;

}

template<typename T>
std::vector<CompressionValue<sizeof(T)>>
Variable<T>::GetDataAsCompressionValues() const
{
    std::vector<CompressionValue<sizeof(T)>> compression_data;
    compression_data.reserve(data_.size());

    for (auto iter = data_.begin(); iter != data_.end(); ++iter)
    {
        compression_data.emplace_back(*iter);
    }
    
    return compression_data;
}

#if 0
/* Create an sc_array_t wrapper of the variable's data */
    //sc_array_t in_data;
    //in_data.elem_size = sizeof(T);
    //in_data.elem_count = data_.size();
    //in_data.byte_alloc = static_cast<ssize_t>(in_data.elem_size * in_data.elem_count);
    //in_data.array = reinterpret_cast<char*>(data_.data());
    //in_data.array = static_cast<char*>(static_cast<void*>(data_.data()));

/* Create a wrapper for the partitioned data */
    //sc_array_t out_data;
    //out_data.elem_size = sizeof(T);
    //out_data.elem_count = data_new_.size();
    //out_data.byte_alloc = static_cast<ssize_t>(out_data.elem_size * out_data.elem_count);
    //out_data.array = reinterpret_cast<char*>(data_new_.data());
    //out_data.array = static_cast<char*>(static_cast<void*>(data_new_.data()));

#endif

template<typename T>
void
Variable<T>::RepartitionData(t8_forest_t adapted_forest, t8_forest_t partitioned_forest)
{
    /* Create an sc_array_t wrapper of the variable's data */
    sc_array_t* in_data = sc_array_new_data (static_cast<void*>(data_.data()), sizeof(T), data_.size());

    cmc_debug_msg("Size of data_: ", data_.size());
    cmc_debug_msg("Num local elems: ", t8_forest_get_local_num_elements(adapted_forest));
    cmc_debug_msg("In Elem size: ", in_data->elem_size);

    /* Allocate memory for the partitioned data */
    const t8_locidx_t new_num_elems = t8_forest_get_local_num_elements(partitioned_forest);
    data_new_ = std::vector<T>(new_num_elems);

    cmc_debug_msg("Size of data_new: ", data_new_.size());
    cmc_debug_msg("Num local elems (partitioned_foerst): ", new_num_elems);

    /* Create a wrapper for the freshly allocated partitioned data */
    sc_array_t* out_data = sc_array_new_data (static_cast<void*>(data_new_.data()), sizeof(T), data_new_.size());

    cmc_debug_msg("Out Elem size: ", out_data->elem_size);

    /* Partition the variables data */
    t8_forest_partition_data(adapted_forest, partitioned_forest, in_data, out_data);

    /* Destroy the array wrappers */
    sc_array_destroy(in_data);
    sc_array_destroy(out_data);

    /* Set the variable's data to the newly partitioned data */
    SwitchToDataNew();
    cmc_debug_msg("Partition of variable data finsihed");
    //cmc_err_msg("We stop here");
    
    /* Repartition the inaccuracy tracker */
    utilities_.RepartitionInaccuracyData(adapted_forest, partitioned_forest);
}


/** VAR MEMBER FUNCTIONS **/
inline
Var& 
Var::operator=(const Var& other)
{
    id_ = other.id_;
    type_ = other.type_;
    var_ = std::move(other.var_);
    return *this;
};

inline
AmrMesh& 
Var::GetAmrMesh()
{
    return std::visit([](auto&& arg) -> AmrMesh& {
        return arg.GetAmrMesh();
    }, var_);
}

inline
const AmrMesh& 
Var::GetAmrMesh() const
{
    return std::visit([](auto&& arg) -> const AmrMesh& {
        return arg.GetAmrMesh();
    }, var_);
}

inline
void 
Var::SetAmrMesh(const AmrMesh& mesh)
{
    std::visit([&](auto&& arg) {
        arg.mesh_ = mesh;
    }, var_);
}

inline
int
Var::GetID() const
{
    return id_;
}

inline
int
Var::GetInternalID() const
{
    return internal_id_;
}

inline
const std::string&
Var::GetName() const
{
    return std::visit([](auto&& var) -> const std::string& {
        return var.GetName();
    }, var_);
}

inline
const GeoDomain&
Var::GetGlobalDomain() const
{
    return std::visit([](auto&& var) -> const GeoDomain& {
        return var.GetGlobalDomain();
    }, var_);
}

inline
DataLayout
Var::GetInitialDataLayout() const 
{
    return std::visit([](auto&& var) -> DataLayout {
        return var.GetInitialDataLayout();
    }, var_);
}

inline
void
Var::SetUpCompressionCriteria(const CompressionSpecifications& variable_specifications)
{
    std::visit([&](auto&& var){
        var.SetUpCompressionCriteria(variable_specifications);
    }, var_);
}

inline
void 
Var::SetUpInaccuracyStorage(const size_t size_hint)
{
    std::visit([=](auto&& var){
        var.SetUpInaccuracyStorage(size_hint);
    }, var_);
}

inline
bool
Var::IsValidForCompression() const
{
    return std::visit([=](auto&& var) -> bool {
        return var.IsValidForCompression();
    }, var_);
}

inline
std::vector<PermittedError>
Var::GetPermittedError(const int num_elements, const t8_element_t* elements[], t8_eclass_scheme_c * ts) const
{
    return std::visit([=](auto&& var) -> std::vector<PermittedError> {
        return var.GetPermittedError(num_elements, elements, ts);
    }, var_);
}

inline
ErrorCompliance
Var::EvaluateCoarsening(const std::vector<PermittedError>& permitted_errors, const t8_forest_t previous_mesh, const t8_locidx_t start_index, const int num_elements)
{
    return std::visit([&](auto&& var) -> ErrorCompliance {
        return var.EvaluateCoarsening(permitted_errors, previous_mesh, start_index, num_elements);
    }, var_);
}

inline
void
Var::ApplyInterpolation(const int lelement_id, const ErrorCompliance& evaluation)
{
    std::visit([&](auto&& var) {
        var.ApplyInterpolationResult(lelement_id, evaluation);
    }, var_);
}

inline
void
Var::PopInterpolation()
{
    std::visit([](auto&& var) {
        var.PopInterpolationResult();
    }, var_);
}

inline
void
Var::LeaveElementUnchanged(const t8_locidx_t index_previous_values)
{
    std::visit([=](auto&& var) {
        var.LeaveElementUnchanged(index_previous_values);
    }, var_);
}

inline
void
Var::UpdateCompressionData()
{
    std::visit([](auto&& var) {
        var.UpdateCompressionData();
    }, var_);
}

inline
void 
Var::UpdateDecompressionData()
{
    std::visit([](auto&& var) {
        var.UpdateDecompressionData();
    }, var_);
}

inline
void
Var::InitializeVariableForCompressionIteration()
{
    std::visit([](auto&& var) {
        var.InitializeVariableForCompressionIteration();
    }, var_);
}

inline
void
Var::AllocateForDecompression(const size_t size_hint)
{
    std::visit([&](auto&& var) {
        var.AllocateForDecompression(size_hint);
    }, var_);
}

inline
void
Var::DecompressElementConstantly(const int index, const int num_insertions)
{
    std::visit([&](auto&& var) {
        var.DecompressElementConstantly(index, num_insertions);
    }, var_);
}

inline
Dimension
Var::HasSplitDimension() const
{
    return std::visit([](auto&& var) -> Dimension {
        return var.HasSplitDimension();
    }, var_);
}

inline
CmcUniversalType
Var::GetMissingValue() const
{
    return std::visit([](auto&& var) -> CmcUniversalType {
        return var.GetMissingValue();
    }, var_);
}
    
inline
DataLayout
Var::GetPreCompressionLayout() const 
{
    return std::visit([](auto&& var) -> DataLayout {
        return var.GetPreCompressionLayout();
    }, var_);
};

inline
CmcUniversalType
Var::GetValueAt(const int index) const
{
    return std::visit([&](auto&& var) -> CmcUniversalType {
        return var.GetValueAt(index);
    }, var_);
}

inline
size_t
Var::size() const
{
    return std::visit([](auto&& var) -> size_t {
        return var.size();
    }, var_);
}

inline
int
Var::GetGlobalContextInformation() const
{
    return std::visit([](auto&& var) -> int {
        return var.GetGlobalContextInformation();
    }, var_);
}

inline
void
Var::SetVariableDataInNCFile(const int ncid, const int var_id) const
{
    #ifdef CMC_WITH_NETCDF
    std::visit([&](auto&& var) {
        var.SetVariableDataInNCFile(ncid, var_id);
    }, var_);
    #endif
}

inline
void
Var::SetMissingValueInNCFile(const int ncid, const int var_id, const int nc_type) const
{
    #ifdef CMC_WITH_NETCDF
    std::visit([&](auto&& var) {
        var.SetMissingValueInNCFile(ncid, var_id, nc_type);
    }, var_);
    #endif
}

inline
std::vector<double>
Var::GetDataAsDoubleVector() const
{
    return std::visit([](const auto& var) -> std::vector<double> {
        return var.GetDataAsDoubleVector();
    }, var_);
}

inline
const CmcVariable&
Var::GetInternalVariable() const
{
    return var_;
}

inline
void
Var::RepartitionData(t8_forest_t adapted_forest, t8_forest_t partitioned_forest)
{
    std::visit([&](auto& var) {
        var.RepartitionData(adapted_forest, partitioned_forest);
    }, var_);
}

}

#endif /* !CMC_T8_DATA_VARIABLES_HXX */
