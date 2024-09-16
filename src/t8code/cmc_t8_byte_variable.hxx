#ifndef CMC_T8_BYTE_VARIABLE_HXX
#define CMC_T8_BYTE_VARIABLE_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "t8code/cmc_t8_data_variables_forward.hxx"
#include "utilities/cmc_variable_transformer_forward.hxx"
#include "utilities/cmc_prefix.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_variables.hxx"
#include "utilities/cmc_span.hxx"
#include "t8code/cmc_t8_prefix_encoding.hxx"
#include <t8_forest/t8_forest_vtk.h>


#include <bitset>
#include <vector>

#ifdef CMC_WITH_NETCDF
#include "netcdf/cmc_netcdf.hxx"
#endif


namespace cmc
{

template<typename T> class ByteVariable;
class ByteVar;

using CmcByteVariable = std::variant<ByteVariable<int8_t>, ByteVariable<char>, ByteVariable<int16_t>, ByteVariable<int32_t>, ByteVariable<float>, ByteVariable<double>,
                                     ByteVariable<uint8_t>, ByteVariable<uint16_t>, ByteVariable<uint32_t>, ByteVariable<int64_t>, ByteVariable<uint64_t>>;


template<typename T>
class ByteVariable
{
public:
    ByteVariable() = default;
    ByteVariable(std::vector<CompressionValue<sizeof(T)>>&& serialized_data)
    : byte_values_(std::move(serialized_data)) {};
    ~ByteVariable() {
        ReleaseInitialMesh();
    }

    ByteVariable(const ByteVariable& other) = default;
    ByteVariable& operator=(const ByteVariable& other) = default;
    ByteVariable(ByteVariable&& other) = default;
    ByteVariable& operator=(ByteVariable&& other) = default;

    const std::string& GetName() const {return name_;}

    void PerformTailTruncation();

    void InitializeCompressionIteration();
    bool HasAtLeastOnePrefixBeenFoundInTheLastIteration() const;
    void IndicatePrefixFound(const CompressionValue<sizeof(T)>& prefix);
    void IndicatePrefixFound(CompressionValue<sizeof(T)>&& prefix);
    void IndicateNoPrefixFound();
    void PopLevelwisePrefixEncoding();
    void LeaveCoarsePrefixUnchanged(const int elem_index);
    void LeaveInitialValueUnchanged(const int elem_index);
    bool EvaluateCommonPrefixFromInitialData(const int elem_start_index, const int num_elements);
    bool EvaluateCommonPrefixFromPreviousPrefix(const int elem_start_index, const int num_elements);
    void IndicatePrefixFoundEGU(const CompressionValue<sizeof(T)>& prefix);
    void IndicateNoPrefixFoundEGU();
    void LeaveCoarsePrefixUnchangedEGU(const int elem_index);
    void EvaluateCommonPrefixFromInitialDataEGU(const int elem_start_index, const int num_elements);
    void EvaluateCommonPrefixFromPreviousPrefixEGU(const int elem_start_index, const int num_elements);
    void FinalizeCompressionIterationEGU();
    t8_forest_t GetMesh() const;
    void SetMesh(t8_forest_t forest);
    AmrMesh& GetAmrMesh() {return mesh_;};
    const AmrMesh& GetAmrMesh() const {return mesh_;};
    DataLayout GetInitialDataLayout() const {return attributes_.GetInitialDataLayout();};
    const GeoDomain& GetGlobalDomain() const {return global_domain_;};
    int GetGlobalContextInformation() const {return attributes_.GetGlobalContextInformation();};
    DataLayout GetPreCompressionLayout() const {return attributes_.GetPreCompressionLayout();};
    void SetMissingValueInNCFile(const int ncid, const int var_id, const int nc_type) const;
    int GetNumberOfSetPrefixIndicationBits() const;
    std::tuple<int, int, int, int, int, int, int, std::vector<uint8_t>> GatherSerializedCompressionData() const;
    std::tuple<int, int, int, int, int, int, std::vector<uint8_t>> GatherSerializedCompressionDataEGU() const;
    void StoreInitialMesh();
    void ReleaseInitialMesh(); 
    friend class TransformerCompressionToByteVariable;

private:
    CompressionValue<sizeof(T)> GetMaximumTailClearedValue(const int, const std::vector<PermittedError>&, const CompressionValue<sizeof(T)>&, const T&) const;
    CompressionValue<sizeof(T)> GetMaximumTailToggledValue(const int, const std::vector<PermittedError>&, const CompressionValue<sizeof(T)>&, const T&) const;
    std::vector<PermittedError> GetPermittedError(const int index) const;
    std::vector<PermittedError> GetRemainingMaxAllowedAbsoluteError(const int index) const;
    
    std::string name_; //!< The name of the variable
    std::vector<T> initial_data_; //!< The actual data of the variable
    AmrMesh mesh_; //!< The mesh on which the variable is defined
    VariableAttributes<T> attributes_;
    VariableUtilities<T> utilities_;
    GeoDomain global_domain_;

    std::vector<CompressionValue<sizeof(T)>> byte_values_;
    std::vector<LevelwisePrefixes<T>> prefixes_;

    std::vector<CompressionValue<sizeof(T)>> previous_prefixes_;
    std::vector<CompressionValue<sizeof(T)>> previous_prefixes_new_;

    bool is_initial_mesh_stored_{true};//TODO:revert to false
    t8_forest_t initial_mesh_{nullptr};
};

class ByteVar
{
public:
    ByteVar() = delete;
    ByteVar(const int id, const CmcType type, CmcByteVariable&& variable)
    : id_{id}, type_{type}, var_{std::move(variable)}{};

    ~ByteVar() = default;

    const std::string& GetName() const;
    int GetID() const {return id_;}
    CmcType GetType() const {return type_;}
    void PerformTailTruncation();
    void InitializeCompressionIteration();
    bool HasAtLeastOnePrefixBeenFoundInTheLastIteration() const;
    void IndicateNoPrefixFound();
    void PopLevelwisePrefixEncoding();
    void LeaveCoarsePrefixUnchanged(const int elem_index);
    void LeaveInitialValueUnchanged(const int elem_index);
    bool EvaluateCommonPrefixFromInitialData(const int elem_start_index, const int num_elements);
    bool EvaluateCommonPrefixFromPreviousPrefix(const int elem_start_index, const int num_elements);
    void EvaluateCommonPrefixFromInitialDataEGU(const int elem_start_index, const int num_elements);
    void EvaluateCommonPrefixFromPreviousPrefixEGU(const int elem_start_index, const int num_elements);
    void LeaveCoarsePrefixUnchangedEGU(const int elem_index);
    void FinalizeCompressionIterationEGU();
    template<typename T> void IndicatePrefixFound(const CompressionValue<sizeof(T)>& prefix);
    template<typename T> void IndicatePrefixFound(CompressionValue<sizeof(T)>&& prefix);
    void SetMesh(t8_forest_t forest);
    t8_forest_t GetMesh() const;
    AmrMesh& GetAmrMesh();
    const AmrMesh& GetAmrMesh() const;
    DataLayout GetInitialDataLayout() const;
    const GeoDomain& GetGlobalDomain() const;
    int GetGlobalContextInformation() const;
    DataLayout GetPreCompressionLayout() const;
    void SetMissingValueInNCFile(const int ncid, const int var_id, const int nc_type) const;
    void StoreInitialMesh();
    int GetNumberOfSetPrefixIndicationBits() const;
    std::tuple<int, int, int, int, int, int, int, std::vector<uint8_t>> GatherSerializedCompressionData() const;
    std::tuple<int, int, int, int, int, int, std::vector<uint8_t>> GatherSerializedCompressionDataEGU() const;
private:
    int id_;
    CmcType type_;
    CmcByteVariable var_;
};

template<int N>
std::vector<float>
GetDataAsType(const std::vector<CompressionValue<N>>& cr_values)
{
    cmc_assert(N == 4);
    std::vector<float> data;
    data.reserve(cr_values.size());

    for (auto iter = cr_values.begin(); iter != cr_values.end(); ++iter)
    {
        //const T current_value = *iter;
        //std::array<uint8_t, sizeof(T)> serialized_value;
        //std::memcpy(serialized_value.data(), &current_value, sizeof(T));
        //compression_data.emplace_back(std::move(serialized_value));
        float value;

        const std::array<uint8_t, N>& pref_mem = iter->GetMemoryForReading();
        std::memcpy(&value, pref_mem.data(), N);
        /* The storage is big endian but, on this machine we have little endian */
        //#if 1
        //std::array<uint8_t, 4> serialized_value{(*iter)[3],(*iter)[2],(*iter)[1],(*iter)[0]};
        //std::memcpy(&value, serialized_value.data(), 4);
        //#else
        //std::array<uint8_t, 4> serialized_value{(*iter)[0],(*iter)[1],(*iter)[2],(*iter)[3]};
        //std::memcpy(&value, serialized_value.data(), 4);
        //#endif
        data.push_back(value);
    }
    
    return data;
}

inline
const std::string&
ByteVar::GetName() const
{
    return std::visit([](auto&& var) -> const std::string& {
        return var.GetName();
    }, var_);
}

inline
DataLayout
ByteVar::GetInitialDataLayout() const 
{
    return std::visit([](auto&& var) -> DataLayout {
        return var.GetInitialDataLayout();
    }, var_);
}

inline
const GeoDomain&
ByteVar::GetGlobalDomain() const
{
    return std::visit([](auto&& var) -> const GeoDomain& {
        return var.GetGlobalDomain();
    }, var_);
}

inline
int
ByteVar::GetGlobalContextInformation() const
{
    return std::visit([](auto&& var) -> int {
        return var.GetGlobalContextInformation();
    }, var_);
}

inline
DataLayout
ByteVar::GetPreCompressionLayout() const 
{
    return std::visit([](auto&& var) -> DataLayout {
        return var.GetPreCompressionLayout();
    }, var_);
};

inline
void
ByteVar::StoreInitialMesh()
{
    std::visit([](auto&& var){
        var.StoreInitialMesh();
    }, var_);
}

inline
void
ByteVar::PerformTailTruncation()
{
    std::visit([](auto&& var){
        var.PerformTailTruncation();
    }, var_);
}

inline
bool
ByteVar::HasAtLeastOnePrefixBeenFoundInTheLastIteration() const
{
    return std::visit([](auto&& var){
        return var.HasAtLeastOnePrefixBeenFoundInTheLastIteration();
    }, var_);
}

inline
int
ByteVar::GetNumberOfSetPrefixIndicationBits() const
{
    return std::visit([](auto&& var){
        return var.GetNumberOfSetPrefixIndicationBits();
    }, var_);
}

inline
void
ByteVar::IndicateNoPrefixFound()
{
    std::visit([](auto&& var){
        var.IndicateNoPrefixFound();
    }, var_);
}

inline
void
ByteVar::LeaveCoarsePrefixUnchanged(const int elem_index)
{
    std::visit([&](auto&& var){
        var.LeaveCoarsePrefixUnchanged(elem_index);
    }, var_);
}

inline
void
ByteVar::LeaveInitialValueUnchanged(const int elem_index)
{
    std::visit([&](auto&& var){
        var.LeaveInitialValueUnchanged(elem_index);
    }, var_);
}

template<typename T>
void
ByteVar::IndicatePrefixFound(const CompressionValue<sizeof(T)>& prefix)
{
    cmc_assert(std::holds_alternative<ByteVariable<T>>(var_));
    std::visit([&](auto&& var){
        var.IndicatePrefixFound(prefix);
    }, var_);
}
template<typename T>
void
ByteVar::IndicatePrefixFound(CompressionValue<sizeof(T)>&& prefix)
{
    cmc_assert(std::holds_alternative<ByteVariable<T>>(var_));
    std::visit([&](auto&& var){
        var.IndicatePrefixFound(std::move(prefix));
    }, var_);
}

inline
void
ByteVar::PopLevelwisePrefixEncoding()
{
    std::visit([](auto&& var){
        var.PopLevelwisePrefixEncoding();
    }, var_);
}

inline
bool
ByteVar::EvaluateCommonPrefixFromInitialData(const int elem_start_index, const int num_elements)
{
    return std::visit([&](auto&& var){
        return var.EvaluateCommonPrefixFromInitialData(elem_start_index, num_elements);
    }, var_);
}

inline
bool
ByteVar::EvaluateCommonPrefixFromPreviousPrefix(const int elem_start_index, const int num_elements)
{
    return std::visit([&](auto&& var){
        return var.EvaluateCommonPrefixFromPreviousPrefix(elem_start_index, num_elements);
    }, var_);
}

inline
void
ByteVar::EvaluateCommonPrefixFromInitialDataEGU(const int elem_start_index, const int num_elements)
{
    std::visit([&](auto&& var){
        var.EvaluateCommonPrefixFromInitialDataEGU(elem_start_index, num_elements);
    }, var_);
}

inline
void
ByteVar::EvaluateCommonPrefixFromPreviousPrefixEGU(const int elem_start_index, const int num_elements)
{
    std::visit([&](auto&& var){
        var.EvaluateCommonPrefixFromPreviousPrefixEGU(elem_start_index, num_elements);
    }, var_);
}

inline
void
ByteVar::LeaveCoarsePrefixUnchangedEGU(const int elem_index)
{
    std::visit([&](auto&& var){
        var.LeaveCoarsePrefixUnchangedEGU(elem_index);
    }, var_);
}

inline
t8_forest_t
ByteVar::GetMesh() const
{
    return std::visit([](auto&& var) -> t8_forest_t {
        return var.GetMesh();
    }, var_);
}

inline
void
ByteVar::SetMesh(t8_forest_t forest)
{
    std::visit([&](auto&& var){
        var.SetMesh(forest);
    }, var_);
}

inline
AmrMesh&
ByteVar::GetAmrMesh()
{
    return std::visit([](auto&& var) -> AmrMesh& {
        return var.GetAmrMesh();
    }, var_);
}

inline
const AmrMesh&
ByteVar::GetAmrMesh() const
{
    return std::visit([](auto&& var) -> const AmrMesh& {
        return var.GetAmrMesh();
    }, var_);
}

inline
void
ByteVar::InitializeCompressionIteration()
{
    std::visit([](auto&& var){
        var.InitializeCompressionIteration();
    }, var_);
}

inline
void
ByteVar::FinalizeCompressionIterationEGU()
{
    std::visit([](auto&& var){
        var.FinalizeCompressionIterationEGU();
    }, var_);
}

inline
std::tuple<int, int, int, int, int, int, int, std::vector<uint8_t>>
//std::vector<uint8_t>
ByteVar::GatherSerializedCompressionData() const
{
    return std::visit([](auto&& var) -> std::tuple<int, int, int, int, int, int, int, std::vector<uint8_t>> {
        return var.GatherSerializedCompressionData();
    }, var_);
}

inline
std::tuple<int, int, int, int, int, int, std::vector<uint8_t>>
ByteVar::GatherSerializedCompressionDataEGU() const
{
    return std::visit([](auto&& var) -> std::tuple<int, int, int, int, int, int, std::vector<uint8_t>> {
        return var.GatherSerializedCompressionDataEGU();
    }, var_);
}

inline
void
ByteVar::SetMissingValueInNCFile(const int ncid, const int var_id, const int nc_type) const
{
    #ifdef CMC_WITH_NETCDF
    std::visit([&](auto&& var) {
        var.SetMissingValueInNCFile(ncid, var_id, nc_type);
    }, var_);
    #endif
}

template<typename T>
t8_forest_t
ByteVariable<T>::GetMesh() const
{
    return mesh_.GetMesh();
}

template<typename T>
void
ByteVariable<T>::SetMesh(t8_forest_t forest)
{
    mesh_.SetMesh(forest);
}

template<typename T>
int
ByteVariable<T>::GetNumberOfSetPrefixIndicationBits() const
{
    int num_set_bits = 0;
    for (auto pf_iter = prefixes_.begin(); pf_iter != prefixes_.end(); ++pf_iter)
    {
    for (auto iter = pf_iter->prefix_indicator_.begin(); iter != pf_iter->prefix_indicator_.end(); ++iter)
    {
        for (int bit_index = 0; bit_index < CHAR_BIT; ++bit_index)
        {
            if ((*iter) & (0b00000001 << bit_index))
            {
                ++num_set_bits;
            }
        }
    }
    }
    //cmc_debug_msg("Num Set Bits: ", num_set_bits);
    return num_set_bits;
}

template<typename T>
CompressionValue<sizeof(T)>
ByteVariable<T>::GetMaximumTailToggledValue(const int index, const std::vector<PermittedError>& permitted_errors, const CompressionValue<sizeof(T)>& initial_serialized_value, const T& missing_value) const
{
    bool is_toogling_progressing = true;
    CompressionValue<sizeof(T)> toggled_value = initial_serialized_value;

    int iteration_count = 0;
    const int max_iteration_count = sizeof(T) * CHAR_BIT;

    while (is_toogling_progressing && iteration_count < max_iteration_count)
    {
        const CompressionValue<sizeof(T)> save_previous_value = toggled_value;

        /* Toggle all ones up until the next unset bit (inclusive) */
        toggled_value.ToggleTailUntilNextUnsetBit();
        const T reinterpreted_value = toggled_value.template ReinterpretDataAs<T>();

        /* Check if it is error compliant */
        const ErrorCompliance error_evaluation = utilities_.IsValueErrorCompliant(permitted_errors, initial_data_[index], reinterpreted_value, missing_value);

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
CompressionValue<sizeof(T)>
ByteVariable<T>::GetMaximumTailClearedValue(const int index, const std::vector<PermittedError>& permitted_errors, const CompressionValue<sizeof(T)>& initial_serialized_value, const T& missing_value) const
{
    bool is_clearing_progressing = true;
    CompressionValue<sizeof(T)> cleared_value = initial_serialized_value;

    int iteration_count = 0;
    const int max_iteration_count = sizeof(T) * CHAR_BIT;

    while (is_clearing_progressing && iteration_count < max_iteration_count)
    {
        const CompressionValue<sizeof(T)> save_previous_value = cleared_value;

        /* Toggle all ones up until the next unset bit (inclusive) */
        cleared_value.ClearNextSetBitFromTail();
        const T reinterpreted_value = cleared_value.template ReinterpretDataAs<T>();

        /* Check if it is error compliant */
        const ErrorCompliance error_evaluation = utilities_.IsValueErrorCompliant(permitted_errors, initial_data_[index], reinterpreted_value, missing_value);

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
ByteVariable<T>::PerformTailTruncation()
{
    const T missing_value = attributes_.GetMissingValue();

    /* Iterate through the serialized values and try to emplace as many zeros at the tail as possible (compliant to the error threshold) */
    int index = 0;
    for (auto val_iter = byte_values_.begin(); val_iter != byte_values_.end(); ++val_iter, ++index)
    {
        if (!ApproxCompare(initial_data_[index], missing_value))
        {
            /* Get the permitted error for the current values */
            //TODO: Revert! Only for now with absolute errors
            const std::vector<PermittedError> permitted_errors = GetRemainingMaxAllowedAbsoluteError(index);
            //const std::vector<PermittedError> permitted_errors = GetPermittedError(index);

            /* Get the value which has been transformed by toggling as many ones from the back while setting the succeeding 'zero' bits to one */
            const CompressionValue<sizeof(T)> toggled_value = GetMaximumTailToggledValue(index, permitted_errors, *val_iter, missing_value);

            /* Get the value which has been transformed by clearing as many of the last set bits as possible */
            const CompressionValue<sizeof(T)> cleared_value = GetMaximumTailClearedValue(index, permitted_errors, *val_iter, missing_value);

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
            val_iter->UpdateTrailBitCount();
            //cmc_debug_msg("Trailing Zeros: Toggled: ", num_toogled_trailing_zeros, ", Cleared: ", num_cleared_trailing_zeros);
        } else
        {
            /* In order to not chenge missing values, we are just able to trim their trailing zeros */
            (*val_iter).UpdateTrailBitCount();
        }
    }
}

template<typename T>
std::vector<PermittedError>
ByteVariable<T>::GetPermittedError(const int index) const
{
    const t8_eclass_t eclass = t8_forest_get_eclass(mesh_.GetMesh(), 0);

    const t8_element_t* element = t8_forest_get_element_in_tree(mesh_.GetMesh(), 0, index);
    
    t8_eclass_scheme_c* ts =  t8_forest_get_eclass_scheme (mesh_.GetMesh(), eclass);

    /* Iterate over all error domains for all elements and find the minimum error */
    bool is_abs_error_present{false};
    bool is_rel_error_present{false};
    double abs_err{std::numeric_limits<double>::max()};
    double rel_err{std::numeric_limits<double>::max()};
    
    /* The first error domain is always the general criterion on the whole domain */
    auto ed_iter = utilities_.GetErrorDomainsBegin();
    if (IsAnyElementWithinGlobalDomain(1, &element, ts, ed_iter->GetDomain(), mesh_.GetInitialRefinementLevel(), attributes_.GetInitialDataLayout()))
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
        if (IsAnyElementWithinGeoDomain(1, &element, ts, ed_iter->GetDomain(), mesh_.GetInitialRefinementLevel(), attributes_.GetInitialDataLayout()))
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
std::vector<PermittedError>
ByteVariable<T>::GetRemainingMaxAllowedAbsoluteError(const int index) const
{
    std::vector<PermittedError> permitted_errors = GetPermittedError(index);
    cmc_assert(permitted_errors.size() == 1 && permitted_errors[0].criterion == CompressionCriterion::AbsoluteErrorThreshold);
    const double current_abs_deviation = utilities_.GetPreviousDeviation(index);

    if (current_abs_deviation > 0.0)
    {
        return std::vector<PermittedError>{PermittedError{CompressionCriterion::AbsoluteErrorThreshold, (permitted_errors[0].error - current_abs_deviation >= 0.0 ? permitted_errors[0].error - current_abs_deviation : 0.0)}};
    } else 
    {
        return permitted_errors;
    }
}

template<typename T>
void
ByteVariable<T>::InitializeCompressionIteration()
{
    prefixes_.emplace_back(mesh_.GetNumberLocalElements());
}

template<typename T>
bool
ByteVariable<T>::HasAtLeastOnePrefixBeenFoundInTheLastIteration() const
{
    for (auto pfi_iter = prefixes_.back().prefix_indicator_.begin(); pfi_iter != prefixes_.back().prefix_indicator_.end(); ++pfi_iter)
    {
        if (*pfi_iter != (uint8_t)0U)
            return true;
    }
    return false;
}

template<typename T>
void
ByteVariable<T>::PopLevelwisePrefixEncoding()
{
    prefixes_.pop_back();
}

template<typename T>
void
ByteVariable<T>::LeaveCoarsePrefixUnchangedEGU(const int elem_index)
{
    if (prefixes_.size() >= 2)
    {
        std::vector<uint8_t>& prev_prefix_indicator = (*std::prev(prefixes_.end(), 2)).prefix_indicator_;
        const int containing_byte_id = static_cast<int>(elem_index / CHAR_BIT);
        const int index_to_clear = elem_index % CHAR_BIT;
        prev_prefix_indicator[containing_byte_id] &= GetBitMaskToClearBit(index_to_clear);

        prefixes_.back().SetPrefixIndicatorBit(true);
        //prefixes_.back().SetPrefixIndicatorBit(false);
        prefixes_.back().SetPrefix((*std::prev(prefixes_.end(), 2)).prefixes_[elem_index]);
        previous_prefixes_new_.push_back((*std::prev(prefixes_.end(), 2)).prefixes_[elem_index]);
        (*std::prev(prefixes_.end(), 2)).prefixes_[elem_index] = CompressionValue<sizeof(T)>();
    } else
    {
        prefixes_.back().SetPrefixIndicatorBit(true);
        //prefixes_.back().SetPrefixIndicatorBit(false);
        prefixes_.back().SetPrefix(byte_values_[elem_index]);
        previous_prefixes_new_.push_back(byte_values_[elem_index]);
        byte_values_[elem_index] = CompressionValue<sizeof(T)>();
    }
}

template<typename T>
void
ByteVariable<T>::IndicateNoPrefixFoundEGU()
{
    prefixes_.back().SetPrefixIndicatorBit(false);
    prefixes_.back().SetPrefix(CompressionValue<sizeof(T)>());
    previous_prefixes_new_.push_back(CompressionValue<sizeof(T)>());
}

template<typename T>
void
ByteVariable<T>::IndicatePrefixFoundEGU(const CompressionValue<sizeof(T)>& prefix)
{
    prefixes_.back().SetPrefixIndicatorBit(true);
    prefixes_.back().SetPrefix(prefix);
    previous_prefixes_new_.push_back(prefix);
}

template<typename T>
void
ByteVariable<T>::FinalizeCompressionIterationEGU()
{
    previous_prefixes_ = previous_prefixes_new_;
    previous_prefixes_new_.clear();
}

template<typename T>
void
ByteVariable<T>::IndicateNoPrefixFound()
{
    prefixes_.back().SetPrefixIndicatorBit(false);
    prefixes_.back().SetPrefix(CompressionValue<sizeof(T)>());//Prefix<sizeof(T)>(), 0, sizeof(T) * CHAR_BIT
    //prefixes_.back().SetPrefix(CompressionValue<sizeof(T)>(Prefix<sizeof(T)>(), kNoPrefixIndicationBit, kNoPrefixIndicationBit));
}

template<typename T>
void
ByteVariable<T>::IndicatePrefixFound(const CompressionValue<sizeof(T)>& prefix)
{
    prefixes_.back().SetPrefixIndicatorBit(true);
    prefixes_.back().SetPrefix(prefix);
}

template<typename T>
void
ByteVariable<T>::IndicatePrefixFound(CompressionValue<sizeof(T)>&& prefix)
{
    prefixes_.back().SetPrefixIndicatorBit(true);
    prefixes_.back().SetPrefix(std::move(prefix));
}

template<typename T>
void
ByteVariable<T>::LeaveCoarsePrefixUnchanged(const int elem_index)
{
    IndicatePrefixFound(*((*std::prev(prefixes_.end(), 2)).prefixes_.data() + elem_index));

    (*std::prev(prefixes_.end(), 2)).prefixes_[elem_index].SetFrontBit(sizeof(T) * CHAR_BIT - (*std::prev(prefixes_.end(), 2)).prefixes_[elem_index].GetTrailBit());

    std::vector<uint8_t>& prev_prefix_indicator = (*std::prev(prefixes_.end(), 2)).prefix_indicator_;
    const int containing_byte_id = static_cast<int>(elem_index / CHAR_BIT);
    const int index_to_clear = elem_index % CHAR_BIT;
    prev_prefix_indicator[containing_byte_id] &= GetBitMaskToClearBit(index_to_clear);
}

template<typename T>
void
ByteVariable<T>::LeaveInitialValueUnchanged(const int elem_index)
{
    prefixes_.back().SetPrefixIndicatorBit(false);
    //(prefixes_.size() >= 2 ? (*std::prev(prefixes_.end(), 2)).prefixes_[elem_index] : byte_values_[elem_index])
    //prefixes_.back().SetPrefix((prefixes_.size() >= 2 ? (*std::prev(prefixes_.end(), 2)).prefixes_[elem_index] : byte_values_[elem_index]));
    if (prefixes_.size() >= 2)
    {
        prefixes_.back().SetPrefix((*std::prev(prefixes_.end(), 2)).prefixes_[elem_index]);
        (*std::prev(prefixes_.end(), 2)).prefixes_[elem_index] = CompressionValue<sizeof(T)>();
    } else
    {
        prefixes_.back().SetPrefix(byte_values_[elem_index]);
        byte_values_[elem_index] = CompressionValue<sizeof(T)>();
    }
    //IndicatePrefixFound(*(byte_values_.data() + elem_index));
    //byte_values_[elem_index].SetFrontBit(sizeof(T) * CHAR_BIT - (*(byte_values_.data() + elem_index)).GetTrailBit());
}

template<typename T>
bool
ByteVariable<T>::EvaluateCommonPrefixFromInitialData(const int elem_start_index, const int num_elements)
{   
    cmc_assert(num_elements >= 2);

    /* Read only byte values */
    VectorView<CompressionValue<sizeof(T)>> byte_values(byte_values_.data() + elem_start_index, num_elements);
    /* If any byte value is empty, we cannot find a prefix */
    for (auto cv_iter = byte_values.begin(); cv_iter != byte_values.end(); ++cv_iter)
    {
        if (cv_iter->IsEmpty())
        {
            /* There cannot be a common prefix */
            IndicateNoPrefixFound();
            return false;
        }
    }

    /* Determine a common prefix of the first two values */
    CompressionValue<sizeof(T)> prefix = GetCommonPrefix(byte_values[0], byte_values[1]);

    if (prefix.GetCountOfSignificantBits() == 0)
    {
        /* There is no common prefix */
        IndicateNoPrefixFound();
        return false;
    }

    CompressionValue<sizeof(T)> prefix_temp = prefix;

    for (int index = 2; index < num_elements; ++index)
    {
        /* Find a common prefix of all variables */
        prefix = GetCommonPrefix(prefix_temp, byte_values[index]);

        //cmc_debug_msg("Index: ", index, ", Trail bit of prefix: ", prefix.GetTrailBit(), ", Length Commpn Prefix: ", prefix.GetCountOfSignificantBits());
        if (prefix.GetCountOfSignificantBits() == 0)
        {
            /* There is no common prefix */
            IndicateNoPrefixFound();
            return false;
        }
        prefix_temp = prefix;
    }
    
    /* Set the prefix as extracted */
    for (int index = 0; index < num_elements; ++index)
    {
        byte_values_[elem_start_index + index].SetFrontBit(sizeof(T) * CHAR_BIT - prefix.GetTrailBit());
    }
    
    /* If we reach this line, a common prefix has been found */
    IndicatePrefixFound(std::move(prefix));
    return true;
}

template<typename T>
bool
ByteVariable<T>::EvaluateCommonPrefixFromPreviousPrefix(const int elem_start_index, const int num_elements)
{
    cmc_assert(num_elements >= 2);

    VectorView<CompressionValue<sizeof(T)>> byte_values((*std::prev(prefixes_.end(), 2)).prefixes_.data() + elem_start_index, num_elements);

    /* Check if all elements are corresponding to a prefix */
    for (auto cv_iter = byte_values.begin(); cv_iter != byte_values.end(); ++cv_iter)
    {
        if (cv_iter->IsEmpty())
        {
            /* There cannot be a common prefix */
            IndicateNoPrefixFound();
            return false;
        }
    }

    /* Determine a common prefix of the first two values */
    CompressionValue<sizeof(T)> prefix = GetCommonPrefix(byte_values[0], byte_values[1]);

    if (prefix.GetCountOfSignificantBits() == 0)
    {
        /* There is no common prefix */
        IndicateNoPrefixFound();
        return false;
    }

    for (int index = 2; index < num_elements; ++index)
    {
        /* Find a common prefix of all variables */
        prefix = GetCommonPrefix(prefix, byte_values[index]);

        if (prefix.GetCountOfSignificantBits() == 0)
        {
            /* There is no common prefix */
            IndicateNoPrefixFound();
            return false;
        }

    }

    //cmc_debug_msg("lelemet_id: ", elem_start_index, " num_elements: ", num_elements, ", Length Commpn Prefix: ", prefix.GetCountOfSignificantBits());
    /* Set the prefix as extracted */
    std::vector<CompressionValue<sizeof(T)>>& prev_prefixes = (*std::prev(prefixes_.end(), 2)).prefixes_;

    for (int index = 0; index < num_elements; ++index)
    {
        //std::cout << "trail: " << prev_prefixes[elem_start_index + index].GetTrailBit() << "Calculated front bit: " << sizeof(T) * CHAR_BIT - prev_prefixes[elem_start_index + index].GetTrailBit() << std::endl;
        
        prev_prefixes[elem_start_index + index].SetFrontBit(sizeof(T) * CHAR_BIT - prefix.GetTrailBit());

        /* In case we empty the current prefix, we need to unset it's prefix bit indicator, since we have forwarded the prefix to a higher (coarser) level */
        if (prev_prefixes[elem_start_index + index].IsEmpty())
        {
            //TODO: Check if correct
            std::vector<uint8_t>& prev_prefix_indicator = (*std::prev(prefixes_.end(), 2)).prefix_indicator_;
            const int containing_byte_id = static_cast<int>((elem_start_index + index) / CHAR_BIT);
            const int index_to_clear = (elem_start_index + index) % CHAR_BIT;
            prev_prefix_indicator[containing_byte_id] &= GetBitMaskToClearBit(index_to_clear);
            //cmc_debug_msg("prefix unset");
        }

    }

    /* If we reach this line, a common prefix has been found */
    IndicatePrefixFound(std::move(prefix));

    return true;
}



template<typename T>
void
ByteVariable<T>::EvaluateCommonPrefixFromInitialDataEGU(const int elem_start_index, const int num_elements)
{   
    cmc_assert(num_elements >= 2);

    /* Read only byte values */
    VectorView<CompressionValue<sizeof(T)>> byte_values(byte_values_.data() + elem_start_index, num_elements);

    /* If any byte value is empty, we cannot find a prefix */
    for (auto cv_iter = byte_values.begin(); cv_iter != byte_values.end(); ++cv_iter)
    {
        if (cv_iter->IsEmpty())
        {
            /* There cannot be a common prefix */
            IndicateNoPrefixFoundEGU();
            return;
        }
    }

    /* Determine a common prefix of the first two values */
    CompressionValue<sizeof(T)> prefix = GetCommonPrefix(byte_values[0], byte_values[1]);

    if (prefix.GetCountOfSignificantBits() == 0)
    {
        /* There is no common prefix */
        IndicateNoPrefixFoundEGU();
        return;
    }

    CompressionValue<sizeof(T)> prefix_temp = prefix;

    for (int index = 2; index < num_elements; ++index)
    {
        /* Find a common prefix of all variables */
        prefix = GetCommonPrefix(prefix_temp, byte_values[index]);

        //cmc_debug_msg("Index: ", index, ", Trail bit of prefix: ", prefix.GetTrailBit(), ", Length Commpn Prefix: ", prefix.GetCountOfSignificantBits());
        if (prefix.GetCountOfSignificantBits() == 0)
        {
            /* There is no common prefix */
            IndicateNoPrefixFoundEGU();
            return;
        }
        prefix_temp = prefix;
    }
    
    /* Set the prefix as extracted and update the front bit position of the suffixes */
    for (int index = 0; index < num_elements; ++index)
    {
        byte_values_[elem_start_index + index].SetFrontBit(sizeof(T) * CHAR_BIT - prefix.GetTrailBit());
    }
    
    /* If we reach this line, a common prefix has been found */
    //IndicatePrefixFound(std::move(prefix));
    IndicatePrefixFoundEGU(prefix);
}


template<typename T>
void
ByteVariable<T>::EvaluateCommonPrefixFromPreviousPrefixEGU(const int elem_start_index, const int num_elements)
{
    cmc_assert(num_elements >= 2);
    //cmc_debug_msg("Elem start index is: ", elem_start_index);
    //VectorView<CompressionValue<sizeof(T)>> byte_values((*std::prev(prefixes_.end(), 2)).prefixes_.data() + elem_start_index, num_elements);
    VectorView<CompressionValue<sizeof(T)>> byte_values(previous_prefixes_.data() + elem_start_index, num_elements);

    /* Check if all elements are corresponding to a prefix */
    for (auto cv_iter = byte_values.begin(); cv_iter != byte_values.end(); ++cv_iter)
    {
        if (cv_iter->IsEmpty())
        {
            /* There cannot be a common prefix */
            IndicateNoPrefixFoundEGU();
            return;
        }
    }

    /* Determine a common prefix of the first two values */
    CompressionValue<sizeof(T)> prefix = GetCommonPrefix(byte_values[0], byte_values[1]);

    if (prefix.GetCountOfSignificantBits() == 0)
    {
        /* There is no common prefix */
        IndicateNoPrefixFoundEGU();
        return;
    }

    CompressionValue<sizeof(T)> prefix_temp = prefix;

    for (int index = 2; index < num_elements; ++index)
    {
        /* Find a common prefix of all variables */
        prefix = GetCommonPrefix(prefix_temp, byte_values[index]);

        if (prefix.GetCountOfSignificantBits() == 0)
        {
            /* There is no common prefix */
            IndicateNoPrefixFoundEGU();
            return;
        }
        prefix_temp = prefix;
    }

    /* Get the previous prefixes */
    std::vector<CompressionValue<sizeof(T)>>& prev_prefixes = (*std::prev(prefixes_.end(), 2)).prefixes_;

    for (int index = 0; index < num_elements; ++index)
    {
        if (sizeof(T) * CHAR_BIT - prefix.GetTrailBit() + prev_prefixes[elem_start_index + index].GetTrailBit() > 32)
        {
            cmc_debug_msg("Hier ist prefix.GetTrailBit() = ", prefix.GetTrailBit(), " und prev_prefixes[elem_start_index + index].GetTrailBit() = ", prev_prefixes[elem_start_index + index].GetTrailBit());
        }
        prev_prefixes[elem_start_index + index].SetFrontBit(sizeof(T) * CHAR_BIT - prefix.GetTrailBit());

        /* In case we empty the current prefix, we need to unset it's prefix bit indicator, since we have forwarded the prefix to a higher (coarser) level */
        if (prev_prefixes[elem_start_index + index].IsEmpty())
        {
            std::vector<uint8_t>& prev_prefix_indicator = (*std::prev(prefixes_.end(), 2)).prefix_indicator_;
            const int containing_byte_id = static_cast<int>((elem_start_index + index) / CHAR_BIT);
            const int index_to_clear = (elem_start_index + index) % CHAR_BIT;
            prev_prefix_indicator[containing_byte_id] &= GetBitMaskToClearBit(index_to_clear);
            /* Explicitly overwirte the prefix (should not make a difference) */
            (*std::prev(prefixes_.end(), 2)).prefixes_[elem_start_index + index] = CompressionValue<sizeof(T)>();
        }
    }

    /* If we reach this line, a common prefix has been found */
    //IndicatePrefixFound(std::move(prefix));
    IndicatePrefixFoundEGU(prefix);
}


int
DetermineForestRefinementBits(std::vector<uint8_t>& serialized_variable, t8_forest_t forest);

template<typename T>
std::tuple<int, int, int, int, int, int, std::vector<uint8_t>>
ByteVariable<T>::GatherSerializedCompressionDataEGU() const
{
    std::vector<uint8_t> serialized_compressed_variable;

    /* We are loosely allocating the the memory for the serialized variable */
    serialized_compressed_variable.reserve(sizeof(T) * initial_data_.size());

    int num_serialized_bytes = 0;

    cmc_debug_msg("\nSize of levelwise prefixes: ", prefixes_.size());

    int num_bytes_prefix_indication = 0;

    /* First copy all level-wise prefix indication bits */
    for (auto lw_pfi_iter = prefixes_.rbegin(); lw_pfi_iter != prefixes_.rend(); ++lw_pfi_iter)
    {
        num_bytes_prefix_indication += lw_pfi_iter->prefix_indicator_.size();
        /* Copy the prefix indication bits of the current level */
        std::copy_n(lw_pfi_iter->prefix_indicator_.begin(), lw_pfi_iter->prefix_indicator_.size(), std::back_inserter(serialized_compressed_variable));
    }

    num_serialized_bytes += num_bytes_prefix_indication;

    /* Copy the encoded prefixes */
    const auto [num_bytes_prefix_lengths, num_bytes_prefix_encodings] = EncodeAndAppendLevelwisePrefixesEGU(serialized_compressed_variable, prefixes_);

    num_serialized_bytes += num_bytes_prefix_lengths + num_bytes_prefix_encodings;

    /* Copy the remaining suffixes */
    const auto [num_bytes_run_length_suffix_indicator, num_bytes_suffix_lengths, num_bytes_suffix_encodings] = EncodeAndAppendSuffixesEGU(serialized_compressed_variable, byte_values_);

    num_serialized_bytes += num_bytes_run_length_suffix_indicator + num_bytes_suffix_lengths + num_bytes_suffix_encodings;

    cmc_debug_msg("The variable was encoded/serialized to ", num_serialized_bytes, " bytes");

    cmc_debug_msg("The serialized variable has size(): ", serialized_compressed_variable.size());

    return std::make_tuple(num_bytes_prefix_indication,
                           num_bytes_prefix_lengths, num_bytes_prefix_encodings,
                           num_bytes_run_length_suffix_indicator,
                           num_bytes_suffix_lengths, num_bytes_suffix_encodings,
                           std::move(serialized_compressed_variable));
}

template<typename T>
void
ByteVariable<T>::SetMissingValueInNCFile(const int ncid, const int var_id, const int nc_type) const
{
    #ifdef CMC_WITH_NETCDF
    attributes_.SetMissingValueInNCFile(ncid, var_id, nc_type);
    #endif
}

template<typename T>
void
ByteVariable<T>::StoreInitialMesh()
{
    initial_mesh_ = mesh_.GetMesh();
    t8_forest_ref(initial_mesh_);
    is_initial_mesh_stored_ = true;
}

template<typename T>
void
ByteVariable<T>::ReleaseInitialMesh()
{
    if (is_initial_mesh_stored_ && initial_mesh_ != nullptr)
    {
        t8_forest_unref(&initial_mesh_);
        is_initial_mesh_stored_ = false;
    }
}


}

#endif
