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
#include "netcdf/cmc_nc_io.hxx"
#include "netcdf/cmc_nc_writer.hxx"
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
    void FinalizeCompressionIteration();
    void PopLevelwisePrefixEncoding();

    void LeaveCoarsePrefixUnchanged(const int elem_index);
    std::pair<bool, CompressionValue<sizeof(T)>> EvaluateCommonPrefix(VectorView<CompressionValue<sizeof(T)>>& byte_values);
    void ExtractCommonPrefix(std::vector<CompressionValue<sizeof(T)>>& prefixes, const int elem_start_index, const int num_elements);
    void ExtractCommonPrefixFromInitialData(const int elem_start_index, const int num_elements);
    void ExtractCommonPrefixFromPreviousPrefixes(const int elem_start_index, const int num_elements);
    void IndicatePrefixFound(const CompressionValue<sizeof(T)>& prefix);
    void IndicateNoPrefixFound();
    NcVariable WriteCompressedData(const int id) const;

    t8_forest_t GetMesh() const;
    void SetMesh(t8_forest_t forest);
    AmrMesh& GetAmrMesh() {return mesh_;};
    const AmrMesh& GetAmrMesh() const {return mesh_;};
    DataLayout GetInitialDataLayout() const {return attributes_.GetInitialDataLayout();};
    const GeoDomain& GetGlobalDomain() const {return global_domain_;};
    int GetGlobalContextInformation() const {return attributes_.GetGlobalContextInformation();};
    DataLayout GetPreCompressionLayout() const {return attributes_.GetPreCompressionLayout();};
    void SetMissingValueInNCFile(const int ncid, const int var_id, const int nc_type) const;
    
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
    std::vector<LevelwisePrefixData<T>> prefixes_;

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
    void PopLevelwisePrefixEncoding();

    void ExtractCommonPrefixFromInitialData(const int elem_start_index, const int num_elements);
    void ExtractCommonPrefixFromPreviousPrefixes(const int elem_start_index, const int num_elements);
    void LeaveCoarsePrefixUnchanged(const int elem_index);
    void FinalizeCompressionIteration();
    NcVariable WriteCompressedData() const;

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
void
ByteVar::PopLevelwisePrefixEncoding()
{
    std::visit([](auto&& var){
        var.PopLevelwisePrefixEncoding();
    }, var_);
}

inline
void
ByteVar::ExtractCommonPrefixFromInitialData(const int elem_start_index, const int num_elements)
{
    std::visit([&](auto&& var){
        var.ExtractCommonPrefixFromInitialData(elem_start_index, num_elements);
    }, var_);
}

inline
void
ByteVar::ExtractCommonPrefixFromPreviousPrefixes(const int elem_start_index, const int num_elements)
{
    std::visit([&](auto&& var){
        var.ExtractCommonPrefixFromPreviousPrefixes(elem_start_index, num_elements);
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
ByteVar::FinalizeCompressionIteration()
{
    std::visit([](auto&& var){
        var.FinalizeCompressionIteration();
    }, var_);
}

inline NcVariable
ByteVar::WriteCompressedData() const
{
    return std::visit([this](auto&& var) -> NcVariable {
        return var.WriteCompressedData(this->GetID());
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
    /* Emplace a new LevelPrefixData-object which will hold the prefix data of the current iteration of prefix extraction */
    prefixes_.emplace_back(mesh_.GetNumberLocalElements());
}

template<typename T>
void
ByteVariable<T>::PopLevelwisePrefixEncoding()
{
    prefixes_.pop_back();
}

/* If a single element is investigated and not a family of elements during the prefix extraction,
 * we are not altering the value/previous-prefix and keep it as a prefx for the next level */
template<typename T>
void
ByteVariable<T>::LeaveCoarsePrefixUnchanged(const int elem_index)
{
    /* In case that the element stays the same (the sibling elements are not yet present in the mesh).
     * We indicate the whole data/prefix as a common prefix and will evaluate in a later adaptation step. */
    prefixes_.back().SetRefinementIndicatorBit(false);
    prefixes_.back().SetPrefixIndicatorBit(true);

    /* Copy the data/prefix to the next 'coarser level' */
    if (prefixes_.size() >= 2)
    {
        /* In any other prefix extraction iteration than the first, we just copy the previous prefix */
        LevelwisePrefixData<T>& previous_prefixes = (*std::prev(prefixes_.end(), 2));
        prefixes_.back().SetPrefix(previous_prefixes.prefixes[elem_index]);
        //previous_prefixes_new_.push_back(previous_prefixes.prefixes[elem_index]);
        previous_prefixes.prefixes[elem_index] = CompressionValue<sizeof(T)>();
    } else
    {
        /* In the first prefix extraction iteration, we just keep the initial byte value and consider it in a 'coarser' iteration */
        prefixes_.back().SetPrefix(byte_values_[elem_index]);
        //previous_prefixes_new_.push_back(byte_values_[elem_index]);
        /* Remove the 'prefix' from the initial data and replace it with an empty 'dummy value' */
        byte_values_[elem_index] = CompressionValue<sizeof(T)>();
    }
}

template<typename T>
void
ByteVariable<T>::IndicateNoPrefixFound()
{
    prefixes_.back().SetRefinementIndicatorBit(true);
    prefixes_.back().SetPrefixIndicatorBit(false);
    prefixes_.back().SetPrefix(CompressionValue<sizeof(T)>());
}

template<typename T>
void
ByteVariable<T>::IndicatePrefixFound(const CompressionValue<sizeof(T)>& prefix)
{
    prefixes_.back().SetRefinementIndicatorBit(true);
    prefixes_.back().SetPrefixIndicatorBit(true);
    prefixes_.back().SetPrefix(prefix);
}

template<typename T>
void
ByteVariable<T>::FinalizeCompressionIteration()
{
    //TODO:Repartition here the data/prefixes?
    return;
}

template<typename T>
void
ByteVariable<T>::ExtractCommonPrefixFromInitialData(const int elem_start_index, const int num_elements)
{
    /* Extract a common prefix from the initial data */
    ExtractCommonPrefix(byte_values_, elem_start_index, num_elements);   
}

template<typename T>
void
ByteVariable<T>::ExtractCommonPrefixFromPreviousPrefixes(const int elem_start_index, const int num_elements)
{
    /* Get the previous prefixes */
    std::vector<CompressionValue<sizeof(T)>>& prev_prefixes = (*std::prev(prefixes_.end(), 2)).prefixes;

    /* Extract a common prefix from the initial data */
    ExtractCommonPrefix(prev_prefixes, elem_start_index, num_elements);  
}

template<typename T>
void
ByteVariable<T>::ExtractCommonPrefix(std::vector<CompressionValue<sizeof(T)>>& prefixes, const int elem_start_index, const int num_elements)
{
    /* Create a read-only view of the considered prefixes */
    VectorView<CompressionValue<sizeof(T)>> prefix_view(prefixes.data() + elem_start_index, num_elements);

    /* Evaluate whether there is a common prefix to extract */
    auto [is_prefix_found, prefix] = EvaluateCommonPrefix(prefix_view);

    if (is_prefix_found)
    {
        IndicatePrefixFound(prefix);

        /* We need to trim the previous prefixes by the extracted common prefix */
        for (int index = 0; index < num_elements; ++index)
        {
            byte_values_[elem_start_index + index].SetFrontBit(sizeof(T) * CHAR_BIT - prefix.GetTrailBit());
        }
    } else
    {
        IndicateNoPrefixFound();
    }
}

template<typename T>
std::pair<bool, CompressionValue<sizeof(T)>>
ByteVariable<T>::EvaluateCommonPrefix(VectorView<CompressionValue<sizeof(T)>>& byte_values)
{
    cmc_assert(byte_values.size() >= 2);

    /* Check if all elements are holding an actual prefix */
    for (auto cv_iter = byte_values.begin(); cv_iter != byte_values.end(); ++cv_iter)
    {
        if (cv_iter->IsEmpty())
        {
            /* Since this prefix value is empty, there cannot be a common prefix */
            return std::make_pair(false, CompressionValue<sizeof(T)>());
        }
    }

    /* Determine a common prefix of the first two values */
    CompressionValue<sizeof(T)> prefix = GetCommonPrefix(byte_values[0], byte_values[1]);

    /* Check if there is common prefix between the first two values */
    if (prefix.GetCountOfSignificantBits() == 0)
    {
        /* There is no common prefix */
        return std::make_pair(false, CompressionValue<sizeof(T)>());
    }

    /* Check if there is a common prefix with the other values within the view */
    for (size_t index = 2; index < byte_values.size() ; ++index)
    {
        /* Find a common prefix of all variables */
        prefix = GetCommonPrefix(prefix, byte_values[index]);

        if (prefix.GetCountOfSignificantBits() == 0)
        {
            /* There is no common prefix */
            return std::make_pair(false, CompressionValue<sizeof(T)>());
        }
    }

    /* If the function arrives here, we do have found a common prefix which can be extracted from the 'previous prefixes' */
    return std::make_pair(true, prefix);
}

int
DetermineForestRefinementBits(std::vector<uint8_t>& serialized_variable, t8_forest_t forest);


template<typename T>
NcVariable
ByteVariable<T>::WriteCompressedData(const int id) const
{
    /* Declare a vector which will hold the level-wise encoded data */
    std::vector<std::vector<uint8_t>> buffered_data;
    buffered_data.reserve(prefixes_.size() + 1);


    size_t num_bytes = 0;
    /* Encode and retrieve the extracted prefixes */
    for (auto lw_prefix_iter = prefixes_.rbegin(); lw_prefix_iter != prefixes_.rend(); ++lw_prefix_iter)
    {
        buffered_data.emplace_back(lw_prefix_iter->EncodeLevelData());
        cmc_debug_msg("Num bytes for level: ", buffered_data.back().size());
        num_bytes += buffered_data.back().size();
    }

    /* Encode and append the leftover suffixes (on the finest level) that could not be truncated or extracted */
    buffered_data.emplace_back(EncodeSuffixes(byte_values_));
    cmc_debug_msg("Num bytes for suffixes: ", buffered_data.back().size());
    num_bytes += buffered_data.back().size();
    cmc_debug_msg("Overall bytes for this variable: ", num_bytes);


    //TODO: Make mechanism for parallel output
    //The complete global byte size needs to gathered and the level data has to be filled in parallel to it

    /* Create a netCDF variable to put out */
    NcSpecificVariable<uint8_t> compressed_variable{GetName(), id};
    compressed_variable.Reserve(num_bytes);

    /* Put the buffered data into the variable to put out */
    for (auto lvl_data_iter = buffered_data.begin(); lvl_data_iter != buffered_data.end(); ++lvl_data_iter)
    {
        compressed_variable.PushBack(*lvl_data_iter);
    }

    /* Assign some attributes to it */
    std::vector<NcAttribute> attributes;
    attributes.emplace_back("initial_refinement_level", mesh_.GetInitialRefinementLevel());
    attributes.emplace_back("missing_value", attributes_.GetMissingValue());
    attributes.emplace_back("initial_layout", static_cast<int>(attributes_.GetInitialDataLayout()));
    attributes.emplace_back("pre_compression_layout", static_cast<int>(attributes_.GetPreCompressionLayout()));
    attributes.emplace_back("global_context", attributes_.GetGlobalContextInformation());
    attributes.emplace_back("datatype", static_cast<int>(ConvertToCmcType<T>()));
    //TODO: Check if scaling and offset have been applied, if not we need to add those attributes

    /* Get the global domain of this variable */
    const GeoDomain& var_domain = GetGlobalDomain();

    /* Write the domain lengths as attributes */
    if (const int lon = var_domain.GetDimensionLength(Dimension::Lon);
        lon > 1)
    {
        attributes.emplace_back("lon", lon);
    }
    if (const int lat = var_domain.GetDimensionLength(Dimension::Lat);
        lat > 1)
    {
        attributes.emplace_back("lat", lat);
    }
    if (const int lev = var_domain.GetDimensionLength(Dimension::Lev);
        lev > 1)
    {
        attributes.emplace_back("lev", lev);
    }
    if (const int time = var_domain.GetDimensionLength(Dimension::Time);
        time > 1)
    {
        attributes.emplace_back("time", time);
    }


    return NcVariable(std::move(compressed_variable), std::move(attributes));
    #if 0
    const std::string file_name = "test_nc_output_file.nc";
    const int netcdf_format = NC_CDF5;

    /* Create a writer which outputs the data */
    NcWriter writer(file_name, netcdf_format);
    cmc_debug_msg("Writer has been created");
    writer.AddVariable(std::move(compressed_variable), std::move(attributes));
    cmc_debug_msg("Variable has been added");
    writer.Write();
    cmc_debug_msg("Write complete");
    #endif
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
