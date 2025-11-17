#ifndef CMC_T8_BYTE_VARIABLE_HXX
#define CMC_T8_BYTE_VARIABLE_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "t8code/cmc_t8_data_variables_forward.hxx"
#include "utilities/cmc_variable_transformer_forward.hxx"
#include "utilities/cmc_prefix.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_variables.hxx"
#include "utilities/cmc_vector_view.hxx"
#include "utilities/cmc_prefix_encoding.hxx"
#include <t8_forest/t8_forest_vtk.h>
#include "utilities/cmc_huffman.hxx"
#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_arithmetic_encoder.hxx"

#include <vector>
#include <cmath>
#include <iterator>
#include <unordered_map>
#include <bitset>

#ifdef CMC_WITH_NETCDF
#include "netcdf/cmc_netcdf.hxx"
#include "netcdf/cmc_nc_io.hxx"
#include "netcdf/cmc_nc_writer.hxx"
#endif

#define _PERFORM_TEST_COMPARISON 1

namespace cmc
{

template<typename T> class ByteVariable;
class ByteVar;

using CmcByteVariable = std::variant<ByteVariable<int8_t>, ByteVariable<char>, ByteVariable<int16_t>, ByteVariable<int32_t>, ByteVariable<float>, ByteVariable<double>,
                                     ByteVariable<uint8_t>, ByteVariable<uint16_t>, ByteVariable<uint32_t>, ByteVariable<int64_t>, ByteVariable<uint64_t>>;

#if _PERFORM_TEST_COMPARISON
struct _TestComparisonEncodings
{
    bit_vector::BitVector encoded_lzc;
    bit_vector::BitVector encoded_residuals;
};
#endif

template<typename T>
class ByteVariable
{
public:
    ByteVariable() = default;
    ByteVariable(std::vector<CompressionValue<sizeof(T)>>&& serialized_data)
    : byte_values_(std::move(serialized_data)) {
        cmc_debug_msg("Variable has obtained data of size: ", byte_values_.size());
    };
    ByteVariable(const std::string& name, const GeoDomain domain, const DataLayout layout, const DataLayout pre_compression_layout,
                 const int global_context_information, const T missing_value)
    : name_{name}, attributes_(missing_value, layout, pre_compression_layout, global_context_information), global_domain_{domain} {};

    ~ByteVariable() {
        ReleaseInitialMesh();
    }

    ByteVariable(const ByteVariable& other) = default;
    ByteVariable& operator=(const ByteVariable& other) = default;
    ByteVariable(ByteVariable&& other) = default;
    ByteVariable& operator=(ByteVariable&& other) = default;

    const std::string& GetName() const {return name_;}

    void PerformTailTruncation();
    void PerformTailTruncationOnInitialData();
    void PerformTailTruncationRegardingUncompressedStates();

    void InitializeCompressionIteration();
    void FinalizeCompressionIteration();
    void PopLevelwisePrefixEncoding();

    void InitializeDecompressionIteration();
    void FinalizeDecompressionIteration();
    void InitializeSuffixDecompression();
    void InitializeDiffDecompression();
    void SetDiffDecompressionRootElementValue(const CmcUniversalType& root_value);
    void WriteDataToVTK_(const int step_id);

    /** Compression Routines **/
    void LeaveCoarsePrefixUnchanged(const int elem_index);
    std::pair<bool, CompressionValue<sizeof(T)>> EvaluateCommonPrefix(VectorView<CompressionValue<sizeof(T)>>& byte_values);
    std::pair<bool, CompressionValue<sizeof(T)>> EvaluateCommonPrefix(std::vector<CompressionValue<sizeof(T)>>& prefixes, const std::vector<int>& element_ids);
    void ExtractCommonPrefix(std::vector<CompressionValue<sizeof(T)>>& prefixes, const int elem_start_index, const int num_elements);
    void ExtractCommonPrefix(std::vector<CompressionValue<sizeof(T)>>& prefixes, const std::vector<int>& element_ids);
    void ExtractCommonPrefixFromInitialData(const int elem_start_index, const int num_elements);
    void ExtractCommonPrefixFromPreviousPrefixes(const int elem_start_index, const int num_elements);
    void ExtractCommonPrefixFromInitialData(const std::vector<int>& element_ids);
    void ExtractCommonPrefixFromPreviousPrefixes(const std::vector<int>& element_ids);
    void IndicatePrefixFound(const CompressionValue<sizeof(T)>& prefix);
    void IndicateNoPrefixFound();
    cmc::nc::Variable WriteCompressedData(const int id, const int time_step, const SuffixEncoding encoding_scheme) const;
    cmc::nc::Variable WriteCompressedData(const int id, const int time_step, SuffixEncodingFunc<sizeof(T)> suffix_encoder) const;

    /** De-Compression Routines **/
    void AssignCompressionValueForDecompressionStart();
    void DecompressionRefine(const int elem_id);
    void DecompressionLeaveElementUnchanged(const int elem_id);
    void RefineAndApplyCommonPrefix(const int elem_id, std::vector<uint8_t> prefix, const int num_prefix_bits);
    void LeaveElementUnchangedAndApplyPrefix(const int elem_id, std::vector<uint8_t> prefix, const int num_prefix_bits);
    void PrintCompressionValues() const;
    size_t GetCountSignificantBits(const int elem_id) const;

    t8_forest_t GetMesh() const;
    void SetMesh(t8_forest_t forest);
    AmrMesh& GetAmrMesh() {return mesh_;};
    const AmrMesh& GetAmrMesh() const {return mesh_;};
    DataLayout GetInitialDataLayout() const {return attributes_.GetInitialDataLayout();};
    const GeoDomain& GetGlobalDomain() const {return global_domain_;};
    int GetGlobalContextInformation() const {return attributes_.GetGlobalContextInformation();};
    DataLayout GetPreCompressionLayout() const {return attributes_.GetPreCompressionLayout();};
    void SetMissingValueInNCFile(const int ncid, const int var_id, const int nc_type) const;
    
    void KeepInitialData(const bool keep_data);
    void StoreInitialMesh();
    void ReleaseInitialMesh(); 
    friend class TransformerCompressionToByteVariable;

    void XORConsecutiveValues();

    void ExtractMeanAndLeaveDifferencesFromInitialData(const int elem_start_index, const int num_elements);
    void ExtractMeanAndLeaveDifferencesFromPreviousMeans(const int elem_start_index, const int num_elements);

    void LeaveValueUnchangedForNextMeanComputation(const int element_index);
    cmc::nc::Variable WriteCompressedDiffData(const int id, const int time_step) const;
    void InitializeResidualAlphabet();

    void ExtractMeanFromInitialData(const int elem_start_index, const int num_elements);
    void ExtractMeanFromPreviousMeans(const int elem_start_index, const int num_elements);
    void LeaveValueUnchanged(const int elem_start_index);
    void TryFittingPyramid();
    void TryLZCEncoding();
    void ApplyResidualAndStoreElement(const int elem_id, const uint32_t encoded_lzc, const std::vector<uint8_t>& residual_bits);
    void StoreElementUnchanged(const int elem_id);

    void WriteDataToFile(const std::string& file_name) const;

    void _TestComparisonPCPLightAMRFromInitialData(const int elem_start_index, const int num_elements);
    void _TestComparisonPCPLightAMRFromPreviousMeans(const int elem_start_index, const int num_elements);
    void _TestComparisonPCPLightAMR(std::vector<CompressionValue<sizeof(T)>>& values, const int elem_start_index, const int num_elements);
    cmc::nc::Variable _WriteCompressedDiffDataPCPLightAMRComparison(const int id, const int time_step) const;    
    void ApplyResidualWithoutImplicitOneBitAndStoreElement(const int elem_id, const uint32_t encoded_lzc, const std::vector<uint8_t>& residual_bits);
private:
    void ExtractMeanAndLeaveDifferences(std::vector<CompressionValue<sizeof(T)>>& values, const int elem_start_index, const int num_elements);
    void ExtractMean(std::vector<CompressionValue<sizeof(T)>>& values, const int elem_start_index, const int num_elements);
    CompressionValue<sizeof(T)> GetMaximumTailClearedValue(const int, const std::vector<PermittedError>&, const CompressionValue<sizeof(T)>&, const T&) const;
    CompressionValue<sizeof(T)> GetMaximumTailToggledValue(const int, const std::vector<PermittedError>&, const CompressionValue<sizeof(T)>&, const T&) const;
    CompressionValue<sizeof(T)> GetMaximumTailClearedValue(const std::vector<PermittedError>& permitted_errors, const double previous_abs_deviation, const CompressionValue<sizeof(T)>& initial_serialized_value, const T& missing_value) const;
    CompressionValue<sizeof(T)> GetMaximumTailToggledValue(const std::vector<PermittedError>& permitted_errors, const double previous_abs_deviation, const CompressionValue<sizeof(T)>& initial_serialized_value, const T& missing_value) const;
    CompressionValue<sizeof(T)> GetMaximumTailClearedValueRegardingUncompressedStates(const std::vector<PermittedError>& permitted_errors, const std::vector<DomainIndex>& initial_data_indices, const CompressionValue<sizeof(T)>& initial_serialized_value, const T& missing_value) const;
    CompressionValue<sizeof(T)> GetMaximumTailToggledValueRegardingUncompressedStates(const std::vector<PermittedError>& permitted_errors, const std::vector<DomainIndex>& initial_data_indices, const CompressionValue<sizeof(T)>& initial_serialized_value, const T& missing_value) const;
    std::vector<PermittedError> GetPermittedError(const int index) const;
    std::vector<PermittedError> GetRemainingMaxAllowedAbsoluteError(const int index) const;
    std::vector<int> GetPrefixLengthFrequency() const;

    std::string name_; //!< The name of the variable
    std::vector<T> initial_data_; //!< The actual data of the variable
    AmrMesh mesh_; //!< The mesh on which the variable is defined
    VariableAttributes<T> attributes_;
    VariableUtilities<T> utilities_;
    GeoDomain global_domain_;

    std::vector<CompressionValue<sizeof(T)>> byte_values_;
    std::vector<LevelwisePrefixData<T>> prefixes_;

    std::vector<bit_map::BitMap> interpolation_indications_;

    //std::vector<arithmetic_encoding::Letter> alphabet_;
    std::unordered_map<uint32_t, uint32_t> alphabet_;
    std::vector<CompressionValue<sizeof(T)>> byte_values_new_;

    bool is_initial_data_kept_{true};
    bool is_initial_mesh_stored_{false};
    t8_forest_t initial_mesh_{nullptr};

    AmrMesh uncompressed_mesh_;
    std::vector<T> uncompressed_data_;
    bool are_uncompressed_states_stored_{false};

    #if _PERFORM_TEST_COMPARISON
    std::vector<_TestComparisonEncodings> test_encodings_;
    #endif
};

class ByteVar
{
public:
    ByteVar() = delete;
    ByteVar(const int id, const CmcType type, CmcByteVariable&& variable)
    : id_{id}, type_{type}, var_{std::move(variable)}{};
    ByteVar(const int id, const CmcType type, const std::string& name, const GeoDomain domain, const DataLayout layout, const DataLayout pre_compression_layout,
            const int global_context_information, const CmcUniversalType missing_value)
    : id_{id}, type_{type} {
        SetUpByteVariable(type, name, domain, layout, pre_compression_layout, global_context_information, missing_value);
    };
    ~ByteVar() = default;

    const std::string& GetName() const;
    int GetID() const {return id_;}
    CmcType GetType() const {return type_;}

    void InitializeCompressionIteration();
    void FinalizeCompressionIteration();
    void PopLevelwisePrefixEncoding();

    void InitializeDecompressionIteration();
    void FinalizeDecompressionIteration();
    void InitializeSuffixDecompression();

    void WriteDataToVTK_(const int step_id);

    /** Compression Routines **/
    void PerformTailTruncation();
    void PerformTailTruncationOnInitialData();
    void PerformTailTruncationRegardingUncompressedStates();
    void ExtractCommonPrefixFromInitialData(const int elem_start_index, const int num_elements);
    void ExtractCommonPrefixFromPreviousPrefixes(const int elem_start_index, const int num_elements);
    void ExtractCommonPrefixFromInitialData(const std::vector<int>& element_ids);
    void ExtractCommonPrefixFromPreviousPrefixes(const std::vector<int>& element_ids);
    void LeaveCoarsePrefixUnchanged(const int elem_index);
    cmc::nc::Variable WriteCompressedData(const int time_step, const SuffixEncoding encoding_scheme) const;
    
    /** De-Compression Routines **/
    void AssignCompressionValueForDecompressionStart();
    void DecompressionRefine(const int elem_id);
    void DecompressionLeaveElementUnchanged(const int elem_id);
    void RefineAndApplyCommonPrefix(const int elem_id, std::vector<uint8_t> prefix, const int num_prefix_bits);
    void LeaveElementUnchangedAndApplyPrefix(const int elem_id, std::vector<uint8_t> prefix, const int num_prefix_bits);
    void PrintCompressionValues() const;
    size_t GetCountSignificantBits(const int elem_id) const;

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
    void KeepInitialData(const bool keep_data);

    void XORConsecutiveValues();
    void ExtractMeanAndLeaveDifferencesFromInitialData(const int elem_start_index, const int num_elements);
    void ExtractMeanAndLeaveDifferencesFromPreviousMeans(const int elem_start_index, const int num_elements);
    void LeaveValueUnchangedForNextMeanComputation(const int element_index);
    cmc::nc::Variable WriteCompressedDiffData(const int time_step) const;

    void ExtractMeanFromInitialData(const int elem_start_index, const int num_elements);
    void ExtractMeanFromPreviousMeans(const int elem_start_index, const int num_elements);
    void LeaveValueUnchanged(const int elem_start_index);

    void InitializeResidualAlphabet();
    void TryLZCEncoding();

    void InitializeDiffDecompression();
    void SetDiffDecompressionRootElementValue(const CmcUniversalType& root_value);
    void ApplyResidualAndStoreElement(const int elem_id, const uint32_t encoded_lzc, const std::vector<uint8_t>& residual_bits);
    void StoreElementUnchanged(const int elem_id);

    void WriteDataToFile(const std::string& file_name) const;
    void _TestComparisonPCPLightAMRFromInitialData(const int elem_start_index, const int num_elements);
    void _TestComparisonPCPLightAMRFromPreviousMeans(const int elem_start_index, const int num_elements);
    cmc::nc::Variable _WriteCompressedDiffDataPCPLightAMRComparison(const int time_step) const;  
    void ApplyResidualWithoutImplicitOneBitAndStoreElement(const int elem_id, const uint32_t encoded_lzc, const std::vector<uint8_t>& residual_bits);
private:
    void SetUpByteVariable(const CmcType type, const std::string& name, const GeoDomain domain, const DataLayout layout, const DataLayout pre_compression_layout,
                      const int global_context_information, const CmcUniversalType missing_value);
    
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
        float value;

        const std::array<uint8_t, N>& pref_mem = iter->GetMemoryForReading();
        std::memcpy(&value, pref_mem.data(), N);

        data.push_back(value);
    }
    
    return data;
}

template <typename T>
void
ByteVariable<T>::PrintCompressionValues() const
{
    if constexpr (sizeof(T) == 4)
    {
        int32_t val{0};
        int i=0;
        for (auto iter = byte_values_.begin(); iter != byte_values_.end(); ++iter)
        {
            const std::array<uint8_t, 4>& pref_mem = iter->GetMemoryForReading();
            std::memcpy(&val, pref_mem.data(), 4);
            cmc_debug_msg(std::bitset<32>(val), " fuer Elem: ", i);
            ++i;
        }
    }
}

inline
void
ByteVar::XORConsecutiveValues()
{
    std::visit([](auto&& var) {
        var.XORConsecutiveValues();
    }, var_);
}

inline
void
ByteVar::PrintCompressionValues() const
{
    std::visit([](auto&& var){
        var.PrintCompressionValues();
    }, var_);
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
ByteVar::KeepInitialData(const bool keep_data)
{
    std::visit([&keep_data](auto&& var){
        var.KeepInitialData(keep_data);
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
ByteVar::PerformTailTruncationOnInitialData()
{
    std::visit([](auto&& var){
        var.PerformTailTruncationOnInitialData();
    }, var_);
}

inline
void
ByteVar::WriteDataToFile(const std::string& file_name) const
{
    std::visit([&](auto&& var){
        var.WriteDataToFile(file_name);
    }, var_);
}

inline
void
ByteVar::PerformTailTruncationRegardingUncompressedStates()
{
    std::visit([](auto&& var){
        var.PerformTailTruncationRegardingUncompressedStates();
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

inline
void
ByteVar::InitializeDecompressionIteration()
{
    std::visit([](auto&& var){
        var.InitializeDecompressionIteration();
    }, var_);
}

inline
void
ByteVar::FinalizeDecompressionIteration()
{
    std::visit([](auto&& var){
        var.FinalizeDecompressionIteration();
    }, var_);
}

inline
void
ByteVar::InitializeDiffDecompression()
{
    std::visit([](auto&& var){
        var.InitializeDiffDecompression();
    }, var_);
}

inline
void
ByteVar::SetDiffDecompressionRootElementValue(const CmcUniversalType& root_value)
{
    std::visit([&](auto&& var){
        var.SetDiffDecompressionRootElementValue(root_value);
    }, var_);
}


inline nc::Variable
ByteVar::WriteCompressedDiffData(const int time_step) const
{
    return std::visit([this, &time_step](auto&& var) -> nc::Variable {
        return var.WriteCompressedDiffData(this->GetID(), time_step);
    }, var_);
}

inline nc::Variable
ByteVar::_WriteCompressedDiffDataPCPLightAMRComparison(const int time_step) const
{
    #if _PERFORM_TEST_COMPARISON
    return std::visit([this, &time_step](auto&& var) -> nc::Variable {
        return var._WriteCompressedDiffDataPCPLightAMRComparison(this->GetID(), time_step);
    }, var_);
    #else
    return nc::Variable();
    #endif
}


inline nc::Variable
ByteVar::WriteCompressedData(const int time_step, const SuffixEncoding encoding_scheme) const
{
    return std::visit([this, &time_step, &encoding_scheme](auto&& var) -> nc::Variable {
        return var.WriteCompressedData(this->GetID(), time_step, encoding_scheme);
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

inline
void
ByteVar::WriteDataToVTK_(const int step_id)
{
    std::visit([&step_id](auto&& var) {
        var.WriteDataToVTK_(step_id);
    }, var_);
}

template<typename T>
void
ByteVariable<T>::XORConsecutiveValues()
{
    cmc_assert(byte_values_.size() >= 1);

    CompressionValue<sizeof(T)> last_val = byte_values_[0];
    CompressionValue<sizeof(T)> new_val;
    for (size_t index = 1; index < byte_values_.size(); ++index)
    {
        new_val = byte_values_[index];
        byte_values_[index] ^= last_val;
        last_val = new_val;
    }
    byte_values_[0] = CompressionValue<sizeof(T)>(T(0));

}

template<typename T>
std::vector<int>
ByteVariable<T>::GetPrefixLengthFrequency() const
{
    std::vector<int> frequency(sizeof(T) * bit_vector::kCharBit + 1, int{0});

    /* Iterate over all levels of extraction */
    for (auto lvl_iter = prefixes_.begin(); lvl_iter != prefixes_.end(); ++lvl_iter)
    {
        /* Iterate over all prefixes on this level */
        for (auto val_iter = lvl_iter->prefixes.begin(); val_iter != lvl_iter->prefixes.end(); ++val_iter)
        {
            /* Get the count of significant bits, i.e. the length of the prefix */
            const size_t pref_length = val_iter->GetCountOfSignificantBits();

            /* Increment the counter for the given prefix length */
            ++frequency[pref_length];
        }
    }

    int i = 0;
    for (auto fiter = frequency.begin(); fiter != frequency.end(); ++fiter, ++i)
    {
        cmc_debug_msg("Frequency of Length: ", i, " is ", *fiter);
    }
    return frequency;
}

template<typename T>
void
ByteVariable<T>::WriteDataToVTK_(const int step_id)
{
    const std::string file_name = GetName() + "_" + std::to_string(GetGlobalContextInformation()) + "_" + std::to_string(step_id);
    std::vector<float> float_data;

    #if 0
    /* For compression write out */
    if (step_id == 0)
    {
        float_data = GetDataAsType(byte_values_);
    } else
    {
        float_data = GetDataAsType(prefixes_.back().prefixes);
    }
    #else
    /* For decompression write eout */
    float_data = GetDataAsType(byte_values_);
    #endif

    t8_forest_t current_mesh = GetMesh();

    /* Create a new vtk field holding the element data arrays */
    t8_vtk_data_field_t *vtk_data = new t8_vtk_data_field_t[1];

        /* Set the type of the data and pointer to the data */
        std::ignore = snprintf(vtk_data[0].description, 128, "%s", "tracer");
        
        vtk_data[0].type = T8_VTK_SCALAR;

        std::vector<double> converted_data;
        converted_data.reserve(float_data.size());
        for (auto fiter = float_data.begin(); fiter != float_data.end(); ++fiter)
        {
            converted_data.push_back(static_cast<double>(*fiter));
        }

        vtk_data[0].data = converted_data.data();

        const int vtk_err = t8_forest_vtk_write_file(current_mesh, file_name.c_str(), 0, 0, 0, 0, 0, 1, vtk_data);
        
        if (vtk_err == 0)
            cmc_err_msg("An error occrued during the creation of the t8code-forest vtk file.");
    
    delete[] vtk_data;
}

inline
void
ByteVar::DecompressionRefine(const int elem_id)
{
    std::visit([&elem_id](auto&& var){
        var.DecompressionRefine(elem_id);
    }, var_);
}

inline
void
ByteVar::DecompressionLeaveElementUnchanged(const int elem_id)
{
    std::visit([&elem_id](auto&& var){
        var.DecompressionLeaveElementUnchanged(elem_id);
    }, var_);
}

inline
void
ByteVar::RefineAndApplyCommonPrefix(const int elem_id, const std::vector<uint8_t> prefix, const int num_prefix_bits)
{
    std::visit([&](auto&& var){
        var.RefineAndApplyCommonPrefix(elem_id, prefix, num_prefix_bits);
    }, var_);
}

inline
void
ByteVar::AssignCompressionValueForDecompressionStart()
{
    std::visit([](auto&& var){
        var.AssignCompressionValueForDecompressionStart();
    }, var_);
}

inline
void
ByteVar::LeaveElementUnchangedAndApplyPrefix(const int elem_id, std::vector<uint8_t> prefix, const int num_prefix_bits)
{
    std::visit([&](auto&& var){
        var.LeaveElementUnchangedAndApplyPrefix(elem_id, prefix, num_prefix_bits);
    }, var_);
}

inline
void
ByteVar::InitializeSuffixDecompression()
{
    std::visit([](auto&& var){
        var.InitializeSuffixDecompression();
    }, var_);
}

inline
void
ByteVar::InitializeResidualAlphabet()
{
    std::visit([](auto& var){
        var.InitializeResidualAlphabet();
    }, var_);
}

inline
void
ByteVar::TryLZCEncoding()
{
    std::visit([](auto& var){
        var.TryLZCEncoding();
    }, var_);
}

inline
void
ByteVar::ExtractMeanAndLeaveDifferencesFromInitialData(const int elem_start_index, const int num_elements)
{
    std::visit([&](auto& var){
        var.ExtractMeanAndLeaveDifferencesFromInitialData(elem_start_index, num_elements);
    }, var_);
}

inline
void
ByteVar::ExtractMeanAndLeaveDifferencesFromPreviousMeans(const int elem_start_index, const int num_elements)
{
    std::visit([&](auto& var){
        var.ExtractMeanAndLeaveDifferencesFromPreviousMeans(elem_start_index, num_elements);
    }, var_);
}

inline
void
ByteVar::_TestComparisonPCPLightAMRFromInitialData(const int elem_start_index, const int num_elements)
{
    #if _PERFORM_TEST_COMPARISON
    std::visit([&](auto& var){
        var._TestComparisonPCPLightAMRFromInitialData(elem_start_index, num_elements);
    }, var_);
    #endif
}

inline
void
ByteVar::_TestComparisonPCPLightAMRFromPreviousMeans(const int elem_start_index, const int num_elements)
{
    #if _PERFORM_TEST_COMPARISON
    std::visit([&](auto& var){
        var._TestComparisonPCPLightAMRFromPreviousMeans(elem_start_index, num_elements);
    }, var_);
    #endif
}

inline
void
ByteVar::LeaveValueUnchangedForNextMeanComputation(const int elem_start_index)
{
    std::visit([&](auto& var){
        var.LeaveValueUnchangedForNextMeanComputation(elem_start_index);
    }, var_);
}

inline
void
ByteVar::LeaveValueUnchanged(const int elem_start_index)
{
    std::visit([&](auto& var){
        var.LeaveValueUnchanged(elem_start_index);
    }, var_);
}

inline
void
ByteVar::ExtractMeanFromInitialData(const int elem_start_index, const int num_elements)
{
    std::visit([&](auto& var){
        var.ExtractMeanFromInitialData(elem_start_index, num_elements);
    }, var_);
}

inline
void
ByteVar::ExtractMeanFromPreviousMeans(const int elem_start_index, const int num_elements)
{
    std::visit([&](auto& var){
        var.ExtractMeanFromPreviousMeans(elem_start_index, num_elements);
    }, var_);
}

inline
size_t
ByteVar::GetCountSignificantBits(const int elem_id) const
{
    return std::visit([&elem_id](const auto& var) -> size_t {
        return var.GetCountSignificantBits(elem_id);
    }, var_);
}

inline
void
ByteVar::StoreElementUnchanged(const int elem_id)
{
    std::visit([&](auto& var){
        var.StoreElementUnchanged(elem_id);
    }, var_);
}

inline
void
ByteVar::ApplyResidualAndStoreElement(const int elem_id, const uint32_t encoded_lzc, const std::vector<uint8_t>& residual_bits)
{
    std::visit([&](auto& var){
        var.ApplyResidualAndStoreElement(elem_id, encoded_lzc, residual_bits);
    }, var_);
}

inline
void
ByteVar::ApplyResidualWithoutImplicitOneBitAndStoreElement(const int elem_id, const uint32_t encoded_lzc, const std::vector<uint8_t>& residual_bits)
{
    std::visit([&](auto& var){
        var.ApplyResidualWithoutImplicitOneBitAndStoreElement(elem_id, encoded_lzc, residual_bits);
    }, var_);
}

template<typename T>
size_t
ByteVariable<T>::GetCountSignificantBits(const int elem_id) const
{
    return byte_values_[elem_id].GetCountOfSignificantBits();
}

template<typename T>
void
ByteVariable<T>::InitializeSuffixDecompression()
{
    byte_values_new_.reserve(byte_values_.size());
}


template<typename T>
void
ByteVariable<T>::LeaveElementUnchangedAndApplyPrefix(const int elem_id, std::vector<uint8_t> prefix, const int num_prefix_bits)
{
    CompressionValue<sizeof(T)> value = byte_values_[elem_id];
    value.ApplyPrefix(prefix, num_prefix_bits);

    /* Copy the element over */
    byte_values_new_.push_back(value);
}

template<typename T>
void
ByteVariable<T>::AssignCompressionValueForDecompressionStart()
{
    /* Push back an empty CompressionValue */
    byte_values_.push_back(CompressionValue<sizeof(T)>());
}

template<typename T>
void
ByteVariable<T>::DecompressionRefine(const int elem_id)
{
    /* Since only a refinement is applied, we copy the value over several times */
    const int num_copies = std::pow(2, GetDimensionalityOfDataLayout(this->GetInitialDataLayout()));

    /* Insert the element as many times */
    byte_values_new_.insert(byte_values_new_.end(), num_copies, byte_values_[elem_id]);
}

template<typename T>
void
ByteVariable<T>::DecompressionLeaveElementUnchanged(const int elem_id)
{
    /* Copy the element over */
    byte_values_new_.push_back(byte_values_[elem_id]);
}

template<typename T>
void
ByteVariable<T>::RefineAndApplyCommonPrefix(const int elem_id, const std::vector<uint8_t> prefix, const int num_prefix_bits)
{
    /* Since only a refinement is applied, we copy the value over several times */
    const int num_copies = std::pow(2, GetDimensionalityOfDataLayout(this->GetInitialDataLayout()));

    CompressionValue<sizeof(T)> value = byte_values_[elem_id];
    value.ApplyPrefix(prefix, num_prefix_bits);

    /* Insert the element as many times */
    byte_values_new_.insert(byte_values_new_.end(), num_copies, value);
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
ByteVariable<T>::GetMaximumTailToggledValue(const std::vector<PermittedError>& permitted_errors, const double previous_abs_deviation, const CompressionValue<sizeof(T)>& initial_serialized_value, const T& missing_value) const
{
    bool is_toogling_progressing = true;
    CompressionValue<sizeof(T)> toggled_value = initial_serialized_value;

    const T initial_value = initial_serialized_value.template ReinterpretDataAs<T>();

    int iteration_count = 0;
    const int max_iteration_count = sizeof(T) * CHAR_BIT;

    while (is_toogling_progressing && iteration_count < max_iteration_count)
    {
        const CompressionValue<sizeof(T)> save_previous_value = toggled_value;

        /* Toggle all ones up until the next unset bit (inclusive) */
        toggled_value.ToggleTailUntilNextUnsetBit();
        const T reinterpreted_value = toggled_value.template ReinterpretDataAs<T>();

        /* Check if it is error compliant */
        const ErrorCompliance error_evaluation = utilities_.IsValueErrorCompliantRegardingPreviousDeviations(permitted_errors, initial_value, reinterpreted_value, previous_abs_deviation, missing_value);

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
ByteVariable<T>::GetMaximumTailClearedValue(const std::vector<PermittedError>& permitted_errors, const double previous_abs_deviation, const CompressionValue<sizeof(T)>& initial_serialized_value, const T& missing_value) const
{
    bool is_clearing_progressing = true;
    CompressionValue<sizeof(T)> cleared_value = initial_serialized_value;

    const T initial_value = initial_serialized_value.template ReinterpretDataAs<T>();

    int iteration_count = 0;
    const int max_iteration_count = sizeof(T) * CHAR_BIT;

    while (is_clearing_progressing && iteration_count < max_iteration_count)
    {
        const CompressionValue<sizeof(T)> save_previous_value = cleared_value;

        /* Clear the next set bit from the tail */
        cleared_value.ClearNextSetBitFromTail();
        const T reinterpreted_value = cleared_value.template ReinterpretDataAs<T>();

        /* Check if it is error compliant */
        const ErrorCompliance error_evaluation = utilities_.IsValueErrorCompliantRegardingPreviousDeviations(permitted_errors, initial_value, reinterpreted_value, previous_abs_deviation, missing_value);

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

        /* Clear the next set bit from the tail */
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
CompressionValue<sizeof(T)>
ByteVariable<T>::GetMaximumTailToggledValueRegardingUncompressedStates(const std::vector<PermittedError>& permitted_errors, const std::vector<DomainIndex>& initial_data_indices, const CompressionValue<sizeof(T)>& initial_serialized_value, const T& missing_value) const
{
    cmc_assert(are_uncompressed_states_stored_);

    bool is_toogling_progressing = true;
    CompressionValue<sizeof(T)> toggled_value = initial_serialized_value;

    int iteration_count = 0;
    const int max_iteration_count = sizeof(T) * CHAR_BIT;

    while (is_toogling_progressing && iteration_count < max_iteration_count)
    {
        /* Temporarily store the previous value */
        const CompressionValue<sizeof(T)> save_previous_value = toggled_value;

        /* Toggle all ones up until the next unset bit (inclusive) */
        toggled_value.ToggleTailUntilNextUnsetBit();
        const T reinterpreted_value = toggled_value.template ReinterpretDataAs<T>();

        /* Check if the toggled value complies to the error thresholds */
        const ErrorCompliance error_evaluation = utilities_.IsCoarseningErrorCompliantRegardingInitialData(permitted_errors, uncompressed_data_, initial_data_indices, reinterpreted_value, missing_value);

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
ByteVariable<T>::GetMaximumTailClearedValueRegardingUncompressedStates(const std::vector<PermittedError>& permitted_errors, const std::vector<DomainIndex>& initial_data_indices, const CompressionValue<sizeof(T)>& initial_serialized_value, const T& missing_value) const
{
    cmc_assert(are_uncompressed_states_stored_);

    bool is_clearing_progressing = true;
    CompressionValue<sizeof(T)> cleared_value = initial_serialized_value;

    int iteration_count = 0;
    const int max_iteration_count = sizeof(T) * CHAR_BIT;

    while (is_clearing_progressing && iteration_count < max_iteration_count)
    {
        /* Temporarily store the previous value */
        const CompressionValue<sizeof(T)> save_previous_value = cleared_value;

        /* Clear the next set bit from the tail */
        cleared_value.ClearNextSetBitFromTail();
        const T reinterpreted_value = cleared_value.template ReinterpretDataAs<T>();

        /* Check if it is error compliant */
        const ErrorCompliance error_evaluation = utilities_.IsCoarseningErrorCompliantRegardingInitialData(permitted_errors, uncompressed_data_, initial_data_indices, reinterpreted_value, missing_value);

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

#if 0
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
                /* If the toggling approach has been more successful */
                *val_iter = toggled_value;
            } else
            {
                /* If the clearing approach has been more successful */
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
#endif


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
            //const std::vector<PermittedError> permitted_errors = GetRemainingMaxAllowedAbsoluteError(index);
            const std::vector<PermittedError> permitted_errors = GetPermittedError(index);
            //const double current_abs_deviation = utilities_.GetPreviousDeviation(index);
            const double current_abs_deviation = 0.0;

            /* Get the value which has been transformed by toggling as many ones from the back while setting the succeeding 'zero' bits to one */
            const CompressionValue<sizeof(T)> toggled_value = GetMaximumTailToggledValue(permitted_errors, current_abs_deviation, *val_iter, missing_value);

            /* Get the value which has been transformed by clearing as many of the last set bits as possible */
            const CompressionValue<sizeof(T)> cleared_value = GetMaximumTailClearedValue(permitted_errors, current_abs_deviation, *val_iter, missing_value);

            /* Check which approach leads to more zero bits at the end */
            const int num_toogled_trailing_zeros = toggled_value.GetNumberTrailingZeros();
            const int num_cleared_trailing_zeros = cleared_value.GetNumberTrailingZeros();

            /* Replace the initial value with the transformed one */
            if (num_toogled_trailing_zeros >= num_cleared_trailing_zeros)
            {
                /* If the toggling approach has been more successful */
                *val_iter = toggled_value;
            } else
            {
                /* If the clearing approach has been more successful */
                *val_iter = cleared_value;
            }

            /* Update the trail bit count for the new value */
            val_iter->UpdateTrailBitCount();
            
            //if (index % 12000 == 0)
                //cmc_debug_msg("Trailing Zeros: Toggled: ", num_toogled_trailing_zeros, ", Cleared: ", num_cleared_trailing_zeros);
        } else
        {
            /* In order to not chenge missing values, we are just able to trim their trailing zeros */
            (*val_iter).UpdateTrailBitCount();
        }
    }
}

//TODO: For load-balancing, we would need a weighted partition of the mesh/data for this functionality 
template<typename T>
void
ByteVariable<T>::PerformTailTruncationRegardingUncompressedStates()
{
    cmc_assert(are_uncompressed_states_stored_ == true);
    cmc_assert(t8_forest_get_local_num_leaf_elements(mesh_.GetMesh()) == byte_values_.size());
    cmc_assert(is_initial_data_kept_ == true);

    /* The misssing value associated with the variable */
    const T missing_value = attributes_.GetMissingValue();
    
    /* Number of local elements/data points*/
    const t8_locidx_t num_elements = t8_forest_get_local_num_leaf_elements(mesh_.GetMesh());
    
    /* We need the eclass of the tree as well as the scheme */
    const t8_eclass_t eclass = t8_forest_get_tree_class (mesh_.GetMesh(), 0);
    const t8_scheme_c* ts =  t8_forest_get_scheme (mesh_.GetMesh());

    cmc_debug_msg("In perform tail truncation before element loop, num elems: ", num_elements);
    cmc_debug_msg("Size of byte values_: ", byte_values_.size());
    /* Iterate over all elements/data points */
    for (t8_locidx_t idx = 0; idx < num_elements; ++idx)
    {
        if (idx % 30000 == 0)
        {
            cmc_debug_msg("current index in iter : ", idx);
        }
        if (!ApproxCompare(byte_values_[idx].template ReinterpretDataAs<T>(), missing_value))
        {
            /* Get the permitted error for the current values */
            const std::vector<PermittedError> permitted_errors = GetPermittedError(idx);

            //cmc_debug_msg("Permitted Error has been accessed");
            /* Get the current element */
            const t8_element_t* element = t8_forest_get_leaf_element_in_tree(mesh_.GetMesh(), 0, idx);
            //cmc_debug_msg("The element with idx ", idx, " has been accessed");

            /* Get all initial elements indices that are covered by this element */
            const std::vector<DomainIndex> initial_value_coverage = GetInitialElementCoverage(uncompressed_mesh_, ts, element);
            if (idx % 30000 == 0)
            {
                cmc_debug_msg("Initial value coverage size is: ", initial_value_coverage.size());
            }
            //cmc_debug_msg("The initial value coverage has been computed: ", initial_value_coverage.size());

            /* Get the value which has been transformed by toggling as many ones from the back while setting the succeeding 'zero' bits to one */
            const CompressionValue<sizeof(T)> toggled_value = GetMaximumTailToggledValueRegardingUncompressedStates(permitted_errors, initial_value_coverage, byte_values_[idx], missing_value);
            //cmc_debug_msg("The toggled value has been computed");
            /* Get the value which has been transformed by clearing as many of the last set bits as possible */
            const CompressionValue<sizeof(T)> cleared_value = GetMaximumTailClearedValueRegardingUncompressedStates(permitted_errors, initial_value_coverage, byte_values_[idx], missing_value);
            //cmc_debug_msg("The cleared value has been computed");
            /* Check which approach leads to more zero bits at the end */
            const int num_toogled_trailing_zeros = toggled_value.GetNumberTrailingZeros();
            const int num_cleared_trailing_zeros = cleared_value.GetNumberTrailingZeros();

            /* Replace the initial value with the transformed one */
            if (num_toogled_trailing_zeros >= num_cleared_trailing_zeros)
            {
                /* If the toggling approach has been more successful */
                byte_values_[idx] = toggled_value;
            } else
            {
                /* If the clearing approach has been more successful */
                byte_values_[idx] = cleared_value;
            }
            //cmc_debug_msg("Best approach has beene valuated and assigned");
            /* Update the tail bit count for the new value */
            byte_values_[idx].UpdateTrailBitCount();
            //cmc_debug_msg("Idx: ", idx, ", initial_value_coverage: ", initial_value_coverage.size(), ", toggled trail zeros: ", num_toogled_trailing_zeros, ", cleared trail zeros: ", num_cleared_trailing_zeros);
        } else
        {
            /* In order to not change missing values, we are just able to trim their trailing zeros */
            byte_values_[idx].UpdateTrailBitCount();
            //cmc_debug_msg("In else, idx: ", idx);
        }
    }



    /* Iterate through the serialized values and try to emplace as many zeros at the tail as possible (compliant to the error threshold) */
    int index = 0;
    for (auto val_iter = byte_values_.begin(); val_iter != byte_values_.end(); ++val_iter, ++index)
    {
        if (!ApproxCompare(initial_data_[index], missing_value))
        {
            /* Get the permitted error for the current values */
            //TODO: Revert! Only for now with absolute errors
            const std::vector<PermittedError> permitted_errors = GetPermittedError(index);

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
            /* In order to not change missing values, we are just able to trim their trailing zeros */
            (*val_iter).UpdateTrailBitCount();
        }
    }
}

template<typename T>
void
ByteVariable<T>::PerformTailTruncationOnInitialData()
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
            const std::vector<PermittedError> permitted_errors = GetPermittedError(index);

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

    const t8_element_t* element = t8_forest_get_leaf_element_in_tree(mesh_.GetMesh(), 0, index);
    
    const t8_scheme_c* ts =  t8_forest_get_scheme (mesh_.GetMesh());

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
    } else
    {
        cmc_debug_msg("No permitted error could be associated with the element (Index: ", index, ").");
        /* In case no error is associated, we state that no error is permitted */
        is_abs_error_present = true;
        abs_err = 0.0;
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
    interpolation_indications_.emplace_back();
    interpolation_indications_.back().Reserve(mesh_.GetNumberLocalElements());

    #if _PERFORM_TEST_COMPARISON
    test_encodings_.emplace_back();
    #endif
}

template<typename T>
void
ByteVariable<T>::InitializeDecompressionIteration()
{
    /* Allocate an upper bound of memory. This resembles the case, when all elements are refined during the decompression */
    byte_values_new_.reserve(mesh_.GetNumberLocalElements() * std::pow(2, GetDimensionalityOfDataLayout(this->GetInitialDataLayout())));
}

template<typename T>
void
ByteVariable<T>::PopLevelwisePrefixEncoding()
{
    prefixes_.pop_back();
}

/* If a single element is investigated and not a family of elements during the prefix extraction,
 * we are not altering the value/previous-prefix and keep it as a prefix for the next level */
template<typename T>
void
ByteVariable<T>::LeaveCoarsePrefixUnchanged(const int elem_index)
{
    /* In case that the element stays the same (the sibling elements are not yet present in the mesh).
     * We indicate the whole data/prefix as a common prefix and will evaluate in a later adaptation step. */
    prefixes_.back().SetRefinementIndicatorBit(false);
    prefixes_.back().SetPrefixIndicatorBit(true); //TODO: Wont this be erased later by the WriteCompressed nevertheless, can we directly set it to false?

    /* Copy the data/prefix to the next 'coarser level' */
    if (prefixes_.size() >= 2)
    {
        /* In any other prefix extraction iteration than the first, we just copy the previous prefix */
        LevelwisePrefixData<T>& previous_prefixes = (*std::prev(prefixes_.end(), 2));
        prefixes_.back().SetPrefix(previous_prefixes.prefixes[elem_index]);
        previous_prefixes.prefixes[elem_index] = CompressionValue<sizeof(T)>();
    } else
    {
        /* In the first prefix extraction iteration, we just keep the initial byte value and consider it in a 'coarser' iteration */
        prefixes_.back().SetPrefix(byte_values_[elem_index]);
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
ByteVariable<T>::FinalizeDecompressionIteration()
{
    byte_values_ = std::move(byte_values_new_);
    byte_values_new_.clear();
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
ByteVariable<T>::InitializeResidualAlphabet()
{
    /* Size of the type of this variable */
    const size_t type_size = sizeof(T);
    /* We want to indicate the position of the first 1 in binary represenatation index, e.g. 4 byte type -> [0,32] plus the signum */
    #if 1
    const size_t num_residual_codes_wo_sign = type_size * bit_map::kCharBit + 1;
    const size_t num_residual_codes = 2 * num_residual_codes_wo_sign;
    #else
    const size_t num_residual_codes_wo_sign = type_size * bit_map::kCharBit + 1;
    const size_t num_residual_codes = num_residual_codes_wo_sign;
    #endif

    /* Allocate memory for the alphabet */
    alphabet_.reserve(num_residual_codes);

    uint32_t ResidualIndicatorFirstOne = 0;

    //TODO: Try to remove letters from the alphabet that do not occur 
    //....

    for (size_t i = 0; i < num_residual_codes_wo_sign; ++i, ++ResidualIndicatorFirstOne)
    {
        /* Emplace back the index of the first one with a positive and negative signum */
        //alphabet_.emplace_back{arithmetic_encoding::Letter{ResidualIndicatorFirstOne, 1}};
        //alphabet_.emplace_back{arithmetic_encoding::Letter{ResidualIndicatorFirstOne + arithmetic_encoding::kMSBBit, 1}};
        #if 1
        //alphabet_.insert({ResidualIndicatorFirstOne, 1});
        //alphabet_.insert({arithmetic_encoding::kMSBBit + ResidualIndicatorFirstOne, 1});
        #else 
        alphabet_.insert({ResidualIndicatorFirstOne, 1});
        #endif
    }
}


template<typename T>
void
ByteVariable<T>::TryLZCEncoding()
{
    if constexpr (std::is_same_v<T, float>)
    {
        const T missing_value = attributes_.GetMissingValue();

        std::vector<T> lzc_byte_vals;
        lzc_byte_vals.reserve(byte_values_.size());

        uint32_t sig_bit_counter = 0;
        for(auto val_iter = byte_values_.begin(); val_iter != byte_values_.end(); ++val_iter)
        {
            uint32_t lzc{sizeof(T) * CHAR_BIT};

            T reinterpreted_val = val_iter->template ReinterpretDataAs<T>();

            if (!ApproxCompare(reinterpreted_val, missing_value))
            {
                lzc = static_cast<uint32_t>(val_iter->GetNumberLeadingZeros());
            }

            float flzc;
            std::memcpy(&flzc, &lzc, 4);

            *val_iter = flzc;

            if (lzc > 0 && lzc < 32)
            {
                ++lzc;
            }
            sig_bit_counter += (sizeof(T) * CHAR_BIT - lzc);

        }

        cmc_debug_msg("\n\n\n\nNum remaining significant bits are: ", sig_bit_counter, " (in bytes: ", sig_bit_counter / 8, ")");
    }
}

template<typename T>
T
GetInterpoaltionMaximizingLZC(std::vector<CompressionValue<sizeof(T)>>& values, const int elem_start_index, const int num_elements, const T missing_value)
{
    if constexpr (std::is_same_v<float, T>)
    {
    const std::vector<T> vals = GetCompressionValuesAs<T>(values, elem_start_index, num_elements);

    /* Get a view on the values */
    //VectorView<uint32_t> value_view(vals.data(), num_elements);
    const VectorView<T> value_view(vals.data(), num_elements);

    int max_lzc = -1;
    int index = 0;

    int pred_idx = 0;
    for (auto pred_iter = value_view.begin(); pred_iter != value_view.end(); ++pred_iter, ++pred_idx)
    {
        T current_predictor = *pred_iter;
        uint32_t cp;
        std::memcpy(&cp, &current_predictor, 4);

        int lzc_count = 0;

        for (auto val_iter = value_view.begin(); val_iter != value_view.end(); ++val_iter)
        {
            T val = *val_iter;
            uint32_t uiv;
            std::memcpy(&uiv, &val, 4);

            uint32_t diff = (uiv >= cp ? uiv - cp : cp - uiv);

            CompressionValue<sizeof(T)> crv(diff);
            lzc_count += crv.GetNumberLeadingZeros();
        }

        if (lzc_count > max_lzc)
        {
            max_lzc = lzc_count;
            index = pred_idx;
        }
    }

    //Try midrange and arithmetic mean as well
    const T mid_range = InterpolateToMidRange(value_view, num_elements, missing_value);
    uint32_t ui_mid_range;
    std::memcpy(&ui_mid_range, &mid_range, 4);

    int lzc_count = 0;

    for (auto val_iter = value_view.begin(); val_iter != value_view.end(); ++val_iter)
    {
        T val = *val_iter;
        uint32_t uiv;
        std::memcpy(&uiv, &val, 4);

        uint32_t diff = (uiv >= ui_mid_range ? uiv - ui_mid_range : ui_mid_range - uiv);

        CompressionValue<sizeof(T)> crv(diff);
        lzc_count += crv.GetNumberLeadingZeros();
    }

    if (lzc_count > max_lzc)
    {
        return mid_range;
    }

    return value_view[index];

    } else
    {
        return T(0);
    }
}

template<typename T>
T
GetInterpoaltionMaximizingLZCwoMissingValue(std::vector<CompressionValue<sizeof(T)>>& values, const int elem_start_index, const int num_elements, const T missing_value)
{
    if constexpr (std::is_same_v<float, T>)
    {
    /* Get all converted values */
    const std::vector<T> vals = GetCompressionValuesAs<T>(values, elem_start_index, num_elements);

    std::vector<T> non_missing_vals;
    non_missing_vals.reserve(vals.size());

    /* Extarct all non-missing values */
    for(auto val_iter = vals.begin(); val_iter != vals.end(); ++val_iter)
    {
        if (not ApproxCompare(*val_iter, missing_value))
        {
            non_missing_vals.push_back(*val_iter);
        }
    }

    /* If only missing values are present, we return a missig value */
    if (non_missing_vals.empty())
    {
        return missing_value;
    }

    /* Get a view on the values */
    const VectorView<T> value_view(vals.data(), num_elements);

    int max_lzc = -1;
    int index = 0;

    int pred_idx = 0;

    /* We choose a non missing value as a predictor in order to have a constructive prediction, but check its effect on all values of the family (i.e. missing values) */
    for (auto pred_iter = non_missing_vals.begin(); pred_iter != non_missing_vals.end(); ++pred_iter, ++pred_idx)
    {
        T current_predictor = *pred_iter;
        uint32_t cp;
        std::memcpy(&cp, &current_predictor, 4);

        int lzc_count = 0;

        for (auto val_iter = value_view.begin(); val_iter != value_view.end(); ++val_iter)
        {
            T val = *val_iter;
            uint32_t uiv;
            std::memcpy(&uiv, &val, 4);

            uint32_t diff = (uiv >= cp ? uiv - cp : cp - uiv);

            CompressionValue<sizeof(T)> crv(diff);
            lzc_count += crv.GetNumberLeadingZeros();
        }

        if (lzc_count > max_lzc)
        {
            max_lzc = lzc_count;
            index = pred_idx;
        }
    }

    //Try midrange and arithmetic mean as well
    const T mid_range = InterpolateToMidRange(VectorView<T>(non_missing_vals), static_cast<int>(non_missing_vals.size()), missing_value);
    uint32_t ui_mid_range;
    std::memcpy(&ui_mid_range, &mid_range, 4);

    int lzc_count = 0;

    for (auto val_iter = value_view.begin(); val_iter != value_view.end(); ++val_iter)
    {
        T val = *val_iter;
        uint32_t uiv;
        std::memcpy(&uiv, &val, 4);

        const uint32_t diff = (uiv >= ui_mid_range ? uiv - ui_mid_range : ui_mid_range - uiv);

        CompressionValue<sizeof(T)> crv(diff);
        lzc_count += crv.GetNumberLeadingZeros();
    }

    if (lzc_count > max_lzc)
    {
        return mid_range;
    }

    return non_missing_vals[index];

    } else
    {
        return T(0);
    }
}



template<typename T>
T
GetInterpoaltionMaximizingLZCInXOR(std::vector<CompressionValue<sizeof(T)>>& values, const int elem_start_index, const int num_elements, const T missing_value)
{
    if constexpr (std::is_same_v<float, T>)
    {
    const std::vector<T> vals = GetCompressionValuesAs<T>(values, elem_start_index, num_elements);

    /* Get a view on the values */
    const VectorView<T> value_view(vals.data(), num_elements);

    int max_lzc = -1;
    int index = 0;

    int pred_idx = 0;
    for (auto pred_iter = value_view.begin(); pred_iter != value_view.end(); ++pred_iter, ++pred_idx)
    {
        T current_predictor = *pred_iter;
        uint32_t cp;
        std::memcpy(&cp, &current_predictor, 4);

        uint32_t xor_result{0};

        for (auto val_iter = value_view.begin(); val_iter != value_view.end(); ++val_iter)
        {
            T val = *val_iter;
            uint32_t uiv;
            std::memcpy(&uiv, &val, 4);

            const uint32_t diff = cp ^ uiv;

            xor_result |= diff;
        }        

        CompressionValue<sizeof(T)> crv(xor_result);
        const int lzc_count = crv.GetNumberLeadingZeros();

        if (lzc_count > max_lzc)
        {
            max_lzc = lzc_count;
            index = pred_idx;
        }
    }

    //Try midrange and arithmetic mean as well
    const T mid_range = InterpolateToMidRange(value_view, num_elements, missing_value);
    uint32_t ui_mid_range;
    std::memcpy(&ui_mid_range, &mid_range, 4);

    uint32_t mr_xor_result{0};

    for (auto val_iter = value_view.begin(); val_iter != value_view.end(); ++val_iter)
    {
        T val = *val_iter;
        uint32_t uiv;
        std::memcpy(&uiv, &val, 4);

        const uint32_t diff = ui_mid_range ^ uiv;

        mr_xor_result |= diff;
    }    

    CompressionValue<sizeof(T)> mr_crv(mr_xor_result);
    const int mr_lzc_count = mr_crv.GetNumberLeadingZeros();

    /* Check the arithmetic mean as well */
    const T mean = InterpolateToArithmeticMean(value_view, num_elements, missing_value);
    uint32_t ui_mean;
    std::memcpy(&ui_mean, &mean, 4);

    uint32_t mean_xor_result{0};

    for (auto val_iter = value_view.begin(); val_iter != value_view.end(); ++val_iter)
    {
        T val = *val_iter;
        uint32_t uiv;
        std::memcpy(&uiv, &val, 4);

        const uint32_t diff = ui_mean ^ uiv;

        mean_xor_result |= diff;
    }    

    CompressionValue<sizeof(T)> mean_crv(mean_xor_result);
    const int mean_lzc_count = mean_crv.GetNumberLeadingZeros();

    /* Check which value to return */
    if (max_lzc >= mean_lzc_count && max_lzc >= mr_lzc_count)
    {
        return value_view[index];
    }  else if (mean_lzc_count >= max_lzc && mean_lzc_count >= mr_lzc_count)
    {
        return mean;
    } else if (mr_lzc_count >= max_lzc && mr_lzc_count >= mean_lzc_count)
    {
        return mid_range;
    } else
    {
        cmc_err_msg("One of the above cases should have been selected");
        return T(0.0);
    }

    } else
    {
        return T(0);
    }
}


template<typename T>
void
ByteVariable<T>::ExtractMeanAndLeaveDifferences(std::vector<CompressionValue<sizeof(T)>>& values, const int elem_start_index, const int num_elements)
{
    if constexpr (std::is_same_v<float, T>)
    {

    const std::vector<T> vals = GetCompressionValuesAs<T>(values, elem_start_index, num_elements);

    /* Get a view on the values */
    VectorView<T> value_view(vals.data(), num_elements);

    const T missing_value = attributes_.GetMissingValue();

    //Choose interpolation fucntion here

    //const T real_interpolation_result = InterpolateToArithmeticMean(value_view, num_elements, missing_value);

    //const T real_interpolation_result = InterpolateToMidRange(value_view, num_elements, missing_value);

    //auto min_iter = std::min_element(vals.begin(), vals.end());
    //const T real_interpolation_result = *min_iter;

    //auto min_iter = std::max_element(vals.begin(), vals.end());
    //const T real_interpolation_result = *min_iter;

    #if 0
    const T real_interpolation_result = GetInterpoaltionMaximizingLZC<T>(values, elem_start_index, num_elements, missing_value);
    #else
    const T real_interpolation_result = GetInterpoaltionMaximizingLZCwoMissingValue<T>(values, elem_start_index, num_elements, missing_value);
    #endif

    uint32_t interpolation_result;
    std::memcpy(&interpolation_result, &real_interpolation_result, sizeof(T));

    #if 0
    uint32_t int_missing_value;
    std::memcpy(&int_missing_value, &missing_value, sizeof(T));

    const uint32_t interpolation_result = InterpolateToArithmeticMean(value_view, GetMesh(), elem_start_index, num_elements, int_missing_value);
    #endif
    
    for (int val_idx = 0; val_idx < num_elements; ++val_idx)
    {
        uint32_t intval;
        std::memcpy(&intval, &vals[val_idx], sizeof(T));

        uint32_t encoded_lzc;

        if (interpolation_result >= intval)
        {
            const uint32_t diff = interpolation_result - intval;

            values[elem_start_index + val_idx] = CompressionValue<sizeof(T)>(diff);
            /* Indicate that the interpolation result is greater than the value */
            interpolation_indications_.back().AppendSetBit();
            /* Get the encoded LZC in order to model the frequency */
            encoded_lzc = values[elem_start_index + val_idx].GetNumberLeadingZeros();
        } else
        {
            const uint32_t diff = intval - interpolation_result;
            values[elem_start_index + val_idx] = CompressionValue<sizeof(T)>(diff);
            /* Indicate that the interpolation result is greater than the value */
            interpolation_indications_.back().AppendUnsetBit();
            /* Get the encoded LZC in order to model the frequency */
            encoded_lzc = arithmetic_encoding::kMSBBit + values[elem_start_index + val_idx].GetNumberLeadingZeros();
        }

        /* Increment the counter for the encoded LZC */
        auto alphaiter = alphabet_.find(encoded_lzc);
        if (alphaiter == alphabet_.end()) {alphabet_.insert({encoded_lzc, 1});}
        else{++(alphabet_[encoded_lzc]);}
    }

    /* After the differences have been calculated and stored in the child elements, we store the interpolation result on the coarser level */
    prefixes_.back().SetPrefix(CompressionValue<sizeof(T)>(real_interpolation_result));

    } else
    {
        cmc_err_msg("Only float compression possible");
    }
}

template<typename T>
void
ByteVariable<T>::_TestComparisonPCPLightAMR(std::vector<CompressionValue<sizeof(T)>>& values, const int elem_start_index, const int num_elements)
{
    #if _PERFORM_TEST_COMPARISON
    if constexpr (std::is_same_v<float, T>)
    {
    const std::vector<T> vals = GetCompressionValuesAs<T>(values, elem_start_index, num_elements);

    /* Get a view on the values */
    VectorView<T> value_view(vals.data(), num_elements);

    const T missing_value = attributes_.GetMissingValue();

    #if 1
    const T real_interpolation_result = GetInterpoaltionMaximizingLZCInXOR<T>(values, elem_start_index, num_elements, missing_value);
    #else
    //const T real_interpolation_result = GetInterpoaltionMaximizingLZCwoMissingValue<T>(values, elem_start_index, num_elements, missing_value);
    #endif

    uint32_t interpolation_result;
    std::memcpy(&interpolation_result, &real_interpolation_result, sizeof(T));
    

    uint32_t xor_result{0};

    for (int val_idx = 0; val_idx < num_elements; ++val_idx)
    {
        uint32_t uiv;
        std::memcpy(&uiv, &vals[val_idx], sizeof(T));

        const uint32_t xor_residual = interpolation_result ^ uiv;

        values[elem_start_index + val_idx] = CompressionValue<sizeof(T)>(xor_residual);

        xor_result |= xor_residual;
    }

    CompressionValue<sizeof(T)> crv(xor_result);
    int lzc = crv.GetNumberLeadingZeros();

    /* Check whether the LZC exceeds 15 */
    if (lzc > 15)
    {
        lzc = 15;
    }

    const uint8_t final_lzc = static_cast<uint8_t>(lzc);

    test_encodings_.back().encoded_lzc.AppendFourBits(final_lzc);

    /* Append all remaining bit sequences from the family of elements */
    for (int val_idx = 0; val_idx < num_elements; ++val_idx)
    {
        CompressionValue<sizeof(T)>& cr_val = values[elem_start_index + val_idx];

        cr_val.SetFrontBit(lzc);

        EncodeAndAppendPrefix(test_encodings_.back().encoded_residuals, cr_val.GetSignificantBitsInBigEndianOrdering(), cr_val.GetCountOfSignificantBits());
    }

    /* After the differences have been calculated and stored in the child elements, we store the interpolation result on the coarser level */
    prefixes_.back().SetPrefix(CompressionValue<sizeof(T)>(real_interpolation_result));

    } else
    {
        cmc_err_msg("Only float compression possible");
    }

    #endif
}


template<typename T>
void
ByteVariable<T>::ExtractMeanAndLeaveDifferencesFromInitialData(const int elem_start_index, const int num_elements)
{
    /* Extract the mean from the initial data */
    ExtractMeanAndLeaveDifferences(byte_values_, elem_start_index, num_elements);
}

template<typename T>
void
ByteVariable<T>::ExtractMeanAndLeaveDifferencesFromPreviousMeans(const int elem_start_index, const int num_elements)
{
    /* Extract a common prefix from the initial data */
    ExtractMeanAndLeaveDifferences((*std::prev(prefixes_.end(), 2)).prefixes, elem_start_index, num_elements);
}

template<typename T>
void
ByteVariable<T>::_TestComparisonPCPLightAMRFromInitialData(const int elem_start_index, const int num_elements)
{
    #if _PERFORM_TEST_COMPARISON
    /* Extract the mean from the initial data */
    _TestComparisonPCPLightAMR(byte_values_, elem_start_index, num_elements);
    #endif
}

template<typename T>
void
ByteVariable<T>::_TestComparisonPCPLightAMRFromPreviousMeans(const int elem_start_index, const int num_elements)
{
    #if _PERFORM_TEST_COMPARISON
    /* Extract a common prefix from the initial data */
    _TestComparisonPCPLightAMR((*std::prev(prefixes_.end(), 2)).prefixes, elem_start_index, num_elements);
    #endif
}

template<typename T>
void
ByteVariable<T>::LeaveValueUnchangedForNextMeanComputation(const int elem_start_index)
{
    if constexpr (std::is_same_v<T, float>)
    {
    /* Copy the data/prefix to the next 'coarser level' */
    if (prefixes_.size() >= 2)
    {
        /* In any other prefix extraction iteration than the first, we just copy the previous prefix */
        LevelwisePrefixData<T>& previous_means = (*std::prev(prefixes_.end(), 2));
        prefixes_.back().SetPrefix(previous_means.prefixes[elem_start_index]);
        previous_means.prefixes[elem_start_index] = CompressionValue<sizeof(T)>(uint32_t{0});
    } else
    {
        /* In the first prefix extraction iteration, we just keep the initial byte value and consider it in a 'coarser' iteration */
        prefixes_.back().SetPrefix(byte_values_[elem_start_index]);
        /* Remove the 'prefix' from the initial data and replace it with an empty 'dummy value' */
        byte_values_[elem_start_index] = CompressionValue<sizeof(T)>(uint32_t{0});
    }
    /* There wont be a residual since, we just drag the value along */
    interpolation_indications_.back().AppendSetBit();
    const uint32_t encoded_lzc = sizeof(T) * bit_map::kCharBit;
    auto alphaiter = alphabet_.find(encoded_lzc);
    if (alphaiter == alphabet_.end()) {alphabet_.insert({encoded_lzc, 1});}
    else{++(alphabet_[encoded_lzc]);}
    }
}

template<typename T>
void
ByteVariable<T>::TryFittingPyramid()
{
    //Copy first three levels
    
}

template<typename T>
void
ByteVariable<T>::StoreElementUnchanged(const int elem_id)
{
    /* Just copy the value over to the new vector */
    byte_values_new_.push_back(byte_values_[elem_id]);
}


template<typename T>
void
ByteVariable<T>::ApplyResidualAndStoreElement(const int elem_id, const uint32_t encoded_lzc, const std::vector<uint8_t>& residual_bits)
{
    /* Get the current mean value and add/subtract the residual to/from it */
    CompressionValue<sizeof(T)> val = byte_values_[elem_id];
    //cmc_debug_msg("In Apply Residual, elem_id: ", elem_id, ", Current byte value: ", val.template ReinterpretDataAs<T>());
    val.AddIntegerResidual(encoded_lzc, residual_bits);
    /* The altered value is stores in the new vetor */
    byte_values_new_.push_back(val);
}

template<typename T>
void
ByteVariable<T>::ApplyResidualWithoutImplicitOneBitAndStoreElement(const int elem_id, const uint32_t encoded_lzc, const std::vector<uint8_t>& residual_bits)
{
    /* Get the current mean value and add/subtract the residual to/from it */
    CompressionValue<sizeof(T)> val = byte_values_[elem_id];
    val.AddXORResidualWithoutImplicitOneBit(encoded_lzc, residual_bits);
    /* The altered value is stores in the new vetor */
    byte_values_new_.push_back(val);
}

template<typename T>
void
ByteVariable<T>::ExtractMean(std::vector<CompressionValue<sizeof(T)>>& values, const int elem_start_index, const int num_elements)
{
    if constexpr (std::is_same_v<float, T>)
    {
    const std::vector<T> vals = GetCompressionValuesAs<T>(values, elem_start_index, num_elements);

    /* Get a view on the values */
    VectorView<T> value_view(vals.data(), num_elements);

    const T missing_value = attributes_.GetMissingValue();

    //const T real_interpolation_result = InterpolateToArithmeticMean(value_view, num_elements, missing_value);
    const T real_interpolation_result = InterpolateToMidRange(value_view, num_elements, missing_value);
    //auto min_iter = std::min_element(vals.begin(), vals.end());
    //const T real_interpolation_result = *min_iter;

    /* After the differences have been calculated and stores in the child elements, we store the interpolation result on the coarser level */
    prefixes_.back().SetPrefix(CompressionValue<sizeof(T)>(real_interpolation_result));

    } else
    {
        cmc_err_msg("Only float compression possible");
    }
}


template<typename T>
void
ByteVariable<T>::ExtractMeanFromInitialData(const int elem_start_index, const int num_elements)
{
    /* Extract the mean from the initial data */
    ExtractMean(byte_values_, elem_start_index, num_elements);
}

template<typename T>
void
ByteVariable<T>::ExtractMeanFromPreviousMeans(const int elem_start_index, const int num_elements)
{
    /* Get the previous means */
    //std::vector<CompressionValue<sizeof(T)>>& prev_means = (*std::prev(prefixes_.end(), 2)).prefixes;
    //cmc_debug_msg("Size of second to last prefixes: ", (*std::prev(prefixes_.end(), 2)).prefixes.size());
    /* Extract a common prefix from the initial data */
    ExtractMean((*std::prev(prefixes_.end(), 2)).prefixes, elem_start_index, num_elements);
}

template<typename T>
void
ByteVariable<T>::LeaveValueUnchanged(const int elem_start_index)
{
    if constexpr (std::is_same_v<T, float>)
    {
    /* Copy the data/prefix to the next 'coarser level' */
    if (prefixes_.size() >= 2)
    {
        /* In any other prefix extraction iteration than the first, we just copy the previous prefix */
        LevelwisePrefixData<T>& previous_means = (*std::prev(prefixes_.end(), 2));
        prefixes_.back().SetPrefix(previous_means.prefixes[elem_start_index]);
    } else
    {
        /* In the first prefix extraction iteration, we just keep the initial byte value and consider it in a 'coarser' iteration */
        prefixes_.back().SetPrefix(byte_values_[elem_start_index]);
    }
    }
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
ByteVariable<T>::ExtractCommonPrefixFromInitialData(const std::vector<int>& elem_ids)
{
    /* Extract a common prefix from the initial data */
    ExtractCommonPrefix(byte_values_, elem_ids);
}

template<typename T>
void
ByteVariable<T>::ExtractCommonPrefixFromPreviousPrefixes(const std::vector<int>& elem_ids)
{
    /* Get the previous prefixes */
    std::vector<CompressionValue<sizeof(T)>>& prev_prefixes = (*std::prev(prefixes_.end(), 2)).prefixes;

    /* Extract a common prefix from the initial data */
    ExtractCommonPrefix(prev_prefixes, elem_ids);  
}

inline
void
ByteVar::ExtractCommonPrefixFromInitialData(const std::vector<int>& elem_ids)
{
    std::visit([&](auto&& var){
        var.ExtractCommonPrefixFromInitialData(elem_ids);
    }, var_);
}

inline
void
ByteVar::ExtractCommonPrefixFromPreviousPrefixes(const std::vector<int>& elem_ids)
{
    std::visit([&](auto&& var){
        var.ExtractCommonPrefixFromPreviousPrefixes(elem_ids);
    }, var_);
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
            prefixes[elem_start_index + index].SetFrontBit(sizeof(T) * CHAR_BIT - prefix.GetTrailBit());
        }
    } else
    {
        IndicateNoPrefixFound();
    }
}


template<typename T>
void
ByteVariable<T>::ExtractCommonPrefix(std::vector<CompressionValue<sizeof(T)>>& prefixes, const std::vector<int>& element_ids)
{
    /* Evaluate whether there is a common prefix to extract */
    auto [is_prefix_found, prefix] = EvaluateCommonPrefix(prefixes, element_ids);

    if (is_prefix_found)
    {
        IndicatePrefixFound(prefix);

        /* We need to trim the previous prefixes by the extracted common prefix */
        for (auto elem_id_iter = element_ids.begin(); elem_id_iter != element_ids.end(); ++elem_id_iter)
        {
            prefixes[*elem_id_iter].SetFrontBit(sizeof(T) * CHAR_BIT - prefix.GetTrailBit());
        }
    } else
    {
        IndicateNoPrefixFound();
    }
}


template<typename T>
std::pair<bool, CompressionValue<sizeof(T)>>
ByteVariable<T>::EvaluateCommonPrefix(std::vector<CompressionValue<sizeof(T)>>& prefixes, const std::vector<int>& element_ids)
{
    cmc_assert(element_ids.size() >= 1);

    /* Initialize the first prefix with the first value */
    CompressionValue<sizeof(T)> prefix = prefixes[element_ids.front()];

    if (prefix.IsEmpty()) {return  std::make_pair(false, CompressionValue<sizeof(T)>());}

    /* Iterate over all leftover values */
    for (auto elem_id_iter = ++(element_ids.begin()); elem_id_iter != element_ids.end(); ++elem_id_iter)
    {
        /* Find a common prefix of all variables */
        prefix = GetCommonPrefix(prefix, prefixes[*elem_id_iter]);

        /* Check if there is a prefix */
        if (prefix.IsEmpty()) {return  std::make_pair(false, CompressionValue<sizeof(T)>());}
    }

    /* If the function arrives here, we do have found a common prefix of all values */
    return std::make_pair(true, prefix);
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
inline nc::Variable
ByteVariable<T>::WriteCompressedData(const int id, const int time_step, const SuffixEncoding encoding_scheme) const
{
    switch (encoding_scheme)
    {
        case SuffixEncoding::Plain:
            //return WriteCompressedData(id, time_step, EncodeSuffixesXORSchemeNew2);
            return WriteCompressedData(id, time_step, EncodePlainSuffixes);
        break;
        case SuffixEncoding::LengthEncoding:
            return WriteCompressedData(id, time_step, EncodeSuffixLengths);
        break;
        case SuffixEncoding::ArithmeticLengthEncoding:
            return WriteCompressedData(id, time_step, EncodeSuffixLengthsArithmeticEncoder);
        break;
        case SuffixEncoding::EncodeFirstOne:
            return WriteCompressedData(id, time_step, EncodeFirstOneInSuffixes);
        break;
        default:
            cmc_err_msg("An unknown suffix encoding scheme has been supplied.");
            return nc::Variable();
    }
}

template<typename T>
nc::Variable
ByteVariable<T>::WriteCompressedData(const int id, const int time_step, SuffixEncodingFunc<sizeof(T)> suffix_encoder) const
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
    buffered_data.emplace_back(suffix_encoder(byte_values_));

    cmc_debug_msg("Num bytes for suffixes: ", buffered_data.back().size());
    num_bytes += buffered_data.back().size();
    cmc_debug_msg("Overall bytes for this variable: ", num_bytes);

    //TODO: Make mechanism for parallel output
    //The complete global byte size needs to gathered and the level data has to be filled in parallel to it

    /* Generate a (potentially new) context information for the variable */
    const int global_context_info = (attributes_.GetGlobalContextInformation() != kNoGlobalContext ? attributes_.GetGlobalContextInformation() : 0);

    /* Create a netCDF variable to put out */
    std::string var_name = GetName() + "_" + std::to_string(global_context_info) + "_" + std::to_string(time_step);
    cmc_debug_msg("WriteComrpessed: Var_name: ", var_name);
    nc::SpecificVariable<uint8_t> compressed_variable{var_name, id};
    compressed_variable.Reserve(num_bytes);

    /* Put the buffered data into the variable to put out */
    for (auto lvl_data_iter = buffered_data.begin(); lvl_data_iter != buffered_data.end(); ++lvl_data_iter)
    {
        compressed_variable.PushBack(*lvl_data_iter);
    }

    /* Assign some attributes to it */
    std::vector<nc::Attribute> attributes;
    attributes.emplace_back("id", id);
    attributes.emplace_back("time_step", time_step);
    attributes.emplace_back("initial_refinement_level", mesh_.GetInitialRefinementLevel());
    attributes.emplace_back("initial_layout", static_cast<int>(attributes_.GetInitialDataLayout()));
    attributes.emplace_back("pre_compression_layout", static_cast<int>(attributes_.GetPreCompressionLayout()));
    attributes.emplace_back("global_context", global_context_info);
    attributes.emplace_back("data_type", static_cast<int>(ConvertToCmcType<T>()));
    attributes.emplace_back("missing_value", attributes_.GetMissingValue());
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


    return nc::Variable(std::move(compressed_variable), std::move(attributes));
}

template<typename T>
void
ByteVariable<T>::WriteDataToFile(const std::string& file_name) const
{
    FILE* file = fopen(file_name.c_str(), "wb");
    std::vector<uint8_t> buffer;
    buffer.reserve(byte_values_.size() * sizeof(T));

    for (auto val_iter = byte_values_.begin(); val_iter != byte_values_.end(); ++val_iter)
    {
        std::copy_n(val_iter->GetMemoryForReading().data(), sizeof(T), std::back_inserter(buffer));
    }

    fwrite(buffer.data(), sizeof(uint8_t), buffer.size(), file);
    fclose(file);
}

template<typename T>
nc::Variable
ByteVariable<T>::WriteCompressedDiffData(const int id, const int time_step) const
{
    std::vector<arithmetic_encoding::Letter> alphabet;
    alphabet.reserve(alphabet_.size());

    for (auto al_iter = alphabet_.begin(); al_iter != alphabet_.end(); ++al_iter)
    {
        alphabet.emplace_back(arithmetic_encoding::Letter{al_iter->first, al_iter->second});
    }

    arithmetic_encoding::StaticFrequencyModel frequency_model(alphabet);
    arithmetic_encoding::Encoder arm_encoder(std::make_unique<cmc::arithmetic_encoding::StaticFrequencyModel>(frequency_model));

    /* Levelwise encoded residuals */
    std::vector<bit_vector::BitVector> encoded_residuals;
    encoded_residuals.reserve(prefixes_.size());

    /* Levelwise encoded leading zero count of residuals */
    std::vector<bit_map::BitMap> encoded_residual_lzc;
    encoded_residual_lzc.reserve(prefixes_.size());

    /* We store the root element's value explicitly and afterwards start the encoding of the residuals */
    const CompressionValue<sizeof(T)> root_elem_value = prefixes_.back().prefixes.front();

    size_t residual_index = prefixes_.size() - 1;

    /* Encode all residuals on the succeeding refinement levels */
    for (auto lw_prefix_iter = std::next(prefixes_.rbegin()); lw_prefix_iter != prefixes_.rend(); ++lw_prefix_iter, --residual_index)
    {
        cmc_debug_msg("Residual index = ", residual_index, "; Size of prefixes: ", lw_prefix_iter->prefixes.size());

        /* Emplace a new Bit_Vector */
        encoded_residuals.emplace_back();
        encoded_residuals.back().Reserve(2 * lw_prefix_iter->prefixes.size());

        cmc_debug_msg("Refinement level: ", prefixes_.size() - residual_index);

        /* Get a view on the interpolation indications, indicating whether the residuals has to be added or subtracted */
        bit_map::BitMapView residual_sign_indications(interpolation_indications_[residual_index]);

        /* Iterate over all residuals and encode them */
        for (auto iter = lw_prefix_iter->prefixes.begin(); iter != lw_prefix_iter->prefixes.end(); ++iter)
        {
            /* Get the current residual */
            CompressionValue<sizeof(T)> val = *iter;

            /* Get the encoded LZC */
            const uint32_t signum = (residual_sign_indications.GetNextBit() == true ? 0 : arithmetic_encoding::kMSBBit);
            const uint32_t first_one_bit = val.GetNumberLeadingZeros();

            /* Encode the LZC */
            arm_encoder.EncodeSymbol(signum + first_one_bit);

            /* The first one is implicitly given by the leading zero count; therefore we set the "front bit" in order to discard the leading zeros and the following one */
            val.SetFrontBit(first_one_bit + 1);
            //val.SetFrontBit(first_one_bit);

            /* If there are remaining bits in the residual, append them to the encoded residuals of this level */
            if (not val.IsEmpty())
            {
                EncodeAndAppendPrefix(encoded_residuals.back(), val.GetSignificantBitsInBigEndianOrdering(), val.GetCountOfSignificantBits());
            }
        }

        /* Once we have finished the encoding of this level, we get the current LZC encodings and emplace a new one for the next level */
        arm_encoder.FinishEncoding();
        encoded_residual_lzc.push_back(arm_encoder.GetEncodedBitStream());

        arm_encoder.ClearBitStream();
    }

    if constexpr (std::is_same_v<T, float>)
    {
        cmc_debug_msg("The root elem value is: ", root_elem_value.template ReinterpretDataAs<float>());
    }
    /* Emplace a new Bit_Vector */
    encoded_residuals.emplace_back();
    encoded_residuals.back().Reserve(2 * byte_values_.size());

    /* Get a view on the indication for the residual summation/subtraction */
    cmc_assert(residual_index == 0);
    bit_map::BitMapView residual_sign_indications(interpolation_indications_[residual_index]);

    /* Iterate over all suffixes and encode them */
    for (auto val_iter = byte_values_.begin(); val_iter != byte_values_.end(); ++val_iter)
    {
        /* Get the current residual */
        CompressionValue<sizeof(T)> val = *val_iter;

        /* Get the encoded LZC */
        const uint32_t signum = (residual_sign_indications.GetNextBit() == true ? 0 : arithmetic_encoding::kMSBBit);
        const uint32_t first_one_bit = val.GetNumberLeadingZeros();

        /* Encode the LZC */
        arm_encoder.EncodeSymbol(signum + first_one_bit);

        /* The first one is implicitly given by the leading zero count; therefore we set the "front bit" in order to discard the leading zeros and the following one */
        //val.SetFrontBit(first_one_bit + 1);
        val.SetFrontBit(first_one_bit + 1);
        /* If there are remaining bits in the residual, append them to the encoded residuals of this level */
        if (not val.IsEmpty())
        {
            EncodeAndAppendPrefix(encoded_residuals.back(), val.GetSignificantBitsInBigEndianOrdering(), val.GetCountOfSignificantBits());
        }
    }

    /* Once we have finished the encoding of this level, we get the current LZC encodings and emplace a new one for the next level */
    arm_encoder.FinishEncoding();
    encoded_residual_lzc.push_back(arm_encoder.GetEncodedBitStream());
    arm_encoder.ClearBitStream();

    cmc_assert(encoded_residuals.size() == encoded_residual_lzc.size());

    /* Cout all bytes together and create a single stream of the data */
    size_t cumulative_byte_count = 0;
    for (size_t idx = 0; idx < encoded_residuals.size(); ++idx)
    {
        cumulative_byte_count += encoded_residuals[idx].size();
        cumulative_byte_count += encoded_residual_lzc[idx].size_bytes();
    }

    /* Encode the alphabet */
    bit_vector::BitVector encoded_arm_alphabet = arm_encoder.EncodeAlphabet();
    cumulative_byte_count += encoded_arm_alphabet.size();

    /* Additionally to the encoded residuals, we store some "size_t"'s indicating the lenghts of the several streams */
    cumulative_byte_count += sizeof(size_t) * (2 * encoded_residuals.size() + 1) + sizeof(T);

    cmc_debug_msg("The overall byte count for the variable amounts to: ", cumulative_byte_count, " bytes.");

    /* We allocate memory for the final output */
    std::vector<uint8_t> encoded_byte_stream;
    encoded_byte_stream.reserve(cumulative_byte_count);

    cmc_debug_msg("After allocation: BitVector size: ", encoded_byte_stream.size());

    /* 1) First, we store the overall bytes */
    std::array<uint8_t, sizeof(size_t)> serialized_value = SerializeValue(cumulative_byte_count, Endian::Big);
    std::copy_n(serialized_value.begin(), serialized_value.size(), std::back_inserter(encoded_byte_stream));
    //encoded_variable_data.AppendBytes<sizeof(size_t)>(serialized_value);
    cmc_debug_msg("After total byte count: BitVector size: ", encoded_byte_stream.size());

    /* 2) Afterwards, we append the encoded alphabet */
    //encoded_variable_data.AppendBits(encoded_arm_alphabet);
    std::copy_n(encoded_arm_alphabet.GetData().begin(), encoded_arm_alphabet.GetData().size(), std::back_inserter(encoded_byte_stream));
    cmc_debug_msg("After encoded alphabet: BitVector size: ", encoded_byte_stream.size());

    /* 3) Next, we store the root element's value */
    //cmc_debug_msg("Root Elem value has num significant bits: ", root_elem_value.GetCountOfSignificantBits());
    const T root_value = root_elem_value.template ReinterpretDataAs<T>();
    std::array<uint8_t, sizeof(T)> serialized_root_value = SerializeValue(root_value, Endian::Big);
    std::copy_n(serialized_root_value.begin(), serialized_root_value.size(), std::back_inserter(encoded_byte_stream));

    //encoded_variable_data.AppendBits(root_elem_value.GetSignificantBitsInBigEndianOrdering());
    //std::copy_n(root_elem_value.GetSignificantBitsInBigEndianOrdering().begin(), sizeof(T), std::back_inserter(encoded_byte_stream));

    cmc_debug_msg("After root value: BitVector size: ", encoded_byte_stream.size());
    //encoded_variable_data.AddPaddingToFullByte();
    /* 4) Now, the residual levels are stored/appended in order from coarse to fine, preceeded by
     * the number of bytes for the encoded LZC and the remaining residual bits */
    for (size_t idx = 0; idx < encoded_residuals.size(); ++idx)
    {
        cmc_debug_msg("Copy of LZC count of level starts at: ", encoded_byte_stream.size());
        /* a) We start with the leading zero bytes count */
        serialized_value = SerializeValue(encoded_residual_lzc[idx].size_bytes(), Endian::Big);
        //encoded_variable_data.AppendBytes<sizeof(size_t)>(serialized_value);
        std::copy_n(serialized_value.begin(), serialized_value.size(), std::back_inserter(encoded_byte_stream));

        cmc_debug_msg("Copy of residual count of level starts at: ", encoded_byte_stream.size());
        /* b) Next follows the remaining residual bytes count */
        serialized_value = SerializeValue(encoded_residuals[idx].size(), Endian::Big);
        //encoded_variable_data.AppendBytes<sizeof(size_t)>(serialized_value);
        std::copy_n(serialized_value.begin(), serialized_value.size(), std::back_inserter(encoded_byte_stream));

        /* Now we copy the actual encoded leading zero bytes */
        //encoded_variable_data.AppendBytes(encoded_residual_lzc[idx].GetByteData());
        std::copy_n(encoded_residual_lzc[idx].GetByteData().begin(), encoded_residual_lzc[idx].GetByteData().size(), std::back_inserter(encoded_byte_stream));

        //encoded_variable_data.AddPaddingToFullByte();

        /* And afterwards the remaining residual bytes */
        //encoded_variable_data.AppendBits(encoded_residuals[idx]);
        std::copy_n(encoded_residuals[idx].begin(), encoded_residuals[idx].size(), std::back_inserter(encoded_byte_stream));

        //encoded_variable_data.AddPaddingToFullByte();
        cmc_debug_msg("After lvl: ", idx, ", BitVector size: ", encoded_byte_stream.size());
    }

    cmc_debug_msg("After everything has been added, size of encoded variable byte stream is: ", encoded_byte_stream.size());


    //TODO: Make mechanism for parallel output
    //The complete global byte size needs to gathered and the level data has to be filled in parallel to it

    /* Generate a (potentially new) context information for the variable */
    const int global_context_info = (attributes_.GetGlobalContextInformation() != kNoGlobalContext ? attributes_.GetGlobalContextInformation() : 0);

    /* Create a netCDF variable to put out */
    std::string var_name = GetName() + "_" + std::to_string(global_context_info) + "_" + std::to_string(time_step);
    cmc_debug_msg("WriteComrpessed: Var_name: ", var_name);

    nc::SpecificVariable<uint8_t> compressed_variable{var_name, id};
    compressed_variable.Reserve(encoded_byte_stream.size());

    //TODO Move instead of copy encoded byte stream
    /* Put the buffered data into the variable to put out */
    compressed_variable.PushBack(encoded_byte_stream);

    /* Assign some attributes to it */
    std::vector<nc::Attribute> attributes;
    attributes.emplace_back("id", id);
    attributes.emplace_back("time_step", time_step);
    attributes.emplace_back("initial_refinement_level", mesh_.GetInitialRefinementLevel());
    attributes.emplace_back("initial_layout", static_cast<int>(attributes_.GetInitialDataLayout()));
    attributes.emplace_back("pre_compression_layout", static_cast<int>(attributes_.GetPreCompressionLayout()));
    attributes.emplace_back("global_context", global_context_info);
    attributes.emplace_back("data_type", static_cast<int>(ConvertToCmcType<T>()));
    attributes.emplace_back("missing_value", attributes_.GetMissingValue());
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


    return nc::Variable(std::move(compressed_variable), std::move(attributes));
}



template<typename T>
nc::Variable
ByteVariable<T>::_WriteCompressedDiffDataPCPLightAMRComparison(const int id, const int time_step) const
{
    #if _PERFORM_TEST_COMPARISON
    size_t cumulative_byte_count = 0;

    for (auto lvl_iter = test_encodings_.rbegin(); lvl_iter != test_encodings_.rend(); ++lvl_iter)
    {
        cmc_debug_msg("Level Sizes of LZC: ",lvl_iter->encoded_lzc.size(), " and residual bits: ", lvl_iter->encoded_residuals.size());
        cumulative_byte_count += lvl_iter->encoded_lzc.size();
        cumulative_byte_count += lvl_iter->encoded_residuals.size();
    }

    /* Additionally to the encoded residuals, we store some "size_t"'s indicating the lenghts of the several streams */
    cumulative_byte_count += sizeof(size_t) * (2 * test_encodings_.size() + 1) + sizeof(T);

    cmc_debug_msg("The overall byte count for the variable amounts to: ", cumulative_byte_count, " bytes.");

    /* We allocate memory for the final output */
    std::vector<uint8_t> encoded_byte_stream;
    encoded_byte_stream.reserve(cumulative_byte_count);

    cmc_debug_msg("After allocation: BitVector size: ", encoded_byte_stream.size());

    /* 1) First, we store the overall bytes */
    std::array<uint8_t, sizeof(size_t)> serialized_value = SerializeValue(cumulative_byte_count, Endian::Big);
    std::copy_n(serialized_value.begin(), serialized_value.size(), std::back_inserter(encoded_byte_stream));
    //encoded_variable_data.AppendBytes<sizeof(size_t)>(serialized_value);
    cmc_debug_msg("After total byte count: BitVector size: ", encoded_byte_stream.size());

    /* 2) Next, we store the root element's value */
    /* We store the root element's value explicitly and afterwards start the encoding of the residuals */
    const CompressionValue<sizeof(T)> root_elem_value = prefixes_.back().prefixes.front();
    const T root_value = root_elem_value.template ReinterpretDataAs<T>();
    std::array<uint8_t, sizeof(T)> serialized_root_value = SerializeValue(root_value, Endian::Big);
    std::copy_n(serialized_root_value.begin(), serialized_root_value.size(), std::back_inserter(encoded_byte_stream));

    cmc_debug_msg("After root value: BitVector size: ", encoded_byte_stream.size());

    /* 3) Now, the residual levels are stored/appended in order from coarse to fine, preceeded by
     * the number of bytes for the encoded LZC and the remaining residual bits */
    for (auto lvl_iter = test_encodings_.rbegin(); lvl_iter != test_encodings_.rend(); ++lvl_iter)
    {
        cmc_debug_msg("Copy of LZC count of level starts at: ", encoded_byte_stream.size());
        /* a) We start with the leading zero bytes count */
        serialized_value = SerializeValue(lvl_iter->encoded_lzc.size(), Endian::Big);

        std::copy_n(serialized_value.begin(), serialized_value.size(), std::back_inserter(encoded_byte_stream));

        cmc_debug_msg("Copy of residual count of level starts at: ", encoded_byte_stream.size());

        /* b) Next follows the remaining residual bytes count */
        serialized_value = SerializeValue(lvl_iter->encoded_residuals.size(), Endian::Big);

        std::copy_n(serialized_value.begin(), serialized_value.size(), std::back_inserter(encoded_byte_stream));

        /* Now we copy the actual encoded leading zero bytes */
        std::copy_n(lvl_iter->encoded_lzc.begin(), lvl_iter->encoded_lzc.size(), std::back_inserter(encoded_byte_stream));

        /* And afterwards the remaining residual bytes */
        std::copy_n(lvl_iter->encoded_residuals.begin(), lvl_iter->encoded_residuals.size(), std::back_inserter(encoded_byte_stream));

        cmc_debug_msg("After next lvl BitVector size: ", encoded_byte_stream.size());
    }

    cmc_debug_msg("After everything has been added, size of encoded variable byte stream is: ", encoded_byte_stream.size());


    //TODO: Make mechanism for parallel output
    //The complete global byte size needs to gathered and the level data has to be filled in parallel to it

    /* Generate a (potentially new) context information for the variable */
    const int global_context_info = (attributes_.GetGlobalContextInformation() != kNoGlobalContext ? attributes_.GetGlobalContextInformation() : 0);

    /* Create a netCDF variable to put out */
    std::string var_name = GetName() + "_" + std::to_string(global_context_info) + "_" + std::to_string(time_step);
    cmc_debug_msg("WriteComrpessed: Var_name: ", var_name);

    nc::SpecificVariable<uint8_t> compressed_variable{var_name, id};
    //compressed_variable.Reserve(encoded_variable_data.size());
    compressed_variable.Reserve(encoded_byte_stream.size());

    //TODO Move instead of copy encoded byte stream
    /* Put the buffered data into the variable to put out */
    compressed_variable.PushBack(encoded_byte_stream);

    /* Assign some attributes to it */
    std::vector<nc::Attribute> attributes;
    attributes.emplace_back("id", id);
    attributes.emplace_back("time_step", time_step);
    attributes.emplace_back("initial_refinement_level", mesh_.GetInitialRefinementLevel());
    attributes.emplace_back("initial_layout", static_cast<int>(attributes_.GetInitialDataLayout()));
    attributes.emplace_back("pre_compression_layout", static_cast<int>(attributes_.GetPreCompressionLayout()));
    attributes.emplace_back("global_context", global_context_info);
    attributes.emplace_back("data_type", static_cast<int>(ConvertToCmcType<T>()));
    attributes.emplace_back("missing_value", attributes_.GetMissingValue());
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


    return nc::Variable(std::move(compressed_variable), std::move(attributes));
    #else
    return nc::Variable();
    #endif
}



template<typename T>
void
ByteVariable<T>::InitializeDiffDecompression()
{
    byte_values_new_.reserve(byte_values_.size() * std::pow(2, GetGlobalDomain().GetDimensionality()));
}

template<typename T>
void
ByteVariable<T>::SetDiffDecompressionRootElementValue(const CmcUniversalType& root_value)
{
    cmc_assert(std::holds_alternative<T>(root_value));

    byte_values_.clear();
    /* Get the root element value and assign it to the vector */
    const CompressionValue<sizeof(T)> root_elem_value = CompressionValue<sizeof(T)>(std::get<T>(root_value));
    byte_values_.push_back(root_elem_value);
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
ByteVariable<T>::KeepInitialData(const bool keep_data)
{
    is_initial_data_kept_ = keep_data;

    if (!is_initial_data_kept_)
    {
        /* Try to release the memory of the stored intial data */
        initial_data_.clear();
        initial_data_.shrink_to_fit();
    }
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
