#ifndef CMC_T8_ADAPT_HXX
#define CMC_T8_ADAPT_HXX
/**
 * @file cmc_t8_adapt.hxx
 */

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_compression_settings.hxx"
#include "t8code/cmc_t8_data_variables.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_adapt_track_inaccuracy.hxx"
#include "t8code/cmc_t8_byte_variable.hxx"
#include "utilities/cmc_prefix_encoding.hxx"
#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_diff_decoder.hxx"
#include "comparison/cmc_test_comparison_pcp_decoder.hxx"

#include <vector>
#include <memory>

namespace cmc {

class AdaptData
{
public:
    AdaptData() = delete;
    AdaptData(const CompressionSettings& settings, std::vector<Var>& compression_variables, const CoarseningSample& adaptation_sample, const CompressionMode mode)
    : compression_settings_{settings}, variables_{compression_variables}, corresponding_variable_id_{adaptation_sample.corresponding_variable_id}, mode_{mode}{};
    ~AdaptData(){};
    
    AdaptData(const AdaptData& other) = default;
    AdaptData& operator=(const AdaptData& other) = delete;
    AdaptData(AdaptData&& other) = default;
    AdaptData& operator=(AdaptData&& other) = delete;

    bool IsCompressionProgressing() const;
    t8_forest_t GetCurrentMesh() const;
    int GetAdaptationStepCount() const;
    int GetInitialRefinementLevelOfMesh() const;
    void SetCurrentMesh(t8_forest_t forest);
    t8_forest_t RepartitionData(t8_forest_t adapted_forest);
    t8_forest_adapt_t GetAdaptationFunction() const;
    
    Var& GetCurrentCompressionVariable();
    const Var& GetCurrentCompressionVariable() const;
    void InitializeCompressionIteration();

    void FinalizeCompressionIteration();

    void UpdateCompressionData();

    void IndicateCoarsening();
    void IndicateElementStaysUnchanged();
    [[nodiscard]] std::vector<bit_map::BitMap> TransferIndicationBits();

    void IndicateToKeepInitialData() {do_variables_keep_initial_data_ = true;};
    void IndicateToStoreInitialMesh() {GetCurrentCompressionVariable().StoreCurrentMeshAsInitialMesh();};
    
private:
    const CompressionSettings& compression_settings_;
    std::vector<Var>& variables_;
    std::vector<bit_map::BitMap> refinement_indications_;
    t8_gloidx_t previous_number_of_elements_{-1};
    t8_gloidx_t new_number_of_elements_{0};
    int corresponding_variable_id_;
    CompressionMode mode_;
    int count_adaptation_step_{0};
    bool do_variables_keep_initial_data_{false};
};


class PrefixAdaptData
{
public:
    PrefixAdaptData() = delete;
    PrefixAdaptData(ByteVar& variable)
    : byte_variable_{variable}{
        new_number_of_elements_ = t8_forest_get_global_num_elements(variable.GetAmrMesh().GetMesh());
    };

    /* The compression runs until the (single) root element of the mesh is reached */
    bool IsCompressionProgressing()
    {
        return (new_number_of_elements_ > 1 ? true : false);
    }

    /* Allocate memory and prepare the variable for the prefix extraction */
    void InitializeCompressionIteration()
    {
        byte_variable_.InitializeCompressionIteration();
        if (count_adaptation_step_ == 0)
        {
            //byte_variable_.WriteDataToVTK_(11100);
            //byte_variable_.PrintCompressionValues();
        }
    }

    void FinalizeCompressionIteration()
    {
        new_number_of_elements_ = t8_forest_get_global_num_elements(GetCurrentMesh());
        byte_variable_.FinalizeCompressionIteration();
        ++count_adaptation_step_;

        //byte_variable_.WriteDataToVTK_(11100 + count_adaptation_step_);
    }

    /* Indicate that a single element has been passed to the adaptation function, which is why a common
     * prefix could not be evaluated/extracted. Therefore, the element/prefix stays unchanged. */
    void LeaveCoarsePrefixUnchanged(const int elem_index)
    {
        byte_variable_.LeaveCoarsePrefixUnchanged(elem_index);
    }

    /* This function extracts a common prefix of the family of elements passed to the adaptation function */
    void ExtractCommonPrefix(const int start_index, const int num_elements)
    {
        /* Extract a common prefix of the family of elements indicated by the start_index and the element count */
        if (count_adaptation_step_ == 0)
        {
            return byte_variable_.ExtractCommonPrefixFromInitialData(start_index, num_elements);
        } else
        {
            return byte_variable_.ExtractCommonPrefixFromPreviousPrefixes(start_index, num_elements);
        }
    }

    /* This function extracts a common prefix of the family of elements passed to the adaptation function */
    void ExtractCommonPrefix(const std::vector<int>& elem_ids)
    {
        /* Extract a common prefix of the family of elements indicated by the start_index and the element count */
        if (count_adaptation_step_ == 0)
        {
            return byte_variable_.ExtractCommonPrefixFromInitialData(elem_ids);
        } else
        {
            return byte_variable_.ExtractCommonPrefixFromPreviousPrefixes(elem_ids);
        }
    }

    /* Get the mesh the data/prefixes of the variable are currently defined on */
    t8_forest_t GetCurrentMesh() const
    {
        return byte_variable_.GetMesh();
    }

    /* Update the mesh of the variable (to be called after an adaptation step) */
    void SetCurrentMesh(t8_forest_t mesh)
    {
        byte_variable_.SetMesh(mesh);
    }

    /* Receive the initial refinement level of the data */
    int GetInitialRefinementLevelOfMesh() const
    {
        return byte_variable_.GetAmrMesh().GetInitialRefinementLevel();
    }

    const GeoDomain& GetGlobalDomain() const
    {
        return byte_variable_.GetGlobalDomain();
    }

    DataLayout GetInitialDataLayout() const
    {
        return byte_variable_.GetInitialDataLayout();
    }

    /* Get the number of iterations (of the prefix extraction) that have been performed */
    int GetAdaptationStepCount() const
    {
        return count_adaptation_step_;
    }

private:
    ByteVar& byte_variable_;
    t8_gloidx_t new_number_of_elements_{2};
    int count_adaptation_step_{0};
};

class FittingPyramidAdaptData
{
public:

};


class DiffAdaptData
{
public:
    DiffAdaptData() = delete;
    DiffAdaptData(ByteVar& variable)
    : byte_variable_{variable}{
        new_number_of_elements_ = t8_forest_get_global_num_elements(variable.GetAmrMesh().GetMesh());
    };

    /* The compression runs until the (single) root element of the mesh is reached */
    bool IsCompressionProgressing()
    {
        return (new_number_of_elements_ > 1 ? true : false);
    }

    /* Allocate memory and prepare the variable for the prefix extraction */
    void InitializeCompressionIteration()
    {
        byte_variable_.InitializeCompressionIteration();
        if (count_adaptation_step_ == 0)
        {
            //byte_variable_.WriteDataToVTK_(11100);
            //byte_variable_.PrintCompressionValues();
        }
    }

    void FinalizeCompressionIteration()
    {
        new_number_of_elements_ = t8_forest_get_global_num_elements(GetCurrentMesh());
        byte_variable_.FinalizeCompressionIteration();
        ++count_adaptation_step_;

        //byte_variable_.WriteDataToVTK_(11100 + count_adaptation_step_);
    }


    void ExtractMeanAndLeaveDifferences(const int elem_index, const int num_elements)
    {
        /* Extract an interpolated value from the family of elements indicated by the start_index and the element count
         * and leave the differences behind */
        if (count_adaptation_step_ > 0)
        {
            byte_variable_.ExtractMeanAndLeaveDifferencesFromPreviousMeans(elem_index, num_elements);
        } else
        {
            byte_variable_.ExtractMeanAndLeaveDifferencesFromInitialData(elem_index, num_elements);
        }
    }

    void _TestComparisonPCPLightAMR(const int elem_index, const int num_elements)
    {
        /* Extract an interpolated value from the family of elements indicated by the start_index and the element count
         * and leave the differences behind */
        if (count_adaptation_step_ > 0)
        {
            byte_variable_._TestComparisonPCPLightAMRFromPreviousMeans(elem_index, num_elements);
        } else
        {
            byte_variable_._TestComparisonPCPLightAMRFromInitialData(elem_index, num_elements);
        }
    }

    void LeaveValueUnchanged(const int elem_index)
    {
        byte_variable_.LeaveValueUnchangedForNextMeanComputation(elem_index);
    }

    /* Get the mesh the data/prefixes of the variable are currently defined on */
    t8_forest_t GetCurrentMesh() const
    {
        return byte_variable_.GetMesh();
    }

    /* Update the mesh of the variable (to be called after an adaptation step) */
    void SetCurrentMesh(t8_forest_t mesh)
    {
        byte_variable_.SetMesh(mesh);
    }

    #if 0
    /* Receive the initial refinement level of the data */
    int GetInitialRefinementLevelOfMesh() const
    {
        return byte_variable_.GetAmrMesh().GetInitialRefinementLevel();
    }

    const GeoDomain& GetGlobalDomain() const
    {
        return byte_variable_.GetGlobalDomain();
    }

    DataLayout GetInitialDataLayout() const
    {
        return byte_variable_.GetInitialDataLayout();
    }
    #endif
    /* Get the number of iterations (of the prefix extraction) that have been performed */
    int GetAdaptationStepCount() const
    {
        return count_adaptation_step_;
    }

private:
    ByteVar& byte_variable_;
    t8_gloidx_t new_number_of_elements_{2};
    int count_adaptation_step_{0};
};




class MeanAdaptData
{
public:
    MeanAdaptData() = delete;
    MeanAdaptData(ByteVar& variable)
    : byte_variable_{variable}{
        new_number_of_elements_ = t8_forest_get_global_num_elements(variable.GetAmrMesh().GetMesh());
    };

    /* The compression runs until the (single) root element of the mesh is reached */
    bool IsCompressionProgressing()
    {
        return (new_number_of_elements_ > 1 ? true : false);
    }

    /* Allocate memory and prepare the variable for the prefix extraction */
    void InitializeCompressionIteration()
    {
        byte_variable_.InitializeCompressionIteration();
        if (count_adaptation_step_ == 0)
        {
            //byte_variable_.WriteDataToVTK_(11100);
            //byte_variable_.PrintCompressionValues();
        }
    }

    void FinalizeCompressionIteration()
    {
        new_number_of_elements_ = t8_forest_get_global_num_elements(GetCurrentMesh());
        byte_variable_.FinalizeCompressionIteration();
        ++count_adaptation_step_;

        //byte_variable_.WriteDataToVTK_(11100 + count_adaptation_step_);
    }


    void ExtractMean(const int elem_index, const int num_elements)
    {
        /* Extract an interpolated value from the family of elements indicated by the start_index and the element count
         * and leave the differences behind */
        //cmc_debug_msg("In Ectract Mean and leave differences: adaptation count: ", count_adaptation_step_);
        if (count_adaptation_step_ > 0)
        {
            //cmc_debug_msg("Why is this not called?????");
            byte_variable_.ExtractMeanFromPreviousMeans(elem_index, num_elements);
        } else
        {
            //cmc_debug_msg("\n\n\n Wird der FAll hier aufgereufen");
            byte_variable_.ExtractMeanFromInitialData(elem_index, num_elements);
        }
    }

    void LeaveValueUnchanged(const int elem_index)
    {
        byte_variable_.LeaveValueUnchangedForNextMeanComputation(elem_index);
    }

    
    /* Get the mesh the data/prefixes of the variable are currently defined on */
    t8_forest_t GetCurrentMesh() const
    {
        return byte_variable_.GetMesh();
    }

    /* Update the mesh of the variable (to be called after an adaptation step) */
    void SetCurrentMesh(t8_forest_t mesh)
    {
        byte_variable_.SetMesh(mesh);
    }

    /* Receive the initial refinement level of the data */
    int GetInitialRefinementLevelOfMesh() const
    {
        return byte_variable_.GetAmrMesh().GetInitialRefinementLevel();
    }

    const GeoDomain& GetGlobalDomain() const
    {
        return byte_variable_.GetGlobalDomain();
    }

    DataLayout GetInitialDataLayout() const
    {
        return byte_variable_.GetInitialDataLayout();
    }

    /* Get the number of iterations (of the prefix extraction) that have been performed */
    int GetAdaptationStepCount() const
    {
        return count_adaptation_step_;
    }

private:
    ByteVar& byte_variable_;
    t8_gloidx_t new_number_of_elements_{2};
    int count_adaptation_step_{0};
};


class DecompressPrefixAdaptData
{
public:
    DecompressPrefixAdaptData() = delete;
    DecompressPrefixAdaptData(ByteVar& variable, const std::vector<uint8_t>& compressed_byte_stream)
    : byte_variable_{variable}, byte_stream_{compressed_byte_stream} {
        decoder_ = std::make_unique<PrefixDecoder>(byte_stream_);
    };
    DecompressPrefixAdaptData(ByteVar& variable, std::vector<uint8_t>&& compressed_byte_stream)
    : byte_variable_{variable}, byte_stream_{std::move(compressed_byte_stream)} {
        decoder_ = std::make_unique<PrefixDecoder>(byte_stream_);
    };
    ~DecompressPrefixAdaptData() = default;

    /* Get the mesh the data/prefixes of the variable are currently defined on */
    t8_forest_t GetCurrentMesh() const
    {
        return byte_variable_.GetMesh();
    }

    /* Update the mesh of the variable (to be called after an adaptation step) */
    void SetCurrentMesh(t8_forest_t mesh)
    {
        byte_variable_.SetMesh(mesh);
    }

    void InitializePlainSuffixDecompression()
    {
        byte_variable_.InitializeSuffixDecompression();
        decoder_->MoveToPlainEncodedSuffixLevel();
    }

    void FinalizePlainSuffixDecompression()
    {
        byte_variable_.FinalizeDecompressionIteration();
        ++count_decompression_step;

        //byte_variable_.WriteDataToVTK_(count_decompression_step);
        //byte_variable_.PrintCompressionValues();
    }

    void InitializeSuffixDecompression()
    {
        byte_variable_.InitializeSuffixDecompression();
        decoder_->MoveToSuffixLevel();
    }

    void FinalizeSuffixDecompression()
    {
        byte_variable_.FinalizeDecompressionIteration();
        ++count_decompression_step;

        //byte_variable_.WriteDataToVTK_(count_decompression_step);
        //byte_variable_.PrintCompressionValues();
    }

    /* Allocate memory and prepare the variable for the prefix extraction */
    void InitializeDecompressionIteration()
    {
        byte_variable_.InitializeDecompressionIteration();

        if (count_decompression_step > 0)
        {
            /* We are moving to the next level within the decoding byte stream.
             * During the first iteration, we do not need to call this function, because
             * the constructor of the PrefixDecoder sets it up in the first iteration */
            decoder_->MoveToNextLevel();
        }

        /* In parallel, we need to get to the process-local start within the decoding byte stream */
        const t8_gloidx_t global_elem_offset = t8_forest_get_first_local_element_id(GetCurrentMesh());
        decoder_->MoveToGlobalPositionWithinLevel(global_elem_offset);
    }

    void FinalizeDecompressionIteration()
    {
        byte_variable_.FinalizeDecompressionIteration();
        ++count_decompression_step;

        //byte_variable_.WriteDataToVTK_(count_decompression_step);
    }

    bool IsPrefixGiven()
    {
        return decoder_->GetNextPrefixIndicatorBit();
    }

    bool IsRefinementGiven()
    {
        return decoder_->GetNextRefinementIndicatorBit();
    }

    std::pair<std::vector<uint8_t>, int>
    GetNextPrefix()
    {
        return decoder_->GetNextPrefix();
    }

    std::vector<uint8_t>
    GetNextPlainSuffix(const size_t num_bits)
    {
        return decoder_->GetNextPlainSuffix(num_bits);
    }

    void Refine(const int elem_id)
    {
        byte_variable_.DecompressionRefine(elem_id);
    }

    void LeaveElementUnchanged(const int elem_id)
    {
        byte_variable_.DecompressionLeaveElementUnchanged(elem_id);
    }

    void ApplyPrefixAndRefine(const int elem_id, const std::vector<uint8_t> prefix, const int num_prefix_bits)
    {
        byte_variable_.RefineAndApplyCommonPrefix(elem_id, prefix, num_prefix_bits);
    }

    void ApplyPrefixAndLeaveElementUnchanged(const int elem_id, const std::vector<uint8_t> prefix, const int num_prefix_bits)
    {
        byte_variable_.LeaveElementUnchangedAndApplyPrefix(elem_id, prefix, num_prefix_bits);
    }

    size_t GetCountSignificantBits(const int elem_id) const
    {
        return byte_variable_.GetCountSignificantBits(elem_id);
    }
    size_t GetBitCountOfDataType() const
    {
        return CmcTypeToBytes(byte_variable_.GetType()) * CHAR_BIT;
    }
private: 
    ByteVar& byte_variable_;
    std::vector<uint8_t> byte_stream_;
    std::unique_ptr<PrefixDecoder> decoder_;
    int count_decompression_step{0};
};




class DecompressDiffAdaptData
{
public:
    DecompressDiffAdaptData() = delete;
    DecompressDiffAdaptData(ByteVar& variable, const std::vector<uint8_t>& compressed_byte_stream)
    : byte_variable_{variable}, byte_stream_{compressed_byte_stream} {
        decoder_ = std::make_unique<diff::DiffDecoder>(byte_stream_, byte_variable_.GetType());
    };
    DecompressDiffAdaptData(ByteVar& variable, std::vector<uint8_t>&& compressed_byte_stream)
    : byte_variable_{variable}, byte_stream_{std::move(compressed_byte_stream)} {
        decoder_ = std::make_unique<diff::DiffDecoder>(byte_stream_, byte_variable_.GetType());
    };
    ~DecompressDiffAdaptData() = default;

    /* Get the mesh the data/prefixes of the variable are currently defined on */
    t8_forest_t GetCurrentMesh() const
    {
        return byte_variable_.GetMesh();
    }

    int GetInitialRefinementLevelOfMesh() const
    {
        return byte_variable_.GetAmrMesh().GetInitialRefinementLevel();
    }

    DataLayout GetInitialDataLayout() const
    {
        return byte_variable_.GetInitialDataLayout();
    }

    ByteVar& GetCurrentCompressionVariable() {return byte_variable_;}
    const ByteVar& GetCurrentCompressionVariable() const {return byte_variable_;}

    /* Update the mesh of the variable (to be called after an adaptation step) */
    void SetCurrentMesh(t8_forest_t mesh)
    {
        byte_variable_.SetMesh(mesh);
    }

    void InitializeDiffDecompression()
    {
        cmc_debug_msg("In adapt init diff decompression");
        if (count_decompression_step == 0)
        {
            cmc_debug_msg("in if");
            byte_variable_.SetDiffDecompressionRootElementValue(decoder_->GetRootElementValue());
            cmc_debug_msg("Root elem value is set");
        }
        cmc_debug_msg("before init byte var diff decompression");
        byte_variable_.InitializeDiffDecompression();
        cmc_debug_msg("Before move to next refinement lvl");
        decoder_->MoveToNextRefinementLevel();
        cmc_debug_msg("at the end of adapt init diff decomp");
    }

    void FinalizeDiffDecompression()
    {
        byte_variable_.FinalizeDecompressionIteration();
        ++count_decompression_step;

        //byte_variable_.WriteDataToVTK_(count_decompression_step);
        //byte_variable_.PrintCompressionValues();
    }

    size_t GetBitCountOfDataType() const
    {
        return CmcTypeToBytes(byte_variable_.GetType()) * CHAR_BIT;
    }

    uint32_t GetNextEncodedResidualLength()
    {
        return decoder_->GetNextEncodedResidualLength();
    }

    std::vector<uint8_t> GetNextResidualBitSequence(const size_t num_bits)
    {
        return decoder_->GetNextResidualBitSequence(num_bits);
    }
    
    void StoreElementUnchanged(const int elem_id)
    {
        byte_variable_.StoreElementUnchanged(elem_id);
    }

    void ApplyResidualAndStoreElement(const int elem_id, const uint32_t encoded_lzc, const std::vector<uint8_t>& residual_bits)
    {
        byte_variable_.ApplyResidualAndStoreElement(elem_id, encoded_lzc, residual_bits);
    }

private: 
    ByteVar& byte_variable_;
    std::vector<uint8_t> byte_stream_;
    std::unique_ptr<diff::DiffDecoder> decoder_;
    int count_decompression_step{0};
};

#if _PERFORM_TEST_COMPARISON

class _TestLightAMRPCPDecompressData
{
public:
    _TestLightAMRPCPDecompressData() = delete;
    _TestLightAMRPCPDecompressData(ByteVar& variable, const std::vector<uint8_t>& compressed_byte_stream)
    : byte_variable_{variable}, byte_stream_{compressed_byte_stream} {
        decoder_ = std::make_unique<test_comparison::light_amr_pcp::_TestPCPDecoder>(byte_stream_, byte_variable_.GetType());
    };
    _TestLightAMRPCPDecompressData(ByteVar& variable, std::vector<uint8_t>&& compressed_byte_stream)
    : byte_variable_{variable}, byte_stream_{std::move(compressed_byte_stream)} {
        decoder_ = std::make_unique<test_comparison::light_amr_pcp::_TestPCPDecoder>(byte_stream_, byte_variable_.GetType());
    };
    ~_TestLightAMRPCPDecompressData() = default;

    /* Get the mesh the data/prefixes of the variable are currently defined on */
    t8_forest_t GetCurrentMesh() const
    {
        return byte_variable_.GetMesh();
    }

    int GetInitialRefinementLevelOfMesh() const
    {
        return byte_variable_.GetAmrMesh().GetInitialRefinementLevel();
    }

    DataLayout GetInitialDataLayout() const
    {
        return byte_variable_.GetInitialDataLayout();
    }

    ByteVar& GetCurrentCompressionVariable() {return byte_variable_;}
    const ByteVar& GetCurrentCompressionVariable() const {return byte_variable_;}

    /* Update the mesh of the variable (to be called after an adaptation step) */
    void SetCurrentMesh(t8_forest_t mesh)
    {
        byte_variable_.SetMesh(mesh);
    }

    void InitializeDiffDecompression()
    {
        cmc_debug_msg("In adapt init diff decompression");
        if (count_decompression_step == 0)
        {
            cmc_debug_msg("in if");
            byte_variable_.SetDiffDecompressionRootElementValue(decoder_->GetRootElementValue());
            cmc_debug_msg("Root elem value is set");
        }
        cmc_debug_msg("before init byte var diff decompression");
        byte_variable_.InitializeDiffDecompression();
        cmc_debug_msg("Before move to next refinement lvl");
        decoder_->MoveToNextRefinementLevel();
        cmc_debug_msg("at the end of adapt init diff decomp");
    }

    void FinalizeDiffDecompression()
    {
        byte_variable_.FinalizeDecompressionIteration();
        ++count_decompression_step;

        //byte_variable_.WriteDataToVTK_(count_decompression_step);
        //byte_variable_.PrintCompressionValues();
    }

    size_t GetBitCountOfDataType() const
    {
        return CmcTypeToBytes(byte_variable_.GetType()) * CHAR_BIT;
    }

    size_t GetNextEncodedLZC()
    {
        return decoder_->GetNextEncodedLZC();
    }

    std::vector<uint8_t> GetNextResidualBitSequence(const size_t num_bits)
    {
        return decoder_->GetNextResidualBitSequence(num_bits);
    }
    
    void StoreElementUnchanged(const int elem_id)
    {
        byte_variable_.StoreElementUnchanged(elem_id);
    }

    void ApplyResidualAndStoreElement(const int elem_id, const uint32_t encoded_lzc, const std::vector<uint8_t>& residual_bits)
    {
        byte_variable_.ApplyResidualAndStoreElement(elem_id, encoded_lzc, residual_bits);
    }

private: 
    ByteVar& byte_variable_;
    std::vector<uint8_t> byte_stream_;
    std::unique_ptr<test_comparison::light_amr_pcp::_TestPCPDecoder> decoder_;
    int count_decompression_step{0};
};


#endif

}

#endif /* !CMC_T8_ADAPT_HXX */
