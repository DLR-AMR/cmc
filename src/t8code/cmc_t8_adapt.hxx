#ifndef CMC_T8_ADAPT_HXX
#define CMC_T8_ADAPT_HXX
/**
 * @file cmc_t8_adapt.hxx
 */

#include "utilities/cmc_utilities.hxx"
#include "lossy/cmc_amr_lossy_compression_settings.hxx" 
#include "t8code/cmc_t8_data_variables.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_adapt_track_inaccuracy.hxx"
#include "t8code/cmc_t8_byte_variable.hxx"

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
    
    Var& GetCurrentCompressionVariable()
    {
        return variables_[corresponding_variable_id_];
    };

    void InitializeCompressionIteration()
    {
        previous_number_of_elements_ = t8_forest_get_global_num_elements(GetCurrentMesh());
        GetCurrentCompressionVariable().InitializeVariableForCompressionIteration();
    };

    void FinalizeCompressionIteration()
    {
        new_number_of_elements_ = t8_forest_get_global_num_elements(GetCurrentMesh());
        ++count_adaptation_step_;
    };

    void UpdateCompressionData()
    {
        Var& compression_variable = GetCurrentCompressionVariable();
        compression_variable.UpdateCompressionData();
    }

private:
    const CompressionSettings& compression_settings_;
    std::vector<Var>& variables_;
    t8_gloidx_t previous_number_of_elements_{-1};
    t8_gloidx_t new_number_of_elements_{0};
    int corresponding_variable_id_;
    CompressionMode mode_;
    int count_adaptation_step_{0};
};


class PrefixAdaptData
{
public:
    PrefixAdaptData() = delete;
    PrefixAdaptData(ByteVar& variable)
    : byte_variable_{variable}{
        new_number_of_elements_ = t8_forest_get_global_num_elements(variable.GetAmrMesh().GetMesh());
    };

    void InitializeCompressionIteration()
    {
        previous_number_of_elements_ = new_number_of_elements_;
        byte_variable_.InitializeCompressionIteration();
    }

    bool IsCompressionProgressing()
    {
        return (previous_number_of_elements_ != new_number_of_elements_ ? true : false);
        //if (count_adaptation_step_ == 0)
        //{
        //    return true;
        //}
        //return byte_variable_.HasAtLeastOnePrefixBeenFoundInTheLastIteration();
    }

    void FinalizeCompressionIteration()
    {
        new_number_of_elements_ = t8_forest_get_global_num_elements(GetCurrentMesh());
        ++count_adaptation_step_;
    }
    t8_forest_t GetCurrentMesh() const
    {
        return byte_variable_.GetMesh();
    }
    void SetCurrentMesh(t8_forest_t mesh)
    {
        byte_variable_.SetMesh(mesh);
    }
    void IndicateNoPrefixFound()
    {
        byte_variable_.IndicateNoPrefixFound();
    }
    void LeaveCoarsePrefixUnchanged(const int elem_index)
    {
        byte_variable_.LeaveInitialValueUnchanged(elem_index);
        //IndicateNoPrefixFound(*(byte_values_.data() + elem_index));
        //byte_variable_.SetPrefix(*(byte_values_.data() + elem_index));
        //byte_values_[elem_index].SetFrontBit(sizeof(T) * CHAR_BIT - (*(byte_values_.data() + elem_index)).GetTrailBit());
        //if (count_adaptation_step_ == 0)
        //{
        //    byte_variable_.LeaveInitialValueUnchanged(elem_index);
        //} else
        //{
        //    byte_variable_.LeaveCoarsePrefixUnchanged(elem_index);
        //}
    }
    bool EvaluateCommonPrefix(const int start_index, const int num_elements)
    {
        if (count_adaptation_step_ == 0)
        {
            return byte_variable_.EvaluateCommonPrefixFromInitialData(start_index, num_elements);
        } else
        {
            return byte_variable_.EvaluateCommonPrefixFromPreviousPrefix(start_index, num_elements);
        }
    }
    void PrintNumPrefixIndicationBits() const
    {
        byte_variable_.GetNumberOfSetPrefixIndicationBits();
    }
    int GetInitialRefinementLevelOfMesh() const
    {
        return byte_variable_.GetAmrMesh().GetInitialRefinementLevel();
    }
    int GetAdaptationStepCount() const
    {
        return count_adaptation_step_;
    }
private:
    ByteVar& byte_variable_;
    t8_gloidx_t previous_number_of_elements_{-1};
    t8_gloidx_t new_number_of_elements_{0};
    int count_adaptation_step_{0};
};


class PrefixAdaptDataEGU
{
public:
    PrefixAdaptDataEGU() = delete;
    PrefixAdaptDataEGU(ByteVar& variable)
    : byte_variable_{variable}{
        new_number_of_elements_ = t8_forest_get_global_num_elements(variable.GetAmrMesh().GetMesh());
    };

    void InitializeCompressionIteration()
    {
        previous_number_of_elements_ = new_number_of_elements_;
        byte_variable_.InitializeCompressionIteration();
    }

    bool IsCompressionProgressing()
    {
        return (new_number_of_elements_ > 1 ? true : false);
    }

    void FinalizeCompressionIteration()
    {
        new_number_of_elements_ = t8_forest_get_global_num_elements(GetCurrentMesh());
        byte_variable_.FinalizeCompressionIterationEGU();
        ++count_adaptation_step_;
    }
    t8_forest_t GetCurrentMesh() const
    {
        return byte_variable_.GetMesh();
    }
    void SetCurrentMesh(t8_forest_t mesh)
    {
        byte_variable_.SetMesh(mesh);
    }

    void LeaveCoarsePrefixUnchangedEGU(const int elem_index)
    {
        byte_variable_.LeaveCoarsePrefixUnchangedEGU(elem_index);
    }

    void EvaluateCommonPrefixEGU(const int start_index, const int num_elements)
    {
        if (count_adaptation_step_ == 0)
        {
            return byte_variable_.EvaluateCommonPrefixFromInitialDataEGU(start_index, num_elements);
        } else
        {
            return byte_variable_.EvaluateCommonPrefixFromPreviousPrefixEGU(start_index, num_elements);
        }
    }

    void PrintNumPrefixIndicationBits() const
    {
        byte_variable_.GetNumberOfSetPrefixIndicationBits();
    }
    int GetInitialRefinementLevelOfMesh() const
    {
        return byte_variable_.GetAmrMesh().GetInitialRefinementLevel();
    }
    int GetAdaptationStepCount() const
    {
        return count_adaptation_step_;
    }
private:
    ByteVar& byte_variable_;
    t8_gloidx_t previous_number_of_elements_{INT16_MAX};
    t8_gloidx_t new_number_of_elements_{INT16_MAX};
    int count_adaptation_step_{0};
};

}

#endif /* !CMC_T8_ADAPT_HXX */
