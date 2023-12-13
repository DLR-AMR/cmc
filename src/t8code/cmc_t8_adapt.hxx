#ifndef CMC_T8_ADAPT_HXX
#define CMC_T8_ADAPT_HXX
/**
 * @file cmc_t8_adapt.hxx
 */

#include "lossy/cmc_amr_lossy_compression_settings.hxx" 
#include "t8code/cmc_t8_data_variables.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_utilities.hxx"
#include "t8code/cmc_t8_adapt_track_inaccuracy.hxx"

#include <vector>
#include <memory>

namespace cmc {

//TODO: This enum has no effect by now
constexpr enum TopologyPreserving {NonePreservation = -1, Minimum, Maximum, Mean, All};

class AdaptData
{
public:
    AdaptData() = delete;
    AdaptData(std::shared_ptr<CompressionSettings> settings, const CoarsenMesh& adaptation_sample, std::shared_ptr<std::vector<Variable>> compression_variables)
    : compression_settings_{settings}, variables_{compression_variables}, corresponding_variable_id_{adaptation_sample.variable_id}{
        //Setup up InaccuracyTracker
    };
    ~AdaptData(){};
    
    AdaptData(const AdaptData& other) = default;
    AdaptData& operator=(const AdaptData& other) = default;
    AdaptData(AdaptData&& other) = default;
    AdaptData& operator=(AdaptData&& other) = default;

    bool IsCompressionProgressing() const;
    t8_forest_t GetCurrentMesh() const;
    t8_forest_t SetCurrentMesh(t8_forest_t forest);
    void InterpolateData(t8_forest_t adapted_forest);
    t8_forest_t RepartitionData(t8_forest_t adapted_forest);
    t8_forest_adapt_t GetAdaptationFunction() const;
    
private:
    std::shared_ptr<CompressionSettings> compression_settings_;
    std::shared_ptr<std::vector<Variable>> variables_;
    int corresponding_variable_id_;
    t8_gloidx_t previous_number_of_elements_{-1};
    t8_gloidx_t new_number_of_elements_{0};

    /* Make a vector out of it, for each variable, there is an inaccuracy tracker with the compliant error criterion */
    InaccuracyTracker deviations_;
};


}

#endif /* !CMC_T8_ADAPT_HXX */
