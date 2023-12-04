#ifndef CMC_T8_ADAPT_HXX
#define CMC_T8_ADAPT_HXX
/**
 * @file cmc_t8_data.hxx
 */

#include "lossy/cmc_amr_lossy_compression_settings.hxx" 
#include "t8code/cmc_t8_data_variables.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_utilities.hxx"

#include <vector>
#include <memory>

namespace cmc {

class InaccuracyTracker
{
public:
    InaccuracyTracker(){};
    ~InaccuracyTracker(){};

private:

};

class AdaptData
{
public:
    AdaptData() = delete;
    AdaptData(std::shared_ptr<CompressionSettings> settings, const CoarsenMesh& adaptation_sample, std::shared_ptr<std::vector<Variable>> compression_variables)
    : compression_settings{settings}, variables{compression_variables}, corresponding_variable_id{adaptation_sample.variable_id};
    ~AdaptData(){};
    
    AdaptData(const AdaptData& other) = default;
    AdaptData& operator=(const AdaptData& other) = default;
    AdaptData(AdaptData&& other) = default;
    AdaptData& operator=(AdaptData&& other) = default;

    GetAdaptationFunction() const;
    void InterpolateData(t8_forest_t adapted_forest);
    t8_forest_t RepartitionData(t8_forest_t adapted_forest);

private:
    std::shared_ptr<CompressionSettings> compression_settings;
    std::shared_ptr<std::vector<Variable>> variables;
    int corresponding_variable_id;

    InaccuracyTracker deviations;
};


}

#endif /* !CMC_T8_ADAPT_HXX */
