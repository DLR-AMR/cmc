#ifndef CMC_IFACE_COMPRESSION_ADAPT_DATA_HXX
#define CMC_IFACE_COMPRESSION_ADAPT_DATA_HXX

#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"

namespace cmc
{

class ICompressionAdaptData
{
public:
    virtual bool IsCompressionProgressing() = 0;
    virtual void InitializeCompressionIteration() = 0;
    virtual void FinalizeCompressionIteration() = 0;
    virtual void CompleteInterpolation(const t8_forest_t previous_forest, const t8_forest_t adapted_forest) = 0;
    virtual void RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest) = 0;
    virtual cmc::t8::AdaptationFn GetAdaptationFunction() = 0;
    virtual int EvaluateCoarsening(const int tree_id, const int elem_id, const int num_elements, const std::vector<PermittedError> permitted_errors) = 0;
    virtual int EvaluateCoarsening(const int tree_id, const int elem_id, const int num_elements, const std::vector<PermittedError> permitted_errors, const CmcUniversalType missing_value) = 0;
    virtual int LeaveElementUnchanged(const int tree_id, const int elem_id) = 0;

    virtual ~ICompressionAdaptData(){};
};

}


#endif /* !CMC_IFACE_COMPRESSION_ADAPT_DATA_HXX */
