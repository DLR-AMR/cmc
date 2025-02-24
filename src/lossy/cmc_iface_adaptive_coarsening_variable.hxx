#ifndef CMC_IFACE_ADAPTIVE_COARSENING_VARIABLE_HXX
#define CMC_IFACE_ADAPTIVE_COARSENING_VARIABLE_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_compression_settings.hxx"

namespace cmc::lossy
{

class AbstractAdaptiveCoarseningVariable
{
public:
    void Setup() final;

protected: 
    AbstractAdaptiveCoarseningVariable() = default;
    void SplitVariables(){};
    void CheckConsistencyOfInputVariables(){};
    void BuildInitialMesh(){};
    void ApplyScalingAndOffset(){};
    void SetupCompressionFeatures(){};
};

inline void
AbstractAdaptiveCoarseningVariable::Setup()
{
    this->SplitVariables();

    this->CheckConsistencyOfInputVariables();

    this->BuildInitialMesh();

    this->DistributeDataOnInitialMesh();

    this->ApplyScalingAndOffset();

    this->SetupCompressionFeatures();
}


}

#endif /* !CMC_IFACE_ADAPTIVE_COARSENING_VARIABLE_HXX */
