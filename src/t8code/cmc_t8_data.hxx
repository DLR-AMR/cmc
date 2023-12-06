#ifndef CMC_T8_DATA_HXX
#define CMC_T8_DATA_HXX
/**
 * @file cmc_t8_data.hxx
 */

#include "t8code/cmc_t8_data_variables.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_utilities.hxx"
#include "t8code/cmc_t8_adapt.hxx"

#ifdef CMC_ENABLE_MPI
#include <mpi.h>
#endif

#include <vector>
#include <memory>

namespace cmc {

class AmrData
{
public:
    AmrData() = delete;
    AmrData(const std::vector<Variable>& compression_variables, std::shared_ptr<CompressionSettings> settings)
    : compression_settings{settings} {
        variables = std::make_shared(std::vector<Variable>(compression_variables));
    };
    AmrData(std::vector<Variable>&& compression_variables, std::shared_ptr<CompressionSettings> settings)
    : compression_settings{settings} {
        variables = std::make_shared(std::vector<Variable>(std::move(compression_variables)));
    };

    AmrData(const AmrData& other) = default;
    AmrData& operator=(const AmrData& other) = default;
    AmrData(AmrData&& other) = default;
    AmrData& operator=(AmrData&& other) = default;

    ~AmrData() = default;

    void BuildInitialMesh();
    void DistributeDataOnInitialMesh();
    void CompressByAdaptiveCoarsening(const CompressionMode compression_mode);

private:
    void SplitVariables();
    void ApplyVariableAttributes();
    bool ValidateSetup();

    AdaptData CreateAdaptationData(const AmrMesh& adaptation_sample) const;
    std::vector<AmrMesh> RetrieveMeshesToBeCoarsened(const CompressionMode compression_mode) const;

    MPI_Comm comm;
    AmrMesh initial_mesh;
    std::shared_ptr<CompressionSettings> compression_settings;
    std::shared_ptr<std::vector<Variable>> variables;
};


}

#endif /* !CMC_T8_DATA_HXX */
