#ifndef CMC_T8_DATA_HXX
#define CMC_T8_DATA_HXX
/**
 * @file cmc_t8_data.hxx
 */

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_compression_settings.hxx"
#include "t8code/cmc_t8_data_variables.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_adapt.hxx"
#include "utilities/cmc_output_variable_forward.hxx"
#include "t8code/cmc_t8_byte_variable.hxx"
#include "utilities/cmc_bit_map.hxx"

#ifdef CMC_ENABLE_MPI
#include <mpi.h>
#endif

#include <vector>
#include <memory>

namespace cmc
{

struct AdaptiveCoarseningIndications
{
    AdaptiveCoarseningIndications(const std::vector<bit_map::BitMap>& indications)
    : ac_indicator_bits{indications} {};
    AdaptiveCoarseningIndications(std::vector<bit_map::BitMap>&& indications)
    : ac_indicator_bits{std::move(indications)} {};

    std::vector<bit_map::BitMap> ac_indicator_bits;
};

class AmrData
{
public:
    using input_var_const_iterator = std::vector<InputVar>::const_iterator;

    AmrData() = delete;
    ~AmrData() = default;
    
    AmrData(const std::vector<InputVar>& compression_variables, const CompressionSettings& settings)
    : input_variables_{compression_variables}, compression_settings_{settings} {};
    AmrData(std::vector<InputVar>&& compression_variables, const CompressionSettings& settings)
    : input_variables_{std::move(compression_variables)}, compression_settings_{settings} {};
    AmrData(const std::vector<InputVar>& compression_variables)
    : input_variables_{compression_variables} {};
    AmrData(std::vector<InputVar>&& compression_variables)
    : input_variables_{std::move(compression_variables)} {};

    AmrData(const AmrData& other) = default;
    AmrData& operator=(const AmrData& other) = default;
    AmrData(AmrData&& other) = default;
    AmrData& operator=(AmrData&& other) = default;

    void SetCompressionSettings(const CompressionSettings& settings);
    void SetCompressionSettings(CompressionSettings&& settings);

    void SplitVariables();
    bool CheckConsistencyOfInputVariables() const;
    void BuildInitialMesh();
    void DistributeDataOnInitialMesh();
    void CompressByAdaptiveCoarsening(const CompressionMode compression_mode);
    void CompressByAdaptiveCoarseningWithRegardToInitialData(); //Automatically one for one compression
    void ApplyScalingAndOffset();
    void TransformInputToCompressionVariables();
    void SetupCompressionCriteria();
    void SetupInaccuracyStorage();
    void SetupVariablesForCompression();
    void SetInitialMesh(const AmrMesh& mesh);
    std::vector<ByteVar> GetByteVariablesForCompression();
    [[nodiscard]] std::vector<AdaptiveCoarseningIndications> TransferIndicationBits();

    void IndicateToKeepInitialMeshAndData();
    
    void DecompressToInitialRefinementLevel(const bool restrict_to_global_domain = true);
    std::vector<OutputVar> SeizeRawDecompressedVariable();
    OutputVar DecompressVariable(const int variable_id);

    input_var_const_iterator GetInputVariablesBegin() const { return input_variables_.begin(); };
    input_var_const_iterator GetInputVariablesEnd() const { return input_variables_.end(); };
    input_var_const_iterator GetInputVariablesCBegin() const { return input_variables_.cbegin(); };
    input_var_const_iterator GetInputVariablesCEnd() const { return input_variables_.cend(); };

    size_t GetNumberOfInputVariables() const;
    size_t GetNumberOfCompressionVariables() const;
    bool IsValidForCompression() const;

    int GetMpiSize() const;
    int GetMpiRank() const;

    void WriteCompressedData(const std::string& file_name) const;
    void WriteVTKFilePerVariable(const std::string& file_name) const;

    const std::vector<Var>& GetCompressionVariables() const {return variables_;};
    void FilterDataAsDifferences();
private:
    void SetCompressionSettings(CompressionSettings& settings) {compression_settings_ = settings;};
    std::vector<IndexReduction> UpdateLinearIndicesToTheInitialMesh();
    void TransferToDecompressionData(OutputVar& output_variable, const Var& compression_variable, const GeoDomain& resulting_domain);
    void CreateInternIDsForRedistribution();
    AdaptData CreateAdaptationData(const CoarseningSample& adaptation_sample, const CompressionMode mode);
    std::vector<CoarseningSample> RetrieveMeshesToBeCoarsened(const CompressionMode compression_mode) const;

    [[nodiscard]] std::pair<std::vector<VariableSendMessage>, std::vector<MPI_Request>> SendInitialData();
    std::vector<VariableRecvMessage> ReceiveInitialData();
    void SortInitialDataIntoVariables(const std::vector<VariableRecvMessage>& messages);
    void SortLocalDataOnInitialMesh();

    std::vector<Var> variables_; //To be used after the Setup call
    
    std::vector<InputVar> input_variables_;
    CompressionSettings compression_settings_;
    AmrMesh initial_mesh_;

    std::vector<AdaptiveCoarseningIndications> ac_indications_;

    MPI_Comm comm_{MPI_COMM_WORLD};
};

}

#endif /* !CMC_T8_DATA_HXX */
