#ifndef CMC_PREFIX_DECOMPRESSION_HXX
#define CMC_PREFIX_DECOMPRESSION_HXX

#include "utilities/cmc_utilities.hxx"
#include "cmc_compression.hxx"

namespace cmc
{

namespace prefix
{

struct CompressedVariableInfo;

class Decompressor : public lossy::Decompressor
{
public:
    Decompressor() = default;
    Decompressor(const std::string& file_name)
    : file_name_{file_name} {};
    virtual ~Decompressor() = default;

    void Setup();
    void Decompress() override;
    void DecompressVariable(const std::string& variable_base_name);

private:
    void PerformLevelWiseDecompression();

    const std::string file_name_;

    std::vector<CompressedVariableInfo> variable_info_;

    std::vector<ByteVar> compression_variables_;

    MPI_Comm comm_{MPI_COMM_WORLD};
};


/* A data carrier holding information about the compressed variables in the file */
struct CompressedVariableInfo
{
    CompressedVariableInfo(const std::string& name_)
    : base_name{name_} {};

    std::string base_name;
    int id;
    int num_time_steps{0};
    int num_split_vars{0};
    int global_context_information{kGlobalContextInformationNotGiven};
    int intial_refinement_level;
    GeoDomain domain;
    DataLayout initial_layout{DataLayout::LayoutUndefined};
    DataLayout pre_compression_layout{DataLayout::LayoutUndefined};
    CmcType type{CmcType::TypeUndefined};
    CmcUniversalType missing_value;
};

}

}
#endif /* !CMC_PREFIX_DECOMPRESSION_HXX */
