#ifndef CMC_COMPRESSION_OUTPUT_HXX
#define CMC_COMPRESSION_OUTPUT_HXX

#include "lossy/cmc_ac_compression_variable.hxx"
#include "lossless/cmc_byte_compression_variable.hxx"
#include "lossless/cmc_embedded_byte_compression_variable.hxx"
#include "compression_io/cmc_compression_attr_names.hxx"

#ifdef CMC_WITH_NETCDF
#include "netcdf/cmc_netcdf.hxx"
#include "netcdf/cmc_nc_io.hxx"
#include "netcdf/cmc_nc_writer.hxx"
#endif

#include "mpi/cmc_mpi.hxx"

#include <string>
#include <vector>
#include <algorithm>
#include <utility>

namespace cmc::compression_io
{

struct VariableLevelOffset;
struct IntraLevelStreamOffset;

class Writer
{
public:
    Writer(const std::string& file_name)
    : file_name_{file_name}, nc_writer_(file_name, NC_NETCDF4) {};

    template<typename T> void SetVariable(cmc::lossless::AbstractByteCompressionVariable<T>* variable);
    template<typename T> void SetVariable(cmc::lossless::embedded::AbstractEmbeddedByteCompressionVariable<T>* variable);

    void AddGlobalAttribute(const nc::Attribute& attribute);
    void AddGlobalAttribute(nc::Attribute&& attribute);

    void Write();

private:
    template<typename T> void SetDataVariable(cmc::lossless::AbstractByteCompressionVariable<T>* variable, const std::vector<VariableLevelOffset>& level_byte_counts, const std::vector<IntraLevelStreamOffset>& intra_level_offsets, const int var_id, const int corresponding_mesh_id);
    template<typename T> void SetMeshVariable(cmc::lossless::AbstractByteCompressionVariable<T>* variable, const std::vector<VariableLevelOffset>& level_byte_counts, const std::vector<IntraLevelStreamOffset>& intra_level_offsets, const int mesh_id);
    template<typename T> void SetDataVariable(cmc::lossless::embedded::AbstractEmbeddedByteCompressionVariable<T>* variable, const std::vector<VariableLevelOffset>& level_byte_counts, const std::vector<IntraLevelStreamOffset>& intra_level_offsets, const int var_id, const int corresponding_mesh_id);
    template<typename T> void SetMeshVariable(cmc::lossless::embedded::AbstractEmbeddedByteCompressionVariable<T>* variable, const std::vector<VariableLevelOffset>& level_byte_counts, const std::vector<IntraLevelStreamOffset>& intra_level_offsets, const int mesh_id);

    std::string file_name_;
    cmc::nc::Writer nc_writer_;
    int var_id_counter_{0};
    int mesh_id_counter_{0};
};

struct IntraLevelStreamOffset
{
    IntraLevelStreamOffset() = default;
    IntraLevelStreamOffset(const uint64_t entropy_codes_offset_count, const uint64_t significant_bits_offset_count, const uint64_t mesh_encoding_offset_count)
    : entropy_codes_offset{entropy_codes_offset_count}, significant_bits_offset{significant_bits_offset_count}, mesh_encoding_offset{mesh_encoding_offset_count} {};

    uint64_t entropy_codes_offset{0};
    uint64_t significant_bits_offset{0};
    uint64_t mesh_encoding_offset{0};
};

struct VariableLevelOffset
{
    VariableLevelOffset() = default;
    VariableLevelOffset(const uint64_t entropy_codes_count_, const uint64_t significant_data_count_, const uint64_t encoded_mesh_count_)
    : entropy_codes_count{entropy_codes_count_}, significant_data_count{significant_data_count_}, encoded_mesh_count{encoded_mesh_count_} {};

    uint64_t entropy_codes_count{0};
    uint64_t significant_data_count{0};
    uint64_t encoded_mesh_count{0};
};
 
inline void
CreateIntraLevelStreamOffsetMPIType(MPI_Datatype* intra_level_stream_offset_type)
{
    /* Define the properties of the custom 'LevelOffset' data type */
    const int num_fields = 3;
    int array_of_blocklengths[] = {1,1,1};
    MPI_Aint array_of_displacements[num_fields];
    array_of_displacements[0] = offsetof(IntraLevelStreamOffset, entropy_codes_offset);
    array_of_displacements[1] = offsetof(IntraLevelStreamOffset, significant_bits_offset);
    array_of_displacements[2] = offsetof(IntraLevelStreamOffset, mesh_encoding_offset);
    MPI_Datatype array_of_types[] = {MPI_UINT64_T, MPI_UINT64_T, MPI_UINT64_T};

    const int ret_val_struct = MPI_Type_create_struct(num_fields, array_of_blocklengths, array_of_displacements,
                                                      array_of_types, intra_level_stream_offset_type); 
    MPICheckError(ret_val_struct);
    const int ret_val_type_commit = MPI_Type_commit(intra_level_stream_offset_type);
    MPICheckError(ret_val_type_commit); 
}

inline void
CreateVariableLevelOffsetMPIType(MPI_Datatype* variable_level_offset_type)
{
    /* Define the properties of the custom 'LevelOffset' data type */
    const int num_fields = 3;
    int array_of_blocklengths[] = {1,1,1};
    MPI_Aint array_of_displacements[num_fields];
    array_of_displacements[0] = offsetof(VariableLevelOffset, entropy_codes_count);
    array_of_displacements[1] = offsetof(VariableLevelOffset, significant_data_count);
    array_of_displacements[2] = offsetof(VariableLevelOffset, encoded_mesh_count);
    MPI_Datatype array_of_types[] = {MPI_UINT64_T, MPI_UINT64_T, MPI_UINT64_T};

    const int ret_val_struct = MPI_Type_create_struct(num_fields, array_of_blocklengths, array_of_displacements,
                                                      array_of_types, variable_level_offset_type); 
    MPICheckError(ret_val_struct);
    const int ret_val_type_commit = MPI_Type_commit(variable_level_offset_type);
    MPICheckError(ret_val_type_commit); 
}

inline void
MpiVariableLevelByteSum(void* input_buffer, void* output_buffer, int* count, [[maybe_unused]] MPI_Datatype* datatype)
{
    /* Cast the pointers to the actual type */
    VariableLevelOffset* input = (VariableLevelOffset*) input_buffer;
    VariableLevelOffset* output = (VariableLevelOffset*) output_buffer;
 
    /* Sum up all members */
    for(int i = 0; i < *count; ++i)
    {
        output[i].entropy_codes_count += input[i].entropy_codes_count;
        output[i].significant_data_count += input[i].significant_data_count;
        output[i].encoded_mesh_count += input[i].encoded_mesh_count;
    }
}

inline void
MpiIntraLevelOffsetSum(void* input_buffer, void* output_buffer, int* count, [[maybe_unused]] MPI_Datatype* datatype)
{
    /* Cast the pointers to the actual type */
    IntraLevelStreamOffset* input = (IntraLevelStreamOffset*) input_buffer;
    IntraLevelStreamOffset* output = (IntraLevelStreamOffset*) output_buffer;
 
    /* Sum up all members */
    for(int i = 0; i < *count; ++i)
    {
        output[i].entropy_codes_offset += input[i].entropy_codes_offset;
        output[i].significant_bits_offset += input[i].significant_bits_offset;
        output[i].mesh_encoding_offset += input[i].mesh_encoding_offset;
    }
}

inline void
Writer::AddGlobalAttribute(const nc::Attribute& attribute)
{
    nc_writer_.AddGlobalAttribute(attribute);
}

inline void
Writer::AddGlobalAttribute(nc::Attribute&& attribute)
{
    nc_writer_.AddGlobalAttribute(std::move(attribute));
}

inline void
Writer::Write()
{
    nc_writer_.Write();
}

}

/* The implementation specific fuctionalities for the AbstractByteCompressionVariables are included */
#include "compression_io/cmc_byte_compression_output.txx"

/* The implementation specific fuctionalities for the AbstractEmbeddedByteCompressionVariables are included */
#include "compression_io/cmc_embedded_byte_compression_output.txx"


#endif /* !CMC_COMPRESSION_OUTPUT_HXX */
