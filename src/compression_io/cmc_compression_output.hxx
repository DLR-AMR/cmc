#ifndef CMC_COMPRESSION_OUTPUT_HXX
#define CMC_COMPRESSION_OUTPUT_HXX

#include "lossy/cmc_ac_compression_variable.hxx"
#include "lossless/cmc_byte_compression_variable.hxx"
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

    void AddGlobalAttribute(const nc::Attribute& attribute);
    void AddGlobalAttribute(nc::Attribute&& attribute);

    void Write();

private:
    template<typename T> void SetDataVariable(cmc::lossless::AbstractByteCompressionVariable<T>* variable, const std::vector<VariableLevelOffset>& level_byte_counts, const std::vector<IntraLevelStreamOffset>& intra_level_offsets, const int var_id, const int corresponding_mesh_id);
    template<typename T> void SetMeshVariable(cmc::lossless::AbstractByteCompressionVariable<T>* variable, const std::vector<VariableLevelOffset>& level_byte_counts, const std::vector<IntraLevelStreamOffset>& intra_level_offsets, const int mesh_id);

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

template<typename T>
std::pair<std::vector<VariableLevelOffset>, std::vector<IntraLevelStreamOffset>>
CommunicateOffsets(cmc::lossless::AbstractByteCompressionVariable<T>* variable, const MPI_Comm comm)
{
    /* Get a reference to the encoded streams */
    const std::vector<std::vector<uint8_t>>& levelwise_entropy_codes = variable->GetEncodedEntropyCodes();
    const std::vector<std::vector<uint8_t>>& levelwise_encoded_data = variable->GetEncodedData();
    const std::vector<std::vector<uint8_t>>& levelwise_encoded_mesh = variable->GetEncodedMesh();
    
    /* All streams must have the same amount of levels */
    cmc_assert(levelwise_entropy_codes.size() == levelwise_encoded_data.size() && levelwise_entropy_codes.size() == levelwise_encoded_mesh.size());

    /* Get the number of compression iterations */
    const size_t num_compression_iterations = levelwise_entropy_codes.size();

    /* Define an VariableLevelOffset MPI data type, for indicating the global start positions of
     * each encoded level within the data variable as well as in the mesh variable */
    MPI_Datatype MPI_VARIABLE_LEVEL_BYTE_COUNT;
    CreateVariableLevelOffsetMPIType(&MPI_VARIABLE_LEVEL_BYTE_COUNT);

    /* Define an LevelOffset MPI data type */
    MPI_Datatype MPI_INTRA_LEVEL_STREAM_OFFSET;
    CreateIntraLevelStreamOffsetMPIType(&MPI_INTRA_LEVEL_STREAM_OFFSET);

    std::vector<VariableLevelOffset> local_level_counts;
    local_level_counts.reserve(num_compression_iterations);

    /* Allocate the vector collecting the global offsets */
    std::vector<VariableLevelOffset> global_level_offsets(num_compression_iterations);

    /* Fill the local level counts per variable */
    for (size_t idx = 0; idx < num_compression_iterations; ++idx)
    {
        /* Get the local count of the encoded level-wise variable streams */
        local_level_counts.emplace_back(static_cast<uint64_t>(levelwise_entropy_codes[idx].size()), static_cast<uint64_t>(levelwise_encoded_data[idx].size()),
                                        static_cast<uint64_t>(levelwise_encoded_mesh[idx].size()));
    }

    /* Define a sum operation for the data type */
    MPI_Op count_bytes_levelwise;
    MPI_Op_create(&MpiVariableLevelByteSum, 1, &count_bytes_levelwise);

    /* Reduce the overall amount of bytes per level per variable */
    int ret_val = MPI_Allreduce(local_level_counts.data(), global_level_offsets.data(), static_cast<int>(num_compression_iterations),
                                MPI_VARIABLE_LEVEL_BYTE_COUNT, count_bytes_levelwise, comm);
    MPICheckError(ret_val);

    /* Exchange the intra-level offsets for the data streams */
    std::vector<IntraLevelStreamOffset> local_intra_level_offset_counts;
    local_intra_level_offset_counts.reserve(num_compression_iterations);

    /* Allocate the vector collecting the global intra-level offsets */
    std::vector<IntraLevelStreamOffset> global_intra_level_offsets(num_compression_iterations);

    /* Fill the local intra-level counts */
    for (size_t idx = 0; idx < num_compression_iterations; ++idx)
    {
        /* Get the local levelwise counts of the entropy codes and the encoded bits */
        local_intra_level_offset_counts.emplace_back(levelwise_entropy_codes[idx].size(), levelwise_encoded_data[idx].size(), levelwise_encoded_mesh[idx].size());
    }

    /* Define a sum operation for the data type */
    MPI_Op count_intra_level_offsets;
    MPI_Op_create(&MpiIntraLevelOffsetSum, 1, &count_intra_level_offsets);

    /* Exchange the intra-level offset with an exclusive scan */
    ret_val = MPI_Exscan(local_intra_level_offset_counts.data(), global_intra_level_offsets.data(), num_compression_iterations,
                         MPI_INTRA_LEVEL_STREAM_OFFSET, count_intra_level_offsets, comm);
    MPICheckError(ret_val);
    
    /* Free the MPI resources */
    ret_val = MPI_Op_free(&count_bytes_levelwise);
    MPICheckError(ret_val);
    ret_val = MPI_Op_free(&count_intra_level_offsets);
    MPICheckError(ret_val);
    ret_val = MPI_Type_free(&MPI_VARIABLE_LEVEL_BYTE_COUNT);
    MPICheckError(ret_val);
    ret_val = MPI_Type_free(&MPI_INTRA_LEVEL_STREAM_OFFSET);
    MPICheckError(ret_val);


    return std::make_pair(std::move(global_level_offsets), std::move(global_intra_level_offsets));
}

template<typename T>
void
Writer::SetDataVariable(cmc::lossless::AbstractByteCompressionVariable<T>* variable, const std::vector<VariableLevelOffset>& level_byte_counts, const std::vector<IntraLevelStreamOffset>& intra_level_offsets, const int var_id, const int corresponding_mesh_id)
{
    /* Get the encoded level-wise data */
    std::vector<std::vector<uint8_t>> levelwise_entropy_codes_;
    variable->MoveEncodedEntropyCodesInto(levelwise_entropy_codes_);

    std::vector<std::vector<uint8_t>> levelwise_data_;
    variable->MoveEncodedDataInto(levelwise_data_);

    cmc_assert(levelwise_data_.size() == levelwise_entropy_codes_.size());

    /* Get the number of compression iterations */
    const int num_compression_iterations = static_cast<int>(levelwise_data_.size());

    /* Gather the global count of bytes for this variable */
    uint64_t global_byte_count{0};

    /* Iterate over all codes and data and gather the overall byte count */
    for (int idx = num_compression_iterations - 1; idx >= 0; --idx)
    {
        global_byte_count += level_byte_counts[idx].entropy_codes_count;
        global_byte_count += level_byte_counts[idx].significant_data_count;
    }

    /* Define a new variable for the mesh  */
    nc::SpecificVariable<uint8_t> spec_variable(variable->GetName(), var_id);

    /* Set the global diemnsion length */
    spec_variable.SetGlobalDimensionLength(global_byte_count);

    /* Set the levelwise data compliant to the offsets.
     * We do this in reverse order in order to reconstruct the mesh fromt the root. */
    uint64_t level_offset = 0;
    for (int idx = num_compression_iterations - 1; idx >= 0; --idx)
    {
        /* Get the intra-level offset for the entropy codes based on the rank of this process */
        uint64_t intra_level_offset = intra_level_offsets[idx].entropy_codes_offset;

        /* Set the entropy cods of the current level with the correct offset */
        spec_variable.SetData(std::move(levelwise_entropy_codes_[idx]), level_offset + intra_level_offset);
        
        /* Move to the data part of this level for this rank */
        intra_level_offset = level_byte_counts[idx].entropy_codes_count + intra_level_offsets[idx].significant_bits_offset;

        /* Set the data of this level for this rank */
        spec_variable.SetData(std::move(levelwise_data_[idx]), level_offset + intra_level_offset);

        /* Update the levelwise offset with the global number of bytes in this level */
        level_offset += level_byte_counts[idx].entropy_codes_count + level_byte_counts[idx].significant_data_count;
    }

    /* Now, we define some attributes that will be added to the variable */
    std::vector<cmc::nc::Attribute> attributes;
    attributes.reserve(num_compression_iterations + 5);

    attributes.emplace_back(id_attr, var_id);
    attributes.emplace_back(data_type_attr, static_cast<int>(ConvertToCmcType<T>()));
    attributes.emplace_back(mesh_id_attr, corresponding_mesh_id);
    cmc::CompressionSchema schema = variable->GetCompressionSchema();
    attributes.emplace_back(compression_schema_attr, static_cast<int>(schema));
    attributes.emplace_back(num_compression_iterations_attr, num_compression_iterations);

     /* Store the level-wise byte count */
    for (int idx = num_compression_iterations - 1; idx >= 0; --idx)
    {
        const std::string attr_name = "lvl_" + std::to_string(num_compression_iterations - 1 - idx) + "_num_bytes";
        attributes.emplace_back(attr_name.c_str(), level_byte_counts[idx].entropy_codes_count + level_byte_counts[idx].significant_data_count);
    }

    /* Add the variable for output */
    nc_writer_.AddVariable(nc::Variable(std::move(spec_variable), std::move(attributes)));
}

template<typename T>
void
Writer::SetMeshVariable(cmc::lossless::AbstractByteCompressionVariable<T>* variable, const std::vector<VariableLevelOffset>& level_byte_counts, const std::vector<IntraLevelStreamOffset>& intra_level_offsets, const int mesh_id)
{
    /* Get the encoded mesh data */
    std::vector<std::vector<uint8_t>> levelwise_encoded_mesh;
    variable->MoveEncodedMeshInto(levelwise_encoded_mesh);

    /* Get the number of compression iterations */
    const int num_compression_iterations = static_cast<int>(levelwise_encoded_mesh.size());

    /* Define an name for the mesh variable */
    const std::string mesh_var_name = GenerateMeshVariableName(mesh_id);

    /* Gather the global count of bytes for this variable */
    uint64_t global_byte_count{0};

    /* We iterate over all compression levels */
    for (auto lvl_iter = level_byte_counts.begin(); lvl_iter != level_byte_counts.end(); ++lvl_iter)
    {
        global_byte_count += lvl_iter->encoded_mesh_count;
    }

    /* Define a new variable for the mesh  */
    nc::SpecificVariable<uint8_t> spec_variable(mesh_var_name, mesh_id);

    /* Set the global diemnsion length */
    spec_variable.SetGlobalDimensionLength(global_byte_count);

    /* Set the levelwise data compliant to the offsets.
     * We do this in reverse order in order to reconstruct the mesh fromt the root. */
    uint64_t level_offset = 0;
    for (int idx = num_compression_iterations - 1; idx >= 0; --idx)
    {
        /* Get the intra-level offset based on the rank of this process */
        const uint64_t intra_level_offset = intra_level_offsets[idx].mesh_encoding_offset;

        /* Set the data of the current level with the correct offset */
        spec_variable.SetData(std::move(levelwise_encoded_mesh[idx]), level_offset + intra_level_offset);

        /* Update the levelwise offset with the global number of bytes in this level */
        level_offset += level_byte_counts[idx].encoded_mesh_count;
    }

    /* Now, we define some attributes that will be added to the variable */
    std::vector<cmc::nc::Attribute> attributes;
    attributes.reserve(num_compression_iterations + 2);

    /* Store some general information */
    attributes.emplace_back(id_attr, mesh_id);
    attributes.emplace_back(num_compression_iterations_attr, num_compression_iterations);

    /* Store the level-wise byte count */
    for (int idx = num_compression_iterations - 1; idx >= 0; --idx)
    {
        const std::string attr_name = "lvl_" + std::to_string(num_compression_iterations - 1 - idx) + "_num_bytes";
        attributes.emplace_back(attr_name.c_str(), level_byte_counts[idx].encoded_mesh_count);
    }

    /* Add the variable for output */
    nc_writer_.AddVariable(nc::Variable(std::move(spec_variable), std::move(attributes)));
}

template<typename T>
void
Writer::SetVariable(cmc::lossless::AbstractByteCompressionVariable<T>* variable)
{
    /* Get the communicator */
    MPI_Comm comm = variable->GetMPIComm();

    /* Communicate all offsets */
    const auto [level_global_byte_counts, intra_level_data_offsets] = CommunicateOffsets<T>(variable, comm);

    /* The encoding of the mesh is also byte-aligned (the bits are getting offsetted during the compression at the end of each iteration) */
    SetMeshVariable<T>(variable, level_global_byte_counts, intra_level_data_offsets, mesh_id_counter_);
    
    /* The variable data is already byte-aligned and can be written to the file */
    SetDataVariable<T>(variable, level_global_byte_counts, intra_level_data_offsets, var_id_counter_, mesh_id_counter_);

    /* Update the variable counters */
    ++mesh_id_counter_;
    ++var_id_counter_;
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


#endif /* !CMC_COMPRESSION_OUTPUT_HXX */
