#ifndef CMC_LOSSY_AMR_PAR_MULTI_RES_DECOMPRESSION_IDW_HXX
#define CMC_LOSSY_AMR_PAR_MULTI_RES_DECOMPRESSION_IDW_HXX

#include "cmc.hxx"
#include "amr/lossy/cmc_par_multi_res_extraction_util.hxx"
#include "amr/lossy/cmc_par_multi_res_error_mesh.hxx"
#include "amr/lossy/cmc_par_multi_res_idw_interpolation_util.hxx"
#include "mpi/cmc_mpi.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "utilities/cmc_compression_schema.hxx"
#include "utilities/cmc_bits_stream_decoder.hxx"
#include "amr/lossy/cmc_par_multi_res_idw_interpolation_util.hxx"
#include <t8_forest/t8_forest_partition.h>

#include <string>
#include <span>
#include <filesystem>

namespace cmc::par::lossy::multi_res::idw
{

/* Switch to write out intermediate data from the decompression steps */
constexpr bool kWriteVTKDecompressionStepData = true;

/* Forward declaration of the general compression variable */
template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
class DecompressionVariableMultiData;

/* Typedef for the general case with one data point per variable */
template<ArithmeticType T, int32_t DIM>
using DecompressionVariable = DecompressionVariableMultiData<T, DIM, int32_t{1}>;

struct CompressionInfoStruct
{
    cmc::bits::vector_view
    GetGlobalElementIndicationsStep(const int step)
    {
        cmc_assert(step >= 0 && step < static_cast<int>(this->mesh_compression_levels) - 1);
        cmc_assert((this->GetOffsetPartitionTable() + this->partition_table_size + this->mesh_compression_huffman_codes + this->intra_elem_compression_huffman_codes) % sizeof(uint64_t) == 0);

        int offset = (this->GetOffsetPartitionTable() + this->partition_table_size + this->mesh_compression_huffman_codes + this->intra_elem_compression_huffman_codes) / sizeof(uint64_t);
        for (int step_idx{0}; step_idx < step; ++step_idx)
        {
            offset += global_num_elem_indications_step[step_idx] / (sizeof(uint64_t) * cmc::bits::kCharBit) + (global_num_elem_indications_step[step_idx] % (sizeof(uint64_t) * cmc::bits::kCharBit) != 0 ? 1 : 0);
        }

        return cmc::bits::vector_view(shared_data_init_ptr + offset);
    }

    std::vector<LevelPartition>
    GetLevelPartition(const int level) const
    {
        cmc_assert(level >= 0);
        std::vector<LevelPartition> level_partition;
        level_partition.reserve(this->compression_comm_size);
        cmc_assert((this->GetOffsetPartitionTable() + level * this->compression_comm_size * 2 * sizeof(SizeType)) % sizeof(uint64_t) == 0);

        const int offset = (this->GetOffsetPartitionTable() + level * this->compression_comm_size * 2 * sizeof(SizeType)) / sizeof(uint64_t);

        for(int rank_idx{0}; rank_idx < static_cast<int>(this->compression_comm_size); ++rank_idx)
        {
            level_partition.emplace_back(cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(this->shared_data_init_ptr + offset + 2 * rank_idx)),
                                         cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(this->shared_data_init_ptr + offset + 2 * rank_idx + 1)));
        }
        return level_partition;
    }

    const uint64_t*
    GetMeshCompressionHuffmanCodes() const
    {
        cmc_assert((this->GetOffsetPartitionTable() + this->partition_table_size) % sizeof(uint64_t) == 0);
        const int offset = (this->GetOffsetPartitionTable() + this->partition_table_size) / sizeof(uint64_t);
        return shared_data_init_ptr + offset;
    }

    const uint64_t*
    GetIntraElementCompressionHuffmanCodes() const
    {
        cmc_assert((this->GetOffsetPartitionTable() + this->partition_table_size + this->mesh_compression_huffman_codes) % sizeof(uint64_t) == 0);
        const int offset = (this->GetOffsetPartitionTable() + this->partition_table_size + this->mesh_compression_huffman_codes) / sizeof(uint64_t);
        return shared_data_init_ptr + offset;
    }

    SizeType GetOffsetPartitionTable() const {return sizeof(SizeType) * (15 + global_level_bytes.size() + global_num_elem_indications_step.size())
                                                     + kNumCharsVariableName * sizeof(char);}

    const char* GetName() const {return name.data();}             

    SizeType global_byte_count{0};
    SizeType offset_start_encoding{0};
    std::array<char, kNumCharsVariableName> name;
    SizeType data_type{0};
    SizeType dimensionality{0};
    SizeType data_per_elem{0};
    SizeType compression_comm_size{0};
    SizeType compression_scheme{0};
    SizeType mesh_compression_levels{0};
    SizeType intra_elem_compression_levels{0};
    SizeType pack_size{0};
    std::vector<SizeType> global_level_bytes;
    std::vector<SizeType> global_num_elem_indications_step; 
    SizeType partition_table_size{0};
    SizeType mesh_compression_huffman_codes{0};
    SizeType intra_elem_compression_huffman_codes{0};
    SizeType num_bytes_global_mesh_encoding{0};
    SizeType num_bytes_global_error_mesh_encoding{0};

    uint64_t* shared_data_init_ptr{nullptr};
    MPI_Win shared_data;
};

inline void 
CheckMPIReadCorrectness(const MPI_Status* status, const MPI_Datatype datatype, const int expected_num_elems)
{
    int elem_count_{0};
    const int rv_check_read = MPI_Get_count(status, datatype, &elem_count_);
    MPICheckError(rv_check_read);
    if (elem_count_ != expected_num_elems)
    {
        cmc_global_msg("The expeceted number of bytes could not be read from the file.");
        MPICheckError(MPI_ERR_COUNT);
    }
}

inline CompressionInfoStruct
ConstructCompressionInfoStruct(const std::vector<uint64_t>& encoded_preamble)
{
    static_assert(sizeof(char) == sizeof(uint8_t));

    /* Define a start pointer to the data */
    const uint64_t* start_ptr = encoded_preamble.data();
    size_t offset{0};

    CompressionInfoStruct info;

    info.global_byte_count = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;

    info.offset_start_encoding = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;

    std::copy_n(reinterpret_cast<const char*>(start_ptr + offset), kNumCharsVariableName, info.name.data());
    offset += 32;

    info.data_type = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;

    info.dimensionality = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;

    info.data_per_elem = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;

    info.compression_comm_size = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;

    info.compression_scheme = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;

    info.mesh_compression_levels = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;

    info.intra_elem_compression_levels = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;

    info.pack_size = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;

    const int num_levels = info.mesh_compression_levels + (info.intra_elem_compression_levels > 0 ? 1 : 0);
    info.global_level_bytes.reserve(num_levels);
    for (int lvl_idx{0}; lvl_idx < num_levels; ++lvl_idx)
    {
        info.global_level_bytes.push_back(cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset)));
        ++offset;
    }

    info.global_num_elem_indications_step.reserve(info.mesh_compression_levels - 1);
    for (int lvl_idx{0}; lvl_idx < static_cast<int>(info.mesh_compression_levels) - 1; ++lvl_idx)
    {
        info.global_num_elem_indications_step.push_back(cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset)));
        ++offset;
    }

    info.partition_table_size = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;

    info.mesh_compression_huffman_codes = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;

    info.intra_elem_compression_huffman_codes = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;

    info.num_bytes_global_mesh_encoding = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;

    info.num_bytes_global_error_mesh_encoding = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));

    return info;
}

inline void
DecodeCompressionInfoStruct(CompressionInfoStruct& encoded_info)
{
    static_assert(sizeof(char) == sizeof(uint8_t));

    size_t offset{0};

    encoded_info.global_byte_count = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(encoded_info.shared_data_init_ptr + offset));
    ++offset;

    encoded_info.offset_start_encoding = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(encoded_info.shared_data_init_ptr + offset));
    ++offset;

    static_assert(32 * sizeof(SizeType) == kNumCharsVariableName);
    for (int idx{0}; idx < 32; ++idx)
    {
        const auto serialized_val = cmc::bits::SerializeValueBE<SizeType>(cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(encoded_info.shared_data_init_ptr + offset + idx)));
        std::copy_n(reinterpret_cast<const char*>(serialized_val.data()), sizeof(SizeType), encoded_info.name.data() + idx * sizeof(SizeType));
    }
    offset += 32;

    encoded_info.data_type = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(encoded_info.shared_data_init_ptr + offset));
    ++offset;

    encoded_info.dimensionality = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(encoded_info.shared_data_init_ptr + offset));
    ++offset;

    encoded_info.data_per_elem = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(encoded_info.shared_data_init_ptr + offset));
    ++offset;

    encoded_info.compression_comm_size = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(encoded_info.shared_data_init_ptr + offset));
    ++offset;

    encoded_info.compression_scheme = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(encoded_info.shared_data_init_ptr + offset));
    ++offset;

    encoded_info.mesh_compression_levels = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(encoded_info.shared_data_init_ptr + offset));
    ++offset;

    encoded_info.intra_elem_compression_levels = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(encoded_info.shared_data_init_ptr + offset));
    ++offset;

    encoded_info.pack_size = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(encoded_info.shared_data_init_ptr + offset));
    ++offset;

    const int num_levels = encoded_info.mesh_compression_levels + (encoded_info.intra_elem_compression_levels > 0 ? 1 : 0);
    encoded_info.global_level_bytes.reserve(num_levels);
    for (int lvl_idx{0}; lvl_idx < num_levels; ++lvl_idx)
    {
        encoded_info.global_level_bytes.push_back(cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(encoded_info.shared_data_init_ptr + offset)));
        ++offset;
    }

    encoded_info.global_num_elem_indications_step.reserve(encoded_info.mesh_compression_levels - 1);
    for (int lvl_idx{0}; lvl_idx < static_cast<int>(encoded_info.mesh_compression_levels) - 1; ++lvl_idx)
    {
        encoded_info.global_num_elem_indications_step.push_back(cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(encoded_info.shared_data_init_ptr + offset)));
        ++offset;
    }

    encoded_info.partition_table_size = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(encoded_info.shared_data_init_ptr + offset));
    ++offset;

    encoded_info.mesh_compression_huffman_codes = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(encoded_info.shared_data_init_ptr + offset));
    ++offset;

    encoded_info.intra_elem_compression_huffman_codes = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(encoded_info.shared_data_init_ptr + offset));
    ++offset;

    encoded_info.num_bytes_global_mesh_encoding = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(encoded_info.shared_data_init_ptr + offset));
    ++offset;

    encoded_info.num_bytes_global_error_mesh_encoding = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(encoded_info.shared_data_init_ptr + offset));
}

inline
void
PrintCompressionInfo(const std::vector<uint64_t>& encoded_preamble, const std::string& file_name)
{
    /* Get the compression info from the encoding */
    const CompressionInfoStruct info = ConstructCompressionInfoStruct(encoded_preamble);

    /* Print the gathered information */
    cmc_global_msg("Compression Information retrieved from: ", file_name);
    cmc_global_msg("\t Overall Byte Count: ", info.global_byte_count, " bytes");
    cmc_global_msg("\t Compression Preamble Byte Count: ", info.offset_start_encoding, " bytes");
    cmc_global_msg("\t Name: ", reinterpret_cast<const char*>(info.name.data()));
    cmc_global_msg("\t DataType: ", info.data_type);
    cmc_global_msg("\t Dimensionality: ", info.dimensionality, "D");
    cmc_global_msg("\t Data Per Elem: ", info.data_per_elem);
    cmc_global_msg("\t Compression Communicator Size: ", info.compression_comm_size);
    cmc_global_msg("\t Compression Scheme: ", info.compression_scheme);
    cmc_global_msg("\t Mesh Compression LVLs: ", info.mesh_compression_levels);
    cmc_global_msg("\t Intra Elem Compression LVLs: ", info.intra_elem_compression_levels);
    cmc_global_msg("\t Intra Elem Compaction Stencil: ", info.pack_size);

    cmc_global_msg("\t Mesh Compression Level Byte Count:");
    for (SizeType lvl_idx{0}; lvl_idx < info.mesh_compression_levels; ++lvl_idx)
    {
        cmc_global_msg("\t\t Level ", lvl_idx, ": ", info.global_level_bytes[lvl_idx], " bytes");
    }

    if (info.intra_elem_compression_levels > 0)
    {
        cmc_global_msg("\t\t Intra Elem Compression Byte Counts: ", info.global_level_bytes.back(), " bytes");
    }

    cmc_global_msg("\t Mesh Compression Global Element Count (size: ", info.global_num_elem_indications_step.size(), "): ");
    for (SizeType lvl_idx{0}; lvl_idx < info.mesh_compression_levels - 1; ++lvl_idx)
    {
        cmc_global_msg("\t\t Level ", lvl_idx, ": ", info.global_num_elem_indications_step[lvl_idx], " elements");
    }
    cmc_global_msg("\t Partition Table Byte Count: ", info.partition_table_size, " bytes");

    for (SizeType lvl_idx{0}; lvl_idx < info.mesh_compression_levels; ++lvl_idx)
    {
        const int offset = 15 + info.global_level_bytes.size() + info.global_num_elem_indications_step.size() +
                          (kNumCharsVariableName * sizeof(char) + lvl_idx * info.compression_comm_size * 2 * sizeof(SizeType)) / sizeof(uint64_t);

        for(int rank_idx{0}; rank_idx < static_cast<int>(info.compression_comm_size); ++rank_idx)
        {
            const SizeType elem_offset = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(encoded_preamble[offset + 2 * rank_idx]);
            const SizeType byte_offset = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(encoded_preamble[offset + 2 * rank_idx + 1]);
            cmc_global_msg("\t\t Partition Level: ", lvl_idx, ", Rank: ", rank_idx, ", Elem_Offset: ", elem_offset, ", Coding_Byte_Offset: ", byte_offset);
        }
    }

    cmc_global_msg("\t Serialized Mesh Compression Huffman Codes: ", info.mesh_compression_huffman_codes, " bytes");
    cmc_global_msg("\t Serialized Intra Elem Compression Huffman Codes: ", info.intra_elem_compression_huffman_codes, " bytes");
    cmc_global_msg("\t Serialized Global Mesh Encoding: ", info.num_bytes_global_mesh_encoding, " bytes");
    cmc_global_msg("\t Serialized Global Error-Mesh Encoding: ", info.num_bytes_global_error_mesh_encoding, " bytes");
}

inline
void
ReadCompressionInfo(const std::string& file_name)
{
    MPI_File fhandle;
    /* Open the file for reading */
    const int opening_mode = MPI_MODE_RDONLY;
    const int rv_open = MPI_File_open(MPI_COMM_SELF, file_name.c_str(), opening_mode, MPI_INFO_NULL, &fhandle);
    MPICheckError(rv_open);
    
    cmc_debug_msg("The file ", file_name, " has been opened.");

    /* Set the native data representation */
    const int rv_file_view = MPI_File_set_view(fhandle, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
    MPICheckError(rv_file_view);

    /* Read the first two values */
    std::array<uint64_t, 2> num_bytes;

    MPI_Status status;

    size_t offset{0};

    /* Move to the start of the file */
    const int rv_seek_start = MPI_File_seek(fhandle, 0, MPI_SEEK_SET);
    MPICheckError(rv_seek_start);

    /* Read the first two SizeTypes */
    const int rv_start_bytes = MPI_File_read(fhandle, num_bytes.data(), 2, MPI_UINT64_T, &status);
    MPICheckError(rv_start_bytes);
    CheckMPIReadCorrectness(&status, MPI_UINT64_T, 2);

    [[maybe_unused]] const SizeType num_global_bytes = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(num_bytes[0]);
    ++offset;
    const SizeType num_preamble_bytes = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(num_bytes[1]);
    ++offset;

    const SizeType num_preamble_values = num_preamble_bytes / sizeof(uint64_t);

    /* Move to the start of the file */
    const int rv_seek_start_again = MPI_File_seek(fhandle, 0, MPI_SEEK_SET);
    MPICheckError(rv_seek_start_again);

    /* Read the complete preamble */
    std::vector<uint64_t> preamble(num_preamble_values);
    const int rv_preamble = MPI_File_read(fhandle, preamble.data(), num_preamble_values, MPI_UINT64_T, &status);
    MPICheckError(rv_preamble);
    CheckMPIReadCorrectness(&status, MPI_UINT64_T, num_preamble_values);

    PrintCompressionInfo(preamble, file_name);

    /* Close the file */
    const int rv_close = MPI_File_close(&fhandle);
    MPICheckError(rv_close);

    cmc_debug_msg("The file ", file_name, " has been closed.");
}



constexpr int kMaxDecompressionLevelUndefined = -1;

/* Actual class definition of the decompression variable */
template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
class DecompressionVariableMultiData
{
public:
    DecompressionVariableMultiData() = delete;
    DecompressionVariableMultiData(const std::string& file_name, const MPI_Comm comm, const t8_cmesh_t cmesh, const t8_scheme* scheme);


    //void Decompress(const int max_decompression_level = kMaxDecompressionLevelUndefined);
    void Decompress();

    std::pair<t8_forest_t, std::vector<T>> GetDecompressedData();

private:
    void InquireCompressionInfo();
    const uint64_t* OpenSharedLevelDataWindow(const int level);
    void CloseSharedLevelDataWindow();
    t8_forest_t RepartitionInitMesh(t8_forest_t adapted_forest);
    void RepartitionData(t8_forest_t adapted_forest, t8_forest_t partitioned_forest);
    void DecodeRootLevelValues();
    void SetStartPositionForStreamDecoder(cmc::bits::StreamDecoder<SymbolType>& stream_decoder, const cmc::bits::vector_view level_encoding, const cmc::bits::vector_view lvl_mesh, const std::vector<LevelPartition>& level_partition_info, const SizeType mesh_offset, const SizeType num_local_elems);
    void DecodeMeshCompressionSteps();
    void SetStartPositionForIntraElemStreamDecoder(cmc::bits::StreamDecoder<SymbolType>& stream_decoder, const cmc::bits::vector_view level_encoding, const std::vector<LevelPartition>& level_partition_info, const SizeType mesh_offset, const SizeType num_local_elems);
    void DecodeIntraElementSteps();
    std::pair<t8_forest_t, std::vector<T>> RepartitionForRefinementIteration(t8_forest_t mesh, std::vector<T>& data);
    void CreateGhostLayerOnMesh();
    void ReconstructErrorMesh();

    /* File containing the compressed data */
    const std::string file_name_;

    /* Partition table in order to find the best-suited starting position */
    std::vector<LevelPartition> levelwise_partition_info_;

    /* The handle to the compressed file */
    MPI_File fhandle_;

    /* Infos extarcted from the preamble */
    CompressionInfoStruct compression_info_;

    /* The current AMR mesh */
    AmrMesh mesh_;

    /* The current data during the decompression */
    std::vector<T> data_;

    /* Indicator of the maximum desired decompression level */
    int max_decompression_level{kMaxDecompressionLevelUndefined};

    /* MPI communciators */
    const MPI_Comm comm_{MPI_COMM_NULL};
    int comm_rank_{0}, comm_size_{1};
    MPI_Comm shm_comm_{MPI_COMM_NULL};
    int shm_rank_{0}, shm_size_{1};

    /* Shared memory window */
    MPI_Win lvl_window_;

    /* A decoder */
    cmc::bits::StreamDecoder<SymbolType> stream_decoder_;

    /* The utilized error mesh */
    std::unique_ptr<ErrorMesh> error_mesh_;

    /* A step counter for the compression */
    int decompression_step_idx_{0};

    /* A flag whether the decompression has been carried out and the file been closed */
    bool is_already_decompressed_{false};
};

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
DecompressionVariableMultiData<T, DIM, N>::DecompressionVariableMultiData(const std::string& file_name, const MPI_Comm comm, const t8_cmesh_t cmesh, const t8_scheme *scheme)
: file_name_(file_name), comm_{comm}
{
    /* Check if the compressed file exists */
    if (const std::filesystem::path input_file_path(this->file_name_); not std::filesystem::exists(input_file_path))
    {
        throw std::invalid_argument("The compressed file does not exist!");
    }

    /* Check if an MPI Communicator is given */
    if (this->comm_ == MPI_COMM_NULL)
    {
        throw std::invalid_argument("The MPI_Communicator is NULL!");
    }

    /* Gather the MPI rank and size of comm */
    const int rv_rank_ = MPI_Comm_rank(this->comm_, &(this->comm_rank_));
    MPICheckError(rv_rank_);
    const int rv_size_ = MPI_Comm_size(this->comm_, &(this->comm_size_));
    MPICheckError(rv_size_);

    /* Split communicator into shared memory groups */
    const int rv_split_comm = MPI_Comm_split_type(this->comm_, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &(this->shm_comm_));
    MPICheckError(rv_split_comm);

    /* Get the rank within and the size of the shared communicator */
    const int rv_shm_rank_ = MPI_Comm_rank(this->shm_comm_, &(this->shm_rank_));
    MPICheckError(rv_shm_rank_);
    const int rv_shm_size_ = MPI_Comm_size(this->shm_comm_, &(this->shm_size_));
    MPICheckError(rv_shm_size_);

    /* Open the file on the shared memory communicator */
    const int rv_open = MPI_File_open(this->shm_comm_, this->file_name_.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &(this->fhandle_));
    MPICheckError(rv_open);

    /* Set the native data representation */
    const int rv_file_view = MPI_File_set_view(this->fhandle_, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
    MPICheckError(rv_file_view);

    /* Create a forest mesh from the given parameters */
    mesh_.SetMesh(t8_forest_new_uniform (cmesh, scheme, 0, 0, this->comm_));

    /* Repartition the mesh */
    t8_forest_t partitioned_forest = this->RepartitionInitMesh(mesh_.GetMesh());
    mesh_.SetMesh(partitioned_forest);
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
void
DecompressionVariableMultiData<T, DIM, N>::InquireCompressionInfo()
{
    /* Read the first two SizeTypes from the file */
    std::array<uint64_t, 2> num_bytes{};

    MPI_Status status;

    /* We read the beginning from the header only on the root rank */
    /* Move to the start of the file */
    const int rv_seek_start = MPI_File_seek(this->fhandle_, 0, MPI_SEEK_SET);
    MPICheckError(rv_seek_start);

    /* Read the first two SizeTypes */
    const int rv_start_bytes = MPI_File_read(this->fhandle_, num_bytes.data(), 2, MPI_UINT64_T, &status);
    MPICheckError(rv_start_bytes);
    CheckMPIReadCorrectness(&status, MPI_UINT64_T, 2);
    
    /* Read the second value from the file which gives the number of preamble bytes */
    const SizeType num_preamble_bytes = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(num_bytes[1]);

    const MPI_Aint stream_length_in_bytes = (this->shm_rank_ != kRootRank ? 0 : num_preamble_bytes);
    cmc_assert(stream_length_in_bytes % sizeof(uint64_t) == 0); //The stream length is a multiple of 64 bit
    const MPI_Aint stream_length = stream_length_in_bytes / sizeof(uint64_t);

    /* Allocate a shared window for the compression preamble */
    /* Since the data is relatively small, we allocate it from the root rank only and read it from this process only */
    const int rv_win_alloc = MPI_Win_allocate_shared(stream_length_in_bytes, sizeof(uint64_t), MPI_INFO_NULL, this->shm_comm_, &(this->compression_info_.shared_data_init_ptr), &(this->compression_info_.shared_data));
    MPICheckError(rv_win_alloc);

    if (this->shm_rank_ == kRootRank)
    {
        /* Move to the start of the file */
        const int rv_seek_start_again = MPI_File_seek(this->fhandle_, 0, MPI_SEEK_SET);
        MPICheckError(rv_seek_start_again);

        /* Read the complete preamble from the root rank of the shared communicator */
        const int rv_read_data = MPI_File_read(this->fhandle_, this->compression_info_.shared_data_init_ptr, stream_length, MPI_UINT64_T, &status);
        MPICheckError(rv_read_data);
        CheckMPIReadCorrectness(&status, MPI_UINT64_T, stream_length);
    }

    /* Next, we synchronize the window */
    const int rv_sync_win = MPI_Win_sync(this->compression_info_.shared_data);
    MPICheckError(rv_sync_win);

    /* Impose an additional barrier on the shared communicator */
    const int rv_shm_barrier = MPI_Barrier(this->shm_comm_);
    MPICheckError(rv_shm_barrier);

    /* Setup the decoding start for this level correctly */
    const SizeType num_preamble_vals = num_preamble_bytes / sizeof(uint64_t);
    this->compression_info_.shared_data_init_ptr -= (this->shm_rank_ != kRootRank ? num_preamble_vals : 0);

    /* Decode and store the compression info struct */
    DecodeCompressionInfoStruct(this->compression_info_);
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
const uint64_t*
DecompressionVariableMultiData<T, DIM, N>::OpenSharedLevelDataWindow(const int level)
{
    cmc_assert(this->compression_info_.global_level_bytes.size() > static_cast<size_t>(level) && level >= 0);

    /* Compute equal distributions for the whole data level */
    const SizeType num_global_lvl_vals = this->compression_info_.global_level_bytes[level] / sizeof(uint64_t);

    const SizeType proc_offset_val_stream = static_cast<SizeType>(((static_cast<double>(this->shm_rank_) * static_cast<long double>(num_global_lvl_vals)) / static_cast<double>(this->shm_size_)));
    const SizeType next_proc_offset_val_stream = static_cast<SizeType>(((static_cast<double>(this->shm_rank_ + 1) * static_cast<long double>(num_global_lvl_vals)) / static_cast<double>(this->shm_size_)));
    cmc_assert(proc_offset_val_stream <= next_proc_offset_val_stream);
    const int val_stream_length = static_cast<int>(next_proc_offset_val_stream - proc_offset_val_stream);
    const int byte_stream_length = val_stream_length * sizeof(uint64_t);

    /* Allocate a window on the shared memory communicators */
    uint64_t* shm_mem{nullptr};
    const int rv_win_alloc = MPI_Win_allocate_shared(static_cast<MPI_Aint>(byte_stream_length), sizeof(uint64_t), MPI_INFO_NULL, this->shm_comm_, &shm_mem, &(this->lvl_window_));
    MPICheckError(rv_win_alloc);

    /* Move to the correct position in the file for the current level on the current process */
    const MPI_Offset var_level_offset = std::accumulate(this->compression_info_.global_level_bytes.begin(), std::next(this->compression_info_.global_level_bytes.begin(), level), 0);
    const MPI_Offset file_offset = this->compression_info_.offset_start_encoding + var_level_offset + proc_offset_val_stream * sizeof(uint64_t);
    
    const int rv_seek_start = MPI_File_seek(this->fhandle_, file_offset, MPI_SEEK_SET);
    MPICheckError(rv_seek_start);

    cmc_debug_msg("OPenSharedWin FileOFf: ", file_offset, ", var_lvl_off: ", var_level_offset, ", proc_offset_val_stream: ",proc_offset_val_stream, ", next_proc_offset_val_stream: ", next_proc_offset_val_stream);

    /* Read the corresponding data */
    MPI_Status status;
    const int rv_read_data = MPI_File_read(this->fhandle_, shm_mem, val_stream_length, MPI_UINT64_T, &status);
    MPICheckError(rv_read_data);
    CheckMPIReadCorrectness(&status, MPI_UINT64_T, val_stream_length);

    /* Next, we synchronize the window */
    const int rv_sync_win = MPI_Win_sync(this->lvl_window_);
    MPICheckError(rv_sync_win);

    /* Impose an additional barrier on the shared communicator */
    const int rv_shm_barrier = MPI_Barrier(this->shm_comm_);
    MPICheckError(rv_shm_barrier);

    /* Setup the decoding start for this level correctly */
    const uint64_t* global_var_start = shm_mem - proc_offset_val_stream;

    return global_var_start;
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
void
DecompressionVariableMultiData<T, DIM, N>::CloseSharedLevelDataWindow()
{
    /* Impose a barrier to finish all outstanding reads */
    const int rv_shm_barrier = MPI_Barrier(this->shm_comm_);
    MPICheckError(rv_shm_barrier);

    /* Free the allocated window */
    const int rv_win_free = MPI_Win_free(&(this->lvl_window_));
    MPICheckError(rv_win_free);
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
t8_forest_t
DecompressionVariableMultiData<T, DIM, N>::RepartitionInitMesh(t8_forest_t adapted_forest)
{
    /* Allocate a forest */
    t8_forest_t partitioned_forest;
    t8_forest_init(&partitioned_forest);

    /* Partition the forest */
    const int partition_for_coarsening = 0;
    t8_forest_set_partition(partitioned_forest, adapted_forest, partition_for_coarsening);
    t8_forest_commit(partitioned_forest);

    return partitioned_forest;
}


template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
void
DecompressionVariableMultiData<T, DIM, N>::RepartitionData(t8_forest_t adapted_forest, t8_forest_t partitioned_forest)
{
    /* Create an sc_array_t wrapper of the variable's data */
    sc_array_t* in_data = sc_array_new_data (static_cast<void*>(this->data_.data()), sizeof(T), this->data_.size());

    /* Allocate an output vector for the partitioned data */
    std::vector<T> partitioned_data(t8_forest_get_local_num_leaf_elements(partitioned_forest));

    /* Create a wrapper for the freshly allocated partitioned data */
    sc_array_t* out_data = sc_array_new_data (static_cast<void*>(partitioned_data.data()), sizeof(T), partitioned_data.size());

    /* Partition the variables data */
    t8_forest_partition_data(adapted_forest, partitioned_forest, in_data, out_data);

    /* Destroy the array wrappers */
    sc_array_destroy(in_data);
    sc_array_destroy(out_data);

    /* Store the partitioned data */
    std::swap(this->data_, partitioned_data);
}

template<OneByteArithmeticType T>
inline std::vector<T>
GetRootLevelValuesFromView(cmc::bits::vector_view lvl_data_start_view, const t8_gloidx_t mesh_offset, const t8_locidx_t num_local_elems, const std::vector<LevelPartition>& level_partition)
{
    cmc_assert(level_partition.size() >= 1);
    constexpr int type_size = sizeof(T);

    /* Find the first relevant parition that is larger than the mesh offset */
    auto start_partition_iter = level_partition.begin();
    for (auto partition_iter = level_partition.begin(); partition_iter != level_partition.end(); ++partition_iter)
    {
        if (partition_iter->elem_offset > static_cast<SizeType>(mesh_offset))
        {
           break;
        } else
        {
            start_partition_iter = partition_iter;
        }
    }

    /* Set the level view accordingly to the offset */
    lvl_data_start_view.MoveToOffsetBitInStream(start_partition_iter->coding_byte_offset * cmc::bits::kCharBit);

    std::vector<T> data;
    data.reserve(num_local_elems);

    int num_reamining_elements{num_local_elems};

    int last_element_offset = start_partition_iter->elem_offset;
    for (auto partition_iter = std::next(start_partition_iter); partition_iter != level_partition.end(); ++partition_iter)
    {
        /* Compute the next contigupus value sequence and extract it */
        const int max_elems_to_extract = std::min<int>(static_cast<int>(partition_iter->elem_offset - last_element_offset), num_reamining_elements);

        for (int elem_idx{0}; elem_idx < max_elems_to_extract; ++elem_idx)
        {
            const OneByteResidualType uvalue = lvl_data_start_view.GetNextBitSequence<OneByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
            data.push_back(std::bit_cast<T>(uvalue));
        }

        last_element_offset = partition_iter->elem_offset;
        num_reamining_elements -= max_elems_to_extract;

        if (num_reamining_elements <= 0) {break;}

        /* Move to the next value sequence start */
        lvl_data_start_view.MoveToNextVectorValueStart();
    }

    /* If there are still elements missing, they lay contiguoulsy in memory at that point and can be extarcted */
    for (int elem_idx{0}; elem_idx < num_reamining_elements; ++elem_idx)
    {
        const OneByteResidualType uvalue = lvl_data_start_view.GetNextBitSequence<OneByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
        data.push_back(std::bit_cast<T>(uvalue));
    }

    return data;
}

template<TwoByteArithmeticType T>
inline std::vector<T>
GetRootLevelValuesFromView(cmc::bits::vector_view lvl_data_start_view, const t8_gloidx_t mesh_offset, const t8_locidx_t num_local_elems, const std::vector<LevelPartition>& level_partition)
{
    cmc_assert(level_partition.size() >= 1);
    constexpr int type_size = sizeof(T);

    /* Find the first relevant parition that is larger than the mesh offset */
    auto start_partition_iter = level_partition.begin();
    for (auto partition_iter = level_partition.begin(); partition_iter != level_partition.end(); ++partition_iter)
    {
        if (partition_iter->elem_offset > static_cast<SizeType>(mesh_offset))
        {
           break;
        } else
        {
            start_partition_iter = partition_iter;
        }
    }

    /* Set the level view accordingly to the offset */
    lvl_data_start_view.MoveToOffsetBitInStream(start_partition_iter->coding_byte_offset * cmc::bits::kCharBit);

    std::vector<T> data;
    data.reserve(num_local_elems);

    int num_reamining_elements{num_local_elems};

    int last_element_offset = start_partition_iter->elem_offset;
    for (auto partition_iter = std::next(start_partition_iter); partition_iter != level_partition.end(); ++partition_iter)
    {
        /* Compute the next contigupus value sequence and extract it */
        const int max_elems_to_extract = std::min<int>(static_cast<int>(partition_iter->elem_offset - last_element_offset), num_reamining_elements);

        for (int elem_idx{0}; elem_idx < max_elems_to_extract; ++elem_idx)
        {
            const TwoByteResidualType uvalue = lvl_data_start_view.GetNextBitSequence<TwoByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
            data.push_back(std::bit_cast<T>(uvalue));
        }

        last_element_offset = partition_iter->elem_offset;
        num_reamining_elements -= max_elems_to_extract;

        if (num_reamining_elements <= 0) {break;}

        /* Move to the next value sequence start */
        lvl_data_start_view.MoveToNextVectorValueStart();
    }

    /* If there are still elements missing, they lay contiguoulsy in memory at that point and can be extarcted */
    for (int elem_idx{0}; elem_idx < num_reamining_elements; ++elem_idx)
    {
        const TwoByteResidualType uvalue = lvl_data_start_view.GetNextBitSequence<TwoByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
        data.push_back(std::bit_cast<T>(uvalue));
    }

    return data;
}

template<FourByteArithmeticType T>
inline std::vector<T>
GetRootLevelValuesFromView(cmc::bits::vector_view lvl_data_start_view, const t8_gloidx_t mesh_offset, const t8_locidx_t num_local_elems, const std::vector<LevelPartition>& level_partition)
{
    cmc_assert(level_partition.size() >= 1);

    /* Find the first relevant parition that is larger than the mesh offset */
    auto start_partition_iter = level_partition.begin();
    for (auto partition_iter = level_partition.begin(); partition_iter != level_partition.end(); ++partition_iter)
    {
        if (partition_iter->elem_offset > static_cast<SizeType>(mesh_offset))
        {
           break;
        } else
        {
            start_partition_iter = partition_iter;
        }
    }

    /* Set the level view accordingly to the offset */
    lvl_data_start_view.MoveToOffsetBitInStream(start_partition_iter->coding_byte_offset * cmc::bits::kCharBit);

    std::vector<T> data;
    data.reserve(num_local_elems);

    int num_reamining_elements{num_local_elems};

    int last_element_offset = start_partition_iter->elem_offset;
    for (auto partition_iter = std::next(start_partition_iter); partition_iter != level_partition.end(); ++partition_iter)
    {
        /* Compute the next contiguous value sequence and extract it */
        const int max_elems_to_extract = std::min<int>(static_cast<int>(partition_iter->elem_offset - last_element_offset), num_reamining_elements);

        for (int elem_idx{0}; elem_idx < max_elems_to_extract; ++elem_idx)
        {
            const FourByteResidualType uvalue = lvl_data_start_view.GetNextBitSequence<FourByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
            data.push_back(std::bit_cast<T>(uvalue));
        }

        last_element_offset = partition_iter->elem_offset;
        num_reamining_elements -= max_elems_to_extract;

        if (num_reamining_elements <= 0) {break;}

        /* Move to the next value sequence start */
        lvl_data_start_view.MoveToNextVectorValueStart();
    }

    /* If there are still elements missing, they lay contiguoulsy in memory at that point and can be extarcted */
    for (int elem_idx{0}; elem_idx < num_reamining_elements; ++elem_idx)
    {
        const FourByteResidualType uvalue = lvl_data_start_view.GetNextBitSequence<FourByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
        data.push_back(std::bit_cast<T>(uvalue));
    }

    return data;
}

template<EightByteArithmeticType T>
inline std::vector<T>
GetRootLevelValuesFromView(cmc::bits::vector_view lvl_data_start_view, const t8_gloidx_t mesh_offset, const t8_locidx_t num_local_elems, const std::vector<LevelPartition>& level_partition)
{
    cmc_assert(level_partition.size() >= 1);
    constexpr int type_size = sizeof(T);

    /* Find the first relevant parition that is smaller than the mesh offset */
    auto start_partition_iter = level_partition.begin();
    for (auto partition_iter = level_partition.begin(); partition_iter != level_partition.end(); ++partition_iter)
    {
        if (partition_iter->elem_offset > static_cast<SizeType>(mesh_offset))
        {
           break;
        } else
        {
            start_partition_iter = partition_iter;
        }
    }

    /* Set the level view accordingly to the offset */
    lvl_data_start_view.MoveToOffsetBitInStream(start_partition_iter->coding_byte_offset * cmc::bits::kCharBit);

    std::vector<T> data;
    data.reserve(num_local_elems);

    int num_reamining_elements{num_local_elems};

    int last_element_offset = start_partition_iter->elem_offset;
    for (auto partition_iter = std::next(start_partition_iter); partition_iter != level_partition.end(); ++partition_iter)
    {
        /* Compute the next contigupus value sequence and extract it */
        const int max_elems_to_extract = std::min<int>(static_cast<int>(partition_iter->elem_offset - last_element_offset), num_reamining_elements);

        for (int elem_idx{0}; elem_idx < max_elems_to_extract; ++elem_idx)
        {
            const EightByteResidualType uvalue = lvl_data_start_view.GetNextBitSequence<EightByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
            data.push_back(std::bit_cast<T>(uvalue));
        }

        last_element_offset = partition_iter->elem_offset;
        num_reamining_elements -= max_elems_to_extract;

        if (num_reamining_elements <= 0) {break;}

        /* Move to the next value sequence start */
        lvl_data_start_view.MoveToNextVectorValueStart();
    }

    /* If there are still elements missing, they lay contiguoulsy in memory at that point and can be extracted */
    for (int elem_idx{0}; elem_idx < num_reamining_elements; ++elem_idx)
    {
        const EightByteResidualType uvalue = lvl_data_start_view.GetNextBitSequence<EightByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
        data.push_back(std::bit_cast<T>(uvalue));
    }

    return data;
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
void
DecompressionVariableMultiData<T, DIM, N>::DecodeRootLevelValues()
{
    /* Create a shared window on the root level */
    const int root_level{0};
    const uint64_t* shared_lvl_start_ptr = this->OpenSharedLevelDataWindow(root_level);

    /* Create a view on the data and extarct it, since there is no special encoding applied to root elvel values */
    cmc::bits::vector_view data_view(shared_lvl_start_ptr);

    //TODO: Partition the mesh here? 

    /* Get the global offset of the first element */
    const t8_gloidx_t mesh_offset = t8_forest_get_first_local_leaf_element_id(mesh_.GetMesh());

    /* Get the number of local elements */
    const t8_locidx_t num_local_elems = t8_forest_get_local_num_leaf_elements(mesh_.GetMesh());

    /* Get the partition table for the root level */
    const std::vector<LevelPartition> root_level_partition = this->compression_info_.GetLevelPartition(root_level);

    /* De-Serialize the values */
    this->data_ = GetRootLevelValuesFromView<T>(data_view, mesh_offset, num_local_elems, root_level_partition);
    cmc_assert(static_cast<t8_locidx_t>(this->data_.size()) == num_local_elems);

    if constexpr (kWriteVTKDecompressionStepData)
    {
        const std::string file_prefix = std::string("cmc_decompression_") + std::string(this->compression_info_.GetName()) + std::string("_0");
        WriteDataToVTK<T>(this->mesh_.GetMesh(), this->data_, std::string(this->compression_info_.GetName()), file_prefix);
    }

    /* Close the shared level window after all values have been extarcted */
    this->CloseSharedLevelDataWindow();

    /* Update the decompression count */
    ++(this->decompression_step_idx_);
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
void
DecompressionVariableMultiData<T, DIM, N>::SetStartPositionForStreamDecoder(cmc::bits::StreamDecoder<SymbolType>& stream_decoder, const cmc::bits::vector_view level_encoding, const cmc::bits::vector_view lvl_mesh_encoding, const std::vector<LevelPartition>& level_partition_info, const SizeType mesh_offset, const SizeType num_local_elems)
{
    /* We need to determine an coding offset for each process in the decompression communicator.
     * We know, the current (equal) partitioning of elements and the partitionining (elements, coding offset) from the compression communicator.
     * Additionally, we know which elements will be refined wihtin this compression step.
    */
    const SizeType next_rank_starting_pos = mesh_offset + num_local_elems;

    /* Iterate until we find the largest partition bound before the partiton start of the next rank */
    int part_idx{0};
    for (size_t idx{0}; idx < level_partition_info.size(); ++idx)
    {
        if (level_partition_info[idx].elem_offset < next_rank_starting_pos)
        {
            part_idx = idx;
        } else
        {
            break;
        }
    }

    /* Copy the view to determine the correct parallel offsets */
    cmc::bits::vector_view lvl_mesh = lvl_mesh_encoding;

    /* Set the mesh encoding correctly to the start of the offset */
    const SizeType mesh_start_offset = (level_partition_info[part_idx].elem_offset <= mesh_offset ? mesh_offset : level_partition_info[part_idx].elem_offset);
    lvl_mesh.MoveToOffsetBitInStream(mesh_start_offset);

    SizeType num_entropy_codes_correction{0};
    bool once_potential_incomplete_tree{true};

    /* The last rank does not need a correction */
    if (this->shm_rank_ < this->shm_size_ - 1)
    {
        /* Get the offset from the Partition Info */
        const SizeType partition_bound = level_partition_info[part_idx].elem_offset;

        /** Iterate the local elements until we arrive at the elem_offset and count from thereon the number of entropy codes up to the process end **/
        /* Get the number of local trees */
        const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees(this->mesh_.GetMesh());
        
        /* Iterate over the local trees */
        int num_elems_skipped{0};
        for (t8_locidx_t tree_idx{0}; tree_idx < num_local_trees; ++tree_idx)
        {
            /* Get the local number of elements in the tree */
            const t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (this->mesh_.GetMesh(), tree_idx);

            /* Check if we need to iterate through this tree or whether we can skip it completely */
            if (mesh_offset + num_elems_skipped + num_elements_in_tree <= partition_bound)
            {
                /* Skip the whole tree */
                num_elems_skipped += num_elements_in_tree;
                continue;
            }

            /* Get the corresponding tree class */
            const t8_eclass_t tree_class = t8_forest_get_tree_class (this->mesh_.GetMesh(), tree_idx);
            
            /* Get the scheme from the mesh */
            const t8_scheme_c* scheme = t8_forest_get_scheme (this->mesh_.GetMesh());

            t8_locidx_t tree_local_start_idx{0};
            if (once_potential_incomplete_tree && mesh_offset <= partition_bound)
            {
                /* Compute the tree_local start idx for the potential for an potential incomplete tree offset */
                tree_local_start_idx = partition_bound - mesh_offset - num_elems_skipped;
                once_potential_incomplete_tree = false;
            }

            /* Check whether the tree refines regularly */
            const bool refines_regular = not scheme->refines_irregular(tree_class);

            if (refines_regular)
            {
                if (num_elements_in_tree < 1) {continue;}

                /* Get the first element in the tree */
                const t8_element_t *element = t8_forest_get_leaf_element_in_tree (this->mesh_.GetMesh(), tree_idx, 0);
                /* Get the number of children elements an element refines to */
                const int num_children = scheme->element_get_num_children(tree_class, element);

                for (t8_locidx_t elem_idx{tree_local_start_idx}; elem_idx < num_elements_in_tree; ++elem_idx)
                {
                    /* Check whether this element will be refined during this iteration */
                    if (lvl_mesh.GetNextBit())
                    {
                        /* Add this amount of entropy codes to the counter */
                        num_entropy_codes_correction += num_children;
                    }
                }
            } else
            {
                /* If the tree refines irregular, we need to check the number of children elements for each element */
                /* We iterate through all elements and count the number of entropy codes that are stored */
                for (t8_locidx_t elem_idx{tree_local_start_idx}; elem_idx < num_elements_in_tree; ++elem_idx)
                {
                    /* Check whether this element will be refined during this iteration */
                    if (lvl_mesh.GetNextBit())
                    {
                        /* Get the element in the tree */
                        const t8_element_t *element = t8_forest_get_leaf_element_in_tree (this->mesh_.GetMesh(), tree_idx, elem_idx);

                        /* Get the number of elements this element refines to */
                        const int elem_num_children = scheme->element_get_num_children(tree_class, element);

                        /* Add this amount of entropy codes to the counter */
                        num_entropy_codes_correction += elem_num_children;
                    }
                }
            }

            /* Update the processed elements */
            num_elems_skipped += num_elements_in_tree;
        }
    }

    /* Now, we need to exchange the corrections */
    std::array<SizeType, 2> exchange_data{mesh_offset, num_entropy_codes_correction};

    /* Allocate an ouput vector locally */
    std::vector<SizeType> offset_array(2 * (this->comm_size_));

    /* Exchange the data */
    const int rv_all_gather_offset = MPI_Allgather(exchange_data.data(), 2, MPI_SIZE_TYPE, offset_array.data(), 2, MPI_SIZE_TYPE, this->comm_);
    MPICheckError(rv_all_gather_offset);

    /* Find first (non-local) partition bound */
    int first_part_idx{0};
    for (size_t idx{0}; idx < level_partition_info.size(); ++idx)
    {
        if (level_partition_info[idx].elem_offset > mesh_offset)
        {
           break;
        } else
        {
            first_part_idx = idx;
        }
    }
    
    /* Get the relevant partition bound offset */
    const SizeType first_partition_bound_offset = level_partition_info[first_part_idx].elem_offset;

    int rank_offset{0};
    /* Now, we need to find the first relevant rank for the entropy count correction */
    for (int rank_idx{0}; rank_idx <= this->comm_rank_; ++rank_idx)
    {
        if (offset_array[2 * rank_idx] > first_partition_bound_offset)
        {
            break;
        } else
        {
            rank_offset = rank_idx;
        }
    }

    /* Count the number of entropy codes we need to correct the compression partition bound */
    int num_entropy_codes_to_correct{0};
    for (int rank_idx{rank_offset}; rank_idx < this->comm_rank_; ++rank_idx)
    {
        num_entropy_codes_to_correct += offset_array[2 * rank_idx + 1];
    }

    /* Set the start to the first relevant partition bound */
    const SizeType level_entropy_offset = level_partition_info[first_part_idx].coding_byte_offset;

    /* Set the level stream decoder correctly */
    cmc::bits::vector_view adjusted_level_data = level_encoding;

    adjusted_level_data.MoveToOffsetBitInStream(level_entropy_offset * cmc::bits::kCharBit);

    /* Set the adjusted in the stream decoder */
    stream_decoder.StartDecoding(adjusted_level_data);

    /* And now, we need to skip the number of entropy codes to correct in order to move to the process-local start position */
    SizeType entropy_codes_skipped{0};

    for (int entropy_skip_idx{0}; entropy_skip_idx < num_entropy_codes_to_correct; ++entropy_skip_idx)
    {
        /* We are reading the next bit, if it is a refinement we discard the entropy codes and move on;
         * in case the element remains unchanged, we just move on with the next bit */

        /* If there will be a refinement, we are skipping entropy codes */
        /* Get the next symbol */
        const SymbolType symbol = std::invoke([&stream_decoder](){
            SymbolType next_symbol = stream_decoder.DecodeNextEntropySymbol();
            if (next_symbol == kProcessEndSymbol) [[unlikely]]
            {
                while (next_symbol == kProcessEndSymbol)
                {
                    stream_decoder.ApplyProcessEndSymbol64Bit();
                    next_symbol = stream_decoder.DecodeNextEntropySymbol();
                }
            }
            return next_symbol;
        });

        /* Check whether the value was predictable during the compression */
        if (symbol == kFlagUnpredictable)
        {
            /* In this case the entropy code is followed by the significant bits of the unpredicted value which we need to skip */
            stream_decoder.SkipNextBits(sizeof(T) * cmc::bits::kCharBit);
        }

        /* Update the discarded entropy codes */
        ++entropy_codes_skipped;
    }
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
struct RefinementIterationData
{
    RefinementIterationData(const std::span<T> current_data, cmc::bits::StreamDecoder<SymbolType>& stream_decoder_, cmc::bits::vector_view process_local_level_mesh_encoding_, const std::unique_ptr<ErrorMesh>& error_mesh_)
    : data(current_data), stream_decoder{stream_decoder_}, process_local_level_mesh_encoding{process_local_level_mesh_encoding_}, error_mesh{error_mesh_}
    {
        fine_level_data.reserve(current_data.size() * 2 * DIM + 1);
    }

    void LeaveElementUnchanged(const int local_idx);
    void PerformRefinement(const std::vector<cmc::par::lossy::idw::util::PointData<T, DIM>>& control_values, const std::vector<cmc::par::lossy::idw::util::Coordinate<DIM>>& eval_coords, const std::vector<float>& permitted_abs_error);
    bool WillNextElementBeRefined() {return process_local_level_mesh_encoding.GetNextBit();}
    T GetData(const int idx) const {return data[idx];}
    inline float GetPermittedAbsError(const t8_scheme_c* scheme, const int global_tree_id, const t8_eclass_t tree_class, const t8_element_t* element) const
    {
        return error_mesh->GetPermittedAbsError(scheme, global_tree_id, tree_class, element);
    }

    void PerformPrediction(const std::vector<T>& control_values, const std::vector<float>& permitted_abs_error);

    const std::span<T> data;
    cmc::bits::StreamDecoder<SymbolType>& stream_decoder;
    cmc::bits::vector_view process_local_level_mesh_encoding;
    const std::unique_ptr<ErrorMesh>& error_mesh;
    std::vector<T> fine_level_data;
};


template<OneByteArithmeticType T>
inline T 
GetUnpredictableValue(cmc::bits::StreamDecoder<SymbolType>& stream_decoder)
{
    const auto value = stream_decoder.GetNextBitSequence<OneByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
    return std::bit_cast<T>(value);
}

template<TwoByteArithmeticType T>
inline T 
GetUnpredictableValue(cmc::bits::StreamDecoder<SymbolType>& stream_decoder)
{
    const auto value = stream_decoder.GetNextBitSequence<TwoByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
    return std::bit_cast<T>(value);
}

template<FourByteArithmeticType T>
inline T 
GetUnpredictableValue(cmc::bits::StreamDecoder<SymbolType>& stream_decoder)
{
    const auto value = stream_decoder.GetNextBitSequence<FourByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
    return std::bit_cast<T>(value);
}

template<EightByteArithmeticType T>
inline T 
GetUnpredictableValue(cmc::bits::StreamDecoder<SymbolType>& stream_decoder)
{
    const auto value = stream_decoder.GetNextBitSequence<EightByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
    return std::bit_cast<T>(value);
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
RefinementIterationData<T, DIM>::PerformRefinement(const std::vector<cmc::par::lossy::idw::util::PointData<T, DIM>>& control_values, const std::vector<cmc::par::lossy::idw::util::Coordinate<DIM>>& eval_coords, const std::vector<float>& permitted_abs_error)
{
    const int num_elems = static_cast<int>(eval_coords.size());

    /* Re-Create the predictions vector for this family of elements */
    /** The first child element is always equal to the first control point **/
    std::vector<T> fam_predictions;
    fam_predictions.reserve(control_values.size());
    fam_predictions.push_back(control_values[0].data);

    /* The first child's prediction is equal to the family control value */
    for (int child_idx{1}; child_idx < num_elems; ++child_idx)
    {
        /* Compute the prediction for the refined children elements */
        fam_predictions.push_back(cmc::par::lossy::idw::util::ComputePrediction<T>(control_values, eval_coords[child_idx]));
    }

    /* Check all the entropy codes for this family of elements and potentially de-quantize them */
    for (int elem_idx{0}; elem_idx < num_elems; ++elem_idx)
    {
        /* Get the next symbol */
        const SymbolType symbol = std::invoke([this](){
            SymbolType next_symbol = stream_decoder.DecodeNextEntropySymbol();
            if (next_symbol == kProcessEndSymbol) [[unlikely]]
            {
                while (next_symbol == kProcessEndSymbol)
                {
                    stream_decoder.ApplyProcessEndSymbol64Bit();
                    next_symbol = stream_decoder.DecodeNextEntropySymbol();
                }
            }
            return next_symbol;
        });

        /* Check, if the value was un-predictable */
        if (symbol == kFlagUnpredictable) [[unlikely]]
        {
            /* Get the unpredicted value from the stream */
            const T unpredicted_value = GetUnpredictableValue<T>(this->stream_decoder);

            /* We replace the predicted value with the actual one */
            fam_predictions[elem_idx] = unpredicted_value;
        } else
        {
            /* Get the quantization bin from the entropy code */
            const auto [was_prediction_greater, quant_bin] = GetQuantizationBinFromEntropySymbol(symbol);

            /* Apply the de-quantization */
            fam_predictions[elem_idx] += ((was_prediction_greater ? -2.0 : +2.0) * permitted_abs_error[elem_idx] * quant_bin);
        }
    }

    /* We need to store the the de-quantized values for the next decompression iteration */
    std::copy_n(fam_predictions.begin(), num_elems, std::back_inserter(this->fine_level_data)); 
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
inline void
RefinementIterationData<T, DIM>::LeaveElementUnchanged(const int local_idx)
{
    this->fine_level_data.push_back(this->data[local_idx]);
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
inline t8_locidx_t
LossyMultiResDecompression (t8_forest_t forest,
                            t8_forest_t forest_from,
                            t8_locidx_t which_tree,
                            const t8_eclass_t tree_class,
                            t8_locidx_t lelement_id,
                            const t8_scheme_c * ts,
                            [[maybe_unused]] const int is_family,
                            [[maybe_unused]] const int num_elements,
                            t8_element_t * elements[])
{
    /* Retrieve the adapt_data */
    RefinementIterationData<T, DIM>* adapt_data = static_cast<RefinementIterationData<T, DIM>*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    /* Compute the start offset in the local contiguous array of the data*/
    const int local_start_index = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;

    /* Check if the current element will be refined/performs a prediction */
    if (adapt_data->WillNextElementBeRefined())
    {
        /** The element will be refined and the children elements need to be predicted **/

        /* Convert the local to a global tree id */
        const int global_tree_idx = t8_forest_global_tree_id (forest_from, which_tree);

        /* Get the number of children elements that need to be predicted */
        const int num_children = ts->element_get_num_children(tree_class, elements[0]);
        cmc_assert(num_children <= kNumMaxChildrenElements<DIM>);

        /* Allocate storage for the children elements */
        std::vector<t8_element_t*> child_elements(num_children, nullptr);
        t8_element_new (ts, tree_class, num_children, child_elements.data());

        /* Create the children */
        t8_element_get_children(ts, tree_class, elements[0], num_children, child_elements.data());

        /* Compute the evaluation coords of the constructed children elements */
        std::vector<cmc::par::lossy::idw::util::Coordinate<DIM>> eval_coords;
        eval_coords.reserve(child_elements.size());

        for (int child_idx{0}; child_idx < num_children; ++child_idx)
        {
            eval_coords.emplace_back(cmc::par::lossy::idw::util::FillElementIDWCoordinates<DIM>(forest_from, which_tree, child_elements[child_idx]));
        }

        /* Get the number of faces of this element */
        const int num_faces = ts->element_get_num_faces(tree_class, elements[0]);

        /* Create a vector for all face_values (and fill it with the default value (this families control value)) */
        std::vector<cmc::par::lossy::idw::util::PointData<T, DIM>> control_points;
        control_points.reserve(1 + num_faces);

        /* Create a point value for the family's control value */
        control_points.emplace_back();
        control_points.back().data = adapt_data->GetData(local_start_index);
        control_points.back().coordinates = eval_coords[0];

        /* Allocate some variables to be filled by the face neighbor calls */
        const t8_element_t** neighbor_leaves;
        int* dual_faces;
        int num_neighbors{0};
        t8_locidx_t* neighbor_element_indices;
        t8_eclass_t neighbor_tree_class;
        t8_gloidx_t gneigh_tree;
        int orientation;

        /* Iterate over all faces and fill the control points */
        for (int face_idx{0}; face_idx < num_faces; ++face_idx)
        {
            /* Gather the face neighbor via this face */
            //t8_forest_leaf_face_neighbors (forest_from, which_tree, elements[0], &neighbor_leaves, face_idx, &dual_faces, &num_neighbors,
            //    &neighbor_element_indices, &neighbor_tree_class);
            t8_forest_leaf_face_neighbors_ext (forest_from, which_tree, elements[0], &neighbor_leaves, face_idx, &dual_faces, &num_neighbors,
                                               &neighbor_element_indices, &neighbor_tree_class, &gneigh_tree, &orientation);

            /* Get the coordinates from the control points */
            for (int neigh_idx{0}; neigh_idx < num_neighbors; ++neigh_idx)
            {
                const t8_locidx_t neigh_tree_idx = t8_forest_get_local_or_ghost_id(forest_from, gneigh_tree);
                cmc_assert(neigh_tree_idx >= 0);

                /* Add a new control value via this face */
                control_points.emplace_back();
                control_points.back().data = adapt_data->GetData(neighbor_element_indices[neigh_idx]);
                control_points.back().coordinates = cmc::par::lossy::idw::util::FillElementIDWCoordinates<DIM>(forest_from, neigh_tree_idx, neighbor_leaves[neigh_idx]);
            }

            /* Deallocate the memory for the face neighbor construction */
            if (num_neighbors > 0) {
                T8_FREE (neighbor_leaves);
                T8_FREE (neighbor_element_indices);
                T8_FREE (dual_faces);
            }
        }

        /* We need to get the permitted errors for the children elements that will be constructed */
        std::vector<float> permitted_abs_errors;
        permitted_abs_errors.reserve(num_children);
        
        /* Get the permitted errors for the child elements */
        for (int child_idx{0}; child_idx < num_children; ++child_idx)
        {
            permitted_abs_errors.emplace_back(adapt_data->GetPermittedAbsError(ts, global_tree_idx, tree_class, child_elements[child_idx]));
        }

        /* Destroy the constructed elements */
        ts->element_destroy(tree_class, num_children, child_elements.data());

        /* Perform the prediction and get the data for encoding */
        adapt_data->PerformRefinement(control_points, eval_coords, permitted_abs_errors);

        return cmc::t8::kRefineElement;
    } else
    {
        /**  The element remains unchanged, therefore no prediction needs to be made **/
        adapt_data->LeaveElementUnchanged(local_start_index);
        
        #ifdef CMC_DEBUG
        /* Compare whether the permitted error still holds */

        #endif
        
        return cmc::t8::kLeaveElementUnchanged;
    }
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
void
DecompressionVariableMultiData<T, DIM, N>::CreateGhostLayerOnMesh()
{
    t8_forest_t mesh = this->mesh_.GetMesh();

    const t8_locidx_t current_local_elems = t8_forest_get_local_num_leaf_elements(mesh);

    /* Allocate a new forest */
    t8_forest_t ghost_mesh;
    t8_forest_init (&ghost_mesh);

    /* Set the base forest mesh */
    t8_forest_set_copy(ghost_mesh, mesh);

    /* Set the forest for creating a face ghost layer */
    t8_forest_set_ghost (ghost_mesh, 1, T8_GHOST_FACES);

    /* Commit the forest, this step will perform the partitioning and ghost layer creation. */
    t8_forest_commit (ghost_mesh);

    /* Get the number of ghost elements of forest. */
    const t8_locidx_t num_ghost_elements = t8_forest_get_num_ghosts (ghost_mesh);

    /* Allocate a vector holding the local data and the ghost data */
    std::vector<T> local_and_ghost_values(current_local_elems + num_ghost_elements);
    std::copy_n(this->data_.begin(), current_local_elems, local_and_ghost_values.begin());

    /* Create a wrapper for the ghost exchange */
    sc_array_t* ghost_exchange_data = sc_array_new_data (local_and_ghost_values.data(), sizeof(T), current_local_elems + num_ghost_elements);

    /* Exchange the ghost data for the newly partitioned data */
    t8_forest_ghost_exchange_data (ghost_mesh, ghost_exchange_data);

    /* Destroy the array wrapper */
    sc_array_destroy(ghost_exchange_data);

    /* Set the mesh and the data */
    this->mesh_.SetMesh(ghost_mesh);
    this->data_ = std::move(local_and_ghost_values);
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
std::pair<t8_forest_t, std::vector<T>>
DecompressionVariableMultiData<T, DIM, N>::RepartitionForRefinementIteration(t8_forest_t mesh, std::vector<T>& data)
{
    const t8_locidx_t current_local_elems = t8_forest_get_local_num_leaf_elements(mesh);
    const t8_locidx_t current_ghost_elems = t8_forest_get_num_ghosts (mesh);

    cmc_assert(current_local_elems + current_ghost_elems == static_cast<t8_locidx_t>(data.size()));

    /* Keep the not-partitioned forest */
    t8_forest_ref(mesh);

    /** Partition the forest correctly and build a halo layer **/
    t8_forest_t partitioned_ghost_mesh;
    t8_forest_init (&partitioned_ghost_mesh);

    /* Set the forest for partitioning */
    t8_forest_set_partition (partitioned_ghost_mesh, mesh, 0);

    /* Set the forest for creating a face ghost layer */
    t8_forest_set_ghost (partitioned_ghost_mesh, 1, T8_GHOST_FACES);

    /* Commit the forest, this step will perform the partitioning and ghost layer creation. */
    t8_forest_commit (partitioned_ghost_mesh);

    /** Exchange the data corectly to set up the halo layer **/
    /* Get the number of local elements */
    const t8_locidx_t num_local_elements = t8_forest_get_local_num_leaf_elements(partitioned_ghost_mesh);
    
    /* Get the number of ghost elements of forest. */
    const t8_locidx_t num_ghost_elements = t8_forest_get_num_ghosts (partitioned_ghost_mesh);

    /* Create an sc_array_t wrapper of the variable's local element data */
    sc_array_t* in_data = sc_array_new_data (static_cast<void*>(data.data()), sizeof(T), current_local_elems);

    /* Allocate an output vector for the partitioned data */
    std::vector<T> partitioned_ghost_data(num_local_elements + num_ghost_elements);

    /* Create a wrapper for the freshly allocated partitioned data */
    sc_array_t* out_data = sc_array_new_data (static_cast<void*>(partitioned_ghost_data.data()), sizeof(T), num_local_elements);

    /* Partition the variables local element data */
    t8_forest_partition_data(mesh, partitioned_ghost_mesh, in_data, out_data);

    /* Destroy the array wrappers */
    sc_array_destroy(in_data);
    sc_array_destroy(out_data);

    /* Free the former forest */
    t8_forest_unref(&mesh);

    /* Create a wrapper for the ghost exchange */
    sc_array_t* ghost_exchange_data = sc_array_new_data (partitioned_ghost_data.data(), sizeof(T), num_local_elements + num_ghost_elements);

    /* Exchange the ghost data for the newly partitioned data */
    t8_forest_ghost_exchange_data (partitioned_ghost_mesh, ghost_exchange_data);

    /* Destroy the array wrapper */
    sc_array_destroy(ghost_exchange_data);

    return std::make_pair(partitioned_ghost_mesh, std::move(partitioned_ghost_data));
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
void
DecompressionVariableMultiData<T, DIM, N>::DecodeMeshCompressionSteps()
{
    /* Set up the Huffman decoder within the stream decoder */
    this->stream_decoder_.StartHuffmanCodesDecoding(this->compression_info_.GetMeshCompressionHuffmanCodes());

    const int num_mesh_compression_levels = this->compression_info_.mesh_compression_levels;

    for (int mesh_lvl{1}; mesh_lvl < num_mesh_compression_levels; ++mesh_lvl)
    {
        cmc_debug_msg("The mesh decompression step ", mesh_lvl," starts...");
        /****** Get the correct view on this level for each process ******/
        /* Get the Level partition information */
        const std::vector<LevelPartition> level_partition = this->compression_info_.GetLevelPartition(mesh_lvl);

        /* Get the view on the mesh encoding */
        cmc::bits::vector_view level_mesh_encoding = this->compression_info_.GetGlobalElementIndicationsStep(mesh_lvl - 1);
        /* Open shared data window on this levels encoding */
        const uint64_t* shared_lvl_start_ptr = this->OpenSharedLevelDataWindow(mesh_lvl);
        /* Create an encoded level data view */
        cmc::bits::vector_view level_data_encoding(shared_lvl_start_ptr);

        /* Get the global offset of the first element */
        const SizeType mesh_offset = (SizeType) t8_forest_get_first_local_leaf_element_id(this->mesh_.GetMesh());
        const SizeType num_local_elems = (SizeType) t8_forest_get_local_num_leaf_elements(this->mesh_.GetMesh());

        /* Move the stream decoder to the correct starting position */
        this->SetStartPositionForStreamDecoder(this->stream_decoder_, level_data_encoding, level_mesh_encoding, level_partition, mesh_offset, num_local_elems);

        /* Set the process local mesh encoding start */
        cmc::bits::vector_view process_local_level_mesh_encoding = level_mesh_encoding;
        process_local_level_mesh_encoding.MoveToOffsetBitInStream(mesh_offset);
    
        /****** Perform the local refinement onto the next finer level as dictated by the mesh encoding ******/
        /* Crerate the adaptation data */
        RefinementIterationData<T, DIM> adapt_data(std::span<T>(this->data_), this->stream_decoder_, process_local_level_mesh_encoding, this->error_mesh_);

        /* Perform a decompression iteration */
        t8_forest_t adapted_forest = t8_forest_new_adapt(mesh_.GetMesh(), LossyMultiResDecompression<T, DIM>, 0, 0, static_cast<void*>(&adapt_data));
        cmc_debug_msg("The mesh adaptation step is finished; resulting in ", t8_forest_get_global_num_leaf_elements(adapted_forest), " global elements");

        /****** Partition the mesh and the data ******/

        /* Repartition the mesh and the data */
        auto [partitioned_mesh, partitioned_data] = this->RepartitionForRefinementIteration(adapted_forest, adapt_data.fine_level_data);

        /* Store the partitioned mesh and the data */
        this->mesh_.SetMesh(partitioned_mesh);
        this->data_ = std::move(partitioned_data);

        if constexpr (kWriteVTKDecompressionStepData)
        {
            const std::string file_prefix = std::string("cmc_decompression_") + std::string(this->compression_info_.GetName()) + std::string("_") + std::to_string(mesh_lvl);
            WriteDataToVTK<T>(this->mesh_.GetMesh(), this->data_, std::string(this->compression_info_.GetName()), file_prefix);
        }

        /****** Clean-Up of the decompression iteration ******/

        /* Close this levels data window */
        this->CloseSharedLevelDataWindow();

        /* Update the decompression count */
        ++(this->decompression_step_idx_);
        cmc_debug_msg("The mesh decompression step ", mesh_lvl," has been completed.");
    }
}


template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
void
DecompressionVariableMultiData<T, DIM, N>::SetStartPositionForIntraElemStreamDecoder(cmc::bits::StreamDecoder<SymbolType>& stream_decoder, const cmc::bits::vector_view level_encoding, const std::vector<LevelPartition>& level_partition_info, const SizeType mesh_offset, [[maybe_unused]] const SizeType num_local_elems)
{   
    /* Find first (potential non-local) partition bound */
    int first_part_idx{0};
    for (size_t idx{0}; idx < level_partition_info.size(); ++idx)
    {
        if (level_partition_info[idx].elem_offset > mesh_offset)
        {
           break;
        } else
        {
            first_part_idx = idx;
        }
    }
    
    /* Get the relevant partition bound offset */
    const SizeType first_partition_bound_offset = level_partition_info[first_part_idx].elem_offset;

    /* Set the start to the first relevant partition bound */
    const SizeType level_entropy_offset = level_partition_info[first_part_idx].coding_byte_offset;

    /* Set the level stream decoder correctly */
    cmc::bits::vector_view adjusted_level_data = level_encoding;

    adjusted_level_data.MoveToOffsetBitInStream(level_entropy_offset * cmc::bits::kCharBit);

    /* Set the adjusted in the stream decoder */
    stream_decoder.StartDecoding(adjusted_level_data);

    /* And now, we need to skip the number of entropy codes for a fixed amount of elements */
    const int num_elems_to_skip = mesh_offset - first_partition_bound_offset;

    for (int elem_idx{0}; elem_idx < num_elems_to_skip; ++elem_idx)
    {
        /* Skip this elements significant bits in the encoded stream */
        SkipToNextCompressedElement<T, DIM, N>(stream_decoder);
    }
}

template<FourByteArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
struct RefinementIntraElementData
{
    RefinementIntraElementData(const std::span<T> current_data, cmc::bits::StreamDecoder<SymbolType>& stream_decoder_)
    : data(current_data), stream_decoder{stream_decoder_}
    {
        static_assert(N >= 1);

        fine_elem_data.reserve(current_data.size() * N);
    }

    void PerformRefinement(const int coarse_value_id);

    const std::span<T> data;
    cmc::bits::StreamDecoder<SymbolType>& stream_decoder;
    std::vector<T> fine_elem_data;
};

template<FourByteArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
inline void
RefinementIntraElementData<T, DIM, N>::PerformRefinement(const int coarse_value_id)
{
    /* Perform the element decoding */
    const std::array<T, N> elem_data = PerformElementDecoding<float, DIM, N>(this->stream_decoder, data[coarse_value_id]);

    /* Copy the data to the fine value output */
    std::copy_n(elem_data.begin(), N, std::back_inserter(fine_elem_data));
}


template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
void
DecompressionVariableMultiData<T, DIM, N>::DecodeIntraElementSteps()
{
    cmc_assert(this->compression_info_.intra_elem_compression_levels == 1);
    cmc_err_msg("Not supported yet");
    /* Set up the Huffman decoder within the stream decoder */
    this->stream_decoder_.StartHuffmanCodesDecoding(this->compression_info_.GetIntraElementCompressionHuffmanCodes());

    const int num_mesh_compression_levels = this->compression_info_.mesh_compression_levels;
    const int num_intra_elem_compression_levels = this->compression_info_.intra_elem_compression_levels > 0 ? 1 : 0;

    for (int intra_lvl{0}; intra_lvl < num_intra_elem_compression_levels; ++intra_lvl)
    {
        cmc_debug_msg("The intra element decompression step ", num_mesh_compression_levels + intra_lvl," starts...");
        /****** Get the correct view on this level for each process ******/
        /* Get the Level partition information */
        const std::vector<LevelPartition> intra_level_data_partition = this->compression_info_.GetLevelPartition(num_mesh_compression_levels + intra_lvl);

        /* Get the global offset of the first element */
        const SizeType mesh_offset = (SizeType) t8_forest_get_first_local_leaf_element_id(mesh_.GetMesh());
        const SizeType num_local_elems = (SizeType) t8_forest_get_local_num_leaf_elements(mesh_.GetMesh());

        /* Open shared data window on this levels encoding */
        const uint64_t* shared_intra_lvl_start_ptr = this->OpenSharedLevelDataWindow(num_mesh_compression_levels + intra_lvl);
        /* Create an encoded level data view */
        cmc::bits::vector_view intra_level_data_encoding(shared_intra_lvl_start_ptr);

        /* Setup intra element decompression view */
        this->SetStartPositionForIntraElemStreamDecoder(this->stream_decoder_, intra_level_data_encoding, intra_level_data_partition, mesh_offset, num_local_elems);
        
        /* Setup the decompressor */
        RefinementIntraElementData<T, DIM, N> intra_level_decompressor(std::span<T>(this->data_), this->stream_decoder_);

        /* Iterate over all local elements and decompress them */
        for (int elem_idx{0}; elem_idx < num_local_elems; ++elem_idx)
        {
            intra_level_decompressor.PerformRefinement(elem_idx);
        }

        /* Store the output data */
        this->data_ = std::move(intra_level_decompressor.fine_elem_data);

        /* CLose this levels data window */
        this->CloseSharedLevelDataWindow();

        /* Update the decompression count */
        ++(this->decompression_step_idx_);
        cmc_debug_msg("The intra element decompression step ", num_mesh_compression_levels + intra_lvl, " has been completed.");
    }
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
void
DecompressionVariableMultiData<T, DIM, N>::ReconstructErrorMesh()
{
    /* The error mesh is stored directly before the encoding start, therfore, we can just move the given bytes back */
    const MPI_Aint file_error_mesh_offset = this->compression_info_.offset_start_encoding - this->compression_info_.num_bytes_global_error_mesh_encoding;
    cmc_assert(file_error_mesh_offset % sizeof(uint64_t) == 0);

    /* Create a shared memory window */
    MPI_Win error_mesh_window;

    /* Compute equal distributions for the encoded error mesh */
    const SizeType num_global_vals = this->compression_info_.num_bytes_global_error_mesh_encoding / sizeof(uint64_t);
    cmc_assert(this->compression_info_.num_bytes_global_error_mesh_encoding % sizeof(uint64_t) == 0);

    const SizeType proc_offset_val_stream = static_cast<SizeType>(((static_cast<double>(this->shm_rank_) * static_cast<long double>(num_global_vals)) / static_cast<double>(this->shm_size_)));
    const SizeType next_proc_offset_val_stream = static_cast<SizeType>(((static_cast<double>(this->shm_rank_ + 1) * static_cast<long double>(num_global_vals)) / static_cast<double>(this->shm_size_)));
    cmc_assert(proc_offset_val_stream <= next_proc_offset_val_stream);
    const int val_stream_length = static_cast<int>(next_proc_offset_val_stream - proc_offset_val_stream);
    const int byte_stream_length = val_stream_length * sizeof(uint64_t);

    /* Allocate a window on the shared memory communicators */
    uint64_t* shm_mem{nullptr};
    const int rv_win_alloc = MPI_Win_allocate_shared(static_cast<MPI_Aint>(byte_stream_length), sizeof(uint64_t), MPI_INFO_NULL, this->shm_comm_, &shm_mem, &error_mesh_window);
    MPICheckError(rv_win_alloc);

    /* Move to the correct position in the file for the current level on the current process */
    const MPI_Offset file_offset = file_error_mesh_offset + proc_offset_val_stream * sizeof(uint64_t);
    
    const int rv_seek_start = MPI_File_seek(this->fhandle_, file_offset, MPI_SEEK_SET);
    MPICheckError(rv_seek_start);

    /* Read the corresponding data */
    MPI_Status status;
    const int rv_read_data = MPI_File_read(this->fhandle_, shm_mem, val_stream_length, MPI_UINT64_T, &status);
    MPICheckError(rv_read_data);
    CheckMPIReadCorrectness(&status, MPI_UINT64_T, val_stream_length);

    /* Next, we synchronize the window */
    const int rv_sync_win = MPI_Win_sync(error_mesh_window);
    MPICheckError(rv_sync_win);

    /* Impose an additional barrier on the shared communicator */
    const int rv_shm_barrier = MPI_Barrier(this->shm_comm_);
    MPICheckError(rv_shm_barrier);

    /* Setup the decoding start for the error mesh correctly */
    const uint64_t* error_mesh_start = shm_mem - proc_offset_val_stream;

    /* Reconstruct the error mesh */
    this->error_mesh_ = std::make_unique<ErrorMesh>(this->mesh_.GetMesh(), error_mesh_start, this->shm_comm_, this->shm_rank_, this->shm_size_);

    /* Delete the window */
    const int rv_win_free = MPI_Win_free(&error_mesh_window);
    MPICheckError(rv_win_free);
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
void
DecompressionVariableMultiData<T, DIM, N>::Decompress()
{
    if (this->is_already_decompressed_) [[unlikely]]
    {
        cmc_err_msg("The variable has already been decompressed!");
    }

    /* Inquire the basic information about the compression */
    this->InquireCompressionInfo();

    /* Reconstruct the error mesh */
    this->ReconstructErrorMesh();

    /* Decode the root level */
    this->DecodeRootLevelValues();

    /* Create a ghost layer on the initial forest with the decompressed root level data */
    this->CreateGhostLayerOnMesh();

    /* Perform the level-wise iterative decompression */
    this->DecodeMeshCompressionSteps();

    /* Perform the intra-element decompression */
    if constexpr (N > 1)
    {
        this->DecodeIntraElementSteps();
    }

    /* Close shared window on compression infos */
    /* Impose a barrier to finish all outstanding reads */
    const int rv_shm_barrier = MPI_Barrier(this->shm_comm_);
    MPICheckError(rv_shm_barrier);

    /* Free the allocated window */
    const int rv_win_free = MPI_Win_free(&(this->compression_info_.shared_data));
    MPICheckError(rv_win_free);

    /* Free the error mesh */
    cmc_assert(this->error_mesh_ != nullptr);
    this->error_mesh_->DestructErrorMeshCollectively();

    /* Close the compressed file */
    const int rv_close = MPI_File_close(&this->fhandle_);
    MPICheckError(rv_close);

    /* Set the flag that the data has already been decompressed */
    this->is_already_decompressed_ = true;

    cmc_debug_msg("Decompression of variable ", reinterpret_cast<const char*>(this->compression_info_.name.data()), " has been completed.");
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
[[nodiscard("The decompressed data will be lost if the return value is discarded!")]]
std::pair<t8_forest_t, std::vector<T>>
DecompressionVariableMultiData<T, DIM, N>::GetDecompressedData()
{
    t8_forest_t mesh = mesh_.GetMesh();
    mesh_.SetNullMesh();
    std::vector<T> decompressed_data{};
    std::swap(decompressed_data, this->data_);
    return std::make_pair(mesh, std::move(decompressed_data));
}

}
#endif /* !CMC_LOSSY_AMR_PAR_MULTI_RES_DECOMPRESSION_IDW_HXX */
