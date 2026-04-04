#ifndef CMC_PAR_MULTI_RES_DECOMPRESSION_HXX
#define CMC_PAR_MULTI_RES_DECOMPRESSION_HXX

#include "cmc.hxx"
#include "amr/lossless/cmc_par_multi_res_extraction_util.hxx"
#include "mpi/cmc_mpi.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "utilities/cmc_compression_schema.hxx"
#include "utilities/cmc_bits_stream_decoder.hxx"

#include <t8_forest/t8_forest_partition.h>

#include <string>
#include <span>
#include <filesystem>


namespace cmc::par::lossless::multi_res
{

    static bool mesh_enc_once = true;

struct CompressionInfoStruct
{
    cmc::bits::vector_view
    GetGlobalElementIndicationsStep(const int step)
    {
        cmc_assert(step >= 0 && step < this->mesh_compression_levels - 1);
        cmc_assert((this->GetOffsetPartitionTable() + this->partition_table_size + this->mesh_compression_huffman_codes + this->intra_elem_compression_huffman_codes) % sizeof(uint64_t) == 0);

        int offset = (this->GetOffsetPartitionTable() + this->partition_table_size + this->mesh_compression_huffman_codes + this->intra_elem_compression_huffman_codes) / sizeof(uint64_t);
        for (int step_idx{0}; step_idx < step; ++step_idx)
        {
            offset += global_num_elem_indications_step[step_idx] / (sizeof(uint64_t) * cmc::bits::kCharBit) + (global_num_elem_indications_step[step_idx] % (sizeof(uint64_t) * cmc::bits::kCharBit) != 0 ? 1 : 0);
        }
        cmc_global_msg("Offset for Elem Indications step ", step, ": ", offset);
        if (mesh_enc_once)
        {
        cmc_global_msg("Complete Global Mesh Encoding:");
        for (int idx{0}; idx < this->num_bytes_global_mesh_encoding  /sizeof(uint64_t); ++idx)
        {
            cmc_global_msg("Bitset: ", std::bitset<64>(*(shared_data_init_ptr + offset + idx)));
        }
        mesh_enc_once = false;
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
        cmc_global_msg("LevelPartition Level ", level, " offset: ", offset);
        for(int rank_idx{0}; rank_idx < this->compression_comm_size; ++rank_idx)
        {
            level_partition.emplace_back(cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(this->shared_data_init_ptr + offset + 2 * rank_idx)),
                                         cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(this->shared_data_init_ptr + offset + 2 * rank_idx + 1)));
            cmc_global_msg("rank_idx: ", rank_idx, ", level_partition.elem_offset: ", level_partition.back().elem_offset, ", level_partition.coding_byte_offset: ", level_partition.back().coding_byte_offset);
        }
        return level_partition;
    }

    const uint64_t*
    GetMeshCompressionHuffmanCodes() const
    {
        cmc_assert((this->GetOffsetPartitionTable() + this->partition_table_size) % sizeof(uint64_t) == 0);
        const int offset = (this->GetOffsetPartitionTable() + this->partition_table_size) / sizeof(uint64_t);
        cmc_global_msg("Offset mesh huffman codes: ", offset);
        return shared_data_init_ptr + offset;
    }

    const uint64_t*
    GetIntraElementCompressionHuffmanCodes() const
    {
        cmc_assert((this->GetOffsetPartitionTable() + this->partition_table_size + this->mesh_compression_huffman_codes) % sizeof(uint64_t) == 0);
        const int offset = (this->GetOffsetPartitionTable() + this->partition_table_size + this->mesh_compression_huffman_codes) / sizeof(uint64_t);
        cmc_global_msg("Offset intra elem huffman codes: ", offset);
        return shared_data_init_ptr + offset;
    }

    SizeType GetOffsetPartitionTable() const {return sizeof(SizeType) * (14 + global_level_bytes.size() + global_num_elem_indications_step.size())
                                                     + kNumCharsVariableName * sizeof(char);}

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

//TODO: Update to uint64_t pointer
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
    cmc_global_msg("Received info.global_byte_count: ", info.global_byte_count);

    info.offset_start_encoding = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;
    cmc_global_msg("Received info.offset_start_encoding: ", info.offset_start_encoding);

    std::copy_n(reinterpret_cast<const char*>(start_ptr + offset), kNumCharsVariableName, info.name.data());
    offset += 32;
    cmc_global_msg("Received info.name: ", std::string(reinterpret_cast<const char*>(start_ptr + offset - 32)));

    info.data_type = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;
    cmc_global_msg("Received info.data_type: ", info.data_type);

    info.dimensionality = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;
    cmc_global_msg("Received info.dimensionality: ", info.dimensionality);

    info.data_per_elem = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;
    cmc_global_msg("Received info.data_per_elem: ", info.data_per_elem);

    info.compression_comm_size = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;
    cmc_global_msg("Received info.compression_comm_size: ", info.compression_comm_size);

    info.compression_scheme = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;
    cmc_global_msg("Received info.compression_scheme: ", info.compression_scheme);

    info.mesh_compression_levels = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;
    cmc_global_msg("Received info.mesh_compression_levels: ", info.mesh_compression_levels);

    info.intra_elem_compression_levels = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;
    cmc_global_msg("Received info.intra_elem_compression_levels: ", info.intra_elem_compression_levels);

    info.pack_size = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;
    cmc_global_msg("Received info.pack_size: ", info.pack_size);

    const int num_levels = info.mesh_compression_levels + (info.intra_elem_compression_levels > 0 ? 1 : 0);
    info.global_level_bytes.reserve(num_levels);
    for (int lvl_idx{0}; lvl_idx < num_levels; ++lvl_idx)
    {
        info.global_level_bytes.push_back(cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset)));
        ++offset;
        cmc_global_msg("Received info.global_level_bytes.back(): ", info.global_level_bytes.back());
    }

    cmc_global_msg("Received info.mesh_compression_levels: ", info.mesh_compression_levels);
    info.global_num_elem_indications_step.reserve(info.mesh_compression_levels - 1);
    for (int lvl_idx{0}; lvl_idx < info.mesh_compression_levels - 1; ++lvl_idx)
    {
        info.global_num_elem_indications_step.push_back(cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset)));
        ++offset;
        cmc_global_msg("Received info.global_num_elem_indications_step.back(): ", info.global_num_elem_indications_step.back());
    }

    info.partition_table_size = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;
    cmc_global_msg("Received info.partition_table_size: ", info.partition_table_size);

    info.mesh_compression_huffman_codes = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;
    cmc_global_msg("Received info.mesh_compression_huffman_codes: ", info.mesh_compression_huffman_codes);

    info.intra_elem_compression_huffman_codes = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    ++offset;
    cmc_global_msg("Received info.intra_elem_compression_huffman_codes: ", info.intra_elem_compression_huffman_codes);

    info.num_bytes_global_mesh_encoding = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(start_ptr + offset));
    //++offset;
    cmc_global_msg("Received info.num_bytes_global_mesh_encoding: ", info.num_bytes_global_mesh_encoding);

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
    for (int lvl_idx{0}; lvl_idx < encoded_info.mesh_compression_levels - 1; ++lvl_idx)
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
    //++offset;
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
    for (int lvl_idx{0}; lvl_idx < info.mesh_compression_levels; ++lvl_idx)
    {
        cmc_global_msg("\t\t Level ", lvl_idx, ": ", info.global_level_bytes[lvl_idx], " bytes");
    }

    if (info.intra_elem_compression_levels > 0)
    {
        cmc_global_msg("\t\t Intra Elem Compression Byte Counts: ", info.global_level_bytes.back(), " bytes");
    }

    cmc_global_msg("\t Mesh Compression Global Element Count (size: ", info.global_num_elem_indications_step.size(), "): ");
    for (int lvl_idx{0}; lvl_idx < info.mesh_compression_levels - 1; ++lvl_idx)
    {
        cmc_global_msg("\t\t Level ", lvl_idx, ": ", info.global_num_elem_indications_step[lvl_idx], " bytes");
    }
    cmc_global_msg("\t Partition Table Byte Count: ", info.partition_table_size, " bytes");
    cmc_global_msg("\t Serialized Mesh Compression Huffman Codes: ", info.mesh_compression_huffman_codes, " bytes");
    cmc_global_msg("\t Serialized Intra Elem Compression Huffman Codes: ", info.intra_elem_compression_huffman_codes, " bytes");
    cmc_global_msg("\t Serialized Global Mesh Encoding: ", info.num_bytes_global_mesh_encoding, " bytes");
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

    const SizeType num_global_bytes = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(num_bytes[0]);
    ++offset;
    const SizeType num_preamble_bytes = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(num_bytes[1]);
    ++offset;

    const SizeType num_preamble_values = num_preamble_bytes / sizeof(uint64_t);
    cmc::cmc_global_msg("Num global bytes: ", num_global_bytes, ", Num preamble bytes: ",num_preamble_bytes, ", num_preamble_values: ", num_preamble_values);

    /* Move to the start of the file */
    const int rv_seek_start_again = MPI_File_seek(fhandle, 0, MPI_SEEK_SET);
    MPICheckError(rv_seek_start_again);

    /* Read the complete preamble */
    std::vector<uint64_t> preamble(num_preamble_values);
    const int rv_preamble = MPI_File_read(fhandle, preamble.data(), num_preamble_values, MPI_UINT64_T, &status);
    MPICheckError(rv_preamble);
    CheckMPIReadCorrectness(&status, MPI_UINT64_T, num_preamble_values);

    int a = 0;
    for (const auto& val : preamble)
    {
        cmc_global_msg("Preamble: ", a, ", Value: ", val, "; NE: ", cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(val), ", Bitset: ", std::bitset<64>(cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(val)));
        ++a;
    }
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
    DecompressionVariableMultiData(const std::string& file_name, const MPI_Comm comm, const t8_cmesh_t cmesh, const t8_scheme *scheme);


    //void Decompress(const int max_decompression_level = kMaxDecompressionLevelUndefined);
    void Decompress();

    std::pair<t8_forest_t, std::vector<T>> GetDecompressedData();

private:
    void InquireCompressionInfo();
    const uint64_t* OpenSharedLevelDataWindow(const int level);
    void CloseSharedLevelDataWindow();
    t8_forest_t RepartitionMesh(t8_forest_t adapted_forest);
    void RepartitionData(t8_forest_t adapted_forest, t8_forest_t partitioned_forest);
    void DecodeRootLevelValues();
    void SetStartPositionForStreamDecoder(cmc::bits::StreamDecoder<SymbolType>& stream_decoder, const cmc::bits::vector_view level_encoding, const cmc::bits::vector_view lvl_mesh, const std::vector<LevelPartition>& level_partition_info, const SizeType mesh_offset, const SizeType num_local_elems);

    void DecodeMeshCompressionSteps();

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

    /* Create a forest mesh from the given parameters */
    mesh_.SetMesh(t8_forest_new_uniform (cmesh, scheme, 0, 0, this->comm_));
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
    this->compression_info_.shared_data_init_ptr -= (this->shm_rank_ != kRootRank ? num_preamble_bytes : 0);

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
    const SizeType num_global_lvl_bytes = this->compression_info_.global_level_bytes[level];
    const SizeType num_global_lvl_vals = this->compression_info_.global_level_bytes[level] / sizeof(uint64_t);

    const SizeType proc_offset_val_stream = static_cast<SizeType>(((static_cast<double>(this->shm_rank_) * static_cast<long double>(num_global_lvl_vals)) / static_cast<double>(this->shm_size_)));
    const SizeType next_proc_offset_val_stream = static_cast<SizeType>(((static_cast<double>(this->shm_rank_ + 1) * static_cast<long double>(num_global_lvl_vals)) / static_cast<double>(this->shm_size_)));
    cmc_assert(proc_offset_val_stream <= next_proc_offset_val_stream);
    const int val_stream_length = static_cast<int>(next_proc_offset_val_stream - proc_offset_val_stream);
    const int byte_stream_length = val_stream_length * sizeof(uint64_t);
    cmc_global_msg("SHM Rank ", this->shm_rank_," starts at value: ", proc_offset_val_stream, " and ends at value: ", next_proc_offset_val_stream, ", val length: ", val_stream_length, " byte stream length: ", byte_stream_length);
    /* Store the offset to the start of the contiguous memeory of the global level data */
    const MPI_Aint global_level_start_offset = 0 - proc_offset_val_stream;

    /* Allocate a window on the shared memory communicators */
    uint64_t* shm_mem{nullptr};
    const int rv_win_alloc = MPI_Win_allocate_shared(static_cast<MPI_Aint>(byte_stream_length), sizeof(uint64_t), MPI_INFO_NULL, this->shm_comm_, &shm_mem, &(this->lvl_window_));
    MPICheckError(rv_win_alloc);

    /* Move to the correct position in the file for the current level on the current process */
    const MPI_Offset var_level_offset = std::accumulate(this->compression_info_.global_level_bytes.begin(), std::next(this->compression_info_.global_level_bytes.begin(), level), 0);
    const MPI_Offset file_offset = this->compression_info_.offset_start_encoding + var_level_offset + proc_offset_val_stream * sizeof(uint64_t);
    
    cmc_global_msg("SHM Rank ", this->shm_rank_," var level offset: ", var_level_offset, ", file_offset: ", file_offset);
    const int rv_seek_start = MPI_File_seek(this->fhandle_, file_offset, MPI_SEEK_SET);
    MPICheckError(rv_seek_start);

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
    const uint64_t* global_var_start = shm_mem - global_level_start_offset;

    cmc_global_msg("First eight encoding levels:");
        for (int idx{0}; idx < val_stream_length; ++idx)
        {
            cmc_global_msg("Bitset: ", std::bitset<64>(cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(global_var_start + idx))), ", value: ", cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(*(global_var_start + idx)));
        }

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
DecompressionVariableMultiData<T, DIM, N>::RepartitionMesh(t8_forest_t adapted_forest)
{
    /* Keep the not-partitioned forest */
    t8_forest_ref(adapted_forest);

    /* Allocate a forest */
    t8_forest_t partitioned_forest;
    t8_forest_init(&partitioned_forest);

    /* Partition the forest */
    const int partition_for_coarsening = 0; //TODO: change to 'one' when partition for coarsening is in t8code
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
GetRootLevelValuesFromView(cmc::bits::vector_view lvl_data_start_view, const t8_gloidx_t mesh_offset, const t8_locidx_t num_local_elems)
{
    /* Move to the start position (the root data is not encoded and therfore directly addressable) */
    lvl_data_start_view.MoveToOffsetBitInStream(mesh_offset * sizeof(T) * cmc::bits::kCharBit);

    std::vector<T> data;
    data.reserve(num_local_elems);

    /* Decode the root values */
    for (t8_locidx_t elem_idx{0}; elem_idx < num_local_elems; ++elem_idx)
    {
        /* De-Serialize this root value and store it in the current data vector */
        const OneByteResidualType uvalue = lvl_data_start_view.GetNextBitSequence<OneByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
        data.push_back(std::bit_cast<T>(uvalue));
        cmc_global_msg("Root value ", elem_idx," is: ", data.back());
    }

    return data;
}

template<TwoByteArithmeticType T>
inline std::vector<T>
GetRootLevelValuesFromView(cmc::bits::vector_view lvl_data_start_view, const t8_gloidx_t mesh_offset, const t8_locidx_t num_local_elems)
{
    /* Move to the start position (the root data is not encoded and therfore directly addressable) */
    lvl_data_start_view.MoveToOffsetBitInStream(mesh_offset * sizeof(T) * cmc::bits::kCharBit);

    std::vector<T> data;
    data.reserve(num_local_elems);

    /* Decode the root values */
    for (t8_locidx_t elem_idx{0}; elem_idx < num_local_elems; ++elem_idx)
    {
        /* De-Serialize this root value and store it in the current data vector */
        const TwoByteResidualType uvalue = lvl_data_start_view.GetNextBitSequence<TwoByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
        data.push_back(std::bit_cast<T>(uvalue));
        cmc_global_msg("Root value ", elem_idx," is: ", data.back());
    }

    return data;
}

template<FourByteArithmeticType T>
inline std::vector<T>
GetRootLevelValuesFromView(cmc::bits::vector_view lvl_data_start_view, const t8_gloidx_t mesh_offset, const t8_locidx_t num_local_elems)
{
    /* Move to the start position (the root data is not encoded and therfore directly addressable) */
    lvl_data_start_view.MoveToOffsetBitInStream(mesh_offset * sizeof(T) * cmc::bits::kCharBit);

    std::vector<T> data;
    data.reserve(num_local_elems);

    /* Decode the root values */
    for (t8_locidx_t elem_idx{0}; elem_idx < num_local_elems; ++elem_idx)
    {
        /* De-Serialize this root value and store it in the current data vector */
        const FourByteResidualType uvalue = lvl_data_start_view.GetNextBitSequence<FourByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
        data.push_back(std::bit_cast<T>(uvalue));
        cmc_global_msg("Root value ", elem_idx," is: ", data.back());
    }

    return data;
}

template<EightByteArithmeticType T>
inline std::vector<T>
GetRootLevelValuesFromView(cmc::bits::vector_view lvl_data_start_view, const t8_gloidx_t mesh_offset, const t8_locidx_t num_local_elems)
{
    /* Move to the start position (the root data is not encoded and therfore directly addressable) */
    lvl_data_start_view.MoveToOffsetBitInStream(mesh_offset * sizeof(T) * cmc::bits::kCharBit);

    std::vector<T> data;
    data.reserve(num_local_elems);

    /* Decode the root values */
    for (t8_locidx_t elem_idx{0}; elem_idx < num_local_elems; ++elem_idx)
    {
        /* De-Serialize this root value and store it in the current data vector */
        const EightByteResidualType uvalue = lvl_data_start_view.GetNextBitSequence<EightByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
        data.push_back(std::bit_cast<T>(uvalue));
        cmc_global_msg("Root value ", elem_idx," is: ", data.back());
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

    /* Partition the mesh */
    t8_forest_t init_mesh = this->mesh_.GetMesh();
    t8_forest_t partitioned_mesh = this->RepartitionMesh(init_mesh);
    t8_forest_unref(&init_mesh);
    this->mesh_.SetMesh(partitioned_mesh);

    /* Get the global offset of the first element */
    const t8_gloidx_t mesh_offset = t8_forest_get_first_local_leaf_element_id(mesh_.GetMesh());

    /* Get the number of local elements */
    const t8_locidx_t num_local_elems = t8_forest_get_local_num_leaf_elements(mesh_.GetMesh());

    /* De-Serialize the values */
    this->data_ = GetRootLevelValuesFromView<T>(data_view, mesh_offset, num_local_elems);

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
    //const auto lower_bound_iter = std::lower_bound(level_partition_info.begin(), level_partition_info.end(), next_rank_starting_pos, [](const LevelPartition& part1, const LevelPartition& part2){return part1.elem_offset < part2.elem_offset;}); 
    cmc_global_msg("Mesh offset: ", mesh_offset, " num local elems:  ", num_local_elems, ", next rank start pos: ", next_rank_starting_pos);
    cmc_global_msg("level_partition_info.size() : ", level_partition_info.size());
    /* Iterate until we find the largest partition bound before the partiton start of the next rank */
    int part_idx{0};
    for (size_t idx{0}; idx < level_partition_info.size(); ++idx)
    {
        if (level_partition_info[idx].elem_offset < next_rank_starting_pos)
        {
            part_idx = idx;
            cmc_global_msg("part_idx: ", part_idx);
        } else
        {
            break;
        }
    }

    /* Copy the view to determine the correct parallel offsets */
    cmc::bits::vector_view lvl_mesh = lvl_mesh_encoding;

    cmc_global_msg("part_idx: ", part_idx, ", (level_partition_info[part_idx].elem_offset : ", level_partition_info[part_idx].elem_offset);
    /* Set the mesh encoding correctly to the start of the offset */
    lvl_mesh.MoveToOffsetBitInStream(mesh_offset + (level_partition_info[part_idx].elem_offset <= mesh_offset ? 0 : level_partition_info[part_idx].elem_offset - mesh_offset));

    SizeType num_entropy_codes_correction{0};
    SizeType num_elements_correction{0};

    bool once_potential_incomplete_tree{true};

    /* The last rank does not need a correction */
    if (part_idx < level_partition_info.size())
    {
        /* Get the offset from the Partition Info */
        const SizeType partition_bound = level_partition_info[part_idx].elem_offset;

        /** Iterate the local elements until we arrive at the elem_offset and count from thereon the number of entropy codes up to the process end **/
        /* Get the number of local trees */
        const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees(this->mesh_.GetMesh());
        
        /* Iterate over the local trees */
        for (t8_locidx_t tree_idx{0}, num_elems_skipped{0}; tree_idx < num_local_trees; ++tree_idx)
        {
            /* Get the corresponding tree class */
            const t8_eclass_t tree_class = t8_forest_get_tree_class (this->mesh_.GetMesh(), tree_idx);
            
            /* Get the local number of elements in the tree */
            const t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (this->mesh_.GetMesh(), tree_idx);
            
            /* Check if we need to iterate through this tree or whether we can skip it completely */
            if (mesh_offset + num_elems_skipped + num_elements_in_tree < partition_bound)
            {
                /* Skip the whole tree */
                num_elems_skipped += num_elements_in_tree;
                continue;
            }

            /* Get the scheme from the mesh */
            const t8_scheme_c* scheme = t8_forest_get_scheme (this->mesh_.GetMesh());

            t8_locidx_t tree_local_start_idx{0};
            if (once_potential_incomplete_tree)
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

            /* Update the number of elements that needs to be corrected */
            num_elements_correction += num_elements_in_tree - tree_local_start_idx;
        }
    }

    /* Now, we need to exchange the corrections */
    std::array<SizeType, 3> exchange_data{mesh_offset, num_entropy_codes_correction, num_elements_correction};

    /* Allocate an ouput vector locally */
    std::vector<SizeType> offset_array(3 * (this->comm_size_));

    cmc_global_msg("Before Allgather");
    /* Exchange the data */
    const int rv_all_gather_offset = MPI_Allgather(exchange_data.data(), 3, MPI_SIZE_TYPE, offset_array.data(), 3, MPI_SIZE_TYPE, this->comm_);
    MPICheckError(rv_all_gather_offset);

    cmc_global_msg("After Allgather");

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
    //while(offset_array[3 * rank_offset] <= first_partition_bound_offset && rank_offset <= this->comm_rank_)
    //{
    //    ++rank_offset;
    //}
    for (int rank_idx{0}; rank_idx <= this->comm_rank_; ++rank_idx)
    {
        if (offset_array[3 * rank_offset] > first_partition_bound_offset)
        {
            break;
        } else
        {
            rank_offset = rank_idx;
        }
    }
    /* Count the number of entropy codes we need to correct the compression partition bound */
    int num_entropy_codes_to_correct{0};
    int num_elements_to_correct{0};
    for (int rank_idx{rank_offset}; rank_idx < this->comm_rank_; ++rank_idx)
    {
        num_entropy_codes_to_correct += offset_array[3 * rank_idx + 1];
        num_elements_to_correct += offset_array[3 * rank_idx + 2];
    }

    cmc::cmc_debug_msg(this->comm_, "rank: ", this->comm_rank_, ", Mesh elem off: ", mesh_offset, ", Offset last part bound: ", level_partition_info[part_idx].elem_offset, ", num_entropy_codes to proc end: ", num_entropy_codes_correction);
    cmc::cmc_debug_msg(this->comm_, "rank: ", this->comm_rank_, ", Mesh elem off: ", mesh_offset, ", First rel part bound: ", first_partition_bound_offset, ", num_entropy_codes to correct: ", num_entropy_codes_to_correct, ", Num elems to correct: ", num_elements_to_correct);
    cmc::cmc_global_msg("level_partition_info[first_part_idx].coding_byte_offset : ", level_partition_info[first_part_idx].coding_byte_offset);
    /* Set the start to the first relevant partition bound */
    const SizeType level_entropy_offset = level_partition_info[first_part_idx].coding_byte_offset;

    /* Set the level stream decoder correctly */
    cmc::bits::vector_view adjusted_level_data = level_encoding;//lvl_mesh_encoding;

    cmc_global_msg("Bit Offset in mesh straem: ", level_entropy_offset * cmc::bits::kCharBit);
    adjusted_level_data.MoveToOffsetBitInStream(level_entropy_offset * cmc::bits::kCharBit);

    /* Set the adjusted in the stream decoder */
    stream_decoder.StartDecoding(adjusted_level_data);

    /* And now, we need to skip the number of entropy codes to correct in order to move to the process-local start position */
    SizeType entropy_codes_skipped{0};
    SizeType elements_skipped{0};
    while (elements_skipped < num_elements_to_correct || entropy_codes_skipped < num_entropy_codes_to_correct)
    {
        cmc_global_msg("Adapt start position in stream decoder");
        /* We are reading the next bit, if it is a refinement we discard the entropy codes and move on;
         * in case the element remains unchanged, we just move on with the next bit */
        if (stream_decoder.GetNextBit())
        {
            /* If there will be a refinement, we are skipping entropy codes */
            /* Get the next symbol */
            const SymbolType symbol = std::invoke([&stream_decoder](){
                SymbolType next_symbol = stream_decoder.DecodeNextEntropySymbol();
                if (next_symbol == kProcessEndSymbol<T>) [[unlikely]]
                {
                    while (next_symbol == kProcessEndSymbol<T>)
                    {
                        stream_decoder.ApplyProcessEndSymbol64Bit();
                        next_symbol = stream_decoder.DecodeNextEntropySymbol();
                    }
                }
                return next_symbol;
            });

            /* Get the LZC */
            const int lzc = GetLZCFromEntropySymbol(symbol);

            /* Check if there are significant residual bits */
            if (lzc < sizeof(T) * cmc::bits::kCharBit - 1) [[likely]]
            {
                /* Compute the residual length */
                const int residual_length = sizeof(T) * cmc::bits::kCharBit - 1 - lzc;

                /* Skip those bits */
                stream_decoder.SkipNextBits(residual_length);
            }

            /* Update the discarded entropy codes */
            ++entropy_codes_skipped;
        }

        /* Update the discarded elements */
        ++elements_skipped;
    }
}


template<FourByteArithmeticType T, int32_t DIM>
requires Dimension<DIM>
struct RefinementIterationData
{
    RefinementIterationData(const std::span<T> current_data, cmc::bits::StreamDecoder<SymbolType>& stream_decoder_)
    : data(current_data), stream_decoder{stream_decoder_}
    {
        fine_level_data.reserve(current_data.size() * 2 *DIM + 1);
    }

    void LeaveElementUnchanged(const int local_idx);
    //void PerformRefinement(const int local_idx, const int num_elements);
    void PerformRefinement(const int coarse_value_id, const int num_elements);
    bool WillNextElementBeRefined() {return stream_decoder.GetNextBit();}

    const std::span<T> data;
    cmc::bits::StreamDecoder<SymbolType>& stream_decoder;
    std::vector<T> fine_level_data;
};

template<FourByteArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
RefinementIterationData<T, DIM>::PerformRefinement(const int local_idx, const int num_elements)
{
    /* Get the utilized predictor */
    const FourByteResidualType predictor = std::bit_cast<FourByteResidualType>(this->data[local_idx]);

    /* Iterate over all finer elements that will be constructed */
    for (int elem_idx{0}; elem_idx < num_elements; ++elem_idx)
    {
        #if 0
        /* Get the next entropy symbol */
        SymbolType symbol = stream_decoder.DecodeNextEntropySymbol();

        /* Check that we do not get an process end symbol */
        if (symbol == kProcessEndSymbol<T>) [[unlikely]]
        {
            while (symbol == kProcessEndSymbol<T>)
            {
                stream_decoder_.ApplyProcessEndSymbol64Bit();
                symbol = stream_decoder.DecodeNextEntropySymbol();
            }
        }
        #endif

        /* Get the next symbol */
        const SymbolType symbol = std::invoke([this](){
            SymbolType next_symbol = stream_decoder.DecodeNextEntropySymbol();
            if (next_symbol == kProcessEndSymbol<T>) [[unlikely]]
            {
                while (next_symbol == kProcessEndSymbol<T>)
                {
                    stream_decoder.ApplyProcessEndSymbol64Bit();
                    next_symbol = stream_decoder.DecodeNextEntropySymbol();
                }
            }
            return next_symbol;
        });

        /* Get the LZC from the symbol */
        const int lzc = GetLZCFromEntropySymbol(symbol);

        if (lzc < sizeof(T) * cmc::bits::kCharBit - 1) [[likely]]
        {
            /* Compute the length of the significant residual bits */
            const int residual_length = sizeof(T) * cmc::bits::kCharBit - 1 - lzc;

            /* We obtain the residual and add the implicit one bit  */
            const FourByteResidualType residual = this->stream_decoder.GetNextBitSequence<FourByteResidualType>(residual_length) | (FourByteResidualType{1} << residual_length);

            /* We create the residual applied value */
            if (IsApproximationGreater(symbol))
            {
                /* Subtract the residual from the prediction */
                const FourByteResidualType value = cmc::bits::IntegerSubtraction(predictor, residual);

                /* Store the residual applied value */
                this->fine_level_data.push_back(std::bit_cast<T>(value));
            } else
            {
                /* Add the residual to the prediction */
                const FourByteResidualType value = cmc::bits::IntegerAddition(predictor, residual);

                /* Store the residual applied value */
                this->fine_level_data.push_back(std::bit_cast<T>(value));
            }
        } else
        {
            /* Compute the residual in case we do not need to extract a bit-sequence */
            const FourByteResidualType residual = (FourByteResidualType{0x70000000} >> lzc);

            /* We create the residual applied value */
            if (IsApproximationGreater(symbol))
            {
                /* Subtract the residual from the prediction */
                const FourByteResidualType value = cmc::bits::IntegerSubtraction(predictor, residual);

                /* Store the residual applied value */
                this->fine_level_data.push_back(std::bit_cast<T>(value));
            } else
            {
                /* Add the residual to the prediction */
                const FourByteResidualType value = cmc::bits::IntegerAddition(predictor, residual);

                /* Store the residual applied value */
                this->fine_level_data.push_back(std::bit_cast<T>(value));
            }
        }
    }
}

template<FourByteArithmeticType T, int32_t DIM>
requires Dimension<DIM>
inline void
RefinementIterationData<T, DIM>::LeaveElementUnchanged(const int local_idx)
{
    this->fine_level_data.push_back(this->data[local_idx]);
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
inline t8_locidx_t
LosslessMultiResDecompression (t8_forest_t forest,
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

    /* Check whether we perform a refinement or not */
    if (adapt_data->WillNextElementBeRefined())
    {
        cmc_global_msg("elem ", lelement_id, " will be refined");
        /* If the element will be refined */
        /* Get the local start index */
        const int local_idx = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;
        /* Get the number of children elements into which this element will be refined */
        const int num_children = ts->element_get_num_children(tree_class, elements[0]);
        /* Perform the refinement */
        adapt_data->PerformRefinement(local_idx, num_children);
        return cmc::t8::kRefineElement;
    } else
    {
        cmc_global_msg("elem ", lelement_id, " remians unchanged");
        /* If the element stays unchanged */
        const int local_idx = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;
        /* Perform the unchangendness */
        adapt_data->LeaveElementUnchanged(local_idx);
        return cmc::t8::kLeaveElementUnchanged;
    }
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
        const SizeType mesh_offset = (SizeType) t8_forest_get_first_local_leaf_element_id(mesh_.GetMesh());
        const SizeType num_local_elems = (SizeType) t8_forest_get_local_num_leaf_elements(mesh_.GetMesh());

        /* Move the stream decoder to the correct starting position */
        this->SetStartPositionForStreamDecoder(this->stream_decoder_, level_data_encoding, level_mesh_encoding, level_partition, mesh_offset, num_local_elems);

        cmc_debug_msg("Size of data before adaptation: ", this->data_.size());
    
        /****** Perform the local refinement onto the next finer level as dictated by the mesh encoding ******/
        /* Crerate the adaptation data */
        RefinementIterationData<T, DIM> adapt_data(std::span<T>(this->data_), this->stream_decoder_);

        /* Perform a decompression iteration */
        t8_forest_t adapted_forest = t8_forest_new_adapt(mesh_.GetMesh(), LosslessMultiResDecompression<T, DIM>, 0, 0, static_cast<void*>(&adapt_data));
        cmc_debug_msg("The mesh adaptation step is finished; resulting in ", t8_forest_get_global_num_leaf_elements(adapted_forest), " global elements");

        cmc_debug_msg("Size of data after adaptation: ", this->data_.size());
        /****** Partition the mesh and the data ******/

        /* Store the refined data */
        this->data_ = std::move(adapt_data.fine_level_data);

        /* Repartition the mesh */
        t8_forest_t partitioned_forest = this->RepartitionMesh(adapted_forest);
        cmc_debug_msg("Num lcoal elems after partition: ", this->data_.size());
        /* repartition thefien data to coincide with the mesh */
        this->RepartitionData(adapted_forest, partitioned_forest);
        cmc_debug_msg("Size of data after partition: ", t8_forest_get_global_num_leaf_elements(partitioned_forest));
        /* Store the partitioned mesh */
        mesh_.SetMesh(partitioned_forest);

        /****** Clean-Up of the decompression iteration ******/
        /* Free the former/coarser forest */
        t8_forest_unref(&adapted_forest);

        /* CLose this levels data window */
        this->CloseSharedLevelDataWindow();

        /* Update the decompression count */
        ++(this->decompression_step_idx_);
        cmc_debug_msg("The mesh decompression step ", mesh_lvl," has been completed.");
    }
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

    cmc_debug_msg("Decompression of variable ", reinterpret_cast<const char*>(this->compression_info_.name.data()), " starts...");

    /* Inquire the basic information about the compression */
    this->InquireCompressionInfo();

    /* Decode the root level */
    this->DecodeRootLevelValues();

    /* Perform the level-wise iterative decompression */
    this->DecodeMeshCompressionSteps();

    #if 0

    /* Perform the intra-element decompression */
    this->DecodeIntraElementSteps();
    #endif

    /* Close shared window on compression infos */
    /* Impose a barrier to finish all outstanding reads */
    const int rv_shm_barrier = MPI_Barrier(this->shm_comm_);
    MPICheckError(rv_shm_barrier);

    /* Free the allocated window */
    const int rv_win_free = MPI_Win_free(&(this->compression_info_.shared_data));
    MPICheckError(rv_win_free);

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
    mesh_.SetMesh(nullptr);
    std::vector<T> decompressed_data{};
    std::swap(decompressed_data, this->data_);
    return std::make_pair(mesh, std::move(decompressed_data));
}


}
#endif /* !CMC_PAR_MULTI_RES_DECOMPRESSION_HXX */
