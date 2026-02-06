#ifndef CMC_COMPRESSION_MPI_IO_UTIL_HXX
#define CMC_COMPRESSION_MPI_IO_UTIL_HXX

#include "cmc.hxx"
#include "utilities/cmc_serialization.hxx"
#include "amr/lossless/cmc_byte_par_compression_variable.hxx"

#include "mpi/cmc_mpi_data.hxx"

#include <array>
#include <vector>
#include <cstdint>
#include <string>
#include <algorithm>
#include <tuple>

namespace cmc::compression_io::mpi
{

using SizeType = uint64_t;
using VarIDType = uint32_t;
using VarDataInfoType = uint32_t;
using FormatVersionType = uint32_t;

constexpr int kRootRank = 0;

constexpr SizeType kNumCharsVariableName = 256;

constexpr size_t kNumBytesFileSignature = 4;
constexpr std::array<uint8_t, kNumBytesFileSignature> kFileSignature{0x43, 0x4D, 0x43, 0x7E};
constexpr size_t kNumBytesFileFormatVersion = 4;
constexpr std::array<uint8_t, kNumBytesFileFormatVersion> kFormatVersion{0x1A, 0x01, 0x00, 0x00};
constexpr SizeType kMaxNumBytesVariableHeader = 512;
constexpr SizeType kMinNumBytesVariableHeader = sizeof(VarDataInfoType) * 4;

inline bool
CheckCorrectFileSignature(const std::vector<uint8_t>& file_header)
{
    if (file_header.size() < kNumBytesFileSignature)
    {
        cmc_err_msg("The supplied file header is to short to evaluate the file signature!");
    }
    if (file_header[0] == 0x43 && file_header[1] == 0x4D && file_header[2] == 0x43 && file_header[3] == 0x7E)
    {
        return true;
    } else
    {
        return false;
    }
}

class IFileHeader
{
public:
    virtual ~IFileHeader() = default;

    virtual const std::vector<uint8_t>& GetSerializedFileHeader() const = 0;
};

template <typename T>
class SerializedVariableInfo
{
public:
    SerializedVariableInfo(const cmc::lossless::par::AbstractByteCompressionVariable<T>& variable, const VarIDType id, const MPI_Comm comm)
    : variable_{variable}, id_{id}, comm_{comm}
    {
        /* Determine the global number of bytes for the data and the mesh */
        this->DetermineNumGlobalBytes();
        /* Generate variable attributes for further information */
        this->GenerateVariableHeaders();
    }
    ~SerializedVariableInfo() = default;

    std::string GetName() const {return variable_.GetName();}
    MPI_Comm GetMPIComm() const {return variable_.GetMPIComm();}
    VarIDType GetVarID() const {return id_;}
    const cmc::lossless::par::AbstractByteCompressionVariable<T>& GetVariable() const {return variable_;}
    SizeType GetGlobalNumBytesData() const {return num_global_bytes_encoding_;}
    SizeType GetGlobalNumBytesMesh() const {return num_global_bytes_mesh_;}
    SizeType GetGlobalNumBytesPartitionTable() const {return num_global_bytes_partition_table_;} //This function is only meaningful on the root rank; on all other ranks zero

    const std::vector<uint8_t>& GetVariableDataHeader() const {return data_header_;}
    const std::vector<uint8_t>& GetVariableMeshHeader() const {return mesh_header_;}
    const std::vector<uint8_t>& GetVariablePartitionTableHeader() const {return partition_header_;}
private:
    void DetermineNumGlobalBytes();
    void GenerateVariableHeaders();

    const cmc::lossless::par::AbstractByteCompressionVariable<T>& variable_;
    VarIDType id_{0};
    const MPI_Comm comm_{MPI_COMM_NULL};
    SizeType num_global_bytes_encoding_{0};
    SizeType num_global_bytes_mesh_{0};
    SizeType num_global_bytes_partition_table_{0}; //!< This value is only correct on the root rank

    std::vector<uint8_t> partition_header_;
    std::vector<uint8_t> mesh_header_;
    std::vector<uint8_t> data_header_;
};

template <typename T>
void
SerializedVariableInfo<T>::DetermineNumGlobalBytes()
{
    std::array<SizeType, 2> local_bytes{0, 0};
    std::array<SizeType, 2> global_bytes{0, 0};

    /* Gather the local bytes per level and accumulate them afterwards */
    const std::vector<std::vector<uint8_t>>& lvl_encoded_data = variable_.GetEncodedData();
    const std::vector<std::vector<uint8_t>>& lvl_encoded_mesh = variable_.GetEncodedMesh();
    const std::vector<std::vector<uint8_t>>& lvl_encoded_partition_table = variable_.GetEncodedPartitionTable();
    const size_t num_lvls = lvl_encoded_data.size();

    cmc_assert(lvl_encoded_data.size() == lvl_encoded_mesh.size() && lvl_encoded_mesh.size() == lvl_encoded_partition_table.size());

    /* Accumulate the levels locally */
    for (size_t lvl_iter{0}; lvl_iter < num_lvls; ++lvl_iter)
    {
        local_bytes[0] += lvl_encoded_data[lvl_iter].size();
        local_bytes[1] += lvl_encoded_mesh[lvl_iter].size();
        num_global_bytes_partition_table_ += lvl_encoded_partition_table[lvl_iter].size(); //!< This sum is only correct on the root rank, all other rank obtain a zero
    }

    //TODO: Add the partition table or Bcast it
    /* Reduce the global byte counts */
    const int rv_allreduce = MPI_Allreduce(local_bytes.data(), global_bytes.data(), 2, ConvertToMPIType<SizeType>(), MPI_SUM, comm_);
    MPICheckError(rv_allreduce);

    /* Store the global byte counts */
    num_global_bytes_encoding_ = global_bytes[0];
    num_global_bytes_mesh_ = global_bytes[1];
}

template <typename T>
void
SerializedVariableInfo<T>::GenerateVariableHeaders()
{
    int rank{0};
    const int rv_rank = MPI_Comm_rank(comm_, &rank);
    MPICheckError(rv_rank);
    int size{1};
    const int rv_size = MPI_Comm_size(comm_, &size);
    MPICheckError(rv_size);

    if (rank == kRootRank)
    {
        const VarDataInfoType var_id = static_cast<VarDataInfoType>(id_);
        
        /* Get the encoded data on the root rank */
        const auto encoded_var_data = variable_.GetEncodedData();

        /*** Define a header for the data variable ***/
        std::vector<uint8_t> data_header;
        data_header.reserve(sizeof(VarDataInfoType) * 4 + sizeof(SizeType) * encoded_var_data.size());

        PushBackValueToByteStream<VarDataInfoType>(data_header, var_id);

        const VarDataInfoType data_type = static_cast<VarDataInfoType>(ConvertToCmcType<T>());
        PushBackValueToByteStream<VarDataInfoType>(data_header, data_type);

        const VarDataInfoType compression_scheme = static_cast<VarDataInfoType>(variable_.GetCompressionSchema());
        PushBackValueToByteStream<VarDataInfoType>(data_header, compression_scheme);

        const VarDataInfoType num_lvls_data = static_cast<VarDataInfoType>(encoded_var_data.size());
        PushBackValueToByteStream<VarDataInfoType>(data_header, num_lvls_data);
        
        /* At the beginning of the root ranks encoding of the global level byte count is stored, we access it and store it as an attribute for the data */
        for (int idx = num_lvls_data - 1; idx >= 0; --idx)
        {
            const SizeType num_bytes = GetValueFromByteStream<SizeType>(encoded_var_data[idx].data());
            PushBackValueToByteStream<SizeType>(data_header, num_bytes);
        }

        /* Check if the header size is within the anticipated range */
        if (data_header.size() > kMaxNumBytesVariableHeader)
        {
            cmc_err_msg("The generated variable header is larger than the maximum anticipated variable header size. The permitted size needs to be adapted.");
        }

        /*** END of data variable header ***/

        /* Get the encoded mesh on the root rank */
        const auto encoded_var_mesh = variable_.GetEncodedMesh();

        /*** Define a header for the encoded mesh ***/
        std::vector<uint8_t> mesh_header;
        mesh_header.reserve(sizeof(VarDataInfoType) * 2 + sizeof(SizeType) * encoded_var_mesh.size());

        PushBackValueToByteStream<VarDataInfoType>(mesh_header, var_id);

        const VarDataInfoType num_lvls_mesh = static_cast<VarDataInfoType>(encoded_var_mesh.size());
        PushBackValueToByteStream<VarDataInfoType>(mesh_header, num_lvls_mesh);
    
        /* At the beginning of the root ranks encoding of the global amount of elements is stored, we access it and store it as an attribute for the mesh */
        for (int idx = num_lvls_mesh - 1; idx >= 0; --idx)
        {
            const SizeType num_elems = GetValueFromByteStream<SizeType>(encoded_var_mesh[idx].data());
            PushBackValueToByteStream<SizeType>(mesh_header, num_elems);
        }
        /*** END of encoded mesh header ***/

        /*** Define a header for the encoded partition table ***/
        /* Get the encoded partition table root rank */
        const auto encoded_partition_table = variable_.GetEncodedPartitionTable();
        std::vector<uint8_t> partition_header;
        partition_header.reserve(sizeof(VarDataInfoType) * 3);

        /* Store the variable id */
        PushBackValueToByteStream<VarDataInfoType>(partition_header, var_id);

        /* Store the size of the communicator */
        const VarDataInfoType comm_size = static_cast<VarDataInfoType>(size);
        PushBackValueToByteStream<VarDataInfoType>(partition_header, comm_size);

        /* Store the number of levels */
        const VarDataInfoType num_lvls_partition = static_cast<VarDataInfoType>(encoded_partition_table.size());
        PushBackValueToByteStream<VarDataInfoType>(partition_header, num_lvls_partition);
        /*** END of encoded partition table header ***/

        /* Store the headers */
        data_header_ = std::move(data_header);
        mesh_header_ = std::move(mesh_header);
        partition_header_ = std::move(partition_header);
    }
}

const auto [mesh_var_id, mesh_lvl_byte_counts, mesh_lvls_encoding]
inline
std::tuple<VarDataInfoType, std::vector<SizeType>, std::vector<uint8_t>>
DecodeMeshStream(const std::vector<uint8_t>& var_mesh_header_w_encoding)
{
    cmc_assert(var_mesh_header_w_encoding.size() >= 2 * sizeof(VarDataInfoType));

    /* Get the var ID from the stream */
    const VarDataInfoType var_id = GetValueFromByteStream<VarDataInfoType>(var_mesh_header_w_encoding.data());
    size_t offset = sizeof(VarDataInfoType);

    /* Get the number of encoded levels */
    const VarDataInfoType num_lvls = GetValueFromByteStream<VarDataInfoType>(var_mesh_header_w_encoding.data() + offset);
    offset += sizeof(VarDataInfoType);

    /* Allocate a vector for the reamining encoded stream */
    std::vector<SizeType> global_level_num_elems;
    global_level_num_elems.reserve(num_lvls);

    /* Get the global bytes of the single levels */
    for (VarDataInfoType lvl_iter{0}; lvl_iter < num_lvls; ++lvl_iter)
    {
        /* Get the number of elements on this level */
        const SizeType num_elems = GetValueFromByteStream<SizeType>(var_mesh_header_w_encoding.data() + offset);
        offset += sizeof(SizeType);

        global_level_num_elems.push_back(num_elems);
    }

    /* Copy the remaining encoding */
    cmc_assert(var_mesh_header_w_encoding.size() >= offset);
    const size_t remaining_bytes_encoding = var_mesh_header_w_encoding.size() - offset;

    /* Allocate a new vector for the level-wise encoding */
    std::vector<uint8_t> level_wise_mesh_encoding;
    level_wise_mesh_encoding.reserve(remaining_bytes_encoding);

    /* Copy the encoding over */
    std::copy_n(var_mesh_header_w_encoding.data() + offset, remaining_bytes_encoding, std::back_inserter(level_wise_mesh_encoding));

    /* Return tuple with the extracted information */
    return std::make_tuple(var_id, std::move(global_level_num_elems), std::move(level_wise_mesh_encoding));
}

inline
std::tuple<VarDataInfoType, CmcType, CompressionSchema, std::vector<SizeType>, SizeType>
DecodeVariableHeaderStream(const uint8_t* var_data_header_ptr)
{
    /* Get the var id */
    const VarDataInfoType var_id = GetValueFromByteStream<VarDataInfoType>(var_data_header_ptr);
    SizeType offset = sizeof(VarDataInfoType);

    /* Get the CmcType of the data */
    const CmcType type = static_cast<CmcType>(GetValueFromByteStream<VarDataInfoType>(var_data_header_ptr + offset));
    offset += sizeof(VarDataInfoType);

    /* Get the compression scheme */
    const CompressionSchema cr_scheme = static_cast<CompressionSchema>(GetValueFromByteStream<VarDataInfoType>(var_data_header_ptr + offset));
    offset += sizeof(VarDataInfoType);

    /* Get the levels of the encoding */
    const VarDataInfoType num_lvls = GetValueFromByteStream<VarDataInfoType>(var_data_header_ptr + offset);
    offset += sizeof(VarDataInfoType);

    /* Allocate a vector holding the global bytes per level */
    std::vector<SizeType> global_level_bytes;
    global_level_bytes.reserve(num_lvls);

    /* Get the global bytes of the single levels */
    for (VarDataInfoType lvl_iter{0}; lvl_iter < num_lvls; ++lvl_iter)
    {
        /* Get the number of elements on this level */
        const SizeType num_bytes = GetValueFromByteStream<SizeType>(var_data_header_ptr + offset);
        offset += sizeof(SizeType);

        global_level_bytes.push_back(num_bytes);
    }
    
    /* Return tuple with the extracted information */
    return std::make_tuple(var_id, type, cr_scheme, std::move(global_level_bytes), offset);
}

template <typename T>
class DefaultFileHeader : public IFileHeader
{
public:
    DefaultFileHeader(const std::vector<SerializedVariableInfo<T>>& variables)
    : IFileHeader()
    {
        this->GenerateFileHeader(variables);
    };

    ~DefaultFileHeader() = default;

    const std::vector<uint8_t>& GetSerializedFileHeader() const override;
    static size_t GetFileHeaderSize(const size_t num_variables) const;
private:
    void GenerateFileHeader(const std::vector<SerializedVariableInfo<T>>& variables);
    std::vector<uint8_t> serialized_file_header_;
};

template <typename T>
static size_t
DefaultFileHeader<T>::GetFileHeaderSize(const size_t num_variables) const
{
    /* The file header consists of the file singature, the num header bytes, the number of variables and for each variable of a id, name, num bytes mesh and num bytes data */
    static_assert(sizeof(uint8_t) == sizeof(char));
    return (kNumBytesFileSignature * sizeof(uint8_t) //File Signature
            + kNumBytesFileFormatVersion * sizeof(uint8_t) //Format Version of data format
            + sizeof(SizeType) //Num Header bytes
            + sizeof(SizeType) //Num Variables
            + num_variables * //For each variable
            (
                sizeof(VarDataInfoType) //Variable ID
                + kNumCharsVariableName * sizeof(char) //Variable Name
                + sizeof(SizeType) //Num bytes global mesh encoding
                + sizeof(SizeType) //Num bytes global data encoding
                + sizeof(SizeType) //Num bytes global partition table
            ));
}

template <typename T>
const std::vector<uint8_t>&
DefaultFileHeader<T>::GetSerializedFileHeader() const
{
    cmc_assert(not serialized_file_header_.empty());
    return serialized_file_header_;
}

//TODO: If all compression information lengths are exchanged, the header could be filled with placeholders for the global data and mesh encoding sizes and replaced afetrwards.
//This could remove the global MPI_Reduce operation when the Variable Information is generated.
template <typename T>
void
DefaultFileHeader<T>::GenerateFileHeader(const std::vector<SerializedVariableInfo<T>>& variables)
{
    /* Compute the size of the file header */
    const size_t file_header_size = this->GetFileHeaderSize(variables.size());

    /* Allocate a header */
    std::vector<uint8_t> file_header;
    file_header.reserve(file_header_size);

    /* Calculate the header information byte count and store it */
    const SizeType num_info_header_bytes = file_header_size //File hedaer
                                           - kNumBytesFileSignature * sizeof(uint8_t) //Minus File Signature
                                           - kNumBytesFileFormatVersion * sizeof(uint8_t);//Minus File Format Version

    /* We start by storing the file signature */
    std::copy_n(kFileSignature.begin(), kFileSignature.size(), std::back_inserter(file_header));
    
    /* We continue with the data format version */
    std::copy_n(kFormatVersion.begin(), kFormatVersion.size(), std::back_inserter(file_header));

    /* We store the bytes for the information header */                                            
    PushBackValueToByteStream<SizeType>(file_header, num_info_header_bytes);

    /* Store the number of variables */
    PushBackValueToByteStream<SizeType>(file_header, static_cast<SizeType>(variables.size()));

    /* Store the id, names and the encoding length of the variables */
    for (auto var_iter = variables.begin(); var_iter != variables.end(); ++var_iter)
    {
        std::vector<char> var_name(kNumCharsVariableName, '\0');
        if (var_iter->GetName().size() > kNumCharsVariableName)
        {
            /* Copy on the the number of bytes that are within the permitted range */
            cmc_global_msg("The variable name '", var_iter->GetName(), "' is too long; it gets trimmed to ", kNumCharsVariableName, " characters.");
            std::copy_n(var_iter->GetName().begin(), kNumCharsVariableName, var_name.begin());
        } else
        {
            /* Copy the name of the variable */
            std::copy_n(var_iter->GetName().begin(), var_iter->GetName().size(), var_name.begin());
        }

        /* We store the ID of the variable */
        const VarDataInfoType var_id = static_cast<VarDataInfoType>(var_iter->GetVarID());
        PushBackValueToByteStream<VarDataInfoType>(file_header, var_id);

        /* Store the variable name and the corresponding encoding length */
        /* Reinterpret the variable name as uint8_t's and store the variable name */
        static_assert(sizeof(uint8_t) == sizeof(char));
        if (sizeof(uint8_t) != sizeof(char)) {cmc_err_msg("The byte count of uint8_t and char does not match and therefore the data cannot be written to disk properly.");}        
        uint8_t* var_name_ptr = reinterpret_cast<uint8_t*>(var_name.data());
        std::copy_n(var_name_ptr, kNumCharsVariableName, std::back_inserter(file_header));
    
        /* Store the variable encoding length of the mesh */
        PushBackValueToByteStream<SizeType>(file_header, var_iter->GetGlobalNumBytesMesh());
        cmc_debug_msg("Num mesh bytes that have been stored: ", var_iter->GetGlobalNumBytesMesh());
        /* Store the variable encoding length of the data  */
        PushBackValueToByteStream<SizeType>(file_header, var_iter->GetGlobalNumBytesData());
        cmc_debug_msg("Num data bytes that have been stored: ", var_iter->GetGlobalNumBytesData());
        /* Store the partition table */
        PushBackValueToByteStream<SizeType>(file_header, var_iter->GetGlobalNumBytesPartitionTable());
        cmc_debug_msg("Num partition table bytes that have been stored: ", var_iter->GetGlobalNumBytesPartitionTable());
    }

    /* Store the header */
    serialized_file_header_ = std::move(file_header);
}

}

#endif /* !CMC_COMPRESSION_MPI_IO_UTIL_HXX */
