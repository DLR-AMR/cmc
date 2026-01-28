#ifndef CMC_DECOMPRESSION_MPI_IO_INPUT_HXX
#define CMC_DECOMPRESSION_MPI_IO_INPUT_HXX

#include "cmc.hxx"
#include "mpi/cmc_mpi.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "compression_io/cmc_compression_mpi_io_util.hxx"
#include "amr/lossless/cmc_byte_par_decompression_variable.hxx"

#include <vector>
#include <string>

namespace cmc::compression_io::mpi
{

struct VariableHull;

class Reader
{
public:
    Reader() = delete;
    ~Reader() = default;

    Reader(const Reader& other) = default;
    Reader& operator=(const Reader& other) = default;
    Reader(Reader&& other) = default;
    Reader& operator=(Reader&& other) = default;

    Reader(const std::string& file_name, const MPI_Comm comm)
    : file_name_{file_name}, comm_{comm} {}

    void ReadVariableHulls();
    std::vector<VariableHull> GetVariableHulls() const;
    template <typename T> std::unique_ptr<cmc::decompression::par::AbstractByteParDecompressionVariable<T>> ReadVariableForDecompression(const std::string& var_name);

private:
    void OpenFile();
    void CloseFile();

    const std::string file_name_;
    const MPI_Comm comm_{MPI_COMM_NULL};
    MPI_File fhandle_;

    SizeType file_num_header_bytes_{0}; 
    SizeType file_num_variables_{0};

    std::vector<VariableHull> variable_hulls_;

    bool has_file_been_opened_{false};
};

struct VariableHull
{
    VariableHull() = delete;
    VariableHull(const VarIDType id_, const std::string name_, const SizeType global_num_mesh_bytes_, const SizeType global_num_data_bytes_, const SizeType global_file_offset_mesh_start_, const SizeType global_file_offset_data_start_)
    : id{id_}, name(name_), global_num_mesh_bytes{global_num_mesh_bytes_}, global_num_data_bytes{global_num_data_bytes_}, global_file_offset_mesh_start{global_file_offset_mesh_start_}, global_file_offset_data_start{global_file_offset_data_start_} {}

    VarIDType id;
    std::string name;
    SizeType global_num_mesh_bytes;
    SizeType global_num_data_bytes;
    SizeType global_file_offset_mesh_start;
    SizeType global_file_offset_data_start;
};

inline void
Reader::OpenFile()
{
    /* Open the file for reading */
    const int opening_mode = MPI_MODE_RDONLY;
    const int rv_open = MPI_File_open(this->comm_, this->file_name_.c_str(), opening_mode, MPI_INFO_NULL, &(this->fhandle_));
    MPICheckError(rv_open);

    cmc_debug_msg(this->comm_, "The file ", this->file_name_, " has been opened.");
}

inline void
Reader::CloseFile()
{
    /* Close the file */
    const int rv_close = MPI_File_close(&(this->fhandle_));
    MPICheckError(rv_close);

    cmc_debug_msg(this->comm_, "The file ", this->file_name_, " has been closed.");
}

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

inline 
std::vector<VariableHull>
Reader::ReadVariableHulls()
{
    /* Open the file */
    if (not has_file_been_opened_)
    {
        this->OpenFile();
    }

    /* Move to the start of the file */
    const int rv_seek_start = MPI_File_seek(this->fhandle_, 0, MPI_SEEK_SET);
    MPICheckError(rv_seek_start);

    /* Read in the first bytes to alraedy gather information about the header */
    const int num_start_bytes{kNumBytesFileSignature + 2 * sizeof(SizeType)};
    std::vector<uint8_t> header_info(num_start_bytes);

    MPI_Status status;

    /* Read the first bytes */
    const int rv_start_bytes = MPI_File_read(this->fhandle_, header_info.data(), num_start_bytes, MPI_UINT8_T, &status);
    MPICheckError(rv_start_bytes);
    CheckMPIReadCorrectness(&status, MPI_UINT8_T, num_start_bytes);

    /* Check the file signature */
    if (not CheckCorrectFileSignature(header_info))
    {
        cmc_err_msg("The file has not been created by cmc and therefore cannot be read!");
    }
    
    size_t offset = kNumBytesFileSignature;
    cmc_debug_msg("The file signature consisted of ", kNumBytesFileSignature, " bytes.");

    /* Next, read the number of header bytes */
    file_num_header_bytes_ = GetValueFromByteStream<SizeType>(header_info.data() + offset);
    offset += sizeof(SizeType);

    cmc_debug_msg("The general file header consists of additional ", file_num_header_bytes_, " bytes.");

    /* Now, we read the number of variables */
    file_num_variables_ = GetValueFromByteStream<SizeType>(header_info.data() + offset);
    offset += sizeof(SizeType);

    if (file_num_variables_ == 0)
    {
        cmc_err_msg("There are no variables in the file to read.");
    }

    cmc_debug_msg("The file contains ", file_num_variables_, " encoded variable(s).");

    /* Compute the number of bytes needed for the remaining header */
    const int num_remaining_bytes = static_cast<int>(file_num_header_bytes_ - 2 * sizeof(SizeType));
    cmc_debug_msg("Number of remaining header bytes to be read are ", num_remaining_bytes, " bytes.");
    cmc_assert(num_remaining_bytes > 0);

    /* Allocate memory for the variable information */
    std::vector<uint8_t> variable_infos(num_remaining_bytes);

    /* Read in the remaining header information */
    const int rv_remaining_header = MPI_File_read(this->fhandle_, variable_infos.data(), num_remaining_bytes, MPI_UINT8_T, &status);
    MPICheckError(rv_remaining_header);
    CheckMPIReadCorrectness(&status, MPI_UINT8_T, num_remaining_bytes);

    static_assert(sizeof(char) == sizeof(uint8_t));

    std::vector<VariableHull> variable_hulls;
    variable_hulls.reserve(file_num_variables_);

    SizeType offset_header{0};
    SizeType global_offset = kNumBytesFileSignature + file_num_header_bytes_;

    /* Iterate over the variables and extract the stored information */
    for (SizeType var_iter{0}; var_iter < file_num_variables_; ++var_iter)
    {
        /* Get the var_id */
        const VarDataInfoType var_id = GetValueFromByteStream<VarDataInfoType>(variable_infos.data() + offset_header);
        offset_header += sizeof(VarDataInfoType);

        /* Get the variable name */
        const char* var_name_ptr = reinterpret_cast<const char*>(variable_infos.data() + offset_header);
        /* Read the name */
        std::string var_name(var_name_ptr);
        /* Update the offset counter */
        offset_header += kNumCharsVariableName * sizeof(uint8_t);

        /* Get the number of global mesh bytes */
        const SizeType global_mesh_bytes = GetValueFromByteStream<SizeType>(variable_infos.data() + offset_header);
        offset_header += sizeof(SizeType);

        /* Get the number of global data bytes */
        const SizeType global_data_bytes = GetValueFromByteStream<SizeType>(variable_infos.data() + offset_header);
        offset_header += sizeof(SizeType);

        //TODO:Add later
        /* Get the number of the variable's offset bytes */
        //const SizeType global_offset_bytes = GetValueFromByteStream<SizeType>(variable_infos.data() + offset_header);
        //offset_header += sizeof(SizeType);

        /* Compute the file offsets for the variable */
        const SizeType global_mesh_offset = global_offset;
        global_offset += global_mesh_bytes;
        const SizeType global_data_offset = global_offset;
        global_offset += global_data_bytes;

        cmc_debug_msg("The variable (ID: ", var_id, ") ", var_name, " consists of ", global_mesh_bytes, " encoded mesh bytes and ", global_data_bytes, " global encoded data bytes. The global mesh offset in the file is ", global_mesh_offset, " bytes and the global data offset in the file is ", global_data_offset, " bytes.");
        
        /* Construct the variable hull */
        variable_hulls.emplace_back(var_id, var_name, global_mesh_bytes, global_data_bytes, global_mesh_offset, global_data_offset);
    }

    /* Close the file */
    if (not has_file_been_opened_)
    {
        this->CloseFile();
    }

    /* Store the variable hulls */
    variable_hulls_ = std::move(variable_hulls);
}

inline 
std::vector<VariableHull>
Reader::GetVariableHulls() const 
{
    return variable_hulls_;
}

template <typename T>
std::unique_ptr<cmc::decompression::par::AbstractByteParDecompressionVariable<T>>
Reader::ReadVariableForDecompression(const std::string& var_name)
{
    /* Open the file */
    this->OpenFile();
    has_file_been_opened_ = true;

    /* Check if the variable hulls have been read from the file already and if not, do so */
    if (variable_hulls_.empty())
    {
        this->ReadVariableHulls();
    }

    /* Check for the variable */
    const auto var_hull_iter = std::find_if(variable_hulls_.begin(), variable_hulls_.end(),
                                            [&var_name](const VariableHull& var_hull){return (var_name.compare(var_hull.name) == 0);});
    if (var_hull_iter == variable_hulls_.end())
    {
        cmc_err_msg("The variable ", var_name, " does not exist within the compressed file!");
    }                                            
    cmc_debug_msg("Reading of the compressed variable ", var_name, " starts.");
    cmc_debug_msg("The mesh encoding of this variable starts at ", var_hull_iter->global_file_offset_mesh_start, " and consists of ", var_hull_iter->global_num_mesh_bytes, " bytes.");

    /* Move to the start of this variable's mesh encoding */
    const int rv_seek_start_me = MPI_File_seek(this->fhandle_, var_hull_iter->global_file_offset_mesh_start, MPI_SEEK_SET);
    MPICheckError(rv_seek_start_me);

    /* Read in the complete mesh encoding */
    std::vector<uint8_t> var_mesh_header(var_hull_iter->global_num_mesh_bytes);

    /* Read in global mesh encoding */
    MPI_Status status;
    const int rv_read_mesh_encoding = MPI_File_read(this->fhandle_, var_mesh_header.data(), var_hull_iter->global_num_mesh_bytes, MPI_UINT8_T, &status);
    MPICheckError(rv_read_mesh_encoding);
    CheckMPIReadCorrectness(&status, MPI_UINT8_T, var_hull_iter->global_num_mesh_bytes);

    /* Read the variable mesh attributes from the encoding */
    const auto [mesh_var_id, mesh_lvl_byte_counts, mesh_lvls_encoding] = DecodeMeshStream(var_mesh_header);

    /* We read a slab of data from the file that definetly contains the variable header */
    std::vector<uint8_t> var_data_header(kMaxNumBytesVariableHeader);

    /* Read in the corresponding data stream */
    const int rv_read_var_header_encoding = MPI_File_read(this->fhandle_, var_data_header.data(), kMaxNumBytesVariableHeader, MPI_UINT8_T, &status);
    MPICheckError(rv_read_mesh_encoding);
    /* Check if the read bytes are in range of the minimum and maximum permitted bytes */
    int elem_count{0};
    const int rv_check_read = MPI_Get_count(status, MPI_UINT8_T, &elem_count);
    MPICheckError(rv_check_read);
    if (elem_count < kMinNumBytesVariableHeader || elem_count > kMaxNumBytesVariableHeader)
    {
        cmc_err_msg("The amount of read bytes does not match with the variable header!");
    }

    /* Decode the variable header */
    const auto [data_var_id, data_type, cr_scheme, data_lvl_byte_counts, data_header_processed_bytes] = DecodeVariableHeaderStream(var_data_header.data());
    
    /* Close the file since we have read all information currently needed */
    has_file_been_opened_ = false;
    this->CloseFile();

    /* Check whether the data type is correct */
    if (static_cast<CmcType>(cr_scheme) != ConvertToCmcType<T>())
    {
        cmc_err_msg("The supplied data type does not match the type of the encoded data of the variable ", var_name);
    }

    /* Check whether the variable id is correct */
    if (mesh_var_id != data_var_id)
    {
        cmc_err_msg("The variable mesh id does not match the variable data id".);
    }

    /* Check whether the level counts of the mesh and the data matches */
    if (mesh_lvl_byte_counts.size() != data_lvl_byte_counts.size())
    {
        cmc_err_msg("The number of encoding levels of the mesh and the data does not match.");
    }

    /* Compute the correct offset in the file for the start of this variable's encoded data stream */
    const SizeType file_var_encoding_offset = var_hull_iter->global_file_offset_data_start + data_header_processed_bytes;

    /* Invoke the correct decompressor */
    switch (cr_scheme)
    {
        case CompressionSchema::ParallelMultiResExtraction:
            return std::make_unique<lossless::par::multi_res::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh), num_compression_iterations);
        break;
        default:
            cmc_err_msg("The compression schema of the compressed variable is not recognized.");
            return std::make_unique<lossless::par::multi_res::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh), num_compression_iterations);
    }
}


}

#endif /* !CMC_DECOMPRESSION_MPI_IO_INPUT_HXX */
