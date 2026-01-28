#ifndef CMC_COMPRESSION_MPI_IO_OUTPUT_HXX
#define CMC_COMPRESSION_MPI_IO_OUTPUT_HXX

#include "cmc.hxx"
#include "mpi/cmc_mpi.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "compression_io/cmc_compression_mpi_io_util.hxx"

#include <vector>
#include <numeric>
#include <utility>

namespace cmc::compression_io::mpi
{

template <typename T>
class Writer
{
public:
    Writer(const std::string& file_name, const MPI_Comm comm)
    : file_name_{file_name}, comm_{comm} {};

    void SetVariable(const cmc::lossless::par::AbstractByteCompressionVariable<T>& variable);
    void Write();

private:
    void OpenFile();
    void CloseFile();
    void WriteFileHeader();
    MPI_Offset GetDefaultFileHeaderOffset() const;
    MPI_Offset ComputeGlobalFileStorage(const MPI_Comm comm, const int rank) const;
    std::pair<std::vector<MPI_Offset>, std::vector<MPI_Offset>> ComputeLocalChunkOffsets(const int comm_rank, const int comm_size, const SizeType num_file_header_bytes, const std::vector<SizeType>& global_information_lengths) const;
    MPI_Offset ComputeGlobalFileSizeFromLocalStreamLengths(const SizeType num_file_header_bytes, const std::vector<SizeType>& global_information_lengths) const;
    std::vector<SizeType> ExchangeLocalInformationStreamLenghts(const MPI_Comm comm, const int rank) const;

    const std::string file_name_;
    const MPI_Comm comm_;
    MPI_File fhandle_;

    std::vector<SerializedVariableInfo<T>> variable_infos_;

    int var_id_counter_{0};
};

template <typename T>
void
Writer<T>::SetVariable(const cmc::lossless::par::AbstractByteCompressionVariable<T>& variable)
{
    SerializedVariableInfo<T> var_info(variable, var_id_counter_, comm_);
    variable_infos_.push_back(var_info);

    /* Increase the variable ID counter */
    ++var_id_counter_;
}

template <typename T>
inline void
Writer<T>::OpenFile()
{
    /* Open the file (i.e. create if necessary) in write only mode */
    const int opening_mode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
    const int rv_open = MPI_File_open(this->comm_, this->file_name_.c_str(), opening_mode, MPI_INFO_NULL, &(this->fhandle_));
    MPICheckError(rv_open);

    cmc_debug_msg(this->comm_, "The file ", this->file_name_, " has been opened.");
}

template <typename T>
inline void
Writer<T>::CloseFile()
{
    /* Open the file (i.e. create if necessary) in write only mode */
    const int rv_close = MPI_File_close(&(this->fhandle_));
    MPICheckError(rv_close);

    cmc_debug_msg(this->comm_, "The file ", this->file_name_, " has been closed.");
}

inline void 
CheckMPIWriteCorrectness(const MPI_Status* status, const MPI_Datatype datatype, const int expected_num_elems)
{
    int elem_count_{0};
    const int rv_check_write = MPI_Get_count(status, datatype, &elem_count_);
    MPICheckError(rv_check_write);
    if (elem_count_ != expected_num_elems)
    {
        cmc_global_msg("The expeceted number of bytes could not be written to the file.");
        MPICheckError(MPI_ERR_COUNT);
    }
}

template <typename T>
void
Writer<T>::WriteFileHeader()
{
#ifdef CMC_ENABLE_DEBUG
    int rank{0};
    const int rv_rank = MPI_Comm_rank(comm_, &rank);
    MPICheckError(rv_rank);

    if (rank != kRootRank)
    {
        cmc_warn_msg("WriteFileHeader should only be called by the root rank.");
    }
#endif

    /* Generate the default header */
    DefaultFileHeader<T> file_header(variable_infos_);
    /* Get the serialized hedaer */
    const std::vector<uint8_t>& serialized_header = file_header.GetSerializedFileHeader();
    cmc_debug_msg("The serialized file header consists of ", serialized_header.size(), " bytes.");

    /* Move to the beginning of the file */
    const int rv_seek_beg = MPI_File_seek(this->fhandle_, 0, MPI_SEEK_SET);
    MPICheckError(rv_seek_beg);

    MPI_Status status;

    /* Write the file signature */
    const int rv_file_sign = MPI_File_write(this->fhandle_, kFileSignature.data(), kNumBytesFileSignature, MPI_UINT8_T, &status);
    MPICheckError(rv_file_sign);
    CheckMPIWriteCorrectness(&status, MPI_UINT8_T, static_cast<int>(kNumBytesFileSignature));

    /* Write the file format */

    /* Write the general file header */
    const int rv_file_header = MPI_File_write(this->fhandle_, serialized_header.data(), serialized_header.size(), MPI_UINT8_T, &status);
    MPICheckError(rv_file_header);
    CheckMPIWriteCorrectness(&status, MPI_UINT8_T, static_cast<int>(serialized_header.size()));
}

#if 0
template <typename T>
inline MPI_Offset
Writer<T>::GetDefaultFileHeaderOffset() const
{
    /* The file header consists of the file singature, the num header bytes, the number of variables and for each variable of a id, name, num bytes mesh and num bytes data */
    static_assert(sizeof(uint8_t) == sizeof(char));
    return (kNumBytesFileSignature * sizeof(uint8_t) //File Signature
            + kNumBytesFileFormatVersion * sizeof(uint8_t) //Format Version of data format
            + sizeof(SizeType) //Num Header bytes
            + sizeof(SizeType) //Num Variables
            + variable_infos_.size() * //For each variable
            (
                sizeof(VarDataInfoType) //Variable ID
                + kNumCharsVariableName * sizeof(char) //Variable Name
                + sizeof(SizeType) //Num bytes global mesh encoding
                + sizeof(SizeType) //Num bytes global data encoding
            ));
}
#endif


template <typename T>
MPI_Offset
Writer<T>::ComputeGlobalFileStorage(const MPI_Comm comm, const int rank) const
{
    /* Compute the global file storage (The computation is only correct on the root rank)*/
    SizeType global_file_storage{0};
    
    if (rank == kRootRank)
    {
        /* Add the general file header */
        global_file_storage += static_cast<SizeType>(DefaultFileHeader<T>::GetFileHeaderSize(variable_infos_.size()));

        /* Iterate over the variable information and count the bytes */
        for (auto var_info_iter = variable_infos_.begin(); var_info_iter != variable_infos_.end(); ++var_info_iter)
        {
            global_file_storage += var_info_iter->GetVariableMeshHeader().size();
            global_file_storage += var_info_iter->GetVariableDataHeader().size();
            global_file_storage += var_info_iter->GetVariablePartitionTableHeader().size();
            global_file_storage += var_info_iter->GetGlobalNumBytesMesh();
            global_file_storage += var_info_iter->GetGlobalNumBytesData();
            global_file_storage += var_info_iter->GetGlobalNumBytesPartitionTable();
        }
    }

    /* Broadcast the global file storage */
    const int rv_bcast_file_size = MPI_Bcast(&global_file_storage, 1, ConvertToMPIType<SizeType>(), kRootRank, comm);
    MPICheckError(rv_bcast_file_size);

    return static_cast<MPI_Offset>(global_file_storage);
}

template <typename T>
std::vector<SizeType>
Writer<T>::ExchangeLocalInformationStreamLenghts(const MPI_Comm comm, const int rank) const
{
    int comm_size{0};
    const int rv_size = MPI_Comm_size(comm, &comm_size);
    MPICheckError(rv_size);

    /* Count the global offsets */
    size_t count_local_offsets{0};
    for (auto var_info_iter = variable_infos_.begin(); var_info_iter != variable_infos_.end(); ++var_info_iter)
    {
        /* For each level, we are holding a contiguous piece of information for the mesh and the data on each process */
        count_local_offsets += var_info_iter->GetVariable().GetEncodedData().size();
        count_local_offsets += var_info_iter->GetVariable().GetEncodedMesh().size();
    }

    /* Allocate a local vector holding the local information lengths */
    std::vector<SizeType> local_information_lengths;
    local_information_lengths.reserve(count_local_offsets);

    /* Fill the local information lengths */
    for (auto var_info_iter = variable_infos_.begin(); var_info_iter != variable_infos_.end(); ++var_info_iter)
    {
        /* Get the encoded mesh */
        const auto mesh_encoding = var_info_iter->GetVariable().GetEncodedMesh();

        /* Store all local mesh encoding lengths*/
        for (auto me_iter = mesh_encoding.rbegin(); me_iter != mesh_encoding.rend(); ++me_iter)
        {
            /* Get the length of this local level mesh encoding */
            SizeType mesh_lvl_encoding_length = me_iter->size();
            /* In case of the first iteration and the root rank, we add up additionally the header for the mesh variable */
            if (me_iter == mesh_encoding.rbegin() && rank == kRootRank)
            {
                mesh_lvl_encoding_length += var_info_iter->GetVariableMeshHeader().size();
            }

            /* Store this levels encoding length */
            local_information_lengths.push_back(mesh_lvl_encoding_length);
        }

        /* Get the encoded data */
        const auto data_encoding = var_info_iter->GetVariable().GetEncodedData();

        /* Afterwards, we store all local data encoding lengths */
        for (auto de_iter = data_encoding.rbegin(); de_iter != data_encoding.rend(); ++de_iter)
        {
            /* Get the length of this local level mesh encoding */
            SizeType data_lvl_encoding_length = de_iter->size();
            /* In case of the first iteration and the root rank, we add up additionally the header for the mesh variable */
            if (de_iter == data_encoding.rbegin() && rank == kRootRank)
            {
                data_lvl_encoding_length += var_info_iter->GetVariableDataHeader().size();
            }

            /* Store this levels encoding length */
            local_information_lengths.push_back(data_lvl_encoding_length);
        }
    }

    cmc_assert(local_information_lengths.size() == count_local_offsets);

    /* Allocate a vector to store the offsets */
    std::vector<SizeType> global_offsets(count_local_offsets * comm_size, 0);

    /* Distribute the information about the local encoding lengths */
    const int rv_allg_enc_lengths = MPI_Allgather(local_information_lengths.data(), count_local_offsets, ConvertToMPIType<SizeType>(), global_offsets.data(), count_local_offsets, ConvertToMPIType<SizeType>(), comm);
    MPICheckError(rv_allg_enc_lengths);

    return global_offsets;
}

template <typename T>
inline MPI_Offset
Writer<T>::ComputeGlobalFileSizeFromLocalStreamLengths(const SizeType num_file_header_bytes, const std::vector<SizeType>& global_information_lengths) const
{
    return static_cast<MPI_Offset>(num_file_header_bytes + std::accumulate(global_information_lengths.begin(), global_information_lengths.end(), 0));
}

template <typename T>
std::pair<std::vector<MPI_Offset>, std::vector<MPI_Offset>>
Writer<T>::ComputeLocalChunkOffsets(const int comm_rank, const int comm_size, const SizeType num_file_header_bytes, const std::vector<SizeType>& global_information_lengths) const
{
    /* Compute the number of entries to skip to the next processes entry in the global length vector */
    size_t local_num_mesh_offsets{0};
    size_t local_num_data_offsets{0};
    for (auto tmp_var_info_iter = variable_infos_.begin(); tmp_var_info_iter != variable_infos_.end(); ++tmp_var_info_iter)
    {
        local_num_mesh_offsets += tmp_var_info_iter->GetVariable().GetEncodedMesh().size();
        local_num_data_offsets += tmp_var_info_iter->GetVariable().GetEncodedData().size();
    }

    const int skip_to_next_process_same_var = local_num_mesh_offsets + local_num_data_offsets;

    /* Allocate an output vector for the partition table offsets */
    std::vector<MPI_Offset> partition_table_offsets; //The output is only relevant for the root rank
    partition_table_offsets.reserve(variable_infos_.size());

    /* Allocate output vectors for the mesh and data offsets */
    std::vector<MPI_Offset> mesh_offsets;
    mesh_offsets.reserve(local_num_mesh_offsets);

    std::vector<MPI_Offset> data_offsets;
    data_offsets.reserve(local_num_data_offsets);

    /* Start global offset directly after the general file header */
    MPI_Offset offset = static_cast<MPI_Offset>(num_file_header_bytes);

    /* Counter for the local offset to the next variable */
    int global_var_offset_accessor_id{0};

    /* Iterate over the variables */
    for (auto var_info_iter = variable_infos_.begin(); var_info_iter != variable_infos_.end(); ++var_info_iter)
    {
        /* Get the encoded mesh */
        const auto mesh_encoding = var_info_iter->GetVariable().GetEncodedMesh();
        const int num_mesh_levels = mesh_encoding.size();

        /* Iterate over all levels */
        for (int lvl{0}; lvl < num_mesh_levels; ++lvl)
        {
            int global_length_arr_offset{global_var_offset_accessor_id};

            /* Iterate over all processes */
            for (int rank{0}; rank < comm_size; ++rank)
            {
                if (rank == comm_rank)
                {
                    /* This processes lvl offset has been reached */
                    mesh_offsets.push_back(offset);
                }

                /* Add the processes information length to the offset */
                offset += global_information_lengths[lvl + global_length_arr_offset];

                /* Increase the access counter to the next processes byte stream length */
                global_length_arr_offset += skip_to_next_process_same_var;
            }
        }

        /* Get the encoded mesh */
        const auto data_encoding = var_info_iter->GetVariable().GetEncodedData();
        const int num_data_levels = data_encoding.size();

        /* Iterate over all levels */
        for (int lvl{0}; lvl < num_data_levels; ++lvl)
        {
            /* The data stream lengths come after the mesh byte streams */
            int global_length_arr_offset{num_mesh_levels + global_var_offset_accessor_id};
            
            /* Iterate over all processes */
            for (int rank{0}; rank < comm_size; ++rank)
            {
                if (rank == comm_rank)
                {
                    /* This processes lvl offset has been reached */
                    data_offsets.push_back(offset);
                }

                /* Add the processes information length to the offset */
                offset += global_information_lengths[lvl + global_length_arr_offset];

                /* Increase the access counter to the next processes byte stream length */
                global_length_arr_offset += skip_to_next_process_same_var;
            }
        }

        /* Move the accessor to the start of the next variable */
        global_var_offset_accessor_id += num_mesh_levels + num_data_levels;
    }

    return std::make_pair(std::move(mesh_offsets), std::move(data_offsets));
}


template <typename T>
void
Writer<T>::Write()
{
    cmc_debug_msg("Start write");

    /* Determine the rank of the process and the size of the communicator */
    int rank{0}, size{0};
    const int rv_rank = MPI_Comm_rank(comm_, &rank);
    MPICheckError(rv_rank);
    const int rv_size = MPI_Comm_size(comm_, &size);
    MPICheckError(rv_size);

    cmc_debug_msg("Rank: ", rank, ", Size: ", size);

    /* We exchange all information byte stream lengths for the local processes to be written out globally in order to compute the file offsets for each process */
    const std::vector<SizeType> all_local_stream_lengths = this->ExchangeLocalInformationStreamLenghts(comm_, rank);

    cmc_debug_msg("After exchange local lengths: Size of all_local_stream_lengths: ", all_local_stream_lengths.size());

    /* Open the file */
    this->OpenFile();

    cmc_debug_msg("The file has been openend.");

    /* Get the length of the file header */
    const MPI_Offset general_file_hedaer_size = this->GetDefaultFileHeaderOffset();

    cmc_debug_msg("General file hedaer size: ", general_file_hedaer_size);

    /* Afterwards, we compute the global file storage to store all of the data */
    const MPI_Offset file_size = this->ComputeGlobalFileSizeFromLocalStreamLengths(general_file_hedaer_size, all_local_stream_lengths);
    
    cmc_debug_msg("Computed global file size: ", file_size);

    /* Pre-Allocate global file storage */
    //TODO: This gives problems when an already existent file is opened in parallel; If the file is deleted beforehand (or in a serial call) it runs without error (Permission denied)
    const int rv_file_prealloc = MPI_File_preallocate(this->fhandle_, file_size);
    MPICheckError(rv_file_prealloc);
    cmc_debug_msg("The file has been pre-allocated.");
    /* Compute the global offsets for the process local infromation byte streams for all levels of all variables */
    const auto [mesh_offsets, data_offsets] = this->ComputeLocalChunkOffsets(rank, size, general_file_hedaer_size, all_local_stream_lengths);
    
    cmc_debug_msg("Num mesh offsets: ", mesh_offsets.size(), ", Num data offsets: ", data_offsets.size());
    for (auto miter = mesh_offsets.begin(); miter != mesh_offsets.end(); ++miter)
    {
        cmc_debug_msg("Mesh Offset: ", *miter);
    }
    for (auto miter = data_offsets.begin(); miter != data_offsets.end(); ++miter)
    {
        cmc_debug_msg("Data Offset: ", *miter);
    }

    /* Write the general file header on the root rank */
    if (rank == kRootRank)
    {
        this->WriteFileHeader();
        cmc_debug_msg("The root rank has written the file header.");
    }

    /* Iterate over the variables and store the process-local information chunks them in the global file */
    int accessor_mesh_chunk{0}, accessor_data_chunk{0};
    for (auto var_info_iter = variable_infos_.begin(); var_info_iter != variable_infos_.end(); ++var_info_iter)
    {
        /* Get the encoded mesh */
        const auto mesh_encoding = var_info_iter->GetVariable().GetEncodedMesh();

        /* Iterate over the process-local chunks of each level and write them to the file */
        for (auto mesh_lvl_iter = mesh_encoding.rbegin(); mesh_lvl_iter != mesh_encoding.rend(); ++mesh_lvl_iter, ++accessor_mesh_chunk)
        {
            /* Check during the first iteration and only on the root rank that the mesh header will be written */
            if (mesh_lvl_iter == mesh_encoding.rbegin() && rank == kRootRank)
            {
                /* During the first iteration on the root rank, we need to additionally store the mesh header */
                const std::vector<uint8_t>& mesh_header = var_info_iter->GetVariableMeshHeader();

                MPI_Status status;

                /* Move the file handle pointer to the correct global offset */
                const int rv_pos_fhandle = MPI_File_seek(this->fhandle_, mesh_offsets[accessor_mesh_chunk], MPI_SEEK_SET);
                MPICheckError(rv_pos_fhandle);

                /* Write the header */
                const int rv_write_mesh_header = MPI_File_write(this->fhandle_, mesh_header.data(), mesh_header.size(), MPI_UINT8_T, &status);
                MPICheckError(rv_write_mesh_header);
                CheckMPIWriteCorrectness(&status, MPI_UINT8_T, static_cast<int>(mesh_header.size()));

                /* Potentially, write out the additional information on this level for the root rank afterwards */
                if (not mesh_lvl_iter->empty())
                {
                    /* Write the root level on the root rank */
                    const int rv_write_mesh_lvl = MPI_File_write(this->fhandle_, mesh_lvl_iter->data(), mesh_lvl_iter->size(), MPI_UINT8_T, &status);
                    MPICheckError(rv_write_mesh_lvl);
                    CheckMPIWriteCorrectness(&status, MPI_UINT8_T, static_cast<int>(mesh_lvl_iter->size()));

                    /* Continue with the next level */
                    continue;
                }
            }

            /* We only write something out if there is actual information on this process */
            if (not mesh_lvl_iter->empty())
            {
                MPI_Status status;

                /* Move the file handle to correct global offset */
                const int rv_pos_fhandle = MPI_File_seek(this->fhandle_, mesh_offsets[accessor_mesh_chunk], MPI_SEEK_SET);
                MPICheckError(rv_pos_fhandle);

                /* Write out the level mesh information chunk */
                const int rv_write_mesh_lvl = MPI_File_write(this->fhandle_, mesh_lvl_iter->data(), mesh_lvl_iter->size(), MPI_UINT8_T, &status);
                MPICheckError(rv_write_mesh_lvl);
                CheckMPIWriteCorrectness(&status, MPI_UINT8_T, static_cast<int>(mesh_lvl_iter->size())); 
            }
        }

        cmc_debug_msg("The mesh encoding of the variable ", var_info_iter->GetName(), " has been written to the file.");

        /* After the mesh variable has been written, we will store the encding of the variable's data */
        /* Get the encoded data */
        const auto data_encoding = var_info_iter->GetVariable().GetEncodedData();

        /* Iterate over the process-local chunks of each level and write them to the file */
        for (auto data_lvl_iter = data_encoding.rbegin(); data_lvl_iter != data_encoding.rend(); ++data_lvl_iter, ++accessor_data_chunk)
        {
            /* Check during the first iteration and only on the root rank that the data header will be written */
            if (data_lvl_iter == data_encoding.rbegin() && rank == kRootRank)
            {
                /* During the first iteration on the root rank, we need to additionally store the data header */
                const std::vector<uint8_t>& data_header = var_info_iter->GetVariableDataHeader();

                MPI_Status status;

                /* Move the file handle pointer to the correct global offset */
                const int rv_pos_fhandle = MPI_File_seek(this->fhandle_, data_offsets[accessor_data_chunk], MPI_SEEK_SET);
                MPICheckError(rv_pos_fhandle);

                /* Write the header */
                const int rv_write_data_header = MPI_File_write(this->fhandle_, data_header.data(), data_header.size(), MPI_UINT8_T, &status);
                MPICheckError(rv_write_data_header);
                CheckMPIWriteCorrectness(&status, MPI_UINT8_T, static_cast<int>(data_header.size()));

                /* Potentially, write out the additional information on this level for the root rank afterwards */
                if (not data_lvl_iter->empty())
                {
                    /* Write the root level on the root rank */
                    const int rv_write_data_lvl = MPI_File_write(this->fhandle_, data_lvl_iter->data(), data_lvl_iter->size(), MPI_UINT8_T, &status);
                    MPICheckError(rv_write_data_lvl);
                    CheckMPIWriteCorrectness(&status, MPI_UINT8_T, static_cast<int>(data_lvl_iter->size()));

                    /* Continue with the next level */
                    continue;
                }
            }

            /* We only write something out if there is actual information on this process */
            if (not data_lvl_iter->empty())
            {
                MPI_Status status;

                /* Move the file handle to correct global offset */
                const int rv_pos_fhandle = MPI_File_seek(this->fhandle_, data_offsets[accessor_data_chunk], MPI_SEEK_SET);
                MPICheckError(rv_pos_fhandle);

                /* Write out the level mesh information chunk */
                const int rv_write_data_lvl = MPI_File_write(this->fhandle_, data_lvl_iter->data(), data_lvl_iter->size(), MPI_UINT8_T, &status);
                MPICheckError(rv_write_data_lvl);
                CheckMPIWriteCorrectness(&status, MPI_UINT8_T, static_cast<int>(data_lvl_iter->size())); 
            }
        }

        cmc_debug_msg("The data encoding of the variable ", var_info_iter->GetName(), " has been written to the file.");

        /* After the data has been written, the root rank will additionally wirte out the partition table*/
        /* Write the general file header on the root rank */
        if (rank == kRootRank)
        {
            /* Get the variable partition table header */
            const std::vector<uint8_t>& mesh_header = var_info_iter->GetVariablePartitionTableHeader();

            /* Move to the start of the partition table encoding in the file */
            const MPI_Offset file_start_partition_table_offset = 
            const int rv_seek_part_table =  MPI_File_seek(this->fhandle_, file_start_partition_table_offset, MPI_SEEK_SET);
            MPICheckError(rv_seek_part_table);

            MPI_Status status;

            /* Move the file handle pointer to the correct global offset */
            const int rv_pos_fhandle = MPI_File_seek(this->fhandle_, mesh_offsets[accessor_mesh_chunk], MPI_SEEK_SET);
            MPICheckError(rv_pos_fhandle);

            cmc_debug_msg("The root rank has written the partition table of variable ", var_info_iter->GetName(), " to the file.");
        }
    }

    cmc_debug_msg("All variables have been written.");

    /* Close the file */
    this->CloseFile();

    cmc_debug_msg("The file has been closed");
}


}

#endif /* !CMC_COMPRESSION_MPI_IO_OUTPUT_HXX */
