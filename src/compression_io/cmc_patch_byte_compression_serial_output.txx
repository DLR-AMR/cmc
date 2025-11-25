#ifndef CMC_PATCH_BYTE_COMPRESSION_SERIAL_OUTPUT_TXX
#define CMC_PATCH_BYTE_COMPRESSION_SERIAL_OUTPUT_TXX

#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "compression_io/cmc_serial_io_util.hxx"

#include <cstdio>

namespace cmc::compression_io::serial
{

template<typename T>
void
Writer::SetVariable(cmc::IPatchCompressionVariable<T>* variable)
{
    /* Set the variable for the output */
    variables_.push_back(this->SetDataVariable<T>(variable));
}

template<typename T>
SerialOutputVariable
Writer::SetDataVariable(cmc::IPatchCompressionVariable<T>* variable)
{
    /* Get the name of the variable */
    std::string name = variable->GetName();

    /* Get the encoded level-wise data */
    std::vector<std::vector<uint8_t>> levelwise_entropy_codes_;
    variable->MoveEncodedEntropyCodesInto(levelwise_entropy_codes_);

    std::vector<std::vector<uint8_t>> levelwise_data_;
    variable->MoveEncodedDataInto(levelwise_data_);

    /* Get the pyramidal level dimension lengths */
    const std::vector<std::vector<size_t>>& dim_lengths_pyramid = variable->GetDimensionLengthPyramid();

    cmc_assert(levelwise_data_.size() == levelwise_entropy_codes_.size());

    /* Get the number of compression iterations */
    const int num_compression_iterations = static_cast<int>(levelwise_data_.size());

    /* Gather the global count of bytes for this variable */
    uint64_t global_byte_count{0};

    /* Iterate over all codes and data and gather the overall byte count */
    for (int idx = num_compression_iterations - 1; idx >= 0; --idx)
    {
        global_byte_count += levelwise_entropy_codes_[idx].size();
        global_byte_count += levelwise_data_[idx].size();
    }

    /* We add additional bytes for extra information that is stored */
    global_byte_count += sizeof(attr_type_t) * (kNumHeaderInfo + dim_lengths_pyramid.size() * variable->GetDimensionality());

    cmc_debug_msg("The output data for the variable ", name, " amounts to ", global_byte_count, " bytes.");

    /* Now, we start to serialize the additional information and append the encoded data */

    std::vector<uint8_t> byte_stream;
    byte_stream.reserve(global_byte_count);

    /* Store the data type */
    PushBackValueToByteStream<attr_type_t>(byte_stream, static_cast<attr_type_t>(ConvertToCmcType<T>()));

    /* Store the compression scheme */
    PushBackValueToByteStream<attr_type_t>(byte_stream, static_cast<attr_type_t>(variable->GetCompressionSchema()));

    /* Store the initial data layout */
    PushBackValueToByteStream<attr_type_t>(byte_stream, static_cast<attr_type_t>(variable->GetInitialDataLayout()));

    /* Store the dimensionality */
    PushBackValueToByteStream<attr_type_t>(byte_stream, static_cast<attr_type_t>(variable->GetDimensionality()));

    /* Store the number of compression iterations */
    PushBackValueToByteStream<attr_type_t>(byte_stream, static_cast<attr_type_t>(num_compression_iterations));

    /* Store the levels of the pyramidal dimension lengths */
    PushBackValueToByteStream<attr_type_t>(byte_stream, static_cast<attr_type_t>(dim_lengths_pyramid.size()));

    /* Push back the dimension lengths of each level */
    for (int idx = dim_lengths_pyramid.size() - 1; idx >= 0; --idx)
    {
        for (int dim_id{0}; dim_id < static_cast<attr_type_t>(variable->GetDimensionality()); ++dim_id)
        {
            PushBackValueToByteStream<attr_type_t>(byte_stream, static_cast<attr_type_t>(static_cast<attr_type_t>(dim_lengths_pyramid[idx][dim_id])));
        }
    }

    cmc_debug_msg("File header bytes:");
    for (auto iter = byte_stream.begin(); iter != byte_stream.end(); ++iter)
    {
        cmc_debug_msg("Byte_Stream: ", static_cast<int>(*iter));
    }

    /* Copy the encoded data to the stream */
    for (int idx = num_compression_iterations - 1; idx >= 0; --idx)
    {
        /* Copy the entropy codes */
        std::copy_n(levelwise_entropy_codes_[idx].begin(), levelwise_entropy_codes_[idx].size(), std::back_inserter(byte_stream));
        
        /* Copy the encoded level data */
        std::copy_n(levelwise_data_[idx].begin(), levelwise_data_[idx].size(), std::back_inserter(byte_stream));
    }

    return SerialOutputVariable(std::move(name), std::move(byte_stream));
}


inline void
Writer::Write()
{
    for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter)
    {
        /* Create a file name for the variable */
        const std::string file_name = CreateFileName(file_prefix_, var_iter->name_);
        
        /* Open the file */
        std::FILE* file_out = std::fopen(file_name.c_str(), "wb");
        std::fwrite(var_iter->byte_stream_.data(), sizeof(uint8_t), var_iter->byte_stream_.size(), file_out);
        std::fclose(file_out);
    }

    variables_.clear();
}

}

#endif /* !CMC_PATCH_BYTE_COMPRESSION_SERIAL_OUTPUT_TXX */
