#ifndef CMC_PATCH_BYTE_COMPRESSION_OUTPUT_TXX
#define CMC_PATCH_BYTE_COMPRESSION_OUTPUT_TXX

namespace cmc::compression_io
{


template<typename T>
void
Writer::SetDataVariable(cmc::IPatchCompressionVariable<T>* variable, const int var_id)
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
        global_byte_count += levelwise_entropy_codes_[idx].size();
        global_byte_count += levelwise_data_[idx].size();
    }
    cmc_debug_msg("Global byte count: ", global_byte_count);

    /* We define some attributes that will be added to the variable */
    std::vector<cmc::nc::Attribute> attributes;
    attributes.reserve(num_compression_iterations + 7 + static_cast<int32_t>(variable->GetDimensionality()) * num_compression_iterations);
    attributes.emplace_back(id_attr, var_id);
    attributes.emplace_back(data_type_attr, static_cast<int32_t>(ConvertToCmcType<T>()));
    attributes.emplace_back(compression_schema_attr, static_cast<int32_t>(variable->GetCompressionSchema()));
    attributes.emplace_back(initial_data_layout_attr, static_cast<int32_t>(variable->GetInitialDataLayout()));
    attributes.emplace_back(dim_attr, static_cast<int32_t>(variable->GetDimensionality()));
    attributes.emplace_back(num_compression_iterations_attr, num_compression_iterations);
    uint64_t level_offset = 0;
    /* Store the level-wise byte count */
    for (int idx = num_compression_iterations - 1; idx >= 0; --idx)
    {
        level_offset += levelwise_entropy_codes_[idx].size() + levelwise_data_[idx].size();
        cmc_debug_msg("IDX: ", idx, ", level offset: ", level_offset);
        const std::string attr_name = "lvl_" + std::to_string(num_compression_iterations - 1 - idx) + "_num_bytes";
        attributes.emplace_back(attr_name.c_str(), static_cast<uint64_t>(level_offset));
    }

    /* Get the pyramidal level dimension lengths */
    const std::vector<std::vector<size_t>>& dim_lengths_pyramid = variable->GetDimensionLengthPyramid();
    attributes.emplace_back(dim_lengths_pyra_lvls_attr, static_cast<int32_t>(dim_lengths_pyramid.size()));

    for (int idx = dim_lengths_pyramid.size() - 1; idx >= 0; --idx)
    {
        for (int dim_id{0}; dim_id < static_cast<int32_t>(variable->GetDimensionality()); ++dim_id)
        {
            const std::string attr_name = "dim" + std::to_string(dim_id + 1) + "_" + std::to_string(dim_lengths_pyramid.size() - 1 - idx);
            attributes.emplace_back(attr_name.c_str(), static_cast<int32_t>(dim_lengths_pyramid[idx][dim_id]));
        }
    }

    /* Define a new variable */
    nc::SpecificVariable<uint8_t> spec_variable(variable->GetName(), var_id);

    /* Set the global diemnsion length */
    spec_variable.SetGlobalDimensionLength(global_byte_count);
    
    /* Set the data for the variable */
    level_offset = 0;
    for (int idx = num_compression_iterations - 1; idx >= 0; --idx)
    {
        const uint64_t entropy_codes_size = levelwise_entropy_codes_[idx].size();

        /* Set the entropy cods of the current level with the correct offset */
        spec_variable.SetData(std::move(levelwise_entropy_codes_[idx]), level_offset);

        /* Update the levelwise offset */
        level_offset += entropy_codes_size;

        const uint64_t encoded_lvl_data_size = levelwise_data_[idx].size();

        /* Set the data of this level */
        spec_variable.SetData(std::move(levelwise_data_[idx]), level_offset);

        /* Update the levelwise offset */
        level_offset += encoded_lvl_data_size;
    }

    /* Add the variable for output */
    nc_writer_.AddVariable(nc::Variable(std::move(spec_variable), std::move(attributes)));
}

template<typename T>
void
Writer::SetVariable(cmc::IPatchCompressionVariable<T>* variable)
{
    /* Plac the variable in the file */
    SetDataVariable<T>(variable, var_id_counter_);

    /* Update the variable counters */
    ++var_id_counter_;
}


}

#endif /* !CMC_PATCH_BYTE_COMPRESSION_OUTPUT_TXX */
