#ifndef CMC_DECOMPRESSION_NC_INPUT_HXX
#define CMC_DECOMPRESSION_NC_INPUT_HXX

#include "cmc.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "compression_io/cmc_compression_attr_names.hxx"

#ifdef CMC_WITH_T8CODE
#include "amr/lossless/cmc_byte_decompression_variable.hxx"
#include "amr/lossless/cmc_prefix_extraction_decompression.hxx"
#include "amr/lossless/cmc_prefix_extraction_decompression_plain_suffixes.hxx"
#include "amr/lossless/cmc_multi_res_extraction_decompression.hxx"
#include "amr/lossless/cmc_test_pcp4_decompression.hxx"

#include "embedded/lossless/cmc_embedded_byte_decompression_variable.hxx"
#include "embedded/lossless/cmc_embedded_prefix_extraction_decompression.hxx"
#include "embedded/lossless/cmc_embedded_prefix_extraction_decompression_plain_suffixes.hxx" 
#include "embedded/lossless/cmc_embedded_multi_res_extraction_decompression.hxx"
#include "embedded/lossless/cmc_test_pcp4_embedded_decompression.hxx"
#include "utilities/cmc_embedded_variable_attributes.hxx"

#include "embedded/lossy/cmc_embedded_prefix_quantization_decompression.hxx"

#endif

#include "patch/lossless/cmc_patch_prefix_extraction_plain_suffixes_decompression.hxx"
#include "patch/lossless/cmc_patch_multi_res_extraction_decompression.hxx"


#ifdef CMC_WITH_NETCDF
#include "netcdf/cmc_netcdf.hxx"
#include "netcdf/cmc_nc_io.hxx"
#include "netcdf/cmc_nc_reader.hxx"
#endif

#include "mpi/cmc_mpi.hxx"

#include <string>
#include <memory>
#include <vector>
#include <cstddef>

namespace cmc::compression_io::nc
{

class Reader
{
public:
    Reader() = delete;
    ~Reader() = default;

    Reader(const Reader& other) = default;
    Reader& operator=(const Reader& other) = default;
    Reader(Reader&& other) = default;
    Reader& operator=(Reader&& other) = default;

    Reader(const std::string& file_name, const MPI_Comm comm = MPI_COMM_SELF)
    : file_name_{file_name}, comm_{comm}, reader(file_name, comm) {};

#ifdef CMC_WITH_T8CODE
    template <typename T> std::unique_ptr<decompression::AbstractByteDecompressionVariable<T>> ReadVariableForDecompression(const std::string& var_name);
    template <typename T> std::unique_ptr<cmc::IEmbeddedByteDecompressionVariable<T>> ReadEmbeddedVariableForDecompression(const std::string& var_name);
#endif
    template <typename T> std::unique_ptr<cmc::patch::IPatchDecompressionVariable<T>> ReadPatchVariableForDecompression(const std::string& var_name);

private:
    const std::string file_name_;
    const MPI_Comm comm_;
    cmc::nc::Reader reader;
};

#ifdef CMC_WITH_T8CODE

template <typename T>
std::unique_ptr<decompression::AbstractByteDecompressionVariable<T>>
Reader::ReadVariableForDecompression(const std::string& var_name)
{
    cmc_debug_msg("try to read: ", var_name);
    /* First, we check the attributes of the variable */
    std::vector<cmc::nc::Attribute> attributes = reader.ReadVariableAttributes(var_name);

    /* Get the data type of the compressed variable */
    auto data_type_iter = cmc::nc::FindAttribute(attributes, data_type_attr);
    if (data_type_iter == attributes.end()) {cmc_err_msg("The ", data_type_attr, " attribute has not been found.");}
    const cmc::nc::Attribute& data_type_attribute = *data_type_iter;
    const CmcType type = static_cast<CmcType>(std::get<int>(data_type_attribute.GetValue()));
    if (ConvertToCmcType<T>() != type) {cmc_err_msg("The template parameter of the Reader functionality does not match the data_type of the compressed variable.");}

    /* Get the utilized compression scheme */
    auto compression_scheme_iter = cmc::nc::FindAttribute(attributes, compression_schema_attr);
    if (compression_scheme_iter == attributes.end()) {cmc_err_msg("The ", compression_schema_attr, " attribute has not been found.");}
    const cmc::nc::Attribute& compression_scheme_attribute = *compression_scheme_iter;
    const CompressionSchema compression_scheme = static_cast<cmc::CompressionSchema>(std::get<int>(compression_scheme_attribute.GetValue()));

    /* Get the corresponding mesh id */
    auto mesh_id_iter = cmc::nc::FindAttribute(attributes, mesh_id_attr);
    if (mesh_id_iter == attributes.end()) {cmc_err_msg("The ", mesh_id_attr, " attribute has not been found.");}
    const cmc::nc::Attribute& mesh_id_attribute = *mesh_id_iter;
    const int mesh_id = std::get<int>(mesh_id_attribute.GetValue());

    /* Get the number of compression iterations */
    auto num_compression_iterations_iter = cmc::nc::FindAttribute(attributes, num_compression_iterations_attr);
    if (num_compression_iterations_iter == attributes.end()) {cmc_err_msg("The ", num_compression_iterations_attr, " attribute has not been found.");}
    const cmc::nc::Attribute& num_compression_iterations_attribute = *num_compression_iterations_iter;
    const int num_compression_iterations = std::get<int>(num_compression_iterations_attribute.GetValue());

    /* Generate the mesh-variable name */
    const std::string mesh_name = GenerateMeshVariableName(mesh_id);

    cmc_debug_msg("try to read mesh: ", mesh_name);

    /* Get the attributes for the mesh variable */
    std::vector<cmc::nc::Attribute> mesh_attributes = reader.ReadVariableAttributes(mesh_name);

    /* Get the number of compression iterations */
    auto mesh_num_compression_iterations_iter = cmc::nc::FindAttribute(mesh_attributes, num_compression_iterations_attr);
    if (mesh_num_compression_iterations_iter == mesh_attributes.end()) {cmc_err_msg("The ", num_compression_iterations_attr, " attribute of the mesh has not been found.");}
    const cmc::nc::Attribute& mesh_num_compression_iterations_attribute = *mesh_num_compression_iterations_iter;
    const int mesh_num_compression_iterations = std::get<int>(mesh_num_compression_iterations_attribute.GetValue());

    if (mesh_num_compression_iterations != num_compression_iterations) {cmc_err_msg("The data variable and the mesh variable do not have the same amount of compression iterations.");}

    /* Read the global compression streams for all levels */
    std::vector<uint8_t> encoded_data = reader.ReadVariableData<uint8_t>(var_name);
    std::vector<uint8_t> encoded_mesh = reader.ReadVariableData<uint8_t>(mesh_name);
    
    /* Invoke the correct decompressor */
    switch (compression_scheme)
    {
        case CompressionSchema::PrefixExtraction:
            return std::make_unique<lossless::prefix::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh), num_compression_iterations);
        break;
        case CompressionSchema::PrefixExtractionPlainSuffixes:
            return std::make_unique<lossless::prefix::plain_suffix::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh), num_compression_iterations);
        break;
        case CompressionSchema::MultiResExtraction:
            return std::make_unique<lossless::multi_res::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh), num_compression_iterations);
        break;
        case CompressionSchema::_TestPCP4Extraction:
            return std::make_unique<lossless::test_pcp4::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh), num_compression_iterations);
        break;
        default:
            cmc_err_msg("The compression schema of the compressed variable is not recognized.");
            return std::make_unique<lossless::prefix::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh), num_compression_iterations);
    }
}


template <typename T>
std::unique_ptr<cmc::IEmbeddedByteDecompressionVariable<T>>
Reader::ReadEmbeddedVariableForDecompression(const std::string& var_name)
{
    /* First, we check the attributes of the variable */
    std::vector<cmc::nc::Attribute> attributes = reader.ReadVariableAttributes(var_name);

    /* Get the data type of the compressed variable */
    auto data_type_iter = cmc::nc::FindAttribute(attributes, data_type_attr);
    if (data_type_iter == attributes.end()) {cmc_err_msg("The ", data_type_attr, " attribute has not been found.");}
    const cmc::nc::Attribute& data_type_attribute = *data_type_iter;
    const CmcType type = static_cast<CmcType>(std::get<int32_t>(data_type_attribute.GetValue()));
    if (ConvertToCmcType<T>() != type) {cmc_err_msg("The template parameter of the Reader functionality does not match the data_type of the compressed variable.");}

    /* Get the corresponding mesh id */
    auto mesh_id_iter = cmc::nc::FindAttribute(attributes, mesh_id_attr);
    if (mesh_id_iter == attributes.end()) {cmc_err_msg("The ", mesh_id_attr, " attribute has not been found.");}
    const cmc::nc::Attribute& mesh_id_attribute = *mesh_id_iter;
    const int mesh_id = std::get<int32_t>(mesh_id_attribute.GetValue());

    /* Get the utilized compression scheme */
    auto compression_scheme_iter = cmc::nc::FindAttribute(attributes, compression_schema_attr);
    if (compression_scheme_iter == attributes.end()) {cmc_err_msg("The ", compression_schema_attr, " attribute has not been found.");}
    const cmc::nc::Attribute& compression_scheme_attribute = *compression_scheme_iter;
    const CompressionSchema compression_scheme = static_cast<cmc::CompressionSchema>(std::get<int32_t>(compression_scheme_attribute.GetValue()));

    /* Get the missing value of the variable */
    auto missing_val_iter = cmc::nc::FindAttribute(attributes, missing_value_attr);
    if (missing_val_iter == attributes.end()) {cmc_err_msg("The ", missing_value_attr, " attribute has not been found.");}
    const cmc::nc::Attribute& missing_val_attribute = *missing_val_iter;
    const T missing_value = std::get<T>(missing_val_attribute.GetValue());

    /* Get the initial data layout of the variable */
    auto initial_data_layout_iter = cmc::nc::FindAttribute(attributes, initial_data_layout_attr);
    if (initial_data_layout_iter == attributes.end()) {cmc_err_msg("The ", initial_data_layout_attr, " attribute has not been found.");}
    const cmc::nc::Attribute& initial_data_layout_attribute = *initial_data_layout_iter;
    const DataLayout initial_layout = static_cast<DataLayout>(std::get<int32_t>(initial_data_layout_attribute.GetValue()));

    /* Get the pre-compression data layout of the variable */
    auto pre_compression_layout_iter = cmc::nc::FindAttribute(attributes, pre_compression_layout_attr);
    if (pre_compression_layout_iter == attributes.end()) {cmc_err_msg("The ", pre_compression_layout_attr, " attribute has not been found.");}
    const cmc::nc::Attribute& pre_compression_layout_attribute = *pre_compression_layout_iter;
    const DataLayout pre_compression_layout = static_cast<DataLayout>(std::get<int32_t>(pre_compression_layout_attribute.GetValue()));

    /* Get the global context information */
    auto global_context_information_iter = cmc::nc::FindAttribute(attributes, global_context_information_attr);
    if (global_context_information_iter == attributes.end()) {cmc_err_msg("The ", global_context_information_attr, " attribute has not been found.");}
    const cmc::nc::Attribute& global_context_information_attribute = *global_context_information_iter;
    const int32_t global_context_information = std::get<int32_t>(global_context_information_attribute.GetValue());

    /* Get the number of compression iterations */
    auto num_compression_iterations_iter = cmc::nc::FindAttribute(attributes, num_compression_iterations_attr);
    if (num_compression_iterations_iter == attributes.end()) {cmc_err_msg("The ", num_compression_iterations_attr, " attribute has not been found.");}
    const cmc::nc::Attribute& num_compression_iterations_attribute = *num_compression_iterations_iter;
    const int num_compression_iterations = std::get<int32_t>(num_compression_iterations_attribute.GetValue());

    /* Generate the mesh-variable name */
    const std::string mesh_name = GenerateMeshVariableName(mesh_id);

    /* Get the attributes for the mesh variable */
    std::vector<cmc::nc::Attribute> mesh_attributes = reader.ReadVariableAttributes(mesh_name);

    /* Get the number of compression iterations */
    auto mesh_num_compression_iterations_iter = cmc::nc::FindAttribute(mesh_attributes, num_compression_iterations_attr);
    if (mesh_num_compression_iterations_iter == mesh_attributes.end()) {cmc_err_msg("The ", num_compression_iterations_attr, " attribute of the mesh has not been found.");}
    const cmc::nc::Attribute& mesh_num_compression_iterations_attribute = *mesh_num_compression_iterations_iter;
    const int mesh_num_compression_iterations = std::get<int32_t>(mesh_num_compression_iterations_attribute.GetValue());

    if (mesh_num_compression_iterations != num_compression_iterations) {cmc_err_msg("The data variable and the mesh variable do not have the same amount of compression iterations.");}

    /* Get the flag whether refinement bits are stored */
    auto ref_bit_storgae_iter = cmc::nc::FindAttribute(mesh_attributes, are_refinement_bits_stored_attr);
    if (ref_bit_storgae_iter == mesh_attributes.end()) {cmc_err_msg("The ", are_refinement_bits_stored_attr, " attribute has not been found.");}
    const cmc::nc::Attribute& ref_bit_storgae_attribute = *ref_bit_storgae_iter;
    const bool are_refinement_bits_stored = std::get<int32_t>(ref_bit_storgae_attribute.GetValue());

    /* Read the global compression streams for all levels */
    std::vector<uint8_t> encoded_data = reader.ReadVariableData<uint8_t>(var_name);
    std::vector<uint8_t> encoded_mesh = reader.ReadVariableData<uint8_t>(mesh_name);
    
    /* We set up the attributes for the variable, the GeoDomain of the variable is encoded within the mesh and will be set later when decoded */
    VariableAttributes<T> decompressed_var_attributes(GeoDomain(), missing_value, initial_layout, pre_compression_layout, global_context_information);

    /* Invoke the correct decompressor */
    switch (compression_scheme)
    {
        case CompressionSchema::EmbeddedPrefixExtraction:
            return std::make_unique<lossless::embedded::prefix::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh), std::move(decompressed_var_attributes), are_refinement_bits_stored, comm_);
        break;
        case CompressionSchema::EmbeddedPrefixExtractionPlainSuffixes:
            cmc_debug_msg("Embedded Prefix Extraction Decompression with Plain Suffixes is instantiated.");
            return std::make_unique<lossless::embedded::prefix::plain_suffix::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh), std::move(decompressed_var_attributes), are_refinement_bits_stored, comm_);
        break;
        case CompressionSchema::EmbeddedMultiResExtraction:
            cmc_debug_msg("Embedded MultiRes Decompression is instantiated.");
            return std::make_unique<lossless::embedded::multi_res::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh), std::move(decompressed_var_attributes), are_refinement_bits_stored, comm_);
        break;
        case CompressionSchema::_TestEmbeddedPCP4Extraction:
            cmc_debug_msg("Test MultiRes PCP4 Decompression is instantiated.");
            return std::make_unique<lossless::embedded::test_pcp4::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh), std::move(decompressed_var_attributes), are_refinement_bits_stored, comm_);
        break;
        case CompressionSchema::EmbeddedQuantizedPrefixExtraction:
            cmc_debug_msg("Embedded Quantized PrefixAMR Decompression is instantiated.");
            return std::make_unique<cmc::lossy::embedded::prefix::quantization::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh), std::move(decompressed_var_attributes), are_refinement_bits_stored, comm_);
        break;
        default:
            cmc_err_msg("The compression schema of the compressed variable is not recognized for an embedded variable.");
            return std::make_unique<lossless::embedded::prefix::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh), std::move(decompressed_var_attributes), are_refinement_bits_stored, comm_);
    }
}

#endif

template <typename T>
std::unique_ptr<cmc::patch::IPatchDecompressionVariable<T>>
Reader::ReadPatchVariableForDecompression(const std::string& var_name)
{
    /* First, we check the attributes of the variable */
    std::vector<cmc::nc::Attribute> attributes = reader.ReadVariableAttributes(var_name);

    /* Get the data type of the compressed variable */
    auto data_type_iter = cmc::nc::FindAttribute(attributes, data_type_attr);
    if (data_type_iter == attributes.end()) {cmc_err_msg("The ", data_type_attr, " attribute has not been found.");}
    const cmc::nc::Attribute& data_type_attribute = *data_type_iter;
    const CmcType type = static_cast<CmcType>(std::get<int32_t>(data_type_attribute.GetValue()));
    if (ConvertToCmcType<T>() != type) {cmc_err_msg("The template parameter of the Reader functionality does not match the data_type of the compressed variable.");}

    /* Get the utilized compression scheme */
    auto compression_scheme_iter = cmc::nc::FindAttribute(attributes, compression_schema_attr);
    if (compression_scheme_iter == attributes.end()) {cmc_err_msg("The ", compression_schema_attr, " attribute has not been found.");}
    const cmc::nc::Attribute& compression_scheme_attribute = *compression_scheme_iter;
    const CompressionSchema compression_scheme = static_cast<cmc::CompressionSchema>(std::get<int32_t>(compression_scheme_attribute.GetValue()));

    /* Get the initial data layout of the variable */
    auto initial_data_layout_iter = cmc::nc::FindAttribute(attributes, initial_data_layout_attr);
    if (initial_data_layout_iter == attributes.end()) {cmc_err_msg("The ", initial_data_layout_attr, " attribute has not been found.");}
    const cmc::nc::Attribute& initial_data_layout_attribute = *initial_data_layout_iter;
    const DataLayout initial_layout = static_cast<DataLayout>(std::get<int32_t>(initial_data_layout_attribute.GetValue()));

    /* Get the number of compression iterations */
    auto num_compression_iterations_iter = cmc::nc::FindAttribute(attributes, num_compression_iterations_attr);
    if (num_compression_iterations_iter == attributes.end()) {cmc_err_msg("The ", num_compression_iterations_attr, " attribute has not been found.");}
    const cmc::nc::Attribute& num_compression_iterations_attribute = *num_compression_iterations_iter;
    const int num_compression_iterations = std::get<int32_t>(num_compression_iterations_attribute.GetValue());

    /* Get the dimensionality of the variable */
    auto dim_iter = cmc::nc::FindAttribute(attributes, dim_attr);
    if (dim_iter == attributes.end()) {cmc_err_msg("The ", dim_attr, " attribute has not been found.");}
    const cmc::nc::Attribute& dim_attribute = *dim_iter;
    const int dimensionality = std::get<int32_t>(dim_attribute.GetValue());

    /* Get the levels of the dimension length pyramid */
    auto dim_lengths_pyra_lvls_iter = cmc::nc::FindAttribute(attributes, dim_lengths_pyra_lvls_attr);
    if (dim_lengths_pyra_lvls_iter == attributes.end()) {cmc_err_msg("The ", dim_lengths_pyra_lvls_attr, " attribute has not been found.");}
    const cmc::nc::Attribute& dim_lengths_pyra_lvls_attribute = *dim_lengths_pyra_lvls_iter;
    const int pyra_dim_lenghts_lvls = std::get<int32_t>(dim_lengths_pyra_lvls_attribute.GetValue());

    /* Allocate the dimension length pyramid */
    std::vector<std::vector<size_t>> dim_lengths_pyramid;
    dim_lengths_pyramid.reserve(pyra_dim_lenghts_lvls);

    /* Get the dimension lengths of the levels */
    for (int lvl_idx{0}; lvl_idx < pyra_dim_lenghts_lvls; ++lvl_idx)
    {
        dim_lengths_pyramid.emplace_back();
        dim_lengths_pyramid.back().reserve(dimensionality);

        for (int dim_idx{1}; dim_idx <= dimensionality; ++dim_idx)
        {
            const std::string lvl_dim_attr_name = "dim" + std::to_string(dim_idx) + "_" + std::to_string(lvl_idx);
            auto dim_lengths_lvl_dim_iter = cmc::nc::FindAttribute(attributes, lvl_dim_attr_name);
            if (dim_lengths_lvl_dim_iter == attributes.end()) {cmc_err_msg("The ", lvl_dim_attr_name, " attribute has not been found.");}
            const cmc::nc::Attribute& dim_lengths_lvl_dim_attribute = *dim_lengths_lvl_dim_iter;
            const int dim_length_lvl = std::get<int32_t>(dim_lengths_lvl_dim_attribute.GetValue());

            dim_lengths_pyramid.back().push_back(static_cast<size_t>(dim_length_lvl));
        }
    }

    GeoDomain init_domain;
    const std::vector<Dimension> dims_vec = GetDimensionVectorFromLayout(static_cast<DataLayout>(initial_layout));
    cmc_assert(dims_vec.size() == dim_lengths_pyramid.back().size());

    int dim_int_idx{0};
    for (auto dim_iter = dims_vec.begin(); dim_iter != dims_vec.end(); ++dim_iter, ++dim_int_idx)
    {
        DimensionInterval dim(*dim_iter, 0, dim_lengths_pyramid.back()[dim_int_idx]);
        init_domain.UpdateDimension(dim);
    }

    /* Read the global compression streams for all levels */
    std::vector<uint8_t> encoded_data = reader.ReadVariableData<uint8_t>(var_name);
    cmc_debug_msg("Size of encoded data stream: ", encoded_data.size());
    
    /* Invoke the correct decompressor */
    switch (compression_scheme)
    {
        case CompressionSchema::PatchPrefixExtractionPlainSuffixes:
            switch (dimensionality)
            {
                case 3:
                    return std::make_unique<patch::decompression::prefix::plain_suffix::DecompressionVariable<T, 3>>(std::move(encoded_data), var_name, std::move(dim_lengths_pyramid), init_domain, initial_layout, pyra_dim_lenghts_lvls);
                    break;
                default:
                    cmc_err_msg("The dimensionality of the compressed variable is not recognized for a patch variable.");
                    return std::make_unique<patch::decompression::prefix::plain_suffix::DecompressionVariable<T, 3>>(std::move(encoded_data), var_name, std::move(dim_lengths_pyramid), init_domain, initial_layout, pyra_dim_lenghts_lvls);
            }
        break;
        case CompressionSchema::PatchMultiResExtraction:
            switch (dimensionality)
            {
                case 3:
                    return std::make_unique<patch::decompression::multi_res::DecompressionVariable<T, 3>>(std::move(encoded_data), var_name, std::move(dim_lengths_pyramid), init_domain, initial_layout, pyra_dim_lenghts_lvls);
                    break;
                default:
                    cmc_err_msg("The dimensionality of the compressed variable is not recognized for a patch variable.");
                    return std::make_unique<patch::decompression::multi_res::DecompressionVariable<T, 3>>(std::move(encoded_data), var_name, std::move(dim_lengths_pyramid), init_domain, initial_layout, pyra_dim_lenghts_lvls);
            }
        break;
        default:
        cmc_err_msg("The compression schema of the compressed variable is not recognized for a patch variable.");
            return std::make_unique<patch::decompression::prefix::plain_suffix::DecompressionVariable<T, 3>>(std::move(encoded_data), var_name, std::move(dim_lengths_pyramid), init_domain, initial_layout, pyra_dim_lenghts_lvls);
    }
}

}


#endif /* !CMC_DECOMPRESSION_NC_INPUT_HXX */
