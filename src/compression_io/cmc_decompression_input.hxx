#ifndef CMC_DECOMPRESSION_INPUT_HXX
#define CMC_DECOMPRESSION_INPUT_HXX

#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_iface_abstract_embedded_byte_decompression_variable.hxx"
#include "utilities/cmc_iface_abstract_byte_decompression_variable.hxx"
#include "compression_io/cmc_compression_attr_names.hxx"

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

namespace cmc::compression_io
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

    template <typename T> std::unique_ptr<decompression::AbstractByteDecompressionVariable<T>> ReadVariableForDecompression(const std::string& var_name);
    template <typename T> std::unique_ptr<cmc::IEmbeddedByteDecompressionVariable<T>> ReadEmbeddedVariableForDecompression(const std::string& var_name);

private:
    const std::string file_name_;
    const MPI_Comm comm_;
    cmc::nc::Reader reader;
};

template <typename T>
std::unique_ptr<decompression::AbstractByteDecompressionVariable<T>>
Reader::ReadVariableForDecompression(const std::string& var_name)
{
    cmc_debug_msg("try to read: ", var_name);
    /* First, we check the attributes of the variable */
    std::vector<nc::Attribute> attributes = reader.ReadVariableAttributes(var_name);

    /* Get the data type of the compressed variable */
    auto data_type_iter = nc::FindAttribute(attributes, data_type_attr);
    if (data_type_iter == attributes.end()) {cmc_err_msg("The ", data_type_attr, " attribute has not been found.");}
    const nc::Attribute& data_type_attribute = *data_type_iter;
    const CmcType type = static_cast<CmcType>(std::get<int>(data_type_attribute.GetValue()));
    if (ConvertToCmcType<T>() != type) {cmc_err_msg("The template parameter of the Reader functionality does not match the data_type of the compressed variable.");}

    /* Get the utilized compression scheme */
    auto compression_scheme_iter = nc::FindAttribute(attributes, compression_schema_attr);
    if (compression_scheme_iter == attributes.end()) {cmc_err_msg("The ", compression_schema_attr, " attribute has not been found.");}
    const nc::Attribute& compression_scheme_attribute = *compression_scheme_iter;
    const CompressionSchema compression_scheme = static_cast<cmc::CompressionSchema>(std::get<int>(compression_scheme_attribute.GetValue()));

    /* Get the corresponding mesh id */
    auto mesh_id_iter = nc::FindAttribute(attributes, mesh_id_attr);
    if (mesh_id_iter == attributes.end()) {cmc_err_msg("The ", mesh_id_attr, " attribute has not been found.");}
    const nc::Attribute& mesh_id_attribute = *mesh_id_iter;
    const int mesh_id = std::get<int>(mesh_id_attribute.GetValue());

    /* Get the number of compression iterations */
    auto num_compression_iterations_iter = nc::FindAttribute(attributes, num_compression_iterations_attr);
    if (num_compression_iterations_iter == attributes.end()) {cmc_err_msg("The ", num_compression_iterations_attr, " attribute has not been found.");}
    const nc::Attribute& num_compression_iterations_attribute = *num_compression_iterations_iter;
    const int num_compression_iterations = std::get<int>(num_compression_iterations_attribute.GetValue());

    /* Generate the mesh-variable name */
    const std::string mesh_name = GenerateMeshVariableName(mesh_id);

    cmc_debug_msg("try to read mesh: ", mesh_name);

    /* Get the attributes for the mesh variable */
    std::vector<nc::Attribute> mesh_attributes = reader.ReadVariableAttributes(mesh_name);

    /* Get the number of compression iterations */
    auto mesh_num_compression_iterations_iter = nc::FindAttribute(mesh_attributes, num_compression_iterations_attr);
    if (mesh_num_compression_iterations_iter == mesh_attributes.end()) {cmc_err_msg("The ", num_compression_iterations_attr, " attribute of the mesh has not been found.");}
    const nc::Attribute& mesh_num_compression_iterations_attribute = *mesh_num_compression_iterations_iter;
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
    std::vector<nc::Attribute> attributes = reader.ReadVariableAttributes(var_name);

    /* Get the data type of the compressed variable */
    auto data_type_iter = nc::FindAttribute(attributes, data_type_attr);
    if (data_type_iter == attributes.end()) {cmc_err_msg("The ", data_type_attr, " attribute has not been found.");}
    const nc::Attribute& data_type_attribute = *data_type_iter;
    const CmcType type = static_cast<CmcType>(std::get<int32_t>(data_type_attribute.GetValue()));
    if (ConvertToCmcType<T>() != type) {cmc_err_msg("The template parameter of the Reader functionality does not match the data_type of the compressed variable.");}

    /* Get the corresponding mesh id */
    auto mesh_id_iter = nc::FindAttribute(attributes, mesh_id_attr);
    if (mesh_id_iter == attributes.end()) {cmc_err_msg("The ", mesh_id_attr, " attribute has not been found.");}
    const nc::Attribute& mesh_id_attribute = *mesh_id_iter;
    const int mesh_id = std::get<int32_t>(mesh_id_attribute.GetValue());

    /* Get the utilized compression scheme */
    auto compression_scheme_iter = nc::FindAttribute(attributes, compression_schema_attr);
    if (compression_scheme_iter == attributes.end()) {cmc_err_msg("The ", compression_schema_attr, " attribute has not been found.");}
    const nc::Attribute& compression_scheme_attribute = *compression_scheme_iter;
    const CompressionSchema compression_scheme = static_cast<cmc::CompressionSchema>(std::get<int32_t>(compression_scheme_attribute.GetValue()));

    /* Get the missing value of the variable */
    auto missing_val_iter = nc::FindAttribute(attributes, missing_value_attr);
    if (missing_val_iter == attributes.end()) {cmc_err_msg("The ", missing_value_attr, " attribute has not been found.");}
    const nc::Attribute& missing_val_attribute = *missing_val_iter;
    const T missing_value = std::get<T>(missing_val_attribute.GetValue());

    /* Get the initial data layout of the variable */
    auto initial_data_layout_iter = nc::FindAttribute(attributes, initial_data_layout_attr);
    if (initial_data_layout_iter == attributes.end()) {cmc_err_msg("The ", initial_data_layout_attr, " attribute has not been found.");}
    const nc::Attribute& initial_data_layout_attribute = *initial_data_layout_iter;
    const DataLayout initial_layout = static_cast<DataLayout>(std::get<int32_t>(initial_data_layout_attribute.GetValue()));

    /* Get the pre-compression data layout of the variable */
    auto pre_compression_layout_iter = nc::FindAttribute(attributes, pre_compression_layout_attr);
    if (pre_compression_layout_iter == attributes.end()) {cmc_err_msg("The ", pre_compression_layout_attr, " attribute has not been found.");}
    const nc::Attribute& pre_compression_layout_attribute = *pre_compression_layout_iter;
    const DataLayout pre_compression_layout = static_cast<DataLayout>(std::get<int32_t>(pre_compression_layout_attribute.GetValue()));

    /* Get the global context information */
    auto global_context_information_iter = nc::FindAttribute(attributes, global_context_information_attr);
    if (global_context_information_iter == attributes.end()) {cmc_err_msg("The ", global_context_information_attr, " attribute has not been found.");}
    const nc::Attribute& global_context_information_attribute = *global_context_information_iter;
    const int32_t global_context_information = std::get<int32_t>(global_context_information_attribute.GetValue());

    /* Get the number of compression iterations */
    auto num_compression_iterations_iter = nc::FindAttribute(attributes, num_compression_iterations_attr);
    if (num_compression_iterations_iter == attributes.end()) {cmc_err_msg("The ", num_compression_iterations_attr, " attribute has not been found.");}
    const nc::Attribute& num_compression_iterations_attribute = *num_compression_iterations_iter;
    const int num_compression_iterations = std::get<int32_t>(num_compression_iterations_attribute.GetValue());

    /* Generate the mesh-variable name */
    const std::string mesh_name = GenerateMeshVariableName(mesh_id);

    /* Get the attributes for the mesh variable */
    std::vector<nc::Attribute> mesh_attributes = reader.ReadVariableAttributes(mesh_name);

    /* Get the number of compression iterations */
    auto mesh_num_compression_iterations_iter = nc::FindAttribute(mesh_attributes, num_compression_iterations_attr);
    if (mesh_num_compression_iterations_iter == mesh_attributes.end()) {cmc_err_msg("The ", num_compression_iterations_attr, " attribute of the mesh has not been found.");}
    const nc::Attribute& mesh_num_compression_iterations_attribute = *mesh_num_compression_iterations_iter;
    const int mesh_num_compression_iterations = std::get<int32_t>(mesh_num_compression_iterations_attribute.GetValue());

    if (mesh_num_compression_iterations != num_compression_iterations) {cmc_err_msg("The data variable and the mesh variable do not have the same amount of compression iterations.");}

    /* Get the flag whether refinement bits are stored */
    auto ref_bit_storgae_iter = nc::FindAttribute(mesh_attributes, are_refinement_bits_stored_attr);
    if (ref_bit_storgae_iter == mesh_attributes.end()) {cmc_err_msg("The ", are_refinement_bits_stored_attr, " attribute has not been found.");}
    const nc::Attribute& ref_bit_storgae_attribute = *ref_bit_storgae_iter;
    const bool are_refinement_bits_stored = std::get<int32_t>(ref_bit_storgae_attribute.GetValue());

    /* Read the global compression streams for all levels */
    std::vector<uint8_t> encoded_data = reader.ReadVariableData<uint8_t>(var_name);
    std::vector<uint8_t> encoded_mesh = reader.ReadVariableData<uint8_t>(mesh_name);
    
    /* We set up the attributes for the variable, the GeoDomain of the variable is encoded within the mesh and will be set later when decoded */
    VariableAttributes<T> decomrpessed_var_attributes(GeoDomain(), missing_value, initial_layout, pre_compression_layout, global_context_information);

    /* Invoke the correct decompressor */
    switch (compression_scheme)
    {
        case CompressionSchema::EmbeddedPrefixExtraction:
            return std::make_unique<lossless::embedded::prefix::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh), std::move(decomrpessed_var_attributes), are_refinement_bits_stored, comm_);
        break;
        case CompressionSchema::EmbeddedPrefixExtractionPlainSuffixes:
            cmc_debug_msg("\n\nPlain Suffix Version is instantiated \n\n");
            return std::make_unique<lossless::embedded::prefix::plain_suffix::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh), std::move(decomrpessed_var_attributes), are_refinement_bits_stored, comm_);
        break;
        case CompressionSchema::EmbeddedMultiResExtraction:
            //return std::make_unique<lossless::embedded::multi_res::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh));
            cmc_debug_msg("Embedded MultiRes Decompression is instantiated.");
            return std::make_unique<lossless::embedded::multi_res::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh), std::move(decomrpessed_var_attributes), are_refinement_bits_stored, comm_);
        break;
        case CompressionSchema::_TestEmbeddedPCP4Extraction:
            cmc_debug_msg("Test MultiRes PCP4 Decompression is instantiated.");
            return std::make_unique<lossless::embedded::test_pcp4::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh), std::move(decomrpessed_var_attributes), are_refinement_bits_stored, comm_);
        break;
        case CompressionSchema::EmbeddedQuantizedPrefixExtraction:
            cmc_debug_msg("Embedded Quantized PrefixAMR Decompression is instantiated.");
            return std::make_unique<cmc::lossy::embedded::prefix::quantization::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh), std::move(decomrpessed_var_attributes), are_refinement_bits_stored, comm_);
        break;
        default:
            cmc_err_msg("The compression schema of the compressed variable is not recognized for an embedded variable.");
            return std::make_unique<lossless::embedded::prefix::DecompressionVariable<T>>(var_name, std::move(encoded_data), std::move(encoded_mesh), std::move(decomrpessed_var_attributes), are_refinement_bits_stored, comm_);
    }
}

}


#endif /* !CMC_DECOMPRESSION_INPUT_HXX */
