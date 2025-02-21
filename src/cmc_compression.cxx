#include "cmc_compression.hxx"
#include "netcdf/cmc_nc_reader.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "decompression/cmc_prefix_decompression.hxx"

#include <vector>

namespace cmc
{

/* Start the compression by setting up the comrpessor */
void
Compressor::Start()
{
    compressor_->Setup();
}

void
Compressor::Compress()
{
    
}

void
Compressor::WriteCompressedData(const std::string& file_name)
{

}

/* Check which decompressor to use */
void
Decompressor::Start()
{
    /* Create a reader for the given file */
    nc::Reader reader{file_name_};

    /* Read the global attributes of the file */
    std::vector<nc::Attribute> global_attributes = reader.ReadGlobalAttrtibutes();

    /* Check for a given attribute indicating the compression which has been used */
    auto comp_scheme_iter = FindAttribute(global_attributes, nc::kCompressionSchemeAttrName);

    /* Check if the attribute has been found */
    if (comp_scheme_iter == global_attributes.end())
    {
        cmc_err_msg("There is no 'compression_scheme' attribute attached to this file. Therefore, no decompression is applicable.");
    } else 
    {
        /* Get the attribute and check it's value */
        const nc::Attribute& comp_scheme_attr = *comp_scheme_iter;

        /* Get the value of the compression scheme */
        const nc::CompressionScheme scheme = static_cast<nc::CompressionScheme>(std::get<nc::CompressionSchemeType>(comp_scheme_attr.GetValue()));

        switch (scheme)
        {
            case nc::CompressionScheme::AdaptiveCoarsening:
                cmc_err_msg("This decompressor is not yet integrated in the general usage functions");
            break;
            case nc::CompressionScheme::PrefixExtraction:
                decompressor_ = std::make_unique<prefix::Decompressor>();
            break;
            case nc::CompressionScheme::NoCompression:
                cmc_global_msg("No compression has been applied. Therefore, no decompression is applicable.");
            break;
            default:
                cmc_err_msg("No compression scheme is recognized. Therefore, a decompression is not possible.");
        }
    }
}

void
Decompressor::Decompress()
{
    /* Decompress the given file */
    decompressor_->Decompress();
}

static
std::vector<OutputVar>
GeneralDecompressionByAdaptiveCoarsening()
{
    return std::vector<OutputVar>();
}

static
std::vector<OutputVar>
GeneralDecompressionByPrefixExtraction()
{
    return std::vector<OutputVar>();
}


std::vector<OutputVar>
Decompress(const std::string& file_name)
{
    /* Create a reader for the given file */
    nc::Reader reader{file_name};

    /* Read the global attributes of the file */
    std::vector<nc::Attribute> global_attributes = reader.ReadGlobalAttrtibutes();

    /* Check for a given attribute indicating the compression which has been used */
    auto comp_scheme_iter = FindAttribute(global_attributes, nc::kCompressionSchemeAttrName);

    /* Check if the attribute has been found */
    if (comp_scheme_iter == global_attributes.end())
    {
        cmc_global_msg("There is no 'compression_scheme' attribute attached to this file. Therefore, no decompression is applicable.");
        return std::vector<OutputVar>();
    } else 
    {
        /* Get the attribute and check it's value */
        const nc::Attribute& comp_scheme_attr = *comp_scheme_iter;

        /* Get the value of the compression scheme */
        const nc::CompressionScheme scheme = static_cast<nc::CompressionScheme>(std::get<nc::CompressionSchemeType>(comp_scheme_attr.GetValue()));

        switch (scheme)
        {
            case nc::CompressionScheme::AdaptiveCoarsening:
                return GeneralDecompressionByAdaptiveCoarsening(); 
            break;
            case nc::CompressionScheme::PrefixExtraction:
                return GeneralDecompressionByAdaptiveCoarsening();
            break;
            case nc::CompressionScheme::AdaptiveCoarseningPlusPrefixExtraction:
                return GeneralDecompressionByAdaptiveCoarsening();
            break;
            case nc::CompressionScheme::NoCompression:
                cmc_global_msg("No compression has been applied. Therefore, no decompression is applicable.");
                return std::vector<OutputVar>();
            break;
            default:
                return std::vector<OutputVar>();
        }
    }
}

}
