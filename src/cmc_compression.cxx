#include "cmc_compression.hxx"
#include "netcdf/cmc_nc_reader.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "lossy/cmc_prefix_decompression.hxx"

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
    NcReader reader{file_name_};

    /* Read the global attributes of the file */
    std::vector<NcAttribute> global_attributes = reader.ReadGlobalAttrtibutes();

    /* Check for a given attribute indicating the compression which has been used */
    auto comp_scheme_iter = FindAttribute(global_attributes, kCompressionSchemeAttrName);

    /* Check if the attribute has been found */
    if (comp_scheme_iter == global_attributes.end())
    {
        cmc_err_msg("There is no 'compression_scheme' attribute attached to this file. Therefore, no decompression is applicable.");
    } else 
    {
        /* Get the attribute and check it's value */
        const NcAttribute& comp_scheme_attr = *comp_scheme_iter;

        /* Get the value of the compression scheme */
        const CompressionScheme scheme = static_cast<CompressionScheme>(std::get<CompressionSchemeType>(comp_scheme_attr.GetValue()));

        switch (scheme)
        {
            case CompressionScheme::AdaptiveCoarsening:
                cmc_err_msg("This decompressor is not yet integrated in the general usage functions");
            break;
            case CompressionScheme::PrefixExtraction:
                decompressor_ = std::make_unique<prefix::Decompressor>();
            break;
            case CompressionScheme::NoCompression:
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
    NcReader reader{file_name};

    /* Read the global attributes of the file */
    std::vector<NcAttribute> global_attributes = reader.ReadGlobalAttrtibutes();

    /* Check for a given attribute indicating the compression which has been used */
    auto comp_scheme_iter = FindAttribute(global_attributes, kCompressionSchemeAttrName);

    /* Check if the attribute has been found */
    if (comp_scheme_iter == global_attributes.end())
    {
        cmc_global_msg("There is no 'compression_scheme' attribute attached to this file. Therefore, no decompression is applicable.");
        return std::vector<OutputVar>();
    } else 
    {
        /* Get the attribute and check it's value */
        const NcAttribute& comp_scheme_attr = *comp_scheme_iter;

        /* Get the value of the compression scheme */
        const CompressionScheme scheme = static_cast<CompressionScheme>(std::get<CompressionSchemeType>(comp_scheme_attr.GetValue()));

        switch (scheme)
        {
            case CompressionScheme::AdaptiveCoarsening:
                return GeneralDecompressionByAdaptiveCoarsening(); 
            break;
            case CompressionScheme::PrefixExtraction:
                return GeneralDecompressionByAdaptiveCoarsening();
            break;
            case CompressionScheme::AdaptiveCoarseningPlusPrefixExtraction:
                return GeneralDecompressionByAdaptiveCoarsening();
            break;
            case CompressionScheme::NoCompression:
                cmc_global_msg("No compression has been applied. Therefore, no decompression is applicable.");
                return std::vector<OutputVar>();
            break;
            default:
                return std::vector<OutputVar>();
        }
    }
}

}
