#include "cmc_compression.hxx"
#include "netcdf/cmc_nc_reader.hxx"

namespace cmc
{

std::vector<OutputVar>
Decompress(const std::string& file_name)
{
    /* Create a reader for the given file */
    NcReader reader{file_name};

    /* Read the global attributes of the file */
    std::vector<NcAttribute> global_attributes = reader.ReadGlobalAttrtibutes();

    /* Check for a given attribute indicating the compression which has been used */
    //auto comp_scheme_iter = FindAttribute(global_attributes, "compression_scheme");


}

}
