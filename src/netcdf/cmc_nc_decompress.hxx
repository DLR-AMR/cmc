#ifndef CMC_NC_DECOMPRESS_HXX
#define CMC_NC_DECOMPRESS_HXX
/**
 * @file cmc_nc_decompress.hxx
 */

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_coordinate_array.hxx"
#include "utilities/cmc_hyperslab.hxx"
#include "utilities/cmc_output_variable.hxx"

#include "netcdf/cmc_netcdf.hxx"

#ifdef CMC_WITH_NETCDF
#include "netcdf.h"
#endif
#ifdef CMC_WITH_NETCDF_PAR
#include "netcdf_par.h"
#endif

#include <vector>
#include <string>

namespace cmc::nc
{

#if 0
enum CompressionFormat {UndefinedCompressionFormat = -1, AMR = 0, PrefAMR = 1, CombinedAMR = 2};

class NcDecompress
{
public:
    NcDecompress() = delete;

    NcDecompress(const std::string& path_to_file, const OpeningMode mode, const MPI_Comm comm = MPI_COMM_WORLD)
    {
        NcOpen(path_to_file, mode, comm);
    };

    ~NcDecompress() = default;

    //OutputVar DecompressVariable(const int variable_id); //zu spezifisch hier
    

    void CloseFileHandle();

private:
    void NcOpen(const std::string& path_to_file, const OpeningMode mode, const MPI_Comm comm);
    
    int ncid_;

    CompressionFormat compressed_format{CompressionFormat::UndefinedCompressionFormat};

    bool _file_has_been_closed_{false};
    
};

#endif

}

#endif /* !CMC_NC_DECOMPRESS_HXX */
