#include "netcdf/cmc_nc_decompress.hxx"
#include "utilities/cmc_log_functions.h"

namespace cmc
{


void
NcDecompress::NcOpen(const std::string& path_to_file, const NcOpeningMode mode, const MPI_Comm comm)
{
    #ifdef CMC_WITH_NETCDF
    #ifdef CMC_WITH_NETCDF_PAR
    if (mode == NcOpeningMode::Parallel)
    {
        ncid_ = NcOpenParallel(path_to_file.c_str(), comm);
    } else

    #endif /* CMC_WITH_NETCDF_PAR */
    {
        if (mode == NcOpeningMode::Parallel)
        {
            cmc_warn_msg("The netCDF file is ought to be opened for parallel access, although the parallel functionality is not accessible.");
            cmc_warn_msg("The file is opened for serial access.");
        }
        ncid_ = NcOpenSerial(path_to_file.c_str());
    }

    #else
    cmc_err_msg("CMC is not linked against netCDF. Please recompile the build with netCDF.\n");
    return CMC_ERR;
    #endif
}

void
NcDecompress::CloseFileHandle()
{
    int err = nc_close(ncid_);
    NcCheckError(err);
    _file_has_been_closed_ = true;
}


#if 0
OutputVar
DecompressVariable(const int variable_id)
{
    
}
#endif

void InquireCompressionFile()
{
    //Inquire num dimensions
    //Inquire num variables
    //Inquire global attributes -> Check which Compression has been used
}


}
