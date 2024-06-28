#include "cmc.h"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_hyperslab.hxx"
#include "netcdf/cmc_netcdf.hxx"

int main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    {
        /* Specify the path of the netCDF file */
        const std::string file = "../../data/era5_reanalysis_data_16_11_23.nc";

        /* Create an object which interatcs with the file and opens it */
        cmc::NcData nc_data(file, cmc::NcOpeningMode::Serial);

        /* Define a hyperlsab corresponding to the data we want to read */
        cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 1440),
                                 cmc::DimensionInterval(cmc::Dimension::Lat, 0, 721));

        /* Inquire the hyperslab of data for the given variables */
        nc_data.InquireVariables(hyperslab, "t2m");

        /* Close the file, since we have gathered the data we wanted */
        nc_data.CloseFileHandle();
    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}
