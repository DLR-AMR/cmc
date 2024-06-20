#include "cmc.h"
#include "test/cmc_test.hxx"
#include "netcdf/cmc_nc_io.hxx"
#include "netcdf/cmc_nc_writer.hxx"

#include <vector>

int
main(void)
{
    /* Initialize cmc */
    cmc_initialize();

    /* Initialize a variable for the netcdf output */
    std::vector<cmc::NcVariable> variables;

    /* Define a variable of a specified data type */
    const int testid = 1;
    cmc::NcSpecificVariable<double> dvar{"test_var", testid};
    
    /* Define a global domain */
    const cmc::DomainIndex lon_length = 3;
    const cmc::DomainIndex lat_length = 2;
    cmc::GeoDomain global_domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, lon_length),
                                 cmc::DimensionInterval(cmc::Dimension::Lat, 0, lat_length));

    /* The lcoa domain of the exemplary variable fully filly the global domain */
    cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, lon_length),
                             cmc::DimensionInterval(cmc::Dimension::Lat, 0, lat_length));

    /* Some examplary data */
    std::vector<double> ddata{1.0,2.0,4.0,8.0,16.0,32.0};

    dvar.SetGlobalDomain(global_domain);
    dvar.PushBack(ddata, hyperslab);
    dvar.SetDataLayout(cmc::DataLayout::Lat_Lon);

    /* Create a variable attribute */
    cmc::NcAttribute fattr{"attr", static_cast<float>(1.5)};

    std::vector<cmc::NcAttribute> dvar_attributes;
    dvar_attributes.push_back(fattr);

    /* Create a global attribute */
    cmc::NcAttribute global_attr{"global_attr", static_cast<int16_t>(255)};

    const std::string file_name = "test_nc_output_file.nc";
    const int netcdf_format = NC_64BIT_OFFSET;

    /* Create an output writer */
    cmc::NcWriter writer{file_name, netcdf_format};

    /* Add the variable and the global attribute */
    writer.AddVariable(dvar, dvar_attributes);

    /* Add the global attribute */
    writer.AddGlobalAttribute(global_attr);

    /* Write the variable and the attribute to the file */
    writer.Write();

    /* Finalize cmc */
    cmc_finalize();

    return cmc::CMC_TEST_SUCCESS;
}