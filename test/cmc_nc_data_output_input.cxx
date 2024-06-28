#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "netcdf/cmc_nc_io.hxx"
#include "netcdf/cmc_nc_writer.hxx"
#include "netcdf/cmc_nc_reader.hxx"

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
    const std::vector<double> ddata{1.0,2.0,4.0,8.0,16.0,32.0};

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

    /* Read the data back in and check for compliance */
    cmc::NcReader reader{file_name};

    /* Define a general hyperslab for the data to be read */
    cmc::GeneralHyperslab ghyperslab;
    ghyperslab.start_values = std::vector<size_t>{0, 0};
    ghyperslab.count_values = std::vector<size_t>{lat_length, lon_length};

    /* Set the variables to be read */
    reader.StashVariableForReading("test_var", ghyperslab);

    /* Get the variable with the defined hyperslab and the global attributes of the file */
    auto [read_variables, read_global_attributes] = reader.ReadVariables();

    /* Check some general information */
    cmc::ExpectEQ(reader.GetFileName(), file_name);
    cmc::ExpectEQ(reader.GetNumberOfVariables(), 1);
    cmc::ExpectEQ(reader.GetNumberOfGlobalAttributes(), 1);
    cmc::ExpectEQ(read_global_attributes.size(), static_cast<size_t>(1));

    /* Check the dimensions of the variable */
    const std::vector<cmc::NcDimension> test_var_dims = read_variables.front().GetDimensions();
    cmc::ExpectTrue(!test_var_dims[0].GetName().compare("lat1"));
    cmc::ExpectTrue(!test_var_dims[1].GetName().compare("lon1"));
    cmc::ExpectEQ(test_var_dims[0].GetLength(), static_cast<size_t>(lat_length));
    cmc::ExpectEQ(test_var_dims[1].GetLength(), static_cast<size_t>(lon_length));
    
    /* Check the variable's attributes */
    const std::vector<cmc::NcAttribute> test_var_atts = read_variables.front().GetAttributes();
    cmc::ExpectEQ(test_var_atts.size(), static_cast<size_t>(2));
    cmc::ExpectTrue(!test_var_atts[0].GetName().compare("attr"));
    cmc::ExpectTrue(!test_var_atts[1].GetName().compare("id"));

    cmc::CmcUniversalType attr_val = test_var_atts[0].GetValue();
    cmc::ExpectEQ(std::get<float>(attr_val), static_cast<float>(1.5));

    cmc::CmcUniversalType id_val = test_var_atts[1].GetValue();
    cmc::ExpectEQ(std::get<int>(id_val), testid);

    /* Check the variable's data */
    cmc::NcSpecificVariable<double> spec_var = read_variables.front().DetachVariable<double>();
    const std::vector<double>& var_data = spec_var.GetData();
    int data_accessor = 0;
    for (auto data_iter = var_data.begin(); data_iter != var_data.end(); ++data_iter, ++data_accessor)
    {
        cmc::ExpectEQ(var_data[data_accessor], ddata[data_accessor]);
    }

    /* Check the global attribute */
    cmc::ExpectEQ(read_global_attributes.size(), static_cast<size_t>(1));
    cmc::ExpectTrue(!read_global_attributes.front().GetName().compare("global_attr"));
    cmc::CmcUniversalType global_attribute = read_global_attributes.front().GetValue();
    cmc::ExpectEQ(std::get<int16_t>(global_attribute), static_cast<int16_t>(255));

    /* Finalize cmc */
    cmc_finalize();

    return cmc::CMC_TEST_SUCCESS;
}