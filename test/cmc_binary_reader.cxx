#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_binary_reader.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_geo_domain.hxx"

#include <cstddef>
#include <vector>
#include <cstdio>
#include <variant>

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();

    {
    /* Create a vector of values */
    std::vector<int32_t> vals{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};

    /* Write them out in a file */
    const std::string file_name("_cmc_binaray_reader_test_file.bin");
    FILE* file = fopen(file_name.c_str(), "wb");
    (void) fwrite(vals.data(), sizeof(int32_t), vals.size(), file);
    (void) fclose(file);

    /* Create a binary reader to the file */
    cmc::bin_reader::Reader reader(file_name);

    /* Specify some features of an inputVar that will be created and hold the binary data */
    const std::string var_name = "test_var";
    const int var_id = 0;
    const size_t num_elements = vals.size();
    const cmc::GeoDomain domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 8),
                                cmc::DimensionInterval(cmc::Dimension::Lat, 0, 2));
    const cmc::CmcUniversalType missing_value = int32_t{-1};
    const cmc::DataLayout layout = cmc::DataLayout::Lat_Lon;

    /* Create the variable with the data from the bianry reader */
    cmc::InputVar input_variable_full_domain = reader.CreateVariableFromBinaryData(cmc::CmcType::Int32_t, var_name, var_id, num_elements, missing_value, layout, domain);

    /* We get the internal data structures and will compare the data */
    const cmc::CmcInputVariable& input_var_variant = input_variable_full_domain.GetInternalVariant();
    const cmc::InputVariable<int32_t>& input_variable = std::get<cmc::InputVariable<int32_t>>(input_var_variant);
    const std::vector<int32_t>& var_data = input_variable.GetDataForReading();

    cmc::ExpectTrue(vals.size() == var_data.size());

    /* Comapre the read in values with the initial values */
    for (size_t idx = 0; idx < var_data.size(); ++idx)
    {
        cmc::ExpectTrue(vals[idx] == var_data[idx]);
    }

    /* The initial data of a certain subdomain that we want to read */
    std::vector<int32_t> sub_domain_data{3,4,5,6,11,12,13,14};

    /* The subdomain corresponding to the data we want to have */
    const cmc::GeoDomain sub_domain(cmc::DimensionInterval(cmc::Dimension::Lon, 3, 7),
                                    cmc::DimensionInterval(cmc::Dimension::Lat, 0, 2));

    /* Create an InputVar with only the data of the subdomain */
    cmc::InputVar input_variable_sub_domain = reader.CreateSubDomainVariableFromBinaryData(cmc::CmcType::Int32_t, var_name, var_id, missing_value, layout, domain, sub_domain);

    /* We get the internal data structures and will compare the data */
    const cmc::CmcInputVariable& input_var_sd_variant = input_variable_sub_domain.GetInternalVariant();
    const cmc::InputVariable<int32_t>& input_variable_sd = std::get<cmc::InputVariable<int32_t>>(input_var_sd_variant);
    const std::vector<int32_t>& var_data_sd = input_variable_sd.GetDataForReading();

    cmc::ExpectTrue(sub_domain_data.size() == var_data_sd.size());

    /* Compare the extarcted values */
    for (size_t idx = 0; idx < var_data_sd.size(); ++idx)
    {
        cmc::ExpectTrue(sub_domain_data[idx] == var_data_sd[idx]);
    }

    /* Delete the binary file */
    (void) std::remove(file_name.c_str());

    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
