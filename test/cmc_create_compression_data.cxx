#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "lossy/cmc_amr_lossy_compression_settings.hxx"
#include "lossy/cmc_amr_lossy_compression.hxx"

#include <numeric>

int
main(void)
{
    /* Initialize cmc */
    cmc_initialize();
    
    /* Create some exemlpary variable with data */
    cmc::InputVariable<double> test_var("ex_double_data", 0, cmc::DataLayout::Lon_Lat);
    const cmc::DomainIndex lon_length = 10;
    const cmc::DomainIndex lat_length = 5;
    cmc::GeoDomain global_domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, lon_length),
                                 cmc::DimensionInterval(cmc::Dimension::Lat, 0, lat_length));
    test_var.SetGlobalDomain(global_domain);
    test_var.SetMissingValue(-2.0);
    std::vector<double> double_data(lon_length * lat_length);
    std::iota(double_data.begin(), double_data.end(), 100.0);
    cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, lon_length),
                             cmc::DimensionInterval(cmc::Dimension::Lat, 0, lat_length));
    test_var.PushBack(double_data, std::move(hyperslab));

    /* Create compression settings */
    cmc::CompressionSettings settings;
    const double abs_max_err = 1.0;
    settings.SetAbsoluteErrorCriterion(abs_max_err, cmc::kErrorCriterionHoldsForAllVariables);
    

    std::vector<cmc::InputVar> variables_to_compress;
    variables_to_compress.push_back(std::move(test_var));

    /* Create the compression data */
    [[maybe_unused]] cmc::CompressionData compression_data(variables_to_compress, settings);

    /* Finalize cmc */
    cmc_finalize();

    return cmc::CMC_TEST_SUCCESS;
}
