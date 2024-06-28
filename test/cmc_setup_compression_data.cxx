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
    cmc::CmcInitialize();
    
    {
    
    /* Create a domain for the examplary data */
    const cmc::DomainIndex lon_length = 10;
    const cmc::DomainIndex lat_length = 5;
    cmc::GeoDomain global_domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, lon_length),
                                 cmc::DimensionInterval(cmc::Dimension::Lat, 0, lat_length));
    
    /* Create a examplary variable with double data */
    cmc::InputVariable<double> test_var("ex_double_data", 0, cmc::DataLayout::Lon_Lat);
    test_var.SetGlobalDomain(global_domain);
    test_var.SetMissingValue(-2.0);
    std::vector<double> double_data(lon_length * lat_length);
    std::iota(double_data.begin(), double_data.end(), 100.0);
    cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, lon_length),
                             cmc::DimensionInterval(cmc::Dimension::Lat, 0, lat_length));
    test_var.PushBack(double_data, std::move(hyperslab));

    /* Create a examplary variable with float data */
    cmc::InputVariable<float> test_var2("ex_float_data", 1, cmc::DataLayout::Lon_Lat);
    const float var2_missing_value = -10;
    test_var2.SetMissingValue(var2_missing_value);
    std::vector<float> float_data(lon_length * lat_length);
    std::iota(float_data.begin(), float_data.end(), static_cast<float>(120.0));
    cmc::Hyperslab hyperslab2(cmc::DimensionInterval(cmc::Dimension::Lon, 0, lon_length),
                              cmc::DimensionInterval(cmc::Dimension::Lat, 0, lat_length));
    test_var2.PushBack(float_data, hyperslab2);

    /* Create compression settings */
    cmc::CompressionSettings settings;
    const double abs_max_err = 1.0;
    settings.SetAbsoluteErrorCriterion(abs_max_err, cmc::kErrorCriterionHoldsForAllVariables);
    /* A variable-specific criterion is 'stronger' than the genreal criterion */
    const double rel_max_err = 0.05;
    settings.SetRelativeErrorCriterion(rel_max_err, 1);

    std::vector<cmc::InputVar> variables_to_compress;
    variables_to_compress.push_back(std::move(test_var));
    variables_to_compress.push_back(std::move(test_var2));

    /* Create the compression data */
    cmc::CompressionData compression_data(std::move(variables_to_compress), std::move(settings));

    const size_t num_input_vars = compression_data.GetNumberOfInputVariables();

    /* Setup the example data for the compression */
    compression_data.Setup();

    const size_t num_compression_vars = compression_data.GetNumberOfCompressionVariables();

    cmc::ExpectTrue(num_input_vars == num_compression_vars);

    cmc::ExpectTrue(compression_data.IsValidForCompression());

    }

    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
