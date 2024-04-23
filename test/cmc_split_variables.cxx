#include "cmc.h"
#include "test/cmc_test.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "lossy/cmc_amr_lossy_compression_settings.hxx"
#include "t8code/cmc_t8_data.hxx"

#include <numeric>

int
main(void)
{
    /* Initialize cmc */
    cmc_initialize();
    
    /* Create some examplary variable with data */
    const int var_id = 0;

    cmc::InputVariable<double> test_var("ex_double_data", var_id, cmc::DataLayout::Lon_Lat_Lev);

    const cmc::DomainIndex lon_length = 10;
    const cmc::DomainIndex lat_length = 5;
    const cmc::DomainIndex lev_length = 5;

    cmc::GeoDomain global_domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, lon_length),
                                 cmc::DimensionInterval(cmc::Dimension::Lat, 0, lat_length),
                                 cmc::DimensionInterval(cmc::Dimension::Lev, 0, lev_length));

    test_var.SetGlobalDomain(global_domain);

    test_var.SetMissingValue(-2.0);

    std::vector<double> double_data(lon_length * lat_length * lev_length);
    std::iota(double_data.begin(), double_data.end(), 100.0);

    cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, lon_length),
                             cmc::DimensionInterval(cmc::Dimension::Lat, 0, lat_length),
                             cmc::DimensionInterval(cmc::Dimension::Lev, 0, lev_length));

    test_var.PushBack(double_data, std::move(hyperslab));

    /* Create compression settings */
    cmc::CompressionSettings settings;

    const double abs_max_err = 10.0;
    settings.SetAbsoluteErrorCriterion(abs_max_err, cmc::kErrorCriterionHoldsForAllVariables);
    
    const cmc::Dimension dimension_to_split = cmc::Dimension::Lev;
    cmc::SplitVariable split(var_id, dimension_to_split);
    settings.SplitVariableByDimension(split);

    std::vector<cmc::InputVar> variables_to_compress;
    variables_to_compress.push_back(std::move(test_var));

    /* Create the compression data */
    cmc::AmrData amr_compression_data(std::move(variables_to_compress), settings);

    amr_compression_data.SplitVariables();

    for (auto iter = amr_compression_data.GetInputVariablesBegin(); iter != amr_compression_data.GetInputVariablesEnd(); ++iter)
    {
        cmc::ExpectEQ(iter->GetInitialDataLayout() == cmc::DataLayout::Lon_Lat);
        
        const cmc::GeoDomain& var_global_domain = iter->GetGlobalDomain();
        cmc::ExpectEQ(var_global_domain.GetDimensionality() == 2);
        cmc::ExpectEQ(var_global_domain.GetDimensionStartIndex(cmc::Dimension::Lon) == 0 &&
                      var_global_domain.GetDimensionEndIndex(cmc::Dimension::Lon) == lon_length);
        cmc::ExpectEQ(var_global_domain.GetDimensionStartIndex(cmc::Dimension::Lat) == 0 &&
                      var_global_domain.GetDimensionEndIndex(cmc::Dimension::Lat) == lat_length);
        cmc::ExpectEQ(var_global_domain.GetDimensionStartIndex(cmc::Dimension::Lev) == 0 &&
                      var_global_domain.GetDimensionEndIndex(cmc::Dimension::Lev) == 0);

    }

    /* Finalize cmc */
    cmc_finalize();

    return cmc::CMC_TEST_SUCCESS;
}
