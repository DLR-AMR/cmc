#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "utilities/cmc_input_variable.hxx"

#include <numeric>

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    
    cmc::InputVariable<int32_t> default_var;

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

    cmc::ExpectTrue(test_var.GetActiveDataFormat() == cmc::DataFormat::HyperslabFormat);
    cmc::ExpectTrue(test_var.GetNumberCoordinates() == static_cast<size_t>(lon_length * lat_length));
    cmc::ExpectTrue(test_var.IsValid());

    //TODO: Test invalid variable

    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
