#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "utilities/cmc_input_variable.hxx"

#include <numeric>

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();

    cmc::InputVariable<int32_t> test_var("data_to_be_transformed", 0, cmc::DataLayout::Lon_Lat);

    const cmc::DomainIndex lon_length = 10;
    const cmc::DomainIndex lat_length = 5;
    cmc::GeoDomain global_domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, lon_length),
                                 cmc::DimensionInterval(cmc::Dimension::Lat, 0, lat_length));
    
    test_var.SetGlobalDomain(global_domain);
    test_var.SetMissingValue(-1);

    std::vector<int32_t> int_data(lon_length * lat_length);
    std::iota(int_data.begin(), int_data.end(), 100);

    cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, lon_length),
                             cmc::DimensionInterval(cmc::Dimension::Lat, 0, lat_length));

    test_var.PushBack(int_data, std::move(hyperslab));

    const cmc::CmcUniversalType scale_factor = static_cast<double>(2.0);
    const cmc::CmcUniversalType add_offset = static_cast<float>(10.0);

    test_var.SetScaleFactor(scale_factor);
    test_var.SetAddOffset(add_offset);

    cmc::InputVar initial_var(std::move(test_var));

    initial_var.ApplyScalingAndOffset();

    cmc::ExpectTrue(initial_var.GetType() == cmc::GetDefaultCmcType());

    /* Use copy constructor */
    cmc::CmcInputVariable copied_initial_variant = initial_var.GetInternalVariant();

    /* Get the underlying InputVariable */
    cmc::InputVariable<cmc::CmcDefaultDataType>& copied_initial_var = std::get<cmc::InputVariable<cmc::CmcDefaultDataType>>(copied_initial_variant);

    /* Move the data from the copied variable */
    std::vector<cmc::CmcDefaultDataType> copied_initial_var_data = copied_initial_var.DetachData();
    copied_initial_var.ClearData();

    /* Check manually whether the right scaling and offset has been applied */
    std::vector<cmc::CmcDefaultDataType> correctly_transformed_data{
        210.0,212.0,214.0,216.0,218.0,
        220.0,222.0,224.0,226.0,228.0,
        230.0,232.0,234.0,236.0,238.0,
        240.0,242.0,244.0,246.0,248.0,
        250.0,252.0,254.0,256.0,258.0,
        260.0,262.0,264.0,266.0,268.0,
        270.0,272.0,274.0,276.0,278.0,
        280.0,282.0,284.0,286.0,288.0,
        290.0,292.0,294.0,296.0,298.0,
        300.0,302.0,304.0,306.0,308.0
    };

    cmc::ExpectTrue(std::equal(copied_initial_var_data.begin(), copied_initial_var_data.end(), correctly_transformed_data.begin()));

    //TODO: Check only scaling and only offset
    
    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
