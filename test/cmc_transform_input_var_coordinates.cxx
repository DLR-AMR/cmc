#include "cmc.h"
#include "test/cmc_test.hxx"
#include "utilities/cmc_input_variable.hxx"

#include <numeric>
#include <vector>

int
main(void)
{
    /* Initialize cmc */
    cmc_initialize();

    {
    
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

    /* Transform Hyperslab coordinates to Morton indices */
    test_var.TransformCoordinatesToMortonIndices();

    /* Move the Morton indices fromt the variable */
    std::vector<cmc::LinearIndex> morton_indices = test_var.DetachMortonIndices();
    test_var.ClearMortonIndices();
    
    /* Check the Morton indices manually */
    std::vector<cmc::LinearIndex> correct_morton_indices{
        0,2,8,10,32,
        1,3,9,11,33,
        4,6,12,14,36,
        5,7,13,15,37,
        16,18,24,26,48,
        17,19,25,27,49,
        20,22,28,30,52,
        21,23,29,31,53,
        64,66,72,74,96,
        65,67,73,75,97
    };

    cmc::ExpectTrue(std::equal(morton_indices.begin(), morton_indices.end(), correct_morton_indices.begin()));

    }
    
    /* Finalize cmc */
    cmc_finalize();


    return cmc::CMC_TEST_SUCCESS;
}
