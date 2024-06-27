#include "cmc.h"
#include "test/cmc_test.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_output_variable.hxx"
#include "lossy/cmc_amr_lossy_compression_settings.hxx"
#include "lossy/cmc_amr_lossy_compression.hxx"

#include <numeric>

/** A system test for the lossy amr compression */

int
main(void)
{
    /* Initialize cmc */
    cmc_initialize();
    
    {
    
    /** CREATE A GLOBAL DOMAIN **/
    /* Create a domain for the examplary data */
    const cmc::DomainIndex lon_length = 10;
    const cmc::DomainIndex lat_length = 5;
    cmc::GeoDomain global_domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, lon_length),
                                 cmc::DimensionInterval(cmc::Dimension::Lat, 0, lat_length));
    

    /** CREATE A DOUBLE VARIABLE **/
    /* Create a examplary variable with double data */
    const int double_variable_id = 0;
    cmc::InputVariable<double> test_var("example_double_data", double_variable_id, cmc::DataLayout::Lon_Lat);
    test_var.SetGlobalDomain(global_domain);
    test_var.SetMissingValue(-2.0);
    std::vector<double> double_data(lon_length * lat_length);
    std::iota(double_data.begin(), double_data.end(), 100.0);
    cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, lon_length),
                             cmc::DimensionInterval(cmc::Dimension::Lat, 0, lat_length));
    
    /* Set the data and it's underlying coordinates */
    test_var.PushBack(double_data, std::move(hyperslab));

    /* Set thr inaccuracy storage and the interpolation funciton to use (those options are the default if not specified) */
    test_var.SetInaccuracyTrackingOption(cmc::TrackingOption::TrackFullInaccuracy);
    test_var.SetInterpolation(cmc::InterpoalteToArithmeticMean);


    /** CREATE A FLOAT VARIABLE **/
    /* Create a examplary variable with double data */
    const int float_variable_id = 33;
    cmc::InputVariable<float> another_test_var("example_float_data", float_variable_id, cmc::DataLayout::Lon_Lat);
    another_test_var.SetGlobalDomain(global_domain);
    another_test_var.SetMissingValue(static_cast<float>(-200.0));
    std::vector<float> float_data(lon_length * lat_length);
    std::iota(float_data.begin(), float_data.end(), 15.0);
    /* The ordering of the dimensions does not matter for the hyperslab (it is the same the above one for the double data) */
    cmc::Hyperslab hyperslab2(cmc::DimensionInterval(cmc::Dimension::Lat, 0, lat_length),
                              cmc::DimensionInterval(cmc::Dimension::Lon, 0, lon_length));
    
    /* Set the data and it's underlying coordinates */
    another_test_var.PushBack(float_data, hyperslab2);


    /** SAVE THE DATA FOR COMPARISON */
    /* Copy the initial data */
    std::vector<double> copy_of_initial_double_data = double_data;
    std::vector<float> copy_of_initial_float_data = float_data;


    /** DEFINE AN ABSOLUTE ERROR CRITERION FOR THE DOUBLE DATA AND A RELATIVE ONE FOR THE FLOAT DATA **/
    /* Create compression settings */
    cmc::CompressionSettings settings;

    const double abs_max_err = 5.0;
    settings.SetAbsoluteErrorCriterion(abs_max_err, double_variable_id);

    const double rel_max_err = 0.1; //10 percent permitted error
    settings.SetRelativeErrorCriterion(rel_max_err, float_variable_id);

    /** GATHER THE VARIABLES OUGHT TO BE COMPRESSED **/
    std::vector<cmc::InputVar> variables_to_compress;
    variables_to_compress.push_back(std::move(test_var));
    variables_to_compress.push_back(another_test_var);


    /* Create the compression data */
    cmc::CompressionData compression_data(std::move(variables_to_compress), std::move(settings));

    /* Setup the example data for the compression */
    compression_data.Setup();

    cmc::ExpectTrue(compression_data.IsValidForCompression());

    /* Compress the supplied variables */
    compression_data.Compress(cmc::CompressionMode::OneForOne);

    /** DECOMPRESS THE DOUBLE VARIABLE **/
    /* Decompress the variable and receive the external wrapper to it */
    cmc::OutputVar decomrpessed_double_var = compression_data.DecompressVariable(double_variable_id);
    /* Obtain the actual decompressed variable based on it's data type */
    cmc::OutputVariable<double> decompressed_double_variable = decomrpessed_double_var.SeizeOutputVariable<double>();
    /* Get it's data */
    std::vector<double> decompressed_double_data = decompressed_double_variable.GetData();

    /** CHECK IF THE ASBSOLUTE ERROR CRITERION IS FULLFILLED FOR THE DOUBLE VARIABLE **/
    /* Check if the error criterion has been satisfied */
    int index = 0;
    for (auto data_iter = decompressed_double_data.begin(); data_iter != decompressed_double_data.end(); ++data_iter, ++index)
    {
        cmc::ExpectTrue(std::abs(*data_iter - copy_of_initial_double_data[index]) <= abs_max_err);
    }

    /** DECOMPRESS THE FLOAT VARIABLE **/
    /* Decompress the variable and receive the external wrapper to it */
    cmc::OutputVar decomrpessed_float_var = compression_data.DecompressVariable(float_variable_id);
    /* Obtain the actual decompressed variable based on it's data type */
    cmc::OutputVariable<float> decompressed_float_variable = decomrpessed_float_var.SeizeOutputVariable<float>();
    /* Get it's data */
    std::vector<float> decompressed_float_data = decompressed_float_variable.GetData();

    /** CHECK IF THE RELATIVE ERROR CRITERIO IS FULLFILLED FOR THE FLOAT VARIABLE **/
    index = 0;
    for (auto data_iter = decompressed_float_data.begin(); data_iter != decompressed_float_data.end(); ++data_iter, ++index)
    {
        cmc::ExpectTrue(std::abs(copy_of_initial_float_data[index] - *data_iter) / std::abs(copy_of_initial_float_data[index]) <= rel_max_err);
    }

    
    }

    /* Finalize cmc */
    cmc_finalize();

    return cmc::CMC_TEST_SUCCESS;
}
