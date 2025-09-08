#include "cmc.hxx"
#include "input/cmc_binary_reader.hxx"
#include "utilities/cmc_log_functions.hxx"

#include <vector>
#include <variant>
#include <cmath>

constexpr int kLevLength = 100;
constexpr int kLatLength = 500;
constexpr int kLonLength = 500;

float
GetValue(const std::vector<float>& data, const int lev, const int lat, const int lon)
{
    cmc_assert(lev * (kLatLength * kLonLength) + lat * kLonLength + lon < data.size());
    return data[lev * (kLatLength * kLonLength) + lat * kLonLength + lon];
}

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    {

    /* Read in data */
    #if 1
    const std::string file = "../../programs/data/100x500x500/Uf48.bin.f32";
    #else
    const std::string file = "../../programs/data/SDRBENCH-EXASKY-NYX-512x512x512/baryon_density.f32";
    #endif

    cmc::CmcType type(cmc::CmcType::Float);
    const std::string name("compr_test_var");
    const int id = 1;

    //Missing Values Hurricane ISABEL Dataset
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.0025));//CLOUD
    //cmc::CmcUniversalType missing_value(static_cast<float>(-2.5));//CLOUD.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(3224.4)); //P
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.00755)); //PRECIP
    //cmc::CmcUniversalType missing_value(static_cast<float>(-2.0)); //PRECIPf48.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.00205));//QCLOUD
    //cmc::CmcUniversalType missing_value(static_cast<float>(-2.5));//QCLOUD.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.007295)); //QGraup
    //cmc::CmcUniversalType missing_value(static_cast<float>(-2.0)); //QGraup.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.00085));//QICE
    //cmc::CmcUniversalType missing_value(static_cast<float>(-3.0));//QICE.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.0065));//QRAINf48
    //cmc::CmcUniversalType missing_value(static_cast<float>(-2.0));//QRAINf48.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.000875)); //QSNOW
    //cmc::CmcUniversalType missing_value(static_cast<float>(-3.0)); //QSNOW.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(30.0)); //QVAPOR
    //cmc::CmcUniversalType missing_value(static_cast<float>(29.65)); //TC
    cmc::CmcUniversalType missing_value(static_cast<float>(40.0)); //UF
    //cmc::CmcUniversalType missing_value(static_cast<float>(48.25)); //VF
    //cmc::CmcUniversalType missing_value(static_cast<float>(13.4)); //WF

    //Missing values NYX data
    //cmc::CmcUniversalType missing_value(static_cast<float>(31868.0)); //velocity x
    //cmc::CmcUniversalType missing_value(static_cast<float>(56507.0)); //velocity y
    //cmc::CmcUniversalType missing_value(static_cast<float>(33387.0)); //velocity z
    //cmc::CmcUniversalType missing_value(static_cast<float>(4784.0)); //temp
    //cmc::CmcUniversalType missing_value(static_cast<float>(13780.0)); //dark matter
    //cmc::CmcUniversalType missing_value(static_cast<float>(115864.0)); //baryon denisty

    #if 1
    
    const size_t num_elements = kLonLength * kLatLength * kLevLength;
    cmc::DataLayout layout(cmc::DataLayout::Lev_Lat_Lon);
    cmc::GeoDomain domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, kLonLength),
                          cmc::DimensionInterval(cmc::Dimension::Lat, 0, kLatLength),
                          cmc::DimensionInterval(cmc::Dimension::Lev, 0, kLevLength)
                          );
    #else
    const size_t num_elements = 512 * 512 * 512;
    cmc::DataLayout layout(cmc::DataLayout::Lev_Lat_Lon);
    cmc::GeoDomain domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 512),
                          cmc::DimensionInterval(cmc::Dimension::Lat, 0, 512),
                          cmc::DimensionInterval(cmc::Dimension::Lev, 0, 512)
                          );
    #endif

    /* Generate input variables from the binary file */
    cmc::input::binary::Reader binary_reader(file);
    cmc::input::Var variable = binary_reader.CreateVariableFromBinaryData(type, name, id, num_elements, missing_value, layout, domain);

    const cmc::input::GeneralVariable& variable_variant = variable.GetInternalVariant();
    const cmc::input::Variable<float>& input_variable = std::get<cmc::input::Variable<float>>(variable_variant);
    const std::vector<float>& const_data = input_variable.GetDataForReading();
    std::vector<float> data = const_data;

    cmc::cmc_debug_msg("Data hat size: ", data.size());

    const double rel_err = 0.000001;
    int count = 0;
    double mse = 0.0;
    int count_zero_residuals = 0;

    /* Perform 3D Lorenzo prediciton and check the residual sizes */
    /* We just predict all values which can be 3D predicted and check how they do */
    for (int lev_iter = 1; lev_iter < kLevLength; ++lev_iter)
    {
        const int lev_offset = lev_iter * (kLatLength * kLonLength);

        for (int lat_iter = 1; lat_iter < kLatLength; ++lat_iter)
        {
            const int lat_offset = lat_iter * kLonLength;

            for (int lon_iter = 1; lon_iter < kLonLength; ++lon_iter)
            {
                const int value_idx = lev_offset + lat_offset + lon_iter;

                /* Compute array accesses */
                const float val3_0 = GetValue(data, lev_iter -1, lat_iter -1, lon_iter -1);

                const float val2_0 = GetValue(data, lev_iter - 1, lat_iter - 1, lon_iter);
                const float val2_1 = GetValue(data, lev_iter - 1, lat_iter, lon_iter - 1);
                const float val2_2 = GetValue(data, lev_iter, lat_iter - 1, lon_iter -1);

                const float val1_0 = GetValue(data, lev_iter -1, lat_iter, lon_iter);
                const float val1_1 = GetValue(data, lev_iter, lat_iter -1, lon_iter);
                const float val1_2 = GetValue(data, lev_iter, lat_iter, lon_iter -1);

                
                /* Perform prediction */
                const double prediction = static_cast<double>(val3_0 + val1_0 + val1_1 + val1_2 - val2_0 - val2_1 - val2_2);

                /* Get initial value */
                const double init_value = static_cast<double>(data[value_idx]);

                /* Check residual */
                if (std::abs((init_value - prediction) / init_value) <= rel_err)
                {
                    ++count_zero_residuals;
                    data[value_idx] = val3_0 + val1_0 + val1_1 + val1_2 - val2_0 - val2_1 - val2_2;
                }
                
                /* Substiute intiial value for approximation */
                //data[value_idx] = val3_0 + val1_0 + val1_1 + val1_2 - val2_0 - val2_1 - val2_2;

                ++count;
                mse += std::abs(init_value - prediction) * std::abs(init_value - prediction);
            }
        }
    }

    cmc::cmc_debug_msg("RMSE: ", std::sqrt(mse / static_cast<double>(count)),", MSE: ", mse, ", Count: ", count);
    cmc::cmc_debug_msg("Num zero residuals: ", count_zero_residuals);

    }
    /* Finalize cmc */
    cmc::CmcFinalize();


    return 0;
}