#include "utilities/cmc_utilities.hxx"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <limits>
#include <cmath>

template <typename T>
std::vector<T>
ReadDataFromStream(const std::string& file_name)
{
    std::vector<T> values;

    std::ifstream istrm(file_name, std::ios::binary);
    if (!istrm.is_open())
        std::cout << "failed to open " << file_name << '\n';
    else
    {
        while ( !istrm.eof() ) {
            T value;
            istrm.read(reinterpret_cast<char*>(&value), sizeof(T));

            values.push_back(value);
        }
        istrm.close();
    }

    return values;
}

template <typename T>
void
PrintValuesFromFile(const std::string& file_name)
{
    std::ifstream istrm(file_name, std::ios::binary);
    if (!istrm.is_open())
        std::cout << "failed to open " << file_name << '\n';
    else
    {
        while ( !istrm.eof() ) {
            T value;
            istrm.read(reinterpret_cast<char*>(&value), sizeof(T));

            std::cout << "Value: " << value << std::endl;
        }
        istrm.close();
    }
}

int main(void)
{
    #if 1
    const std::vector<float> initial_data = ReadDataFromStream<float>("../../programs/data/100x500x500/TCf48.bin.f32");
    #else
    const std::vector<float> initial_data = ReadDataFromStream<float>("../../programs/data/SDRBENCH-EXASKY-NYX-512x512x512/baryon_density.f32");
    #endif

    const std::vector<float> decompressed_data = ReadDataFromStream<float>("decompressed_data.cmc");

    const double permitted_rel_error = 0.001;
    const bool print_rel_error_deviations = false;
    const double permitted_abs_error = 0.05;
    const bool print_abs_error_deviations = false;

    std::cout << "Initial data size: " << initial_data.size() << std::endl;
    std::cout << "Decompressed data size: " << decompressed_data.size() << std::endl;

    double max_rel_error{std::numeric_limits<double>::lowest()};
    double max_abs_error{std::numeric_limits<double>::lowest()};
    double mse{0.0};
    size_t mse_count{0};
    double abs_max_val{std::numeric_limits<double>::lowest()};

    float max_rel_err_val_decompressed;
    float max_rel_err_val_initial;
    int64_t max_rel_err_idx = -1;
    size_t count_rel_err_not_compliant = 0;

    for (size_t idx = 0; idx < initial_data.size(); ++idx)
    {
        double rel_err{0.0};
        double abs_err{0.0};

        if (FP_ZERO != std::fpclassify(initial_data[idx]))
        {
            rel_err = (std::abs(initial_data[idx] - decompressed_data[idx])) / std::abs(initial_data[idx]);
        }
        abs_err = std::abs(initial_data[idx] - decompressed_data[idx]);

        if (print_rel_error_deviations && rel_err > permitted_rel_error)
        {
            std::cout << "Rel Error-Bound violation at " << idx << ", permitted: " << permitted_rel_error << ", actual: " << rel_err << ", init value: " << initial_data[idx] << ", decompressed: " << decompressed_data[idx] << std::endl;
        }
        if (print_abs_error_deviations && abs_err > permitted_abs_error)
        {
            std::cout << "Abs Error-Bound violation at " << idx << ", permitted: " << permitted_abs_error << ", actual: " << abs_err << ", init value: " << initial_data[idx] << ", decompressed: " << decompressed_data[idx] << std::endl;
        }
        
        if (rel_err > max_rel_error)
        {
            max_rel_error = rel_err; 
            max_rel_err_val_decompressed = decompressed_data[idx];
            max_rel_err_val_initial = initial_data[idx];
            max_rel_err_idx = idx;
            ++count_rel_err_not_compliant;
        }
        if (abs_err > max_abs_error)
        {
            max_abs_error = abs_err; 
        }

        if (abs_max_val < std::abs(initial_data[idx]))
        {
            abs_max_val = std::abs(initial_data[idx]);
        }

        ++mse_count;
        mse += std::abs(initial_data[idx] - decompressed_data[idx]) * std::abs(initial_data[idx] - decompressed_data[idx]);
    }

    std::cout << "Maximum absolute error: " << max_abs_error << ", Maximum relative error: " << max_rel_error << std::endl << "RMSE: " << std::sqrt(mse / static_cast<double>(mse_count)) << std::endl << "PSNR: " << (10 * std::log10(abs_max_val / mse)) << std::endl;
    std::cout << "Max Rel Error at Index " << max_rel_err_idx << ": Initial: " << max_rel_err_val_initial << ", Decompressed: " << max_rel_err_val_decompressed << std::endl;
    std::cout << "Num of non compliant rel errors: " << count_rel_err_not_compliant << std::endl;
    return 0;
}
