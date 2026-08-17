#include "cmc.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"

#include <unistd.h>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <filesystem>
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

    /* Remove the end of file byte */
    values.pop_back();
    
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

inline void
DisplayHelpMessage()
{
    cmc::cmc_msg("This is cmc (version: ", CMC_VERSION, ").");
    cmc::cmc_msg("This application compares the two supplied binary files.");
    cmc::cmc_msg("It may be used for to compare the intial data to the decompressed data (i.e. <file1> should the initial data and <file2> the decompressed data).");
    cmc::cmc_msg("Usage:\t./cmc_comapre_bin_data_output <file1> <file2>");
    cmc::cmc_msg("or \t./cmc_comapre_bin_data_output -t <data_type> <file1> <file2>");
    cmc::cmc_msg("\t-t \tThe data type of the values stored in the binary files (optional argument); (it is assumed the endiannes of both files matches!)");
    cmc::cmc_msg("\t\tpossible options are float, double, int16_t, uint16_t, int32_t, uint32_t, int64_t or uint64_t");
    cmc::cmc_msg("\t-h \tUse this option to display the help message");
}

cmc::CmcType
EvaluateDataType(const std::string& type)
{
    if (type.compare("float") == 0)
    {
        return cmc::CmcType::Float;
    } else if (type.compare("double") == 0)
    {
        return cmc::CmcType::Double;
    } else if (type.compare("int16_t") == 0)
    {
        return cmc::CmcType::Int16_t;
    } else if (type.compare("uint16_t") == 0)
    {
        return cmc::CmcType::Uint16_t;
    } else if (type.compare("int32_t") == 0)
    {
        return cmc::CmcType::Int32_t;
    } else if (type.compare("uint32_t") == 0)
    {
        return cmc::CmcType::Uint32_t;
    } else if (type.compare("int64_t") == 0)
    {
        return cmc::CmcType::Int64_t;
    } else if (type.compare("uint64_t") == 0)
    {
        return cmc::CmcType::Uint64_t;
    } else
    {
        DisplayHelpMessage();
        cmc::cmc_err_msg("The specified data type (", type, ") is not recognized.");
        return cmc::CmcType::TypeUndefined;
    }
}


template <typename T>
auto ComputeDecompressionMetrics(const std::vector<T>& data1, const std::vector<T>& data2)
 -> std::enable_if_t<std::is_arithmetic_v<T>, void>
{
    if (data1.size() != data2.size())
    {
        cmc::cmc_err_msg("The data vectors are not of the same size!");
    }

    const size_t num_elems = data1.size();

    size_t count_of_equal_values{0};

    for (size_t elem_id{0}; elem_id < num_elems; ++elem_id)
    {
        if (data1[elem_id] == data2[elem_id])
        {
            ++count_of_equal_values;
        }
    }

    const bool are_values_equal = (count_of_equal_values == num_elems);

    if (are_values_equal)
    {
        cmc::cmc_global_msg("Both data sets are bitwise equal!");
        return;
    } else
    {
        /* In case both data sets are not equal */
        cmc::cmc_global_msg("Both data sets are NOT bitwise equal.");
        cmc::cmc_global_msg(count_of_equal_values, " values are equal from overall ", num_elems, " values in the data set.");
    }

    /* Compute some metrics */
    double max_absolute_err{std::numeric_limits<double>::lowest()};
    double min_absolute_err{std::numeric_limits<double>::max()};
    T minimum_value{std::numeric_limits<T>::max()};
    T maximum_value{std::numeric_limits<T>::lowest()};
    
    double max_relative_error{std::numeric_limits<double>::lowest()};
    double min_relative_error{std::numeric_limits<double>::max()};

    double mse{0.0};

    for (size_t idx{0}; idx < num_elems; ++idx)
    {
        if (FP_ZERO != std::fpclassify(data1[idx]))
        {
            const double rel_err = (std::abs(static_cast<double>(data1[idx]) - static_cast<double>(data2[idx]))) / std::abs(static_cast<double>(data1[idx]));
            if (max_relative_error < rel_err)
            {
                max_relative_error = rel_err;
            }
            if (min_relative_error > rel_err)
            {
                min_relative_error = rel_err;
            }
        }

        const double abs_err = std::abs(static_cast<double>(data1[idx]) - static_cast<double>(data2[idx]));
        if (max_absolute_err < abs_err)
        {
            max_absolute_err = abs_err;
        }
        if (min_absolute_err > abs_err)
        {
            min_absolute_err = abs_err;
        }

        if (static_cast<double>(data1[idx]) < minimum_value)
        {
            minimum_value = static_cast<double>(data1[idx]);
        }
        if (static_cast<double>(data1[idx]) > maximum_value)
        {
            maximum_value = static_cast<double>(data1[idx]);
        }

        mse += std::abs(static_cast<double>(data1[idx]) - static_cast<double>(data2[idx])) * std::abs(static_cast<double>(data1[idx]) - static_cast<double>(data2[idx]));
    }


    /* Compute the mean squared error */
    mse = mse / static_cast<double>(num_elems);

    /* Compute the roor MSE*/
    const double rmse = std::sqrt(mse);

    /* Compute the peak signal to noise ratio */
    const double psnr = 20 * std::log10(std::abs(static_cast<double>(maximum_value) - static_cast<double>(minimum_value))) - 10 * std::log10(mse);

    cmc::cmc_global_msg("Maximum initial data value: ", maximum_value, "; Minimum initial data value: ", minimum_value);
    cmc::cmc_global_msg("Maximum absolute error: ", max_absolute_err, "; Minimum absolute error: ", min_absolute_err);
    cmc::cmc_global_msg("Maximum relative error: ", max_relative_error, "; Minimum relative error: ", min_relative_error);
    cmc::cmc_global_msg("Mean Squared Error: ", mse);
    cmc::cmc_global_msg("Root Mean Squared Error: ", rmse);
    cmc::cmc_global_msg("Peak Signal to Noise Ratio: ", psnr);
}

void
CheckEqualityForConvertedValues(const cmc::CmcType type, const std::string& file1, const std::string& file2)
{
    switch (type)
    {
        case cmc::CmcType::Float:
        {
            const std::vector<float> input_data_1 = ReadDataFromStream<float>(file1);
            const std::vector<float> input_data_2 = ReadDataFromStream<float>(file2);
            ComputeDecompressionMetrics<float>(input_data_1, input_data_2);
        }
        break;
        case cmc::CmcType::Double:
        {
            const std::vector<double> input_data_1 = ReadDataFromStream<double>(file1);
            const std::vector<double> input_data_2 = ReadDataFromStream<double>(file2);
            ComputeDecompressionMetrics<double>(input_data_1, input_data_2);
        }
        break;
        case cmc::CmcType::Int16_t:
        {
            const std::vector<int16_t> input_data_1 = ReadDataFromStream<int16_t>(file1);
            const std::vector<int16_t> input_data_2 = ReadDataFromStream<int16_t>(file2);
            ComputeDecompressionMetrics<int16_t>(input_data_1, input_data_2);
        }
        break;
        case cmc::CmcType::Uint16_t:
        {
            const std::vector<uint16_t> input_data_1 = ReadDataFromStream<uint16_t>(file1);
            const std::vector<uint16_t> input_data_2 = ReadDataFromStream<uint16_t>(file2);
            ComputeDecompressionMetrics<uint16_t>(input_data_1, input_data_2);
        }
        break;
        case cmc::CmcType::Int32_t:
        {
            const std::vector<int32_t> input_data_1 = ReadDataFromStream<int32_t>(file1);
            const std::vector<int32_t> input_data_2 = ReadDataFromStream<int32_t>(file2);
            ComputeDecompressionMetrics<int32_t>(input_data_1, input_data_2);
        }
        break;
        case cmc::CmcType::Uint32_t:
        {
            const std::vector<uint32_t> input_data_1 = ReadDataFromStream<uint32_t>(file1);
            const std::vector<uint32_t> input_data_2 = ReadDataFromStream<uint32_t>(file2);
            ComputeDecompressionMetrics<uint32_t>(input_data_1, input_data_2);
        }
        break;
        case cmc::CmcType::Int64_t:
        {
            const std::vector<int64_t> input_data_1 = ReadDataFromStream<int64_t>(file1);
            const std::vector<int64_t> input_data_2 = ReadDataFromStream<int64_t>(file2);
            ComputeDecompressionMetrics<int64_t>(input_data_1, input_data_2);
        }
        break;
        case cmc::CmcType::Uint64_t:
        {
            const std::vector<uint64_t> input_data_1 = ReadDataFromStream<uint64_t>(file1);
            const std::vector<uint64_t> input_data_2 = ReadDataFromStream<uint64_t>(file2);
            ComputeDecompressionMetrics<uint64_t>(input_data_1, input_data_2);
        }
        break;
        default:
            DisplayHelpMessage();
            cmc::cmc_err_msg("A not supported/recognized data type has been specified.");
        break; 
    }
}

constexpr int kNumExpectedArgs = 3;
constexpr int kNumExpectedArgsWithDataTypeOpt = 5;

int
main(int argc, char *argv[])
{
    cmc::CmcInitialize(cmc::kMinimumInitialization);

    int opt;
    bool is_type_given = false;
    cmc::CmcType data_type;

    while ((opt = getopt(argc, argv, "t:h")) != -1)
    {
        switch (opt)
        {
            case 't':
                is_type_given = true;
                data_type = EvaluateDataType(std::string(optarg));
                break;
            case 'h':
                DisplayHelpMessage();
                return 0;
                break;
            default:
            break;
        }
    }

    if (const int num_expected_args = (is_type_given ? kNumExpectedArgsWithDataTypeOpt : kNumExpectedArgs);
        argc != num_expected_args)
    {
        cmc::cmc_msg("The number of supplied arguments does not match the supported calling scheme.");
        DisplayHelpMessage();
        return 1;
    }

    /* Get the file names */
    std::string input_file_1;
    std::string input_file_2;

    if (is_type_given)
    {
        input_file_1 = std::string(argv[3]);
        input_file_2 = std::string(argv[4]);
    } else
    {
        input_file_1 = std::string(argv[1]);
        input_file_2 = std::string(argv[2]);
    }

    cmc::cmc_msg("Input File 1: ", input_file_1);
    cmc::cmc_msg("Input File 2: ", input_file_2);

    /* Check if the input file 1 exists */
    if (const std::filesystem::path input_file_path(input_file_1); not std::filesystem::exists(input_file_path))
    {
        DisplayHelpMessage();
        cmc::cmc_err_msg("The file ", input_file_1, " does not exist.");
    }

    /* Check if the input file 2 exists */
    if (const std::filesystem::path input_file_path(input_file_2); not std::filesystem::exists(input_file_path))
    {
        DisplayHelpMessage();
        cmc::cmc_err_msg("The file ", input_file_2, " does not exist.");
    }

    int cmp;
    {
        const std::vector<char> input_data_1_bytes = ReadDataFromStream<char>(input_file_1);
        const std::vector<char> input_data_2_bytes = ReadDataFromStream<char>(input_file_2);
        cmp = std::memcmp(input_data_1_bytes.data(), input_data_2_bytes.data(), input_data_2_bytes.size() * sizeof(char));
    }
    cmc::cmc_msg("Comparison of both arrays with memcmp yields return value: ", cmp);

    if (is_type_given)
    {
        CheckEqualityForConvertedValues(data_type, input_file_1, input_file_2);
    } else
    {
        if (cmp == 0)
        {
            cmc::cmc_msg("The files ", input_file_1, " and ", input_file_2, " are equal.");
        } else
        {
            cmc::cmc_msg("The files ", input_file_1, " and ", input_file_2, " are not equal.");
        }
    }

    cmc::CmcFinalize();

   return 0;
}
