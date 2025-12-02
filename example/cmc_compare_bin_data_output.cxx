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

inline void
DisplayHelpMessage()
{
    cmc::cmc_msg("This is cmc (version: ", CMC_VERSION, ").");
    cmc::cmc_msg("This application compares the two supplied binary files for equality.");
    cmc::cmc_msg("It may be used for example to compare the intial data to the losslessly decompressed data.");
    cmc::cmc_msg("Usage:\t./cmc_comapre_bin_data_output <file1> <file2>");
    cmc::cmc_msg("or \t./cmc_comapre_bin_data_output -t <data_type> <file1> <file2>");
    cmc::cmc_msg("\t-t \tThe data type of the values stored in the binary files (optional argument); ");
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

bool
CheckEqualityForConvertedValues(const cmc::CmcType type, const std::string& file1, const std::string& file2)
{
    switch (type)
    {
        case cmc::CmcType::Float:
        {
            const std::vector<float> input_data_1 = ReadDataFromStream<float>(file1);
            const std::vector<float> input_data_2 = ReadDataFromStream<float>(file2);
            bool is_input_equal = true;
            if (input_data_1.size() != input_data_2.size()) {cmc::cmc_err_msg("The array sizes do not match (input file 1: ", input_data_1.size(), ", input file 2:", input_data_2.size(), ")");}
            for (size_t idx = 0; idx < input_data_1.size(); ++idx)
            {
                if (input_data_1[idx] != input_data_2[idx])
                {
                    cmc::cmc_debug_msg("No equality at ", idx, ", InputFile1 Value: ", input_data_1[idx], ", InputFile2 Value: ", input_data_2[idx]);
                    is_input_equal = false;
                }
            }
            return is_input_equal;
        }
        break;
        case cmc::CmcType::Double:
        {
            const std::vector<double> input_data_1 = ReadDataFromStream<double>(file1);
            const std::vector<double> input_data_2 = ReadDataFromStream<double>(file2);
            bool is_input_equal = true;
            if (input_data_1.size() != input_data_2.size()) {cmc::cmc_err_msg("The array sizes do not match (input file 1: ", input_data_1.size(), ", input file 2:", input_data_2.size(), ")");}
            for (size_t idx = 0; idx < input_data_1.size(); ++idx)
            {
                if (input_data_1[idx] != input_data_2[idx])
                {
                    cmc::cmc_debug_msg("No equality at ", idx, ", InputFile1 Value: ", input_data_1[idx], ", InputFile2 Value: ", input_data_2[idx]);
                    is_input_equal = false;
                }
            }
            return is_input_equal;
        }
        break;
        case cmc::CmcType::Int16_t:
        {
            const std::vector<int16_t> input_data_1 = ReadDataFromStream<int16_t>(file1);
            const std::vector<int16_t> input_data_2 = ReadDataFromStream<int16_t>(file2);
            bool is_input_equal = true;
            if (input_data_1.size() != input_data_2.size()) {cmc::cmc_err_msg("The array sizes do not match (input file 1: ", input_data_1.size(), ", input file 2:", input_data_2.size(), ")");}
            for (size_t idx = 0; idx < input_data_1.size(); ++idx)
            {
                if (input_data_1[idx] != input_data_2[idx])
                {
                    cmc::cmc_debug_msg("No equality at ", idx, ", InputFile1 Value: ", input_data_1[idx], ", InputFile2 Value: ", input_data_2[idx]);
                    is_input_equal = false;
                }
            }
            return is_input_equal;
        }
        break;
        case cmc::CmcType::Uint16_t:
        {
            const std::vector<uint16_t> input_data_1 = ReadDataFromStream<uint16_t>(file1);
            const std::vector<uint16_t> input_data_2 = ReadDataFromStream<uint16_t>(file2);
            bool is_input_equal = true;
            if (input_data_1.size() != input_data_2.size()) {cmc::cmc_err_msg("The array sizes do not match (input file 1: ", input_data_1.size(), ", input file 2:", input_data_2.size(), ")");}
            for (size_t idx = 0; idx < input_data_1.size(); ++idx)
            {
                if (input_data_1[idx] != input_data_2[idx])
                {
                    cmc::cmc_debug_msg("No equality at ", idx, ", InputFile1 Value: ", input_data_1[idx], ", InputFile2 Value: ", input_data_2[idx]);
                    is_input_equal = false;
                }
            }
            return is_input_equal;
        }
        break;
        case cmc::CmcType::Int32_t:
        {
            const std::vector<int32_t> input_data_1 = ReadDataFromStream<int32_t>(file1);
            const std::vector<int32_t> input_data_2 = ReadDataFromStream<int32_t>(file2);
            bool is_input_equal = true;
            if (input_data_1.size() != input_data_2.size()) {cmc::cmc_err_msg("The array sizes do not match (input file 1: ", input_data_1.size(), ", input file 2:", input_data_2.size(), ")");}
            for (size_t idx = 0; idx < input_data_1.size(); ++idx)
            {
                if (input_data_1[idx] != input_data_2[idx])
                {
                    cmc::cmc_debug_msg("No equality at ", idx, ", InputFile1 Value: ", input_data_1[idx], ", InputFile2 Value: ", input_data_2[idx]);
                    is_input_equal = false;
                }
            }
            return is_input_equal;
        }
        break;
        case cmc::CmcType::Uint32_t:
        {
            const std::vector<uint32_t> input_data_1 = ReadDataFromStream<uint32_t>(file1);
            const std::vector<uint32_t> input_data_2 = ReadDataFromStream<uint32_t>(file2);
            bool is_input_equal = true;
            if (input_data_1.size() != input_data_2.size()) {cmc::cmc_err_msg("The array sizes do not match (input file 1: ", input_data_1.size(), ", input file 2:", input_data_2.size(), ")");}
            for (size_t idx = 0; idx < input_data_1.size(); ++idx)
            {
                if (input_data_1[idx] != input_data_2[idx])
                {
                    cmc::cmc_debug_msg("No equality at ", idx, ", InputFile1 Value: ", input_data_1[idx], ", InputFile2 Value: ", input_data_2[idx]);
                    is_input_equal = false;
                }
            }
            return is_input_equal;
        }
        break;
        case cmc::CmcType::Int64_t:
        {
            const std::vector<int64_t> input_data_1 = ReadDataFromStream<int64_t>(file1);
            const std::vector<int64_t> input_data_2 = ReadDataFromStream<int64_t>(file2);
            bool is_input_equal = true;
            if (input_data_1.size() != input_data_2.size()) {cmc::cmc_err_msg("The array sizes do not match (input file 1: ", input_data_1.size(), ", input file 2:", input_data_2.size(), ")");}
            for (size_t idx = 0; idx < input_data_1.size(); ++idx)
            {
                if (input_data_1[idx] != input_data_2[idx])
                {
                    cmc::cmc_debug_msg("No equality at ", idx, ", InputFile1 Value: ", input_data_1[idx], ", InputFile2 Value: ", input_data_2[idx]);
                    is_input_equal = false;
                }
            }
            return is_input_equal;
        }
        break;
        case cmc::CmcType::Uint64_t:
        {
            const std::vector<uint64_t> input_data_1 = ReadDataFromStream<uint64_t>(file1);
            const std::vector<uint64_t> input_data_2 = ReadDataFromStream<uint64_t>(file2);
            bool is_input_equal = true;
            if (input_data_1.size() != input_data_2.size()) {cmc::cmc_err_msg("The array sizes do not match (input file 1: ", input_data_1.size(), ", input file 2:", input_data_2.size(), ")");}
            for (size_t idx = 0; idx < input_data_1.size(); ++idx)
            {
                if (input_data_1[idx] != input_data_2[idx])
                {
                    cmc::cmc_debug_msg("No equality at ", idx, ", InputFile1 Value: ", input_data_1[idx], ", InputFile2 Value: ", input_data_2[idx]);
                    is_input_equal = false;
                }
            }
            return is_input_equal;
        }
        break;
        default:
            DisplayHelpMessage();
            cmc::cmc_err_msg("A not supported/recognized data type has been specified.");
            return false;
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
        const bool is_input_equal = CheckEqualityForConvertedValues(data_type, input_file_1, input_file_2);
        if (is_input_equal)
        {
            cmc::cmc_msg("The files ", input_file_1, " and ", input_file_2, " are equal.");
        } else
        {
            cmc::cmc_msg("The files ", input_file_1, " and ", input_file_2, " are not equal.");
        }
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
