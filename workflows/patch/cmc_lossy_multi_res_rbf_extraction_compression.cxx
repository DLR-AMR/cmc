#include "cmc.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "patch/lossy/cmc_multi_res_extraction_rbf.hxx"
#include "input/cmc_binary_file_reader.hxx"

#include <unistd.h>
#include <cstdio>
#include <numeric>
#include <algorithm>
#include <vector>
#include <filesystem>
#include <string>

inline void
DisplayHelpMessage()
{
    cmc::cmc_msg("This is cmc (version: ", CMC_VERSION, ").");
    cmc::cmc_msg("The following options are required to run the patch-based lossless multi resolution extraction compression with serial output!");
    cmc::cmc_msg("");
    cmc::cmc_msg("\t-t \tThe data type of the data to compress, possible options are float or double");
    cmc::cmc_msg("\t-i \tThe path to the input file storing the array of binary data to compress");
    cmc::cmc_msg("\t-o \tThe name of the output file storing the compressed data");
    cmc::cmc_msg("\t-a \tThe absolute permitted error, e.g. 0.01");
    cmc::cmc_msg("\t-d \tThe dimensionality of the array storing the binary data, e.g. 3 for 3D data. Currently, 2D and 3D");
    cmc::cmc_msg("\t-s \tA string representing the shape of the data, e.g. '100 200 500' for a 3D array of size 100x200x500,"
                 " such that the slowest varying dimension is the first (e.g. 100) and the fastest varying dimension is the last parameter (e.g. 500)");
    cmc::cmc_msg("\t-h \tTo display this help message");
}

void
Compress(const cmc::CmcType data_type, const std::string& input_file, const std::string& output_file, const int dim, const std::vector<size_t>& dim_lengths, const float abs_error, const std::endian endianness_of_data_in_file = std::endian::native)
{
    cmc::cmc_debug_msg("Performing Patch-Based Lossless MultiRes Extraction Compression.");

    switch (dim)
    {
        /** 2D Compression **/
        case 2:
            switch (data_type)
            {
                case cmc::CmcType::Float:
                {
                    const std::array<int32_t, 2> dims{static_cast<int32_t>(dim_lengths[0]), static_cast<int32_t>(dim_lengths[1])};
                    cmc::input::binary_file::Reader<float, 2> binary_reader(input_file, dims, endianness_of_data_in_file);
                    const std::vector<float> data = binary_reader.ReadData();
                    cmc::patch::lossy::multi_res::rbf::CompressionVariable<float, 2> compression_variable(std::span<const float>(data), dims, abs_error);
                    compression_variable.Compress();
                    compression_variable.WriteCompressedData(output_file);
                }
                break;
                case cmc::CmcType::Double:
                {
                    const std::array<int32_t, 2> dims{static_cast<int32_t>(dim_lengths[0]), static_cast<int32_t>(dim_lengths[1])};
                    cmc::input::binary_file::Reader<double, 2> binary_reader(input_file, dims, endianness_of_data_in_file);
                    const std::vector<double> data = binary_reader.ReadData();
                    cmc::patch::lossy::multi_res::rbf::CompressionVariable<double, 2> compression_variable(std::span<const double>(data), dims, abs_error);
                    compression_variable.Compress();
                    compression_variable.WriteCompressedData(output_file);
                }
                break;
                default:
                    DisplayHelpMessage();
                    cmc::cmc_err_msg("A not supported/recognized data type has been specified.");
                break;
            }
        break;
        /** 3D Compression **/
        case 3:
            switch (data_type)
            {
                case cmc::CmcType::Float:
                {
                    const std::array<int32_t, 3> dims{static_cast<int32_t>(dim_lengths[0]), static_cast<int32_t>(dim_lengths[1]), static_cast<int32_t>(dim_lengths[2])};
                    cmc::input::binary_file::Reader<float, 3> binary_reader(input_file, dims, endianness_of_data_in_file);
                    const std::vector<float> data = binary_reader.ReadData();
                    cmc::patch::lossy::multi_res::rbf::CompressionVariable<float, 3> compression_variable(std::span<const float>(data), dims, abs_error);
                    compression_variable.Compress();
                    compression_variable.WriteCompressedData(output_file);
                }
                break;
                case cmc::CmcType::Double:
                {
                    const std::array<int32_t, 3> dims{static_cast<int32_t>(dim_lengths[0]), static_cast<int32_t>(dim_lengths[1]), static_cast<int32_t>(dim_lengths[2])};
                    cmc::input::binary_file::Reader<double, 3> binary_reader(input_file, dims, endianness_of_data_in_file);
                    const std::vector<double> data = binary_reader.ReadData();
                    cmc::patch::lossy::multi_res::rbf::CompressionVariable<double, 3> compression_variable(std::span<const double>(data), dims, abs_error);
                    compression_variable.Compress();
                    compression_variable.WriteCompressedData(output_file);
                }
                break;
                default:
                    DisplayHelpMessage();
                    cmc::cmc_err_msg("A not supported/recognized data type has been specified.");
                break;
            }
        break;
        default:
            cmc::cmc_err_msg("Only 2D and 3D compression is currently supported.");
    }
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
    } else
    {
        DisplayHelpMessage();
        cmc::cmc_err_msg("The specified data type (", type, ") is not recognized.");
        return cmc::CmcType::TypeUndefined;
    }
}

void
DisplayCompressionResult(cmc::CmcType data_type, const std::vector<size_t>& dim_lengths, const std::string& output_file)
{
    /* Check if the output file exists */
    const std::filesystem::path output_file_path(output_file);
    if (not std::filesystem::exists(output_file_path))
    {
        cmc::cmc_err_msg("The compression output file ", output_file, " does not exist.");
    }

    /* Check the file size */
    const auto file_size = std::filesystem::file_size(output_file_path);

    /* Compute the intial data amount */
    const size_t num_elements = std::reduce(dim_lengths.begin(), dim_lengths.end(), 1, std::multiplies<size_t>());
    const size_t data_amount = num_elements * cmc::CmcTypeToBytes(data_type);

    cmc::cmc_msg("The patch-based lossless multi resolution compression extraction has been applied successfully.");
    cmc::cmc_msg("The compressed data has been written to ", output_file, ".");
    cmc::cmc_msg("The intial data amount accumulated to: \t", data_amount, " bytes.");
    cmc::cmc_msg("The compressed data storage amounts to: \t", file_size, " bytes.");
}

constexpr int kNumArgCRequired = 7;

int
main(int argc, char *argv[])
{
    /* Initialize cmc */
    cmc::CmcInitialize(cmc::kMinimumInitialization);
    {

    cmc::CmcType data_type;
    std::string input_file;
    std::string output_file;
    float abs_error;
    int dim;
    std::string read_dim_lengths;
    std::vector<size_t> dim_lengths;
    
    int opt;
    while ((opt = getopt(argc, argv, "t:i:o:a:d:s:h")) != -1)
    {
        switch (opt)
        {
            case 't':
                cmc::cmc_debug_msg("Data type File: ", std::string(optarg));
                data_type = EvaluateDataType(std::string(optarg));
                break;
            case 'i':
                cmc::cmc_debug_msg("Input File: ", std::string(optarg));
                input_file = std::string(optarg);
                break;
            case 'o':
                cmc::cmc_debug_msg("Output File Prefix: ", std::string(optarg));
                output_file = std::string(optarg);
                break;
            case 'a':
                cmc::cmc_debug_msg("Absolute Permitted Error: ", std::string(optarg));
                abs_error = std::stof(optarg);
                break;
            case 'd':
                cmc::cmc_debug_msg("Dimensionality: ", std::string(optarg));
                dim = atoi(optarg);
                break;
            case 's':
                cmc::cmc_debug_msg("Dimensionality Shape String: ", std::string(optarg));
                read_dim_lengths = std::string(optarg);
                break;
            case 'h':
                DisplayHelpMessage();
                return 0;
                break;
            case '?':
                if (optopt == 't' || optopt == 'i' || optopt == 'o' || optopt == 'a' || optopt == 'd' || optopt == 's') {
                    cmc::cmc_err_msg("The option -", static_cast<char>(optopt), " requires an argument.");
                } else {
                    cmc::cmc_err_msg("An unknown option: ", static_cast<char>(optopt), " has been passed.");
                }
                return 1;
            default:
                return 1;
        }
    }

    /* Check if all arguments have been passed */
    if (argc < kNumArgCRequired)
    {
        DisplayHelpMessage();
        cmc::cmc_err_msg("The application requires arguments to execute the patch-based lossless multi resolution extraction compression.");
    }

    /* Check the input file */
    if (input_file.empty())
    {
        cmc::cmc_err_msg("An input file must be specified.");
    }

    /* Check if the input file exists */
    if (const std::filesystem::path input_file_path(input_file); not std::filesystem::exists(input_file_path))
    {
        DisplayHelpMessage();
        cmc::cmc_err_msg("The file ", input_file, " does not exist.");
    }

    /* Check the output file */
    if (output_file.empty())
    {
        cmc::cmc_err_msg("A output file prefix must be specified.");
    }

    /* Check if the dimensionality is supported */
    if (dim < 2 || dim > 3)
    {
        DisplayHelpMessage();
        cmc::cmc_err_msg("The specified dimensionality (", dim, ") is not supported. Currently, only 2D and 3D compression is supported.");
    }

    int offset{0};
    /* Check the shape string of the dimensionality */
    for (int dim_iter{0}; dim_iter < dim; ++dim_iter)
    {
        std::size_t pos = read_dim_lengths.find(" ", offset);
        if (pos == read_dim_lengths.npos)
        {
            pos = read_dim_lengths.length();
            
            if (dim_iter != dim -1)
            {
                DisplayHelpMessage();
                cmc::cmc_err_msg("The dimension shape string does not hold enough dimension lengths for the specified dimensionality (", dim, ").");
            }
        }

        const std::string dim_length = read_dim_lengths.substr(offset, pos);
        const size_t dimension_length = static_cast<size_t>(std::stoi(dim_length));
        dim_lengths.push_back(dimension_length);
        cmc::cmc_debug_msg("Extracted Dimension length: ", dimension_length);
        offset = pos + 1;
    }

    /* Perform the compresison */
    Compress(data_type, input_file, output_file, dim, dim_lengths, abs_error);

    DisplayCompressionResult(data_type, dim_lengths, output_file);
    
    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}
