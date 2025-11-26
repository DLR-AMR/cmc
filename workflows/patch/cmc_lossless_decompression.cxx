#include "cmc.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "patch/lossless/cmc_patch_prefix_extraction_plain_suffixes_decompression.hxx"
#include "patch/lossless/cmc_patch_multi_res_extraction_decompression.hxx"   
#include "compression_io/cmc_decompression_serial_input.hxx"

#include <unistd.h>
#include <cstdio>
#include <numeric>
#include <algorithm>
#include <memory>
#include <vector>
#include <filesystem>
#include <type_traits>

inline void
DisplayHelpMessage()
{
    cmc::cmc_msg("This is cmc (version: ", CMC_VERSION, ").");
    cmc::cmc_msg("The following options are required to run the patch-based lossless decompression with serial output!");
    cmc::cmc_msg("");
    cmc::cmc_msg("\t-t \tThe data type of the data to decompress, possible options are float, double, int16_t, uint16_t, int32_t, uint32_t, int64_t or uint64_t");
    cmc::cmc_msg("\t-i \tThe path prefix to the input file storing the compressed data");
    cmc::cmc_msg("\t-o \tThe path to the output file which will store the decompressed data");
    cmc::cmc_msg("\t-v \tThe name of the variable to compress");
    cmc::cmc_msg("\t-h \tTo display this help message");
}

template <typename T>
auto
WriteOutDecompressedData(const std::string& output_file, const std::vector<T>& data)
-> std::enable_if_t<std::is_fundamental_v<T>, void>
{
    /* Write this decompressed data out to disk */
    std::FILE* file_out = std::fopen(output_file.c_str(), "wb");
    std::fwrite(data.data(), sizeof(T), data.size(), file_out);
    std::fclose(file_out);
}

void
Decompress(const cmc::CmcType data_type, const std::string& input_file_prefix, const std::string& output_file, const std::string& var_name)
{
    switch (data_type)
    {
        case cmc::CmcType::Float:
        {
            cmc::compression_io::serial::Reader reader(input_file_prefix);
            std::unique_ptr<cmc::patch::IPatchDecompressionVariable<float>> decompression_var = reader.ReadPatchVariableForDecompression<float>(var_name);
            decompression_var->Decompress();
            const std::vector<float> decompressed_data = decompression_var->GetDecompressedData();
            WriteOutDecompressedData<float>(output_file, decompressed_data);
        }
        break;
        case cmc::CmcType::Double:
        {
            cmc::compression_io::serial::Reader reader(input_file_prefix);
            std::unique_ptr<cmc::patch::IPatchDecompressionVariable<double>> decompression_var = reader.ReadPatchVariableForDecompression<double>(var_name);
            decompression_var->Decompress();
            const std::vector<double> decompressed_data = decompression_var->GetDecompressedData();
            WriteOutDecompressedData<double>(output_file, decompressed_data);
        }
        break;
        case cmc::CmcType::Int16_t:
        {
            cmc::compression_io::serial::Reader reader(input_file_prefix);
            std::unique_ptr<cmc::patch::IPatchDecompressionVariable<int16_t>> decompression_var = reader.ReadPatchVariableForDecompression<int16_t>(var_name);
            decompression_var->Decompress();
            const std::vector<int16_t> decompressed_data = decompression_var->GetDecompressedData();
            WriteOutDecompressedData<int16_t>(output_file, decompressed_data);
        }
        break;
        case cmc::CmcType::Uint16_t:
        {
            cmc::compression_io::serial::Reader reader(input_file_prefix);
            std::unique_ptr<cmc::patch::IPatchDecompressionVariable<uint16_t>> decompression_var = reader.ReadPatchVariableForDecompression<uint16_t>(var_name);
            decompression_var->Decompress();
            const std::vector<uint16_t> decompressed_data = decompression_var->GetDecompressedData();
            WriteOutDecompressedData<uint16_t>(output_file, decompressed_data);
        }
        break;
        case cmc::CmcType::Int32_t:
        {
            cmc::compression_io::serial::Reader reader(input_file_prefix);
            std::unique_ptr<cmc::patch::IPatchDecompressionVariable<int32_t>> decompression_var = reader.ReadPatchVariableForDecompression<int32_t>(var_name);
            decompression_var->Decompress();
            const std::vector<int32_t> decompressed_data = decompression_var->GetDecompressedData();
            WriteOutDecompressedData<int32_t>(output_file, decompressed_data);
        }
        break;
        case cmc::CmcType::Uint32_t:
        {
            cmc::compression_io::serial::Reader reader(input_file_prefix);
            std::unique_ptr<cmc::patch::IPatchDecompressionVariable<uint32_t>> decompression_var = reader.ReadPatchVariableForDecompression<uint32_t>(var_name);
            decompression_var->Decompress();
            const std::vector<uint32_t> decompressed_data = decompression_var->GetDecompressedData();
            WriteOutDecompressedData<uint32_t>(output_file, decompressed_data);
        }
        break;
        case cmc::CmcType::Int64_t:
        {
            cmc::compression_io::serial::Reader reader(input_file_prefix);
            std::unique_ptr<cmc::patch::IPatchDecompressionVariable<int64_t>> decompression_var = reader.ReadPatchVariableForDecompression<int64_t>(var_name);
            decompression_var->Decompress();
            const std::vector<int64_t> decompressed_data = decompression_var->GetDecompressedData();
            WriteOutDecompressedData<int64_t>(output_file, decompressed_data);
        }
        break;
        case cmc::CmcType::Uint64_t:
        {
            cmc::compression_io::serial::Reader reader(input_file_prefix);
            std::unique_ptr<cmc::patch::IPatchDecompressionVariable<uint64_t>> decompression_var = reader.ReadPatchVariableForDecompression<uint64_t>(var_name);
            decompression_var->Decompress();
            const std::vector<uint64_t> decompressed_data = decompression_var->GetDecompressedData();
            WriteOutDecompressedData<uint64_t>(output_file, decompressed_data);
        }
        break;
        default:
            DisplayHelpMessage();
            cmc::cmc_err_msg("A not supported/recognized data type has been specified.");
        break;
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

constexpr int kNumArgCRequired = 9;

int
main(int argc, char *argv[])
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    {
    
    cmc::CmcType data_type;
    std::string input_file;
    std::string output_file;
    std::string var_name;

    int opt;
    while ((opt = getopt(argc, argv, "t:i:o:v:h")) != -1)
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
            case 'v':
                cmc::cmc_debug_msg("Variable Name: ", std::string(optarg));
                var_name = std::string(optarg);
                break;
            case 'h':
                DisplayHelpMessage();
                return 0;
                break;
            case '?':
                if (optopt == 't' || optopt == 'i' || optopt == 'o' || optopt == 'v' ) {
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
        cmc::cmc_err_msg("The application requires arguments to execute the patch-based lossless prefix extraction compression.");
    }

    /* Check the input file */
    if (input_file.empty())
    {
        cmc::cmc_err_msg("An input file must be specified.");
    }

    /* Check if the input file exists */
    if (const std::filesystem::path input_file_path(cmc::compression_io::serial::CreateFileName(input_file, var_name));
        not std::filesystem::exists(input_file_path))
    {
        DisplayHelpMessage();
        cmc::cmc_err_msg("The file ", input_file, " does not exist.");
    }

    /* Check the output file */
    if (output_file.empty())
    {
        cmc::cmc_err_msg("A output file prefix must be specified.");
    }

    /* Check the variable name */
    if (var_name.empty())
    {
        cmc::cmc_err_msg("A variable name must be specified.");
    }

    /* Perform the decompresison */
    Decompress(data_type, input_file, output_file, var_name);

    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}


