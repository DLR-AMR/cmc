#include "cmc.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "patch/lossless/cmc_patch_prefix_extraction_plain_suffixes_compression.hxx"
#include "input/cmc_binary_reader.hxx"
#include "compression_io/cmc_compression_serial_output.hxx"

#include <unistd.h>
#include <cstdio>
#include <numeric>
#include <algorithm>
#include <vector>
#include <filesystem>

inline void
DisplayHelpMessage()
{
    cmc::cmc_msg("This is cmc (version: ", CMC_VERSION, ").");
    cmc::cmc_msg("The following options are required to run the patch-based lossless prefix extraction compression with serial output!");
    cmc::cmc_msg("");
    cmc::cmc_msg("\t-t \tThe data type of the data to compress, possible options are float, double, int16_t, uint16_t, int32_t, uint32_t, int64_t or uint64_t");
    cmc::cmc_msg("\t-i \tThe path to the input file storing the array of binary data to compress");
    cmc::cmc_msg("\t-o \tThe prefix of the output file storing the compressed data (the resulting file_name will be $file_prefix + '_' + $variable_name + '.cmc')");
    cmc::cmc_msg("\t-v \tThe name of the variable to compress");
    cmc::cmc_msg("\t-d \tThe dimensionality of the array storing the binary data, e.g. 3 for 3D data");
    cmc::cmc_msg("\t-s \tA string representing the shape of the data, e.g. '100 200 500' for a 3D array of size 100x200x500,"
                 " such that the slowest varying dimension is the first (e.g. 100) and the fastest varying dimension is the last parameter (e.g. 500)");
    cmc::cmc_msg("\t-h \tTo display this help message");
}

void
Compress(const cmc::CmcType data_type, const std::string& input_file, const std::string& output_file, const std::string& var_name, const int dim, const std::vector<size_t>& dim_lengths)
{
    cmc::cmc_debug_msg("Performing Patch-Based Lossless Prefix Extraction Compression.");

    /* Generate input variables from the binary file */
    cmc::input::binary::Reader binary_reader(input_file);

    const size_t num_elements = std::reduce(dim_lengths.begin(), dim_lengths.end(), 1, std::multiplies<size_t>());
    cmc::cmc_debug_msg("Number of elements in the variable: ", num_elements);
    const int arbitrary_var_id = 0;
    const cmc::DataLayout layout = cmc::GetDefaultDataLayout(dim);
    const cmc::GeoDomain domain = cmc::GetDefaultDomain(layout, dim_lengths);
    cmc::input::Var variable = binary_reader.CreateVariableFromBinaryData(data_type, var_name, arbitrary_var_id, num_elements, layout, domain);
    std::vector<cmc::input::Var> input_variables{std::move(variable)};

    switch (dim)
    {
        case 3:
            switch (data_type)
            {
                case cmc::CmcType::Float:
                {
                    cmc::patch::lossless::prefix::plain_suffix::PatchCompressionVariable<float, 3> var(input_variables.front());
                    var.Compress();
                    cmc::compression_io::serial::Writer writer(output_file);
                    writer.SetVariable(&var);
                    writer.Write();
                }
                break;
                case cmc::CmcType::Double:
                {
                    cmc::patch::lossless::prefix::plain_suffix::PatchCompressionVariable<double, 3> var(input_variables.front());
                    var.Compress();
                    cmc::compression_io::serial::Writer writer(output_file);
                    writer.SetVariable(&var);
                    writer.Write();
                }
                break;
                case cmc::CmcType::Int16_t:
                {
                    cmc::patch::lossless::prefix::plain_suffix::PatchCompressionVariable<int16_t, 3> var(input_variables.front());
                    var.Compress();
                    cmc::compression_io::serial::Writer writer(output_file);
                    writer.SetVariable(&var);
                    writer.Write();
                }
                break;
                case cmc::CmcType::Uint16_t:
                {
                    cmc::patch::lossless::prefix::plain_suffix::PatchCompressionVariable<uint16_t, 3> var(input_variables.front());
                    var.Compress();
                    cmc::compression_io::serial::Writer writer(output_file);
                    writer.SetVariable(&var);
                    writer.Write();
                }
                break;
                case cmc::CmcType::Int32_t:
                {
                    cmc::patch::lossless::prefix::plain_suffix::PatchCompressionVariable<int32_t, 3> var(input_variables.front());
                    var.Compress();
                    cmc::compression_io::serial::Writer writer(output_file);
                    writer.SetVariable(&var);
                    writer.Write();
                }
                break;
                case cmc::CmcType::Uint32_t:
                {
                    cmc::patch::lossless::prefix::plain_suffix::PatchCompressionVariable<uint32_t, 3> var(input_variables.front());
                    var.Compress();
                    cmc::compression_io::serial::Writer writer(output_file);
                    writer.SetVariable(&var);
                    writer.Write();
                }
                break;
                case cmc::CmcType::Int64_t:
                {
                    cmc::patch::lossless::prefix::plain_suffix::PatchCompressionVariable<int64_t, 3> var(input_variables.front());
                    var.Compress();
                    cmc::compression_io::serial::Writer writer(output_file);
                    writer.SetVariable(&var);
                    writer.Write();
                }
                break;
                case cmc::CmcType::Uint64_t:
                {
                    cmc::patch::lossless::prefix::plain_suffix::PatchCompressionVariable<uint64_t, 3> var(input_variables.front());
                    var.Compress();
                    cmc::compression_io::serial::Writer writer(output_file);
                    writer.SetVariable(&var);
                    writer.Write();
                }
                break;
                default:
                    DisplayHelpMessage();
                    cmc::cmc_err_msg("A not supported/recognized data type has been specified.");
                break;
            }
        break;
        default:
            cmc::cmc_err_msg("Only 3D compression is currently supported.");
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

constexpr int kNumArgCRequired = 13;

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
    int dim;
    std::string read_dim_lengths;
    std::vector<size_t> dim_lengths;

    int opt;
    while ((opt = getopt(argc, argv, "t:i:o:v:d:s:h")) != -1)
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
                if (optopt == 't' || optopt == 'i' || optopt == 'o' || optopt == 'v' || optopt == 'd' || optopt == 's') {
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

    /* Check the variable name */
    if (var_name.empty())
    {
        cmc::cmc_err_msg("A variable name must be specified.");
    }

    /* Check if the dimensionality is supported */
    if (dim != 3)
    {
        DisplayHelpMessage();
        cmc::cmc_err_msg("The specified dimensionality (", dim, ") is not supported. Currently, only 3D compression is supported.");
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
    Compress(data_type, input_file, output_file, var_name, dim, dim_lengths);

    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}






#if 0





int
main()
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    {

    /* Read in data */
    #if 1
    const std::string file = "../data/100x500x500/TCf48.bin.f32";
    //const std::string file = "../data/era5_reanalysis_t2m_tc03_13_12_23.nc";
    #else
    const std::string file = "../data/SDRBENCH-EXASKY-NYX-512x512x512/baryon_density.f32";
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
    //cmc::CmcUniversalType missing_value(static_cast<float>(40.0)); //UF
    //cmc::CmcUniversalType missing_value(static_cast<float>(48.25)); //VF
    //cmc::CmcUniversalType missing_value(static_cast<float>(13.4)); //WF
    
    //Missing values NYX data
    //cmc::CmcUniversalType missing_value(static_cast<float>(31868.0)); //velocity x
    //cmc::CmcUniversalType missing_value(static_cast<float>(56507.0)); //velocity y
    //cmc::CmcUniversalType missing_value(static_cast<float>(33387.0)); //velocity z
    //cmc::CmcUniversalType missing_value(static_cast<float>(4784.0)); //temp
    //cmc::CmcUniversalType missing_value(static_cast<float>(13780.0)); //dark matter
    //cmc::CmcUniversalType missing_value(static_cast<float>(115864.0)); //baryon denisty


    const size_t num_elements = 500 * 500 * 100;
    cmc::DataLayout layout(cmc::DataLayout::Lev_Lat_Lon);
    cmc::GeoDomain domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 500),
                          cmc::DimensionInterval(cmc::Dimension::Lat, 0, 500),
                          cmc::DimensionInterval(cmc::Dimension::Lev, 0, 100)
                          );

    cmc::CmcUniversalType arbitrary_missing_value(static_cast<float>(-1000.0));

    /* Generate input variables from the binary file */
    cmc::input::binary::Reader binary_reader(file);
    cmc::input::Var variable = binary_reader.CreateVariableFromBinaryData(type, name, id, num_elements, arbitrary_missing_value, layout, domain);
    variable.SetMPIComm(MPI_COMM_SELF);
    std::vector<cmc::input::Var> input_variables{std::move(variable)};


    {
    /* Setup an embedded PrefixAMR (with plain suffix encoding) compression variable from the input variables */           
    //cmc::patch::lossless::prefix::plain_suffix::PatchCompressionVariable<float, 3> var(input_variables.front());
    cmc::patch::lossless::multi_res::PatchCompressionVariable<float, 3> var(input_variables.front());

    /* Perform the compression */
    var.Compress();

#if 1

    #if 1
    {
        /* Write out the compressed data to disk */
        cmc::compression_io::serial::Writer writer("prefix_example_lossless_compression_output");
        writer.SetVariable(&var);
        writer.Write();
    }
    #endif
    }

    #if 1

    /* Create a reader for the compressed output that has been stored */
    cmc::compression_io::serial::Reader reader("prefix_example_lossless_compression_output");
    
    /* Create an embedded decompressor from the compressed data */
    std::unique_ptr<cmc::patch::IPatchDecompressionVariable<float>> decompression_var = reader.ReadPatchVariableForDecompression<float>(name);

    /* Decompress the encoded data */
    decompression_var->Decompress();

    #if 1
    /* De-Mortonize the data and obtain the initial ordering of the input data */
    const std::vector<float> decompressed_data = decompression_var->GetDecompressedData();
    
    /* Write this decompressed data out to disk, in order to be able to compare it to the intiial data */
    FILE* file_out = fopen("decompressed_data.cmc", "wb");
    fwrite(decompressed_data.data(), sizeof(float), decompressed_data.size(), file_out);
    fclose(file_out);

    cmc::cmc_debug_msg("Size of decompressed data: ", decompressed_data.size());
    #endif

    #endif

#endif

    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}
#endif