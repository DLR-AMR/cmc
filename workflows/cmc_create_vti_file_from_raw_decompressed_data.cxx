#include "cmc.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "patch/lossless/cmc_fast_multi_res_extraction.hxx"
#include "input/cmc_binary_file_reader.hxx"

#include <unistd.h>
#include <cstdio>
#include <numeric>
#include <algorithm>
#include <vector>
#include <string>
#include <filesystem>
#include <bit>
#include <fstream>
#include <cstdint>

inline void
DisplayHelpMessage()
{
    cmc::cmc_msg("This is cmc (version: ", CMC_VERSION, ").");
    cmc::cmc_msg("The following options are required to run the patch-based lossless multi resolution extraction compression with serial output!");
    cmc::cmc_msg("");
    cmc::cmc_msg("\t-t \tThe data type of the data to compress, possible options are float or double");
    cmc::cmc_msg("\t-i \tThe path to the input file storing the array of binary data to write into an .vti-file");
    cmc::cmc_msg("\t-o \tThe name of the output .vti-file");
    cmc::cmc_msg("\t-d \tThe dimensionality of the array storing the binary data, e.g. 3 for 3D data");
    cmc::cmc_msg("\t-s \tA string representing the shape of the data, e.g. '100 200 500' for a 3D array of size 100x200x500,"
                 " such that the slowest varying dimension is the first (e.g. 100) and the fastest varying dimension is the last parameter (e.g. 500)");
    cmc::cmc_msg("\t-e \tThe endianness of the data, supported options are little or big");
    cmc::cmc_msg("\t-c \tThe scaling of the spacing between points, must be a positive integer");
    cmc::cmc_msg("\t-h \tTo display this help message");
}

void
WriteVTI(const cmc::CmcType data_type, const std::string& input_file, const std::string& output_file, const int dim, const std::vector<size_t>& dim_lengths, const std::endian endianness_of_data_in_file, const int scale_spacing)
{
    cmc::cmc_debug_msg("Creating .vti file from raw data.");

    /* Create the file */
    std::ofstream file_out(output_file, std::ios::binary);

    /* Compute the number of raw bytes */
    const size_t num_elems = std::reduce(dim_lengths.begin(), dim_lengths.end(), 1, std::multiplies<size_t>());
    const uint64_t num_bytes = num_elems * CmcTypeToBytes(data_type);

    /* Write the header */
    file_out << R"(<?xml version="1.0"?>)" << "\n";
    if (endianness_of_data_in_file == std::endian::little)
    {
        file_out << R"(<VTKFile type="ImageData"
version="1.0"
byte_order="LittleEndian"
header_type="UInt64">)" << "\n";
    } else
    {
        file_out << R"(<VTKFile type="ImageData"
version="1.0"
byte_order="BigEndian"
header_type="UInt64">)" << "\n";
    }

    const double spacing_x = 1.0 * scale_spacing;
    const double spacing_y = 1.0 * scale_spacing;
    const double spacing_z = 1.0 * scale_spacing;

    /* Handle the dimensionality and the data type */
    switch (dim)
    {
        /** 1D Compression **/
        case 1:
        {
            /* Write the data extent */
            file_out << "<ImageData WholeExtent=\""
                     << "0 " << dim_lengths[0] -1 << " "
                     << "0 " << 0 << " "
                     << "0 " << 0 << "\" "
                     << "Origin=\"0 0 0\" "
                     << "Spacing=\""
                     << spacing_x << " "
                     << spacing_y << " "
                     << spacing_z <<"\">\n";
            file_out << "<Piece Extent=\""
                     << "0 " << dim_lengths[0] -1 << " "
                     << "0 " << 0 << " "
                     << "0 " << 0 << "\">\n";



            /* Write the data header */
            file_out << R"(<PointData Scalar="values">)" << "\n";


            if (data_type == cmc::CmcType::Float)
            {
                file_out << R"(<DataArray type="Float32"
Name="values"
NumberOfComponents="1"
format="appended"
offset="0"/>)" << "\n";
            } else
            {
                file_out << R"(<DataArray type="Float64"
Name="values"
NumberOfComponents="1"
format="appended"
offset="0"/>)" << "\n";
            }

            file_out << "</PointData>\n";
            file_out << "<CellData/>\n";
            file_out << "</Piece>\n";
            file_out << "</ImageData>\n";

            file_out << R"(<AppendedData encoding="raw">)" << "\n";
            file_out.put('_');
            file_out.write(reinterpret_cast<const char*>(&num_bytes), sizeof(num_bytes));

            if (data_type == cmc::CmcType::Float)
            {
                /* Create the binary reader */
                const std::array<int32_t, 1> dims{static_cast<int32_t>(dim_lengths[0])};
                cmc::input::binary_file::Reader<float, 1> binary_reader(input_file, dims, endianness_of_data_in_file);

                /* Read the data */
                const std::vector<float> data = binary_reader.ReadData();

                /* Write raw data */
                file_out.write(reinterpret_cast<const char*>(data.data()), static_cast<std::streamsize>(num_bytes));
            } else
            {
                /* Create the binary reader */
                const std::array<int32_t, 1> dims{static_cast<int32_t>(dim_lengths[0])};
                cmc::input::binary_file::Reader<double, 1> binary_reader(input_file, dims, endianness_of_data_in_file);

                /* Read the data */
                const std::vector<double> data = binary_reader.ReadData();

                /* Write raw data */
                file_out.write(reinterpret_cast<const char*>(data.data()), static_cast<std::streamsize>(num_bytes));
            }
        }
        break;
        /** 2D Compression **/
        case 2:
        {    
            /* Write the data extent */
            file_out << "<ImageData WholeExtent=\""
                     << "0 " << dim_lengths[1] -1 << " "
                     << "0 " << dim_lengths[0] -1 << " "
                     << "0 " << 0 << "\" "
                     << "Origin=\"0 0 0\" "
                     << "Spacing=\""
                     << spacing_x << " "
                     << spacing_y << " "
                     << spacing_z <<"\">\n";
            file_out << "<Piece Extent=\""
                     << "0 " << dim_lengths[1] -1 << " "
                     << "0 " << dim_lengths[0] -1 << " "
                     << "0 " << 0 << "\">\n";



            /* Write the data header */
            file_out << R"(<PointData Scalar="values">)" << "\n";


            if (data_type == cmc::CmcType::Float)
            {
                file_out << R"(<DataArray type="Float32"
Name="values"
NumberOfComponents="1"
format="appended"
offset="0"/>)" << "\n";
            } else
            {
                file_out << R"(<DataArray type="Float64"
Name="values"
NumberOfComponents="1"
format="appended"
offset="0"/>)" << "\n";
            }

            file_out << "</PointData>\n";
            file_out << "<CellData/>\n";
            file_out << "</Piece>\n";
            file_out << "</ImageData>\n";

            file_out << R"(<AppendedData encoding="raw">)" << "\n";
            file_out.put('_');
            file_out.write(reinterpret_cast<const char*>(&num_bytes), sizeof(num_bytes));

            if (data_type == cmc::CmcType::Float)
            {
                /* Create the binary reader */
                const std::array<int32_t, 2> dims{static_cast<int32_t>(dim_lengths[0]), static_cast<int32_t>(dim_lengths[1])};
                cmc::input::binary_file::Reader<float, 2> binary_reader(input_file, dims, endianness_of_data_in_file);

                /* Read the data */
                const std::vector<float> data = binary_reader.ReadData();

                /* Write raw data */
                file_out.write(reinterpret_cast<const char*>(data.data()), static_cast<std::streamsize>(num_bytes));
            } else
            {
                /* Create the binary reader */
                const std::array<int32_t, 2> dims{static_cast<int32_t>(dim_lengths[0]), static_cast<int32_t>(dim_lengths[1])};
                cmc::input::binary_file::Reader<double, 2> binary_reader(input_file, dims, endianness_of_data_in_file);

                /* Read the data */
                const std::vector<double> data = binary_reader.ReadData();

                /* Write raw data */
                file_out.write(reinterpret_cast<const char*>(data.data()), static_cast<std::streamsize>(num_bytes));
            }
            }
        break;
        /** 3D Compression **/
        case 3:
        {
            /* Write the data extent */
            file_out << "<ImageData WholeExtent=\""
                     << "0 " << dim_lengths[2] -1 << " "
                     << "0 " << dim_lengths[1] -1 << " "
                     << "0 " << dim_lengths[0] -1  << "\" "
                     << "Origin=\"0 0 0\" "
                     << "Spacing=\""
                     << spacing_x << " "
                     << spacing_y << " "
                     << spacing_z <<"\">\n";
            file_out << "<Piece Extent=\""
                     << "0 " << dim_lengths[2] -1 << " "
                     << "0 " << dim_lengths[1] -1 << " "
                     << "0 " << dim_lengths[0] -1  << "\">\n";



            /* Write the data header */
            file_out << R"(<PointData Scalar="values">)" << "\n";


            if (data_type == cmc::CmcType::Float)
            {
                file_out << R"(<DataArray type="Float32"
Name="values"
NumberOfComponents="1"
format="appended"
offset="0"/>)" << "\n";
            } else
            {
                file_out << R"(<DataArray type="Float64"
Name="values"
NumberOfComponents="1"
format="appended"
offset="0"/>)" << "\n";
            }

            file_out << "</PointData>\n";
            file_out << "<CellData/>\n";
            file_out << "</Piece>\n";
            file_out << "</ImageData>\n";

            file_out << R"(<AppendedData encoding="raw">)" << "\n";
            file_out.put('_');
            file_out.write(reinterpret_cast<const char*>(&num_bytes), sizeof(num_bytes));

            if (data_type == cmc::CmcType::Float)
            {
                /* Create the binary reader */
                const std::array<int32_t, 3> dims{static_cast<int32_t>(dim_lengths[0]), static_cast<int32_t>(dim_lengths[1]), static_cast<int32_t>(dim_lengths[2])};
                cmc::input::binary_file::Reader<float, 3> binary_reader(input_file, dims, endianness_of_data_in_file);

                /* Read the data */
                const std::vector<float> data = binary_reader.ReadData();

                /* Write raw data */
                file_out.write(reinterpret_cast<const char*>(data.data()), static_cast<std::streamsize>(num_bytes));
            } else
            {
                /* Create the binary reader */
                const std::array<int32_t, 3> dims{static_cast<int32_t>(dim_lengths[0]), static_cast<int32_t>(dim_lengths[1]), static_cast<int32_t>(dim_lengths[2])};
                cmc::input::binary_file::Reader<double, 3> binary_reader(input_file, dims, endianness_of_data_in_file);

                /* Read the data */
                const std::vector<double> data = binary_reader.ReadData();

                /* Write raw data */
                file_out.write(reinterpret_cast<const char*>(data.data()), static_cast<std::streamsize>(num_bytes));
            }
        }
        break;
        /** 4D Compression **/
        case 4:
            DisplayHelpMessage();
            cmc::cmc_err_msg("4D Conversion of raw floating-point data to .vti is not possible!");
        break;
        default:
            cmc::cmc_err_msg("Only 1D, 2D and 3D compression is currently supported.");
    }

    /* End writing */
    file_out << "\n";
    file_out << "</AppendedData>\n";
    file_out <<  "</VTKFile>\n";
    file_out.close();

    cmc::cmc_debug_msg("The file ", output_file, " has been written.");
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
        cmc::cmc_err_msg("The specified data type (", type, ") is not recognized. Only float and double data types are currently supported.");
        return cmc::CmcType::TypeUndefined;
    }
}

std::endian
EvaluateEndianness(const std::string& type)
{
    if (type.compare("little") == 0)
    {
        return std::endian::little;
    } else if (type.compare("big") == 0)
    {
        return std::endian::big;
    } else
    {
        DisplayHelpMessage();
        cmc::cmc_err_msg("The specified endianness (", type, ") is not recognized. Only little and big endian orderings are supported.");
        return std::endian::native;
    }
}

constexpr int kNumArgCRequired = 14;

int
main(int argc, char *argv[])
{
    /* Initialize cmc */
    cmc::CmcInitialize(cmc::kMinimumInitialization);
    {

    cmc::CmcType data_type;
    std::string input_file;
    std::string output_file;
    int dim;
    std::string read_dim_lengths;
    std::vector<size_t> dim_lengths;
    std::string endianness_string;
    int scale_spacing{1};

    int opt;
    while ((opt = getopt(argc, argv, "t:i:o:d:s:e:c:h")) != -1)
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
            case 'd':
                cmc::cmc_debug_msg("Dimensionality: ", std::string(optarg));
                dim = atoi(optarg);
                break;
            case 's':
                cmc::cmc_debug_msg("Dimensionality Shape String: ", std::string(optarg));
                read_dim_lengths = std::string(optarg);
                break;
            case 'e':
                cmc::cmc_debug_msg("Endianness: ", std::string(optarg));
                endianness_string = std::string(optarg);
                break;
            case 'c':
                cmc::cmc_debug_msg("Scale Spacing: ", std::string(optarg));
                scale_spacing = atoi(optarg);
                break;
            case 'h':
                DisplayHelpMessage();
                return 0;
                break;
            case '?':
                if (optopt == 't' || optopt == 'i' || optopt == 'o' || optopt == 'd' || optopt == 's' || optopt == 'e') {
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
        cmc::cmc_err_msg("The application requires arguments to perform the conversion from raw output to a .vti-file.");
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
    if (dim < 1 || dim > 3)
    {
        DisplayHelpMessage();
        cmc::cmc_err_msg("The specified dimensionality (", dim, ") is not supported. Currently, only 1D, 2D, and 3D conversion to .vti is supported.");
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

    /* Check sclaing of the spacing */
    if (scale_spacing < 1)
    {
        DisplayHelpMessage();
        cmc::cmc_err_msg("The scaling of the spacing must be a positive interger!");
    }

    /* Check the endianness of the data */
    const std::endian endianness = EvaluateEndianness(endianness_string);

    /* Write the .vti file */
    WriteVTI(data_type, input_file, output_file, dim, dim_lengths, endianness, scale_spacing);
    
    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}
