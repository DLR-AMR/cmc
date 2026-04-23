#ifndef CMC_INPUT_CMC_BINARY_FILE_READER_HXX
#define CMC_INPUT_CMC_BINARY_FILE_READER_HXX

#include "cmc.hxx"
#include "utilities/cmc_log_functions.hxx"

#include <concepts>
#include <string>
#include <vector>
#include <fstream>
#include <bit>
#include <array>
#include <algorithm>
#include <cstdint>
#include <execution>

namespace cmc::input::binary_file
{

template<typename T>
concept OneByteType = (sizeof(T) == 1);

template<typename T>
concept TwoByteType = (sizeof(T) == 2);

template<typename T>
concept FourByteType = (sizeof(T) == 4);

template<typename T>
concept EightByteType = (sizeof(T) == 8);

template<typename T>
concept IntegerType = (std::is_integral_v<T>);

template<typename T>
concept UnsignedIntegerType = (std::is_unsigned_v<T> && std::is_integral_v<T>);

template<typename T>
concept FundamentalArithmeticType = (std::is_arithmetic_v<T> && std::is_fundamental_v<T>);

template <FundamentalArithmeticType T, int DIMENSION>
requires (DIMENSION > 0)
class Reader
{
public:
    Reader() = delete;

    Reader(const std::string& file_name, const std::array<int32_t, DIMENSION>& dimension_lengths, const std::endian endianness_of_data_in_file = std::endian::native)
    : file_name_{file_name}, dim_lengths_(dimension_lengths), endianness_{endianness_of_data_in_file}  {}

    std::vector<T> ReadData() const;

private:
    const std::string file_name_;
    const std::array<int, DIMENSION> dim_lengths_;
    const std::endian endianness_;
};

template <OneByteType T>
inline void
ReorderBytesInPlace([[maybe_unused]] std::vector<T>& data)
{
    /* Nothing to be done here; just for completeness */
}

template <TwoByteType T>
inline void
ReorderBytesInPlace(std::vector<T>& data)
{
    /* Reorder the bytes */
    std::transform(std::execution::par_unseq, data.cbegin(), data.cend(), data.begin(), [](const T& value) -> T {
        const uint16_t uint_value = std::bit_cast<uint16_t>(value);
        const uint16_t swapped_uint_value = std::byteswap<uint16_t>(uint_value);
        return std::bit_cast<T>(swapped_uint_value);
    });
}

template <FourByteType T>
inline void
ReorderBytesInPlace(std::vector<T>& data)
{
    /* Reorder the bytes */
    std::transform(std::execution::par_unseq, data.cbegin(), data.cend(), data.begin(), [](const T& value) -> T {
        const uint32_t uint_value = std::bit_cast<uint32_t>(value);
        const uint32_t swapped_uint_value = std::byteswap<uint32_t>(uint_value);
        return std::bit_cast<T>(swapped_uint_value);
    });
}

template <EightByteType T>
inline void
ReorderBytesInPlace(std::vector<T>& data)
{
    /* Reorder the bytes */
    std::transform(std::execution::par_unseq, data.cbegin(), data.cend(), data.begin(), [](const T& value) -> T {
        const uint64_t uint_value = std::bit_cast<uint64_t>(value);
        const uint64_t swapped_uint_value = std::byteswap<uint64_t>(uint_value);
        return std::bit_cast<T>(swapped_uint_value);
    });
}

template <FundamentalArithmeticType T, int DIMENSION>
requires (DIMENSION > 0)
std::vector<T>
Reader<T, DIMENSION>::ReadData() const
{
    /* Compute the nuber of data elements */
    const int num_elements = std::reduce(this->dim_lengths_.begin(), this->dim_lengths_.end(), 1, std::multiplies<int>());

    /* Allocate a vector */
    std::vector<T> data(num_elements);

    /* Open the file, but delete it first, if it already exists */
    const std::filesystem::path input_file_path(this->file_name_);
    if (not std::filesystem::exists(input_file_path))
    {
        cmc_err_msg("The file ", this->file_name_, " does not exist!");
    }
    std::FILE* file_in = std::fopen(this->file_name_.c_str(), "rb");
    if (not file_in)
    {
        cmc_err_msg("Opening the file ", this->file_name_, " failed!");
    }

    /* Read the data */
    const std::size_t num_read_elems = std::fread(data.data(), sizeof(T), num_elements, file_in);
    if (num_elements != static_cast<int>(num_read_elems))
    {
        cmc_err_msg("An error occured: The number of elements read does not match the number of specified elements in the file!");
    }

    /* Close the file */
    const int rv_close = std::fclose(file_in);
    if (rv_close != 0)
    {
        cmc_err_msg("Closing the file ", this->file_name_, " failed!");
    }

    /* Potentially, reorder the data due to the endianness */
    if (this->endianness_ != std::endian::native)
    {
        static_assert(std::endian::native == std::endian::big || std::endian::native == std::endian::little,
                  "Only little-endian and big-endian systems are suppoprted!");
        
        ReorderBytesInPlace(data);
    }

    return data;
}

}

#endif /* !CMC_INPUT_CMC_BINARY_FILE_READER_HXX */
