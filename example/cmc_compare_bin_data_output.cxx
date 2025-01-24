#include "utilities/cmc_utilities.hxx"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cstring>

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
    
    const std::vector<float> initial_data = ReadDataFromStream<float>("initial_cmc_sfc_data.bin");
    const std::vector<float> decompressed_data = ReadDataFromStream<float>("decompressed_cmc_sfc_data.bin");

    std::cout << "Initial data size: " << initial_data.size() << std::endl;
    std::cout << "Decompressed data size: " << decompressed_data.size() << std::endl;

    const int cmp = std::memcmp(initial_data.data(), decompressed_data.data(), decompressed_data.size() * sizeof(float));
    std::cout << "Comparison of both arrays with memcmp yields return value: " << cmp << std::endl;

    for (size_t idx = 0; idx < initial_data.size(); ++idx)
    {
        if (initial_data[idx] != decompressed_data[idx])
        {
            std::cout << "No equality at " << idx << ", Initial Value: " << initial_data[idx] << ", Decompressed Value: " << decompressed_data[idx] << std::endl;
        }
    }

    //for (auto i = 0; i < decompressed_data.size(); ++i)
    //{
    //    std::cout << "Index i: " << i << ", initial_value: " << initial_data[i] << ", decompressed value: " << decompressed_data[i] << std::endl;
    //}
    

   return 0;
}