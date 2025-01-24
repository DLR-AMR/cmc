#ifndef CMC_BINARY_READER_HXX
#define CMC_BINARY_READER_HXX

#include "utilities/cmc_binary_reader_forward.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_dimension_interval.hxx"
#include "utilities/cmc_geo_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"

#include <fstream>
#include <vector>
#include <variant>

namespace cmc
{

namespace bin_reader
{

struct ReadData;

class Reader
{
public:
    Reader() = delete;

    Reader(const std::string& file_name)
    : file_name_{file_name} {};

    ~Reader() = default;

    InputVar CreateVariableFromBinaryData(const CmcType type, const std::string& name, const int id, const size_t num_elements, const CmcUniversalType missing_value,
                                          const DataLayout layout, const GeoDomain& domain) const;

private:
    const std::string file_name_;
};

//TODO: Endianness needs to be handled
template<typename T>
inline std::vector<T>
ReadBinaryData(const size_t num_elements, const std::string& file_name)
{
    std::vector<T> values(num_elements);

    //TODO: Check for existence 
    
    /* Open the file for input */
    std::ifstream in(file_name, std::ios::binary);
    in.seekg(0, std::ios::beg);

    /* Read in the data */
    in.read(reinterpret_cast<char*>(&values[0]), num_elements * sizeof(T));

    for (size_t index = 0; index < 10; ++index)
    {
        cmc_debug_msg("Binary Data: Index: ", index, " has value: ", values[index]);
    }

    //TODO: Endianness reordering potentially

    return values;
}

struct ReadData
{
public:
    ReadData() = delete;
    ReadData(const size_t num_elements, const std::string& file_name)
    : num_elements_{num_elements}, file_name_{file_name} {};

    void operator()(InputVariable<int8_t>& var) {
        //var.SetData(AccessKey(), ReadBinaryData<int8_t>(num_elements_, file_name_));
        std::vector<int8_t> data = ReadBinaryData<int8_t>(num_elements_, file_name_);
        Hyperslab domain = TransformGeoDomainToHyperslab(var.GetGlobalDomain());
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(domain)});
    }
    void operator()(InputVariable<char>& var) {
        //var.SetData(AccessKey(), ReadBinaryData<char>(num_elements_, file_name_));
    }
    void operator()(InputVariable<int16_t>& var) {
        //var.SetData(AccessKey(), ReadBinaryData<int16_t>(num_elements_, file_name_));
    }
    void operator()(InputVariable<int32_t>& var) {
        //var.SetData(AccessKey(), ReadBinaryData<int32_t>(num_elements_, file_name_));
    }
    void operator()(InputVariable<float>& var) {
        //var.SetData(AccessKey(), ReadBinaryData<float>(num_elements_, file_name_));
        std::vector<float> data = ReadBinaryData<float>(num_elements_, file_name_);
        Hyperslab domain = TransformGeoDomainToHyperslab(var.GetGlobalDomain());
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(domain)});
    }
    void operator()(InputVariable<double>& var) {
        //var.SetData(AccessKey(), ReadBinaryData<double>(num_elements_, file_name_));
    }
    void operator()(InputVariable<uint8_t>& var) {
        //var.SetData(AccessKey(), ReadBinaryData<uint8_t>(num_elements_, file_name_));
    }
    void operator()(InputVariable<uint16_t>& var) {
        //var.SetData(AccessKey(), ReadBinaryData<uint16_t>(num_elements_, file_name_));
    }
    void operator()(InputVariable<uint32_t>& var) {
        //var.SetData(AccessKey(), ReadBinaryData<uint32_t>(num_elements_, file_name_));
    }
    void operator()(InputVariable<int64_t>& var) {
        //var.SetData(AccessKey(), ReadBinaryData<int64_t>(num_elements_, file_name_));
    }
    void operator()(InputVariable<uint64_t>& var) {
        //var.SetData(AccessKey(), ReadBinaryData<uint64_t>(num_elements_, file_name_));
    }
    //template<typename T> void operator()(T& var)
    //{
    //    cmc_err_msg("This type of InputVariable is not suppported.");
    //}
private:
    const size_t num_elements_;
    const std::string& file_name_;
};


inline InputVar
Reader::CreateVariableFromBinaryData(const CmcType type, const std::string& name, const int id, const size_t num_elements, const CmcUniversalType missing_value,
                                     const DataLayout layout, const GeoDomain& domain) const
{
    /* Create a variable with the given specifications */
    InputVar variable(type, name, id, num_elements, missing_value, layout, domain);

    cmc_debug_msg("The variable holds: ", variable.GetInternalVariant().index());
    /* Read the data into the variable */
    //CmcInputVariable& var = variable.GetInternalVariant();
    std::visit(ReadData(num_elements, this->file_name_), variable.GetInternalVariant(AccessKey()));

    return variable;
}

}

}

#endif /* !CMC_BINARY_READER_HXX */
