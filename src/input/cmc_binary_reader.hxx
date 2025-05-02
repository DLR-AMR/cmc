#ifndef INPUT_CMC_BINARY_READER_HXX
#define INPUT_CMC_BINARY_READER_HXX

#include "utilities/cmc_binary_reader_forward.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_dimension_interval.hxx"
#include "utilities/cmc_geo_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "input/cmc_input_variable.hxx"

#include <fstream>
#include <vector>
#include <variant>

namespace cmc::input::binary
{

struct ReadData;

class Reader
{
public:
    Reader() = delete;

    Reader(const std::string& file_name)
    : file_name_{file_name} {};

    ~Reader() = default;

    input::Var CreateVariableFromBinaryData(const CmcType type, const std::string& name, const int id, const size_t num_elements, const CmcUniversalType missing_value,
                                          const DataLayout layout, const GeoDomain& domain) const;

    input::Var CreateSubDomainVariableFromBinaryData(const CmcType type, const std::string& name, const int id, const CmcUniversalType missing_value,
                                                   const DataLayout layout, const GeoDomain& global_domain, const GeoDomain& sub_domain) const;
private:
    const std::string file_name_;
};

//TODO: Endianness needs to be handled
template<typename T>
inline std::vector<T>
ReadBinaryData(const size_t num_elements, const std::string& file_name)
{
    std::vector<T> values(num_elements);

    /* Open the file for input */
    std::ifstream in(file_name, std::ios::binary);

    if (!in.is_open())
    {
        cmc_err_msg("The file ", file_name, " containing the binary data could not be opened.");
    }

    /* Move to the beginning of the file */
    in.seekg(0, std::ios::beg);

    /* Read in the data */
    in.read(reinterpret_cast<char*>(&values[0]), num_elements * sizeof(T));

    cmc_debug_msg("The first ten values of the data will be displayed:");
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

    void operator()(input::Variable<int8_t>& var) {
        std::vector<int8_t> data = ReadBinaryData<int8_t>(num_elements_, file_name_);
        Hyperslab domain = TransformGeoDomainToHyperslab(var.GetGlobalDomain());
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(domain)});
    }
    void operator()(input::Variable<char>& var) {
        std::vector<char> data = ReadBinaryData<char>(num_elements_, file_name_);
        Hyperslab domain = TransformGeoDomainToHyperslab(var.GetGlobalDomain());
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(domain)});
    }
    void operator()(input::Variable<int16_t>& var) {
        std::vector<int16_t> data = ReadBinaryData<int16_t>(num_elements_, file_name_);
        Hyperslab domain = TransformGeoDomainToHyperslab(var.GetGlobalDomain());
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(domain)});
    }
    void operator()(input::Variable<int32_t>& var) {
        std::vector<int32_t> data = ReadBinaryData<int32_t>(num_elements_, file_name_);
        Hyperslab domain = TransformGeoDomainToHyperslab(var.GetGlobalDomain());
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(domain)});
    }
    void operator()(input::Variable<float>& var) {
        std::vector<float> data = ReadBinaryData<float>(num_elements_, file_name_);
        Hyperslab domain = TransformGeoDomainToHyperslab(var.GetGlobalDomain());
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(domain)});
    }
    void operator()(input::Variable<double>& var) {
        std::vector<double> data = ReadBinaryData<double>(num_elements_, file_name_);
        Hyperslab domain = TransformGeoDomainToHyperslab(var.GetGlobalDomain());
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(domain)});
    }
    void operator()(input::Variable<uint8_t>& var) {
        std::vector<uint8_t> data = ReadBinaryData<uint8_t>(num_elements_, file_name_);
        Hyperslab domain = TransformGeoDomainToHyperslab(var.GetGlobalDomain());
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(domain)});
    }
    void operator()(input::Variable<uint16_t>& var) {
        std::vector<uint16_t> data = ReadBinaryData<uint16_t>(num_elements_, file_name_);
        Hyperslab domain = TransformGeoDomainToHyperslab(var.GetGlobalDomain());
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(domain)});
    }
    void operator()(input::Variable<uint32_t>& var) {
        std::vector<uint32_t> data = ReadBinaryData<uint32_t>(num_elements_, file_name_);
        Hyperslab domain = TransformGeoDomainToHyperslab(var.GetGlobalDomain());
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(domain)});
    }
    void operator()(input::Variable<int64_t>& var) {
        std::vector<int64_t> data = ReadBinaryData<int64_t>(num_elements_, file_name_);
        Hyperslab domain = TransformGeoDomainToHyperslab(var.GetGlobalDomain());
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(domain)});
    }
    void operator()(input::Variable<uint64_t>& var) {
        std::vector<uint64_t> data = ReadBinaryData<uint64_t>(num_elements_, file_name_);
        Hyperslab domain = TransformGeoDomainToHyperslab(var.GetGlobalDomain());
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(domain)});
    }
private:
    const size_t num_elements_;
    const std::string& file_name_;
};


template<typename T>
inline std::vector<T>
ReadBinaryDataFromSubDomain(const std::string& file_name, const DataLayout data_layout, const Hyperslab& global_domain, const Hyperslab& sub_domain, const T missing_value)
{
    /* Get the number of elements that will b read */
    const HyperslabIndex num_elements = sub_domain.GetNumberCoordinates();

    /* Get a function with determines indices in the linearized global data to extract in order to obtain the data from the sub-domain */
    LinearIndicesExtractionFn sub_domain_indices_extraction_fn = GetIndicesExtractionFunction(data_layout);

    /* Get the indices to extract the correct data */
    std::vector<HyperslabIndex> sub_domain_indices = sub_domain_indices_extraction_fn(global_domain, sub_domain);

    cmc_assert(static_cast<size_t>(num_elements) == sub_domain_indices.size());

    /* Allocate memory for the data */
    std::vector<T> values;
    values.reserve(num_elements);

    /* Open the file for input */
    std::ifstream in(file_name, std::ios::binary);

    if (!in.is_open())
    {
        cmc_err_msg("The file ", file_name, " containing the binary data could not be opened.");
    }

    /* Get the size of the data */
    const HyperslabIndex type_size = static_cast<HyperslabIndex>(sizeof(T));

    T value;

    /* Get the data from the file */
    for (HyperslabIndex idx = 0; idx < num_elements; ++idx)
    {
        /* Check if there is a valid index, if not we set a missing value */
        if (sub_domain_indices[idx] == kOutsideOfHyperslabDomain)
        {
            values.push_back(missing_value);
            continue;
        }

        /* Move to the position of the index within the file */
        in.seekg(sub_domain_indices[idx] * type_size, std::ios::beg);

        /* Get the data value at this position */
        in.read(reinterpret_cast<char*>(&value), sizeof(T));

        /* Store the value */
        values.push_back(value);
    }

    //TODO: Endianness reordering potentially

    return values;
}

struct ReadSubDomainData
{
public:
    ReadSubDomainData() = delete;
    ReadSubDomainData(const std::string& file_name, const DataLayout layout, const GeoDomain& global_domain, const GeoDomain& sub_domain)
    : file_name_{file_name}, layout_{layout}, global_domain_{global_domain}, sub_domain_{sub_domain} {};

    void operator()(input::Variable<int8_t>& var) {
        const Hyperslab global_domain = TransformGeoDomainToHyperslab(global_domain_);
        Hyperslab sub_domain = TransformGeoDomainToHyperslab(sub_domain_);
        std::vector<int8_t> data = ReadBinaryDataFromSubDomain<int8_t>(file_name_, layout_, global_domain, sub_domain, var.GetMissingValue());
        sub_domain.NullifyStartIndices();
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(sub_domain)});
        var.SetGlobalDomain(sub_domain_.GetZeroOffsetDomain());
    }
    void operator()(input::Variable<char>& var) {
        const Hyperslab global_domain = TransformGeoDomainToHyperslab(global_domain_);
        Hyperslab sub_domain = TransformGeoDomainToHyperslab(sub_domain_);
        std::vector<char> data = ReadBinaryDataFromSubDomain<char>(file_name_, layout_, global_domain, sub_domain, var.GetMissingValue());
        sub_domain.NullifyStartIndices();
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(sub_domain)});
        var.SetGlobalDomain(sub_domain_.GetZeroOffsetDomain());
    }
    void operator()(input::Variable<int16_t>& var) {
        const Hyperslab global_domain = TransformGeoDomainToHyperslab(global_domain_);
        Hyperslab sub_domain = TransformGeoDomainToHyperslab(sub_domain_);
        std::vector<int16_t> data = ReadBinaryDataFromSubDomain<int16_t>(file_name_, layout_, global_domain, sub_domain, var.GetMissingValue());
        sub_domain.NullifyStartIndices();
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(sub_domain)});
        var.SetGlobalDomain(sub_domain_.GetZeroOffsetDomain());
    }
    void operator()(input::Variable<int32_t>& var) {
        const Hyperslab global_domain = TransformGeoDomainToHyperslab(global_domain_);
        Hyperslab sub_domain = TransformGeoDomainToHyperslab(sub_domain_);
        std::vector<int32_t> data = ReadBinaryDataFromSubDomain<int32_t>(file_name_, layout_, global_domain, sub_domain, var.GetMissingValue());
        sub_domain.NullifyStartIndices();
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(sub_domain)});
        var.SetGlobalDomain(sub_domain_.GetZeroOffsetDomain());
    }
    void operator()(input::Variable<float>& var) {
        const Hyperslab global_domain = TransformGeoDomainToHyperslab(global_domain_);
        Hyperslab sub_domain = TransformGeoDomainToHyperslab(sub_domain_);
        std::vector<float> data = ReadBinaryDataFromSubDomain<float>(file_name_, layout_, global_domain, sub_domain, var.GetMissingValue());
        sub_domain.NullifyStartIndices();
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(sub_domain)});
        var.SetGlobalDomain(sub_domain_.GetZeroOffsetDomain());
    }
    void operator()(input::Variable<double>& var) {
        const Hyperslab global_domain = TransformGeoDomainToHyperslab(global_domain_);
        Hyperslab sub_domain = TransformGeoDomainToHyperslab(sub_domain_);
        std::vector<double> data = ReadBinaryDataFromSubDomain<double>(file_name_, layout_, global_domain, sub_domain, var.GetMissingValue());
        sub_domain.NullifyStartIndices();
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(sub_domain)});
        var.SetGlobalDomain(sub_domain_.GetZeroOffsetDomain());
    }
    void operator()(input::Variable<uint8_t>& var) {
        const Hyperslab global_domain = TransformGeoDomainToHyperslab(global_domain_);
        Hyperslab sub_domain = TransformGeoDomainToHyperslab(sub_domain_);
        std::vector<uint8_t> data = ReadBinaryDataFromSubDomain<uint8_t>(file_name_, layout_, global_domain, sub_domain, var.GetMissingValue());
        sub_domain.NullifyStartIndices();
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(sub_domain)});
        var.SetGlobalDomain(sub_domain_.GetZeroOffsetDomain());
    }
    void operator()(input::Variable<uint16_t>& var) {
        const Hyperslab global_domain = TransformGeoDomainToHyperslab(global_domain_);
        Hyperslab sub_domain = TransformGeoDomainToHyperslab(sub_domain_);
        std::vector<uint16_t> data = ReadBinaryDataFromSubDomain<uint16_t>(file_name_, layout_, global_domain, sub_domain, var.GetMissingValue());
        sub_domain.NullifyStartIndices();
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(sub_domain)});
        var.SetGlobalDomain(sub_domain_.GetZeroOffsetDomain());
    }
    void operator()(input::Variable<uint32_t>& var) {
        const Hyperslab global_domain = TransformGeoDomainToHyperslab(global_domain_);
        Hyperslab sub_domain = TransformGeoDomainToHyperslab(sub_domain_);
        std::vector<uint32_t> data = ReadBinaryDataFromSubDomain<uint32_t>(file_name_, layout_, global_domain, sub_domain, var.GetMissingValue());
        sub_domain.NullifyStartIndices();
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(sub_domain)});
        var.SetGlobalDomain(sub_domain_.GetZeroOffsetDomain());
    }
    void operator()(input::Variable<int64_t>& var) {
        const Hyperslab global_domain = TransformGeoDomainToHyperslab(global_domain_);
        Hyperslab sub_domain = TransformGeoDomainToHyperslab(sub_domain_);
        std::vector<int64_t> data = ReadBinaryDataFromSubDomain<int64_t>(file_name_, layout_, global_domain, sub_domain, var.GetMissingValue());
        sub_domain.NullifyStartIndices();
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(sub_domain)});
        var.SetGlobalDomain(sub_domain_.GetZeroOffsetDomain());
    }
    void operator()(input::Variable<uint64_t>& var) {
        const Hyperslab global_domain = TransformGeoDomainToHyperslab(global_domain_);
        Hyperslab sub_domain = TransformGeoDomainToHyperslab(sub_domain_);
        std::vector<uint64_t> data = ReadBinaryDataFromSubDomain<uint64_t>(file_name_, layout_, global_domain, sub_domain, var.GetMissingValue());
        sub_domain.NullifyStartIndices();
        var.SetDataAndCoordinates(std::move(data), std::vector<Hyperslab>{std::move(sub_domain)});
        var.SetGlobalDomain(sub_domain_.GetZeroOffsetDomain());
    }
private:
    const std::string& file_name_;
    const DataLayout layout_;
    const GeoDomain& global_domain_;
    const GeoDomain& sub_domain_;
};

inline input::Var
Reader::CreateVariableFromBinaryData(const CmcType type, const std::string& name, const int id, const size_t num_elements, const CmcUniversalType missing_value,
                                     const DataLayout layout, const GeoDomain& domain) const
{
    if (static_cast<size_t>(domain.GetNumberReferenceCoordsCovered()) != num_elements)
    {
        cmc_err_msg("The amount of elements to be read does not match the number of (reference) elements in the global domain.");
    }

    /* Create a variable with the given specifications */
    input::Var variable(type, name, id, num_elements, missing_value, layout, domain);

    /* Read the data into the variable */
    std::visit(ReadData(num_elements, this->file_name_), variable.GetInternalVariant(AccessKey()));

    cmc_debug_msg("The binary reader created the input::Variable ", name, ".");

    return variable;
}

inline input::Var
Reader::CreateSubDomainVariableFromBinaryData(const CmcType type, const std::string& name, const int id, const CmcUniversalType missing_value,
                                              const DataLayout layout, const GeoDomain& global_domain, const GeoDomain& sub_domain) const
{
    /* Number of elements the variable will hold */
    const size_t num_elements = sub_domain.GetNumberReferenceCoordsCovered();

    /* Create a variable with the given specifications */
    input::Var variable(type, name, id, num_elements, missing_value, layout, sub_domain);

    /* Read the data into the variable */
    std::visit(ReadSubDomainData(this->file_name_, layout, global_domain, sub_domain), variable.GetInternalVariant(AccessKey()));

    cmc_debug_msg("The binary reader created the input::Variable ", name, ".");

    return variable;
}

}

#endif /* !INPUT_CMC_BINARY_READER_HXX */
