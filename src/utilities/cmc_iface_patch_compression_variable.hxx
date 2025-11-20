#ifndef CMC_IFACE_PATCH_COMPRESSION_VARIABLE_HXX
#define CMC_IFACE_PATCH_COMPRESSION_VARIABLE_HXX

#include "cmc.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_iface_compression_variable.hxx"
#include "utilities/cmc_compression_schema.hxx"

#include <cstdint>
#include <vector>

namespace cmc
{

template <typename T>
class IPatchCompressionVariable : public ICompressionVariable<T>
{
public:
    virtual void MoveEncodedEntropyCodesInto(std::vector<std::vector<uint8_t>>& vec_to_hold_encoded_levelwise_entropy_codes) = 0;
    virtual void MoveEncodedDataInto(std::vector<std::vector<uint8_t>>& vec_to_hold_encoded_levelwise_data) = 0;

    virtual const std::vector<std::vector<uint8_t>>& GetEncodedEntropyCodes() const = 0;
    virtual const std::vector<std::vector<uint8_t>>& GetEncodedData() const = 0;

    virtual const std::string& GetName() const = 0;
    virtual size_t Size() const = 0;
    virtual size_t GetDimensionality() const = 0;
    virtual std::vector<size_t> GetInitialDimensionLengths() const = 0;
    virtual std::vector<std::vector<size_t>> GetDimensionLengthPyramid() const = 0;
    virtual CompressionSchema GetCompressionSchema() const = 0;
    virtual DataLayout GetInitialDataLayout() const = 0;
    virtual GeoDomain GetInitialDomain() const = 0;
};

}

#endif /* !CMC_IFACE_PATCH_COMPRESSION_VARIABLE_HXX */