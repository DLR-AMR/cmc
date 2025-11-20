#ifndef CMC_IFACE_PATCH_DECOMPRESSION_VARIABLE_HXX
#define CMC_IFACE_PATCH_DECOMPRESSION_VARIABLE_HXX

#include "utilities/cmc_iface_decompression_variable.hxx"

namespace cmc::patch
{

template <typename T>
class IPatchDecompressionVariable : public IDecompressionVariable<T>
{
public:
    virtual const std::vector<T>& GetDecompressedData() const = 0;

    virtual const std::string& GetName() const = 0;
    virtual size_t Size() const = 0;
    virtual size_t GetDimensionality() const = 0;
    virtual const std::vector<size_t>& GetInitialDimensionLengths() const = 0;
    virtual const std::vector<std::vector<size_t>>& GetDimensionLengthPyramid() const = 0;
    virtual DataLayout GetInitialDataLayout() const = 0;
    virtual GeoDomain GetInitialDomain() const = 0;
};

}


#endif /* !CMC_IFACE_PATCH_DECOMPRESSION_VARIABLE_HXX */
