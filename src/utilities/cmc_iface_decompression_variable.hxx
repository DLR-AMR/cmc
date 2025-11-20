#ifndef CMC_IFACE_DECOMPRESSION_VARIABLE_HXX
#define CMC_IFACE_DECOMPRESSION_VARIABLE_HXX

#include "cmc.hxx"

namespace cmc
{

template <typename T>
class IDecompressionVariable
{
public:
    virtual void Decompress() = 0;
};

}

#endif /* !CMC_IFACE_DECOMPRESSION_VARIABLE_HXX */
