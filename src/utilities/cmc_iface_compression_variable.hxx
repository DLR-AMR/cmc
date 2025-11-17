#ifndef CMC_IFACE_COMPRESSION_VARIABLE_HXX
#define CMC_IFACE_COMPRESSION_VARIABLE_HXX

#include "cmc.hxx"

namespace cmc
{

template <typename T>
class ICompressionVariable
{
public:
    virtual void Compress() = 0;
};

}

#endif /* !CMC_IFACE_COMPRESSION_VARIABLE_HXX */
