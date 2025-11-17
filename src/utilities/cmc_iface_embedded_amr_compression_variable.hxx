#ifndef CMC_IFACE_EMBEDDED_AMR_COMPRESSION_VARIABLE_HXX
#define CMC_IFACE_EMBEDDED_AMR_COMPRESSION_VARIABLE_HXX

#include "cmc.hxx"
#include "utilities/cmc_iface_amr_compression_variable.hxx"
#include "utilities/cmc_embedded_variable_attributes.hxx"

namespace cmc
{

template <typename T>
class IEmbeddedAMRCompressionVariable : public IAMRCompressionVariable<T>
{
public:
    virtual const VariableAttributes<T>& GetVariableAttributes() const = 0;
    virtual bool AreMeshRefinementBitsStored() const = 0;
};

}

#endif /* !CMC_IFACE_EMBEDDED_AMR_COMPRESSION_VARIABLE_HXX */
