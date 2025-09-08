#ifndef CMC_IFACE_ABSTRACT_EMBEDDED_BYTE_DECOMPRESSION_VARIABLE_HXX
#define CMC_IFACE_ABSTRACT_EMBEDDED_BYTE_DECOMPRESSION_VARIABLE_HXX

#include "cmc.hxx"
#include "mpi/cmc_mpi.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "utilities/cmc_compression_schema.hxx"
#include "utilities/cmc_byte_compression_values.hxx"

namespace cmc
{

template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

template <typename T>
class IEmbeddedByteDecompressionVariable
{
public:
    virtual void Decompress() = 0;
    virtual const std::string& GetName() const = 0;
    virtual size_t Size() const = 0;
    virtual const AmrMesh& GetAmrMesh() const = 0;
    virtual const std::vector<CompressionValue<T>>& GetDecompressedData() const = 0;
    virtual MPI_Comm GetMPIComm() const = 0;
    virtual std::vector<T> DeMortonizeData() const = 0;
};

    
}


#endif /* !CMC_IFACE_ABSTRACT_EMBEDDED_BYTE_DECOMPRESSION_VARIABLE_HXX */
