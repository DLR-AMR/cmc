#ifndef CMC_IFACE_ABSTRACT_BYTE_COMPRESSION_VARIABLE_HXX
#define CMC_IFACE_ABSTRACT_BYTE_COMPRESSION_VARIABLE_HXX

#include "cmc.hxx"
#include "mpi/cmc_mpi.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "utilities/cmc_compression_schema.hxx"

#include <cstdint>
#include <vector>

namespace cmc
{

template <typename T>
class IByteCompressionVariable
{
public:
    virtual void MoveEncodedEntropyCodesInto(std::vector<std::vector<uint8_t>>& vec_to_hold_encoded_levelwise_entropy_codes) = 0;
    virtual void MoveEncodedDataInto(std::vector<std::vector<uint8_t>>& vec_to_hold_encoded_levelwise_data) = 0;
    virtual void MoveEncodedMeshInto(std::vector<std::vector<uint8_t>>& vec_to_hold_encoded_levelwise_mesh) = 0;

    virtual const std::vector<std::vector<uint8_t>>& GetEncodedEntropyCodes() const = 0;
    virtual const std::vector<std::vector<uint8_t>>& GetEncodedData() const = 0;
    virtual const std::vector<std::vector<uint8_t>>& GetEncodedMesh() const = 0;

    virtual const std::string& GetName() const = 0;
    virtual size_t Size() const = 0;
    virtual const AmrMesh& GetAmrMesh() const = 0;
    virtual MPI_Comm GetMPIComm() const = 0;
    virtual CompressionSchema GetCompressionSchema() const = 0;

    virtual void Compress() = 0;
};

}

#endif /* !CMC_IFACE_ABSTRACT_BYTE_COMPRESSION_VARIABLE_HXX */
