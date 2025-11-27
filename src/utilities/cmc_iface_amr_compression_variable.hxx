#ifndef CMC_IFACE_AMR_COMPRESSION_VARIABLE_HXX
#define CMC_IFACE_AMR_COMPRESSION_VARIABLE_HXX

#include "cmc.hxx"
#include "utilities/cmc_iface_compression_variable.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "utilities/cmc_compression_schema.hxx"

#ifdef CMC_ENABLE_MPI
#include "mpi/cmc_mpi.hxx"
#endif

#include <cstdint>
#include <vector>

namespace cmc
{

template <typename T>
class IAMRCompressionVariable : public ICompressionVariable<T>
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
    virtual CompressionSchema GetCompressionSchema() const = 0;
#ifdef CMC_ENABLE_MPI
    virtual MPI_Comm GetMPIComm() const = 0;
#endif
};

}

#endif /* !CMC_IFACE_AMR_COMPRESSION_VARIABLE_HXX */