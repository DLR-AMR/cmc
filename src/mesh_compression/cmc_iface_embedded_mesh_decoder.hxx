#ifndef CMC_IFACE_EMBEDDED_MESH_DECODER_HXX
#define CMC_IFACE_EMBEDDED_MESH_DECODER_HXX

#include "mesh_compression/cmc_iface_abstract_mesh_decoder.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "utilities/cmc_geo_domain.hxx"

#include "mpi/cmc_mpi.hxx"

#include <vector>
#include <tuple>
#include <utility>

namespace cmc::mesh_compression
{

class IEmbeddedMeshDecoder : public IAbstractMeshDecoder
{
public:

    std::pair<AmrMesh, GeoDomain> DecodeRootMesh(MPI_Comm comm);

    virtual ~IEmbeddedMeshDecoder(){};

protected:
    IEmbeddedMeshDecoder() = delete;
    IEmbeddedMeshDecoder(const std::vector<uint8_t>& encoded_mesh_byte_stream)
    : IAbstractMeshDecoder(encoded_mesh_byte_stream) {};

    /* Return the decoded root level mesh as well as the processes bytes needed for decoding the mesh */
    virtual std::tuple<AmrMesh, GeoDomain, size_t> DecodeRootLevelMesh(const std::vector<uint8_t>& encoded_mesh_stream_, MPI_Comm comm) = 0;

private:

};

inline std::pair<AmrMesh, GeoDomain>
IEmbeddedMeshDecoder::DecodeRootMesh(MPI_Comm comm)
{
    /* Decode the root level mesh */
    auto [root_mesh, global_domain, num_processed_bytes] = DecodeRootLevelMesh(this->GetEncodedMeshStreamData(), comm);

    /* Store the offset for the already processed bytes of the root mesh */
    this->AddToProcessedBytes(num_processed_bytes);

    /* Return the decoded root mesh and the global domain */
    return std::make_pair(std::move(root_mesh), std::move(global_domain));
}

}


#endif /* !CMC_IFACE_EMBEDDED_MESH_DECODER_HXX */
