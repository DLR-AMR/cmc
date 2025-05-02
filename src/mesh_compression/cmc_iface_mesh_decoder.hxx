#ifndef CMC_IFACE_MESH_DECODER_HXX
#define CMC_IFACE_MESH_DECODER_HXX

#include "mesh_compression/cmc_iface_abstract_mesh_decoder.hxx"
#include "t8code/cmc_t8_mesh.hxx"

#include <vector>
#include <utility>

namespace cmc::mesh_compression
{

class IMeshDecoder : public IAbstractMeshDecoder
{
public:

    t8_forest_t DecodeRootMesh();

    virtual ~IMeshDecoder(){};

protected:
    IMeshDecoder() = delete;
    IMeshDecoder(const std::vector<uint8_t>& encoded_mesh_byte_stream)
    : IAbstractMeshDecoder(encoded_mesh_byte_stream) {};

    /* Return the decoded root level mesh as well as the processes bytes needed for decoding the mesh */
    virtual std::pair<t8_forest_t, size_t> DecodeRootLevelMesh(const std::vector<uint8_t>& encoded_mesh_stream_) = 0;

private:

};

inline t8_forest_t
IMeshDecoder::DecodeRootMesh()
{
    /* Decode the root level mesh */
    auto [root_mesh, num_processed_bytes] = DecodeRootLevelMesh(this->GetEncodedMeshStreamData());

    /* Store the offset for the already processed bytes of the root mesh */
    this->AddToProcessedBytes(num_processed_bytes);

    /* Return the decoded root mesh */
    return root_mesh;
}

}


#endif /* !CMC_IFACE_MESH_DECODER_HXX */
