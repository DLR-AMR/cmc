#ifndef CMC_MESH_ENCODER_HXX
#define CMC_MESH_ENCODER_HXX

#include "mesh_compression/cmc_iface_mesh_encoder.hxx"
#include "utilities/cmc_serialization.hxx"

#include <vector>

namespace cmc::mesh_compression
{

class MeshEncoder : public IMeshEncoder
{
public: 
    MeshEncoder() = default;

    std::vector<uint8_t> EncodeRootLevelMesh(t8_forest_t root_level_mesh) override;    

};


inline std::vector<uint8_t>
MeshEncoder::EncodeRootLevelMesh(t8_forest_t root_level_mesh)
{
    size_t num_bytes_root_mesh_encoding_ = 1;

    std::vector<uint8_t> encoded_root_level;
    encoded_root_level.reserve(sizeof(size_t) + 1);

    /* Push the number of root level bytes to the stream */
    PushBackValueToByteStream(encoded_root_level, num_bytes_root_mesh_encoding_);
    encoded_root_level.push_back(0);

    return encoded_root_level;
}

}


#endif /* !CMC_MESH_ENCODER_HXX */