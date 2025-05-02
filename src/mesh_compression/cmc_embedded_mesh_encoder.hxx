#ifndef CMC_EMBEDDED_MESH_ENCODER_HXX
#define CMC_EMBEDDED_MESH_ENCODER_HXX

#include "mesh_compression/cmc_iface_embedded_mesh_encoder.hxx"
#include "utilities/cmc_serialization.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "utilities/cmc_geo_domain.hxx"

#include <vector>

namespace cmc::mesh_compression
{

class EmbeddedMeshEncoder : public IEmbeddedMeshEncoder
{
public: 
    EmbeddedMeshEncoder() = default;

    std::vector<uint8_t> EncodeRootLevelMesh(AmrMesh root_level_mesh, const GeoDomain& global_domain) override;    
};

inline std::vector<uint8_t>
EmbeddedMeshEncoder::EncodeRootLevelMesh(AmrMesh root_level_mesh, const GeoDomain& global_domain)
{
    cmc_assert(root_level_mesh.GetDimensionality() == global_domain.GetDimensionality());

    uint32_t num_bytes_root_mesh_encoding_ = 3 * sizeof(uint32_t) + 2 * Dimension::NumCoordinates * sizeof(uint64_t);

    std::vector<uint8_t> encoded_root_level;
    encoded_root_level.reserve(num_bytes_root_mesh_encoding_);

    /* Push the number of root level bytes to the stream */
    PushBackValueToByteStream<uint32_t>(encoded_root_level, num_bytes_root_mesh_encoding_);

    /* Push back the dimensionality */
    const uint32_t dimensionality = static_cast<uint32_t>(global_domain.GetDimensionality());
    PushBackValueToByteStream<uint32_t>(encoded_root_level, dimensionality);

    /* Push back the initial refinement level */
    const uint32_t initial_refinement_level = static_cast<uint32_t>(root_level_mesh.GetInitialRefinementLevel()); 
    PushBackValueToByteStream<uint32_t>(encoded_root_level, initial_refinement_level);

    /* Now, we push back the start indices of the global domain */
    for (auto start_idx_iter = global_domain.BeginStartIndices(); start_idx_iter != global_domain.EndStartIndices(); ++start_idx_iter)
    {
        const uint64_t val = static_cast<uint64_t>(*start_idx_iter);
        PushBackValueToByteStream<uint64_t>(encoded_root_level, val);
    }

    /* Afterwards, we push back the end indices of the global domain */
    for (auto end_idx_iter = global_domain.BeginEndIndices(); end_idx_iter != global_domain.EndEndIndices(); ++end_idx_iter)
    {
        const uint64_t val = static_cast<uint64_t>(*end_idx_iter);
        PushBackValueToByteStream<uint64_t>(encoded_root_level, val);
    }

    return encoded_root_level;
}

}


#endif /* !CMC_EMBEDDED_MESH_ENCODER_HXX */