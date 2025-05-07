#ifndef CMC_EMBEDDED_MESH_DECODER_HXX
#define CMC_EMBEDDED_MESH_DECODER_HXX

#include "mesh_compression/cmc_iface_embedded_mesh_decoder.hxx"
#include "utilities/cmc_serialization.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_embedded_mesh_utilities.hxx"

#ifdef CMC_WITH_T8CODE
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#endif

#include <array>
#include <vector>
#include <tuple>

namespace cmc::mesh_compression
{

class EmbeddedMeshDecoder : public IEmbeddedMeshDecoder
{
public: 
    EmbeddedMeshDecoder() = delete;
    EmbeddedMeshDecoder(const std::vector<uint8_t>& encoded_mesh)
    : IEmbeddedMeshDecoder(encoded_mesh) {};

protected:
    std::tuple<AmrMesh, GeoDomain, size_t> DecodeRootLevelMesh(const std::vector<uint8_t>& encoded_mesh_stream_, MPI_Comm comm) override;    
};

inline std::tuple<AmrMesh, GeoDomain, size_t>
EmbeddedMeshDecoder::DecodeRootLevelMesh(const std::vector<uint8_t>& encoded_mesh_stream, MPI_Comm comm)
{
    cmc_debug_msg("The embedded root level mesh will be decoded.");
    /* The root level mesh is encoded at the beginning of the encoded stream */
    size_t processed_bytes = 0;

    /* The number bytes describing the mesh is given at the front */
    [[maybe_unused]] const uint32_t num_bytes_mesh_root_mesh_encoding = GetValueFromByteStream<uint32_t>(encoded_mesh_stream.data());
    processed_bytes += sizeof(uint32_t);

    cmc_debug_msg("The embedded root level mesh is encoded within ", num_bytes_mesh_root_mesh_encoding, " bytes.");
    
    /* Get the diemnsionality */
    const uint32_t dimensionality = GetValueFromByteStream<uint32_t>(encoded_mesh_stream.data() + processed_bytes);
    processed_bytes += sizeof(uint32_t);

    /* Get the initial refinement level */
    const uint32_t initial_refinement_level = GetValueFromByteStream<uint32_t>(encoded_mesh_stream.data() + processed_bytes);
    processed_bytes += sizeof(uint32_t);

    /* Fill the start and end indices arrays from the global domain */
    std::array<DomainIndex, Dimension::NumCoordinates> start_indices;
    std::array<DomainIndex, Dimension::NumCoordinates> end_indices;

    /* Get the start indices */
    for (auto idx = 0; idx < Dimension::NumCoordinates; ++idx)
    {
        const uint64_t value = GetValueFromByteStream<uint64_t>(encoded_mesh_stream.data() + processed_bytes);
        processed_bytes += sizeof(uint64_t);
        start_indices[idx] = static_cast<DomainIndex>(value);
    }

    /* Get the end indices */
    for (auto idx = 0; idx < Dimension::NumCoordinates; ++idx)
    {
        const uint64_t value = GetValueFromByteStream<uint64_t>(encoded_mesh_stream.data() + processed_bytes);
        processed_bytes += sizeof(uint64_t);
        end_indices[idx] = static_cast<DomainIndex>(value);
    }

    /* Reconstruct the global domain from the start and end indices */
    const GeoDomain global_domain = ReconstructGeoDomainFromStartAndEndIndices(start_indices, end_indices);

    /* Create the embedded root level mesh */
    t8_cmesh_t cmesh = t8_cmesh_new_hypercube(EmbeddedMeshDimensionToElementClass(dimensionality), comm, 0, 0, 1);

    /* Construct a forest from the cmesh */
    t8_forest_t base_mesh;
    t8_forest_init(&base_mesh);
    t8_forest_set_cmesh(base_mesh, cmesh, comm);
    t8_forest_set_scheme(base_mesh, t8_scheme_new_default());
    t8_forest_set_level(base_mesh, 0);
    t8_forest_commit(base_mesh);

    /* Create an AmrMesh */
    const AmrMesh root_level_mesh(base_mesh, static_cast<int>(initial_refinement_level), static_cast<int>(dimensionality));

    cmc_debug_msg("The embedded root level mesh has been decoded and reconstructed.");

    return std::make_tuple(root_level_mesh, global_domain, processed_bytes);
}

}


#endif /* !CMC_EMBEDDED_MESH_DECODER_HXX */
