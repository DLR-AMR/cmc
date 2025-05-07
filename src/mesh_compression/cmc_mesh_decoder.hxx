#ifndef CMC_MESH_DECODER_HXX
#define CMC_MESH_DECODER_HXX

#include "mesh_compression/cmc_iface_mesh_decoder.hxx"
#include "utilities/cmc_serialization.hxx"

#include <vector>
#include <utility>

#ifdef CMC_WITH_T8CODE
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#endif

namespace cmc::mesh_compression
{

class MeshDecoder : public IMeshDecoder
{
public: 
    MeshDecoder() = delete;
    MeshDecoder(const std::vector<uint8_t>& encoded_mesh)
    : IMeshDecoder(encoded_mesh) {};

protected:
    std::pair<t8_forest_t, size_t> DecodeRootLevelMesh(const std::vector<uint8_t>& encoded_mesh_stream_) override;

};


inline std::pair<t8_forest_t, size_t>
MeshDecoder::DecodeRootLevelMesh(const std::vector<uint8_t>& encoded_mesh_stream)
{
    //TODO: Make cmesh encoding/decoding

    /* The root level mesh is encoded at the beginning of the encoded stream */
    size_t processed_bytes = 0;

    /* The number bytes describing the mesh is given at the front */
    const size_t num_bytes_mesh_root_mesh_encoding = GetValueFromByteStream<size_t>(encoded_mesh_stream.data());
    processed_bytes += sizeof(size_t);

    /* Reconstruct the base mesh */
    //TODO: WE need an implementation for the cmesh storage
    const int dimensionality = 2;
    /* Get a base cmesh corresponding to the dimensionality */
    //t8_cmesh_t cmesh = t8_cmesh_new_hypercube(T8_ECLASS_TRIANGLE, MPI_COMM_WORLD, 0, 0, 0);
    t8_cmesh_t cmesh = t8_cmesh_new_hypercube(T8_ECLASS_TRIANGLE, MPI_COMM_WORLD, 0, 0, 0);

    /* Construct a forest from the cmesh */
    t8_forest_t base_mesh;
    t8_forest_init(&base_mesh);
    t8_forest_set_cmesh(base_mesh, cmesh, MPI_COMM_WORLD);
    t8_forest_set_scheme(base_mesh, t8_scheme_new_default());
    t8_forest_set_level(base_mesh, 0);
    t8_forest_commit(base_mesh);

    /* Compute the overall byte count in which the base/root mesh has been encoded */
    const size_t num_bytes_encoded_root_level_mesh = num_bytes_mesh_root_mesh_encoding + sizeof(size_t);
    
    return std::make_pair(base_mesh, num_bytes_encoded_root_level_mesh);
}

}


#endif /* !CMC_MESH_DECODER_HXX */
