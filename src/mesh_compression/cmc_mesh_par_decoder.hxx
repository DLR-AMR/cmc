#ifndef CMC_MESH_PAR_DECODER_HXX
#define CMC_MESH_PAR_DECODER_HXX

#include "mesh_compression/cmc_iface_abstract_mesh_par_decoder.hxx"

#include <cstdint>

namespace cmc::mesh_compression
{

class MeshParDecoder : public IAbstractMeshParDecoder
{
public: 
    MeshParDecoder() = delete;
    MeshParDecoder(const uint8_t* encoded_mesh_lvl_byte_stream)
    : IAbstractMeshParDecoder(encoded_mesh) {}

private:
    std::tuple<t8_forest_t, size_t, int> DecodeRootLevel(const uint8_t* encoded_mesh_stream_, const t8_cmesh_t cmesh, const t8_scheme *scheme) override;

};

inline std::tuple<t8_forest_t, size_t, int>
MeshParDecoder::DecodeRootLevel(const uint8_t* encoded_mesh_stream_, const t8_cmesh_t cmesh, const t8_scheme *scheme)
{
    size_t processed_bytes = 0;

    /* The number bytes describing the mesh is given at the front */
    const uint64_t num_root_elements = GetValueFromByteStream<uint64_t>(encoded_mesh_stream.data());
    processed_bytes += sizeof(uint64_t);

    cmc_assert(num_root_elements == static_cast<uint64_t>(t8_cmesh_get_num_trees(cmesh)));

    /* Get the dimension of the cmesh */
    const int dim = t8_cmesh_get_dimension (cmesh);

    /* There is nothin more encoded on the root level */
    return std::make_tuple(base_mesh, processed_bytes, dim);
}

}


#endif /* !CMC_MESH_PAR_DECODER_HXX */
